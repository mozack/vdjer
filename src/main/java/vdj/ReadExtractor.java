package vdj;

import static vdj.VdjRegion.loadRegions;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;



public class ReadExtractor {
	
	private List<VdjRegion> regions;
	private List<VdjRegion> nonConstantRegions;
//	private Map<String, ReadPair> readPairs = new HashMap<String, ReadPair>();
	Set<String> bcrReadIds = new HashSet<String>();
	Set<String> nonConstantBcrReadIds = new HashSet<String>();
	private BufferedWriter bcrOutput;
	private BufferedWriter unalignedOutput;
	
	private ReverseComplementor rc = new ReverseComplementor();
	private int numReads = 0;

	public void extract(String inFile, String bcrOut, String unmappedOut, String regionFile, String nonConstantRegionFile) throws IOException {
		regions = loadRegions(regionFile);
		nonConstantRegions = loadRegions(nonConstantRegionFile);
		
		SAMFileReader reader = new SAMFileReader(new File(inFile));
		reader.setValidationStringency(ValidationStringency.SILENT);
		
		bcrOutput = new BufferedWriter(new FileWriter(bcrOut, false));
		unalignedOutput = new BufferedWriter(new FileWriter(unmappedOut, false));
		
		// Identify reads that overlap the constant regions
		for (VdjRegion region : regions) {
			CloseableIterator<SAMRecord> iter = reader.queryOverlapping(region.chromosome, region.start, region.stop);
			while (iter.hasNext()) {
				SAMRecord read = iter.next();
				bcrReadIds.add(read.getReadName());
			}
			iter.close();
		}

		// Identify reads that overlap the non-constant regions
		for (VdjRegion region : nonConstantRegions) {
			CloseableIterator<SAMRecord> iter = reader.queryOverlapping(region.chromosome, region.start, region.stop);
			while (iter.hasNext()) {
				SAMRecord read = iter.next();
				nonConstantBcrReadIds.add(read.getReadName());
			}
			iter.close();
		}
		
		int vdjReadCount = 0;
		int unmappedReadCount = 0;
			
		for (SAMRecord read : reader) {
			// Only consider primary alignments
			if (!read.getNotPrimaryAlignmentFlag() && !read.getDuplicateReadFlag() && (read.getFlags() & 0x800) == 0) {
				if (bcrReadIds.contains(read.getReadName())) {
					output(read, bcrOutput);
					vdjReadCount += 1;
				} else if (nonConstantBcrReadIds.contains(read.getReadName()) || read.getReadUnmappedFlag() || read.getMateUnmappedFlag()) {
					output(read, unalignedOutput);
					unmappedReadCount += 1;
				}
			}
			
			if (numReads++ % 10000000 == 0) {
				System.out.println("Processed: " + numReads + " reads");
			}
		}
		
		System.out.println("vdj read count: " + vdjReadCount);
		System.out.println("unmapped read count: " + unmappedReadCount);
		
		reader.close();
		bcrOutput.close();
		unalignedOutput.close();
	
	}
		
	private String getOutputString(boolean negativeStrandFlag, String bases, String qualities) {
		StringBuffer readBuffer = new StringBuffer();
		readBuffer.append(negativeStrandFlag ? "1" : "0");
		
		readBuffer.append(bases);
		readBuffer.append(qualities);
		
		return readBuffer.toString();
	}

	
	private void output(SAMRecord read, BufferedWriter output) throws IOException {
		output.append(getOutputString(false, read.getReadString(), read.getBaseQualityString()));
		
		// We don't know the orientation of the read, so reverse complement it as well.
		output.append(getOutputString(false, rc.reverseComplement(read.getReadString()), rc.reverse(read.getBaseQualityString())));
	}
		
	
	public static void main(String[] args) throws Exception {
		String input = args[0];
		String bcrOutput = args[1];
		String unmappedOutput = args[2];
		String regions = args[3];
		String nonConstantRegions = args[4];
		
		new ReadExtractor().extract(input, bcrOutput, unmappedOutput, regions, nonConstantRegions);
	}
}
