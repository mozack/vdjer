package vdj;


import static vdj.VdjRegion.loadRegions;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.List;
import java.util.Set;


/**
 * Extracts all reads mapping to specified regions or that are unmapped or have mate unmapped
 * along with their mates.
 * @author lmose
 *
 */
public class VdjReadExtractor3 {
	
	private Set<String> readIds = new HashSet<String>();
	private Set<String> reads = new HashSet<String>();
	private Set<String> constantReads = new HashSet<String>();
	private int numReads = 0;
	private ReverseComplementor rc = new ReverseComplementor();	
	
	private void getReads(String regionFile, SAMFileReader reader, Set<String> reads) throws IOException {
		List<VdjRegion> regions = null;
		
		if (regionFile.equals("unmapped")) {
			CloseableIterator<SAMRecord> iter = reader.queryUnmapped();
			processReadQuery(iter, reads, regionFile);
		} else if (!regionFile.equals("none")) {
			regions = loadRegions(regionFile);
			
			for (VdjRegion region : regions) {
				CloseableIterator<SAMRecord> iter = reader.queryOverlapping(region.chromosome, region.start, region.stop);
				processReadQuery(iter, reads, regionFile);
			}
		}
	}

	private void processReadQuery(CloseableIterator<SAMRecord> iter, Set<String> reads, String regionFile) {
		while (iter.hasNext()) {
			SAMRecord read = iter.next();
			if (regionFile.equals("unmapped")) {
				// Only include unmapped reads if both read and mate are unmapped
				if (read.getReadUnmappedFlag() && read.getMateUnmappedFlag()) {
					reads.add(read.getReadName());;
				}
			} else {
				reads.add(read.getReadName());
			}
		}
		iter.close();
	}
	
	public void extract(String inFile, String out, String constantOut, String vdjRegionFile, String constantRegionFile) throws IOException {
				
		SAMFileReader reader = new SAMFileReader(new File(inFile));
		reader.setValidationStringency(ValidationStringency.SILENT);
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(out, false));
		BufferedWriter constantWriter = new BufferedWriter(new FileWriter(constantOut, false));
	
		System.out.println("Retrieving constant reads.");
		getReads(constantRegionFile, reader, constantReads);
		System.out.println("Retrieving variable reads");
		getReads(vdjRegionFile, reader, reads);
		System.out.println("Retrieving unmapped reads");
		getReads("unmapped", reader, reads);
		
		System.out.println("Looping over all reads...");

		for (SAMRecord read : reader) {
			
			if (constantReads.contains(read.getReadName())) {
				String id = read.getReadName() + "___" + (read.getFirstOfPairFlag() ? "1" : "2");
				
				if (!readIds.contains(id)) {
					readIds.add(id);
					numReads += 1;
					output(read, constantWriter);
				}
			}
			else if (reads.contains(read.getReadName())) {
				String id = read.getReadName() + "___" + (read.getFirstOfPairFlag() ? "1" : "2");
				
				if (!readIds.contains(id)) {
					readIds.add(id);
					numReads += 1;
					output(read, writer);
				}
			}
		}
		
		reader.close();
		writer.close();
		constantWriter.close();
		
		System.out.println("Processed: " + numReads + " reads");
	}
	
	private void output(SAMRecord read, BufferedWriter writer) throws IOException {
		output(read.getReadString(), read.getBaseQualityString(), writer);
	}
	
	private void output(String bases, String quals, BufferedWriter output) throws IOException {
		output.append(getOutputString(false, bases, quals));
		
		// We don't know the orientation of the read, so reverse complement it as well.
		output.append(getOutputString(false, rc.reverseComplement(bases), rc.reverse(quals)));
	}
	
	private String getOutputString(boolean negativeStrandFlag, String bases, String qualities) {
		StringBuffer readBuffer = new StringBuffer();
		readBuffer.append(negativeStrandFlag ? "1" : "0");
		
		readBuffer.append(bases);
		readBuffer.append(qualities);
		
		return readBuffer.toString();
	}
	
	public static void main(String[] args) throws Exception {
		String input = args[0];
		String output = args[1];
		String constantOut = args[2];
		String vdjRegionsFile = args[3];
		String constantRegionsFile = args[4];
		
		new VdjReadExtractor3().extract(input, output, constantOut, vdjRegionsFile, constantRegionsFile);
	}
}
