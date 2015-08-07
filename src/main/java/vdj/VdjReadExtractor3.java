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
	private Set<String> secondaryReads = new HashSet<String>();
	private int numReads = 0;
	private ReverseComplementor rc = new ReverseComplementor();	
	private int kmerSize;
	private Set<String> vdjKmers;
	
	private void getReads(String regionFile, SAMFileReader reader, Set<String> reads, Set<String> secondaryReads) throws IOException {
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
	
	private boolean hasVdjKmer(String bases, int kmerSize) {
		if (kmerSize == 0) {
			return true;
		}
		
		for (int i=0; i<=bases.length()-kmerSize; i++) {
			String kmer = bases.substring(i, i+kmerSize);
			if (vdjKmers.contains(kmer)) {
				return true;
			}
		}
		
		return false;
	}
	
	private Set<String> loadKmers(String fasta, int kmerSize) throws IOException {
		Set<String> kmers = new HashSet<String>();
		if (kmerSize > 0) {
			BufferedReader reader = new BufferedReader(new FileReader(fasta));
			
			String id = reader.readLine();
			String bases = null;
			if (id != null) {
				bases = reader.readLine();
			}
			
			while (id != null && bases != null) {
				
				for (int i=0; i<=bases.length()-kmerSize; i++) {
					String kmer = bases.substring(i, i+kmerSize); 
					kmers.add(kmer);
					kmers.add(rc.reverseComplement(kmer));
				}
				
				id = reader.readLine();
				bases = null;
				if (id != null) {
					bases = reader.readLine();
				}
			}
			
			reader.close();
		}
		
		return kmers;
	}

	private void processReadQuery(CloseableIterator<SAMRecord> iter, Set<String> reads, String regionFile) {
		while (iter.hasNext()) {
			SAMRecord read = iter.next();
			if (regionFile.equals("unmapped")) {
				// Only include unmapped reads if both read and mate are unmapped
				if (read.getReadUnmappedFlag() && read.getMateUnmappedFlag()) {
					// If an unmapped read has a VDJ kmer, add it to the primary read set
					if (hasVdjKmer(read.getReadString(), kmerSize)) {
						reads.add(read.getReadName());
					} else {
						secondaryReads.add(read.getReadName());
					}
				}
			} else {
				reads.add(read.getReadName());
			}
		}
		iter.close();
	}
	
	public void extract(String inFile, String out, String constantOut, String vdjRegionFile, String constantRegionFile,
			String vdjFasta, int kmerSize) throws IOException {
				
		loadKmers(vdjFasta, kmerSize);
		
		SAMFileReader reader = new SAMFileReader(new File(inFile));
		reader.setValidationStringency(ValidationStringency.SILENT);
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(out, false));
		BufferedWriter constantWriter = new BufferedWriter(new FileWriter(constantOut, false));
	
		// Constant reads always go to secondary read set
		System.out.println("Retrieving constant reads.");
		getReads(constantRegionFile, reader, secondaryReads, secondaryReads);
		
		// Variable reads always go to primary read set
		System.out.println("Retrieving variable reads");
		getReads(vdjRegionFile, reader, reads, reads);
		
		// Unampped reads can go to either read set depending upon content
		System.out.println("Retrieving unmapped reads");
		getReads("unmapped", reader, reads, secondaryReads);
		
		System.out.println("Looping over all reads...");

		for (SAMRecord read : reader) {
			
			if (secondaryReads.contains(read.getReadName())) {
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
		String vdjFasta = args[5];
		int kmerSize = Integer.parseInt(args[6]);
		
		new VdjReadExtractor3().extract(input, output, constantOut, vdjRegionsFile, constantRegionsFile, vdjFasta, kmerSize);
	}
}
