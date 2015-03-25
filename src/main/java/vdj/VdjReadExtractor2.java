package vdj;


import static vdj.VdjRegion.loadRegions;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
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
public class VdjReadExtractor2 {
	
	private Set<String> readIds = new HashSet<String>();
	private Set<String> reads = new HashSet<String>();
	private int numReads = 0;
	private ReverseComplementor rc = new ReverseComplementor();
	private Set<String> vdjKmers;
	
//	private static final int CONSTANT_START = 113256204;
//	private static final int CONSTANT_END = 113422730;
//	private static final int VDJ_START = CONSTANT_END+1;
//	private static final int VDJ_END = 116009954;

	public void extract(String inFile, String outBam, String vdjRegionFile, String constantRegionFile, String vdjFasta, int kmerSize) throws IOException {
		
		List<VdjRegion> vdjRegions = loadRegions(vdjRegionFile);
		List<VdjRegion> constantRegions = null;
		
		if (constantRegionFile != "none") {
			constantRegions = loadRegions(constantRegionFile);
		}
		
		vdjKmers = loadKmers(vdjFasta, kmerSize);
		
		SAMFileReader reader = new SAMFileReader(new File(inFile));
		reader.setValidationStringency(ValidationStringency.SILENT);
	
		SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
		
		SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(
				reader.getFileHeader(), false, new File(outBam)); 
		
		for (VdjRegion region : vdjRegions) {
			CloseableIterator<SAMRecord> iter = reader.queryOverlapping(region.chromosome, region.start, region.stop);
			while (iter.hasNext()) {
				SAMRecord read = iter.next();
				String id = read.getReadName() + "___" + (read.getFirstOfPairFlag() ? "1" : "2");
				if (!readIds.contains(id)) {
					reads.add(read.getReadName());
					readIds.add(id);
					writer.addAlignment(read);
					numReads += 1;
				}
			}
			iter.close();
		}
		
		if (constantRegions != null) {
			for (VdjRegion region : constantRegions) {
				CloseableIterator<SAMRecord> iter = reader.queryOverlapping(region.chromosome, region.start, region.stop);
				while (iter.hasNext()) {
					SAMRecord read = iter.next();
					String id = read.getReadName() + "___" + (read.getFirstOfPairFlag() ? "1" : "2");
					if (!readIds.contains(id)) {
						reads.add(read.getReadName());
						readIds.add(id);
						writer.addAlignment(read);
						numReads += 1;
					}
				}
				iter.close();
			}
		}
		
		System.out.println("Looping over all reads...");
		
		Set<String> bothUnmapped = new HashSet<String>();

		for (SAMRecord read : reader) {
			// Look for mates of previously identified reads
			if (reads.contains(read.getReadName())) {
				String id = read.getReadName() + "___" + (read.getFirstOfPairFlag() ? "1" : "2");
				if (!readIds.contains(id)) {
					readIds.add(id);
					writer.addAlignment(read);
					numReads += 1;
				}
			} else if (read.getReadUnmappedFlag() && read.getMateUnmappedFlag()) {
				// Mark reads containing a vdj kmer
				if (hasVdjKmer(read.getReadString(), kmerSize)) {
					bothUnmapped.add(read.getReadName());
				}
			}
		}
		
		reader.close();
		
		System.out.println("Last pass through unmapped reads");
		reader = new SAMFileReader(new File(inFile));
		reader.setValidationStringency(ValidationStringency.SILENT);
		
		// One last pass through unmapped reads to get unmapped pairs
		CloseableIterator<SAMRecord> iter = reader.queryUnmapped();
		while (iter.hasNext()) {
			SAMRecord read = iter.next();
			if (bothUnmapped.contains(read.getReadName())) {
				writer.addAlignment(read);
				numReads += 1;
			}			
		}
		iter.close();
		
		reader.close();
		writer.close();
		
		System.out.println("Processed: " + numReads + " reads");
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
	
	public static void main(String[] args) throws Exception {
		String input = args[0];
		String output = args[1];
		String vdjRegionsFile = args[2];
		String constantRegionsFile = args[3];
		String vdjFasta = args[4];
		int kmerSize = Integer.parseInt(args[5]);
		
		new VdjReadExtractor2().extract(input, output, vdjRegionsFile, constantRegionsFile, vdjFasta, kmerSize);
	}
}
