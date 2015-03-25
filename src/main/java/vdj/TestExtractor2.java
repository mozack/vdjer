package vdj;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class TestExtractor2 {
	
	private int KMER;
	private Set<String> seqKmers = new HashSet<String>();
	private Set<String> readIds = new HashSet<String>();
	private Map<String, Integer> idKmerMap = new HashMap<String, Integer>();
	private Set<String> bothIds = new HashSet<String>();

	public void run(String inFile, String refSeq, int kmer, String outFile) throws IOException {
		
		KMER = kmer;
		
		loadRef(refSeq);
		
		SAMFileReader reader = new SAMFileReader(new File(inFile));
		reader.setValidationStringency(ValidationStringency.SILENT);
		int count = 0;
		int readCount = 0;
		for (SAMRecord read : reader) {
			if (read.getReadUnmappedFlag() && read.getMateUnmappedFlag()) {
				if (hasKmer(read.getReadString())) {
					count += 1;
					readIds.add(read.getReadName());
					
					if (!idKmerMap.containsKey(read.getReadName())) {
						idKmerMap.put(read.getReadName(), 1);
					} else {
						idKmerMap.put(read.getReadName(), 2);
						bothIds.add(read.getReadName());
					}
				}
			}
			
			if ((readCount++ % 1000000) == 0) {
				System.out.println("Processed: " + readCount);
			}
		}
		
		SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
		
		SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(
				reader.getFileHeader(), false, new File(outFile));
		
		reader.close();
		reader = new SAMFileReader(new File(inFile));
		reader.setValidationStringency(ValidationStringency.SILENT);
		
		for (SAMRecord read : reader) {
			if (!read.getReadUnmappedFlag()) {
				writer.addAlignment(read);
			} else if (readIds.contains(read.getReadName())) {
				writer.addAlignment(read);
			}
		}
		
		reader.close();
		writer.close();
		
		System.out.println("Kmer: " + KMER);
		System.out.println("# reads with kmer: " + count);
		System.out.println("# read ids: " + readIds.size());
		System.out.println("# in both reads: " + bothIds.size());
	}
	
	private boolean hasKmer(String seq) {
		boolean isFound = false;
		
		String[] kmers = getKmers(seq);
		for (String kmer : kmers) {
			if (seqKmers.contains(kmer)) {
				isFound = true;
				break;
			}
		}
		/*

		if (!isFound) {
			kmers = getKmers(new ReverseComplementor().reverseComplement(seq));
			for (String kmer : kmers) {
				if (seqKmers.contains(kmer)) {
					isFound = true;
					break;
				}
			}
		}
		*/
		
		return isFound;
	}
	
	private void loadRef(String refSeq) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(refSeq));
		
		String line = reader.readLine();
		StringBuffer buf = new StringBuffer();
		while (line != null) {
			
			if (line.startsWith(">")) {
				if (buf.length() >= KMER) {
					addKmers(buf.toString());
				}
				buf = new StringBuffer();
			}
			
			//TODO: Include masked sequence?
//			buf.append(line.toUpperCase());
			buf.append(line);
			line = reader.readLine();
		}
		
		// Include the last sequence
		if (buf.length() >= KMER) {
			addKmers(buf.toString());
		}
		
		reader.close();
	}
	
	private void addKmers(String seq) {
		String[] kmers = getKmers(seq);
		for (String kmer : kmers) {
			if (!kmer.contains("N")) {
				seqKmers.add(kmer);
			}
		}
		
		kmers = getKmers(new ReverseComplementor().reverseComplement(seq));
		for (String kmer : kmers) {
			if (!kmer.contains("N")) {
				seqKmers.add(kmer);
			}
		}
	}
	
	private String[] getKmers(String seq) {
		String[] kmers = new String[seq.length()-KMER+1];
		for (int i=0; i<=seq.length()-KMER; i++) {
			String kmer = seq.substring(i, i+KMER);
			kmers[i] = kmer;
		}
		
		return kmers;
	}
	
	public static void main(String[] args) throws Exception {
		String refSeq = args[0];
		String reads = args[1];
		String output = args[2];
		int kmer = Integer.parseInt(args[3]);
//		String refSeq = "/home/lmose/dev/vdj/extraction/igh.fa";
//		String reads = "/home/lmose/dev/vdj/extraction/selected_subset_pass2.sort.bam";
//		String output = "/home/lmose/dev/vdj/extraction/selected_subset_pass2.k15.sort.bam";
//		int kmer = 15;
		
		new TestExtractor2().run(reads, refSeq, kmer, output);
	}
}
