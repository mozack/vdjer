package vdj.hmm;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;


public class ContigSelector {
	
	public static int READ_LENGTH = 0;
	
	private List<HmmThread> threads = new ArrayList<HmmThread>();
	
	private static int MAX_THREADS = 16;

	private List<SAMRecord> loadReads(String input) {
		System.out.println("Loading reads");
		SAMFileReader reader = new SAMFileReader(new File(input));
		reader.setValidationStringency(ValidationStringency.SILENT);

		List<SAMRecord> reads = new ArrayList<SAMRecord>();
		for (SAMRecord read : reader) {
			if (READ_LENGTH == 0) {
				READ_LENGTH = read.getReadLength();
			}
			reads.add(read);
//			// Only consider those reads mapped beyond the constant regions or containing 
//			if (read.getReferenceName().equals("chr12") && read.getAlignmentStart() >= 113256204 && read.getAlignmentEnd() <= 113422730) {
//				// Maps to constant region. filter....
//			} else {
//				reads.add(read);
//			}
		}
		
		reader.close();
		
		Collections.sort(reads, new ReadNameComparator());
		
		System.out.println("Reads loaded: " + reads.size());
		
		return reads;
	}

	private void spawnThread(int readLength, List<SAMRecord> reads, BufferedWriter writer, String line) {
		HmmThread thread = new HmmThread(readLength, reads, 
				line, writer, threads);
		
		while (threads.size() >= MAX_THREADS) {
			try {
				Thread.sleep(50);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		
		synchronized(threads) {
			threads.add(thread);
		}
		
//		System.out.println("Spawning thread...");
		Thread t = new Thread(thread);
		t.start();
	}
	
	public void scoreHaplotypes(String readsBam, String haplotypeFasta, String output) throws IOException {
		int count = 0;
		
		List<SAMRecord> reads = loadReads(readsBam);
		
		BufferedReader reader = new BufferedReader(new FileReader(haplotypeFasta));
		BufferedWriter writer = new BufferedWriter(new FileWriter(output, false));
		
		double max = Double.NEGATIVE_INFINITY;
		String maxId = "None";
		String maxBases = "";
		
		String line = reader.readLine();
		
		long s = System.currentTimeMillis();
		while (line != null) {
						
			spawnThread(READ_LENGTH, reads, writer, line);
			
			count += 1;
			if ((count % 1000) == 0) {
				System.out.println("Processed: " + count + " haplotypes.");
				synchronized(writer) {
					writer.flush();
				}
			}
			
			line = reader.readLine();
		}
		
		System.out.println("Waiting for all threads to complete");
		while (threads.size() > 0) {
			try {
				Thread.sleep(50);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
		}
		
		long e = System.currentTimeMillis();
		System.out.println("count: " + count + "\t" + (e-s)/1000);
		
		reader.close();
		writer.close();
		
//		System.out.println(maxId + "\tscore:\t" + max + "\n" + maxBases);
	}

	static class ReadNameComparator implements Comparator<SAMRecord> {

		@Override
		public int compare(SAMRecord o1, SAMRecord o2) {
			return o1.getReadName().compareTo(o2.getReadName());
		}
		
	}
	
	public static void main(String[] args) throws Exception {
		
		String readsBam = args[0];
		String input = args[1];
		String output = args[2];
		MAX_THREADS = Integer.parseInt(args[3]);
		
		
/*
		String readsBam = "/home/lmose/dev/vdj/scoring/human/cdr3.reads.bam";
		String input = "/home/lmose/dev/vdj/scoring/human/cdr3.txt";
		String output = "/home/lmose/dev/vdj/scoring/human/cdr3.txt.scores";
		MAX_THREADS = 1;
*/		
		ContigSelector cs = new ContigSelector();
		cs.scoreHaplotypes(readsBam, input, output);
	}
}
