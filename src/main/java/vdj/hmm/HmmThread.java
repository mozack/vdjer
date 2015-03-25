package vdj.hmm;

import htsjdk.samtools.SAMRecord;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


public class HmmThread implements Runnable {
	
	private ReadLikelihoodCalc calc = new ReadLikelihoodCalc(ContigSelector.READ_LENGTH);
	
	private String line;
	private List<SAMRecord> reads;
	private BufferedWriter writer;
	private List<HmmThread> threads;
	
	public HmmThread(int readLength, List<SAMRecord> reads, 
			String line, BufferedWriter writer, List<HmmThread> threads) {
		calc = new ReadLikelihoodCalc(readLength);
		this.reads = reads;
		this.writer = writer;
		this.threads = threads;
		this.line = line;
	}
	
//	public void score(String id, String haplotype) {
//		this.id = id;
//		this.haplotype = haplotype;
//	}

	@Override
	public void run() {
		
		String[] fields = line.split("\t");
		
		List<String> haplotypes = new ArrayList<String>();
		for (int i=1; i<fields.length; i++) {
			haplotypes.add(fields[i]);
		}
		
		//double score = calc.getHaplotypeScore2(haplotype, reads);
		double score = calc.getGroupedHaplotypesScore(haplotypes, reads);
//		double score = calc.getHaplotypeScore(haplotype, reads);
		synchronized(writer) {
			try {
				writer.write("" + score + "\t" + line + '\n');
			} catch (IOException e) {
				e.printStackTrace();
				System.exit(-1);
			}
		}
		
//		System.out.println("Releasing thread...");
		synchronized(threads) {
			threads.remove(this);
		}
	}
}
