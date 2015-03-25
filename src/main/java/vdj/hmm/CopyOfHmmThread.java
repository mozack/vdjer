package vdj.hmm;

import htsjdk.samtools.SAMRecord;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;


public class CopyOfHmmThread implements Runnable {
	
	private ReadLikelihoodCalc calc = new ReadLikelihoodCalc(48);
	
	private HmmThreadPoolBak pool;
	private String id;
	private String haplotype;
	private List<SAMRecord> reads;
	private BufferedWriter writer;
	
	public CopyOfHmmThread(HmmThreadPoolBak pool, int readLength, List<SAMRecord> reads, BufferedWriter writer) {
		this.pool = pool;
		calc = new ReadLikelihoodCalc(readLength);
		this.reads = reads;
		this.writer = writer;
	}
	
	public void score(String id, String haplotype) {
		this.id = id;
		this.haplotype = haplotype;
	}

	@Override
	public void run() {
		
		boolean shutdown = false;
		
		while (!shutdown) {
			synchronized(this) {
				if (id == null || haplotype == null) {
					
				} else {
					double score = calc.getHaplotypeScore(haplotype, reads);
					synchronized(writer) {
						try {
							writer.write(id + "\tscore:\t" + score + "\n" + haplotype);
						} catch (IOException e) {
							e.printStackTrace();
							System.exit(-1);
						}
					}
					id = null;
					haplotype = null;
				}
			}
			
//			pool.release(this);
		}
	}
}
