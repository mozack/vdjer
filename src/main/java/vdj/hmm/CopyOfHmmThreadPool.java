package vdj.hmm;

import htsjdk.samtools.SAMRecord;

import java.io.BufferedWriter;
import java.util.List;


public class CopyOfHmmThreadPool {

	private List<CopyOfHmmThread> available;
	
	public void init(int size, int readLength, List<SAMRecord> reads, BufferedWriter writer) {
		for (int i=0; i<size; i++) {
//			available.add(new CopyOfHmmThread(this, readLength, reads, writer));
		}
		
		for (CopyOfHmmThread t : available) {
			Thread thread = new Thread(t);
			thread.start();
		}
	}
	
	public HmmThread acquire() {
		HmmThread thread = null;
		boolean isAcquired = false;
		while (!isAcquired) {
			synchronized(this) {
				if (!available.isEmpty()) {
//					thread = available.get(0);
					isAcquired = true;
				}
			}
			if (!isAcquired) {
				try {
					Thread.sleep(50);
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
		}
		
		return thread;
	}
	
	public void release(HmmThread thread) {
		synchronized(this) {
//			available.add(thread);
		}
	}
}
