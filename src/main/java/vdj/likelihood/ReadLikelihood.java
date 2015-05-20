package vdj.likelihood;

import java.util.Arrays;
import java.util.List;

import vdj.ReverseComplementor;

import htsjdk.samtools.SAMRecord;

public class ReadLikelihood {
	
	public static final double READ_LIKELIHOOD_FLOOR = -100;
	private static double PROB_MATRIX_INIT = Double.MAX_VALUE;
	private static double LOG10_PROB_MATRIX_INIT = Math.log10(PROB_MATRIX_INIT); 
	
	private static double[] phredProb = initPhredProb();
	
	private static final int MAX_BASE_QUAL = 254;
	
	private ReverseComplementor rc = new ReverseComplementor();
	
	private static synchronized double[] initPhredProb() {
		double[] phredProb = new double[MAX_BASE_QUAL+1];
		for (int qual=0; qual<=254; qual++) {
			phredProb[qual] = Math.pow(10.0, qual / -10.0);
		}
		
		return phredProb;
	}
	
	public double calcReadLikelihood(SAMRecord read, String contig) {
		double score1 = score(read.getReadString(), read.getBaseQualities(), contig);
		double score2 = score(rc.reverseComplement(read.getReadString()), rc.reverse(read.getBaseQualities()), contig);
		
		return Math.max(score1, score2);
	}
	
	public double getGroupedContigScore(List<String> contigs, List<SAMRecord> reads) {
		double score = 0;
		
		for (SAMRecord read : reads) {
			double maxReadLikelihood = READ_LIKELIHOOD_FLOOR;
			for (String contig : contigs) {
				double readLikelihood = calcReadLikelihood(read, contig);
				maxReadLikelihood = Math.max(readLikelihood, maxReadLikelihood);
			}
			
			score += maxReadLikelihood;
		}
		
		return score;
	}
	
	public double score(SAMRecord read, String contig) {
		return score(read.getReadString(), read.getBaseQualities(), contig);
	}
	
	public double score(String readBases, byte[] readQuals, String contig) {
		
		// Track 2 matrices.  First for standard SW.  Second for prob tracking.
		int[][] swMatrix = new int[readBases.length()+1][contig.length()+1];
		double[][] probMatrix = new double[readBases.length()+1][contig.length()+1];
		
		
		// Init matrices to large values to avoid numeric underflow
		for (int c=0; c<=contig.length(); c++) {
			probMatrix[0][c] = PROB_MATRIX_INIT;
		}
		
		for (int r=0; r<=readBases.length(); r++) {
			probMatrix[r][0] = PROB_MATRIX_INIT;
		}
		
		for (int r=0; r<readBases.length(); r++) {
			for (int c=0; c<contig.length(); c++) {
				double baseLikelihood = calcQualLikelihood(r, c, readBases, readQuals, contig);
				probMatrix[r+1][c+1] = probMatrix[r][c] * baseLikelihood;
				
				swMatrix[r+1][c+1] = swMatrix[r][c] + (readBases.charAt(r) == contig.charAt(c) ? 1 : 0);
			}
		}
		
		// Get likelihood of "best" alignment in last row (end of read)
		int bestSW = -1;
		double bestProb = Double.NEGATIVE_INFINITY;
		
		// Read is shorter or equal to length of contig.  Look in last
		for (int c=1; c<=contig.length(); c++) {
			if (swMatrix[readBases.length()][c] > bestSW) {
				// New best score
				bestSW   = swMatrix[readBases.length()][c];
				bestProb = probMatrix[readBases.length()][c];
				
			} else if (swMatrix[readBases.length()][c] == bestSW) {
				// Tie for best SW score.  Select the best likelihood
				if (probMatrix[readBases.length()][c] > bestProb) {
					bestProb = probMatrix[readBases.length()][c];
				}
			}
		}
		
		for (int r=1; r<=readBases.length(); r++) {
			if (swMatrix[r][contig.length()] > bestSW) {
				// New best score
				bestSW   = swMatrix[r][contig.length()];
				bestProb = probMatrix[r][contig.length()]; 
			} else if (swMatrix[r][contig.length()] == bestSW) {
				// Tie for best SW score.  Select the best likelihood
				if (probMatrix[r][contig.length()] > bestProb) {
					bestProb = probMatrix[r][contig.length()];
				}
			}
		}
		
//		dumpMatrices(swMatrix, probMatrix);
		
		// Return log scaled prob adjusted by the initial matrix value.
		return Math.log10(bestProb) - LOG10_PROB_MATRIX_INIT;
	}
	
	private void dumpMatrices(int[][] sw, double[][] prob) {
		for (int i=0; i<sw.length; i++) {
			System.out.println(Arrays.toString(sw[i]));	
		}

		for (int i=0; i<sw.length; i++) {
			System.out.println(Arrays.toString(prob[i]));	
		}
	}
	
	private double calcQualLikelihood(int readIndex, int contigIndex, String readBases, byte[] readQuals, String contig) {
		double prior;
		
		if (readBases.charAt(readIndex) == contig.charAt(contigIndex)) {
			prior = 1 - phredProb[readQuals[readIndex]];
		} else {
			prior = phredProb[readQuals[readIndex]];
		}
		
		return prior;
	}
	
	static class MatrixLoc {
		private int row;
		private int col;
		
		public MatrixLoc(int row, int col) {
			this.row = row;
			this.col = col;
		}

		public int getRow() {
			return row;
		}

		public int getCol() {
			return col;
		}
	}
	
	public static void main(String[] args) {
		
		String contig = "ATCGA";
		String read = "TCG";
		
//		byte[] readQuals = new byte[] { 40, 40, 30, 25, 39 };
		byte[] readQuals = new byte[read.length()];
		Arrays.fill(readQuals, (byte) 40);

		ReadLikelihood l = new ReadLikelihood();
		double p;
		
//		p = l.score(read, readQuals, contig);
//		System.out.println("p: " + p);
//		
//		readQuals = new byte[] { 40, 40, 30 };
//		p = l.score(read, readQuals, contig);
//		System.out.println("p2: " + p);
//		
//		read = "CCG";
//		
//		Arrays.fill(readQuals, (byte) 40);
//		p = l.score(read, readQuals, contig);
//		System.out.println("p3: " + p);
//		
//		readQuals = new byte[] { 2, 40, 40 };
//		p = l.score(read, readQuals, contig);
//		System.out.println("p4: " + p);
//		
//		contig = "ACA";
//		read = "CGCAG";
//		readQuals = new byte[] { 40, 10, 40, 40, 40 };
//		p = l.score(read,  readQuals, contig);
//		System.out.println("p5: " + p);
//		
//		contig = "ACACACACAC";
//		read = "ACTC";
//		readQuals = new byte[] { 40,30,20,10 };
//		p = l.score(read,  readQuals, contig);
//		System.out.println("p6: " + p);
//
//		read = "ACACACACA";
//		contig = "ACAC";
//		readQuals = new byte[] { 40,40, 40, 40, 40, 40, 40, 40, 40 }; 
//		p = l.score(read,  readQuals, contig);
//		System.out.println("p7: " + p);
		
//		contig = "GGGGACA";
//		read = "ACAC";
//		readQuals = new byte[] { 40,30,20,10 };
//		p = l.score(read,  readQuals, contig);
//		System.out.println("p6: " + p);
		
		contig = "GGGGACA";
		read = "AGGG";
		readQuals = new byte[] { 40,30,20,10 };
		p = l.score(read,  readQuals, contig);
		
		System.out.println("p6: " + p);
		
		System.out.println("log(p6): " + Math.log10(p));
	}
}
