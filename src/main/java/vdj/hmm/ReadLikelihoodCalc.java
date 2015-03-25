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
import java.util.Arrays;
import java.util.List;


import org.broadinstitute.sting.utils.pairhmm.LoglessPairHMM;

import vdj.ReverseComplementor;

public class ReadLikelihoodCalc {
	
	public static final double READ_LIKELIHOOD_FLOOR = -100;
	
	private byte[] insertionGOP;
	private byte[] deletionGOP;
	private byte[] overallGCP;
	
	private VdjPairHmm hmm;
	
	private ReverseComplementor rc = new  ReverseComplementor();
	
	public ReadLikelihoodCalc(int readLength) {
		insertionGOP = new byte[readLength];
		deletionGOP = new byte[readLength];
		overallGCP = new byte[readLength];
		Arrays.fill(insertionGOP, (byte) 40);
		Arrays.fill(deletionGOP, (byte) 40);
		Arrays.fill(overallGCP, (byte) 10);
		
		hmm = new VdjPairHmm();
		
		hmm.initialize(50, 5000);

	}
	
	private List<SAMRecord> loadReads(String input) {
		System.out.println("Loading reads");
		SAMFileReader reader = new SAMFileReader(new File(input));
		reader.setValidationStringency(ValidationStringency.SILENT);

		List<SAMRecord> reads = new ArrayList<SAMRecord>();
		for (SAMRecord read : reader) {
			// Only consider those reads mapped beyond the constant regions or containing 
			if (read.getReferenceName().equals("chr12") && read.getAlignmentStart() >= 113256204 && read.getAlignmentEnd() <= 113422730) {
				// Maps to constant region. filter....
			} else {
				reads.add(read);
			}
		}
		
		reader.close();
		System.out.println("Reads loaded");
		
		return reads;
	}
	
	public double getHaplotypeScore(String haplotypeStr, List<SAMRecord> reads) {
		double score = 0;
		byte[] haplotype = haplotypeStr.getBytes();
		
		for (SAMRecord read : reads) {
			double readLikelihood = getReadLikelihood(read, haplotype);
			
			score += Math.max(readLikelihood, READ_LIKELIHOOD_FLOOR);
		}
		
		return score;
	}
	
	public double getGroupedHaplotypesScore(List<String> haplotypes, List<SAMRecord> reads) {
		double score = 0;
		
		for (SAMRecord read : reads) {
			double maxReadLikelihood = READ_LIKELIHOOD_FLOOR;
			for (String haplotype : haplotypes) {
				double readLikelihood = getReadLikelihood(read, haplotype.getBytes());
				maxReadLikelihood = Math.max(readLikelihood, maxReadLikelihood);
			}
			
			score += maxReadLikelihood;
		}
		
		return score;
	}
	
	public double getHaplotypeScore2(String haplotypeStr, List<SAMRecord> reads) {
		double score = 0;
		byte[] haplotype = haplotypeStr.getBytes();
		System.out.println("----------------- SCORING -----------------");
		SAMRecord last = null;
		
		for (SAMRecord read : reads) {
			
			double readLikelihood = 0;
			
			if (last != null) {
				if (read.getReadName().equals(last.getReadName())) {
					if (read.getFirstOfPairFlag() != last.getFirstOfPairFlag()) {
						// We have a read pair.  Score them together.
						// TODO: Score with proper orientations 
						readLikelihood = Math.min(getReadLikelihood(read, haplotype), getReadLikelihood(last, haplotype));
						last = null;
					} else {
						System.err.println("*** WARNING *** : Duplicate read: " + read.getReadName());
						last = read;
					}
				} else {
					System.err.println("*** Scoring singleton: " + last.getReadName());
					readLikelihood = getReadLikelihood(last, haplotype);
					last = read;
				}
			} else {
				last = read;
			}
			
			score += Math.max(readLikelihood, READ_LIKELIHOOD_FLOOR);
		}
		
		return score;
	}
	
	public void scoreHaplotypes(String readsBam, String haplotypeFasta, String output) throws IOException {
		int count = 0;
		
		List<SAMRecord> reads = loadReads(readsBam);
		
		BufferedReader reader = new BufferedReader(new FileReader(haplotypeFasta));
		BufferedWriter writer = new BufferedWriter(new FileWriter(output, false));
		
		double max = Double.NEGATIVE_INFINITY;
		String maxId = "None";
		String maxBases = "";
		
		String id = reader.readLine();
		String bases = null;
		if (id != null) {
			bases = reader.readLine();
		}
		
		long s = System.currentTimeMillis();
		while (id != null && bases != null) {
//			double score = getHaplotypeScore(bases, reads);
			double score = getHaplotypeScore2(bases, reads);
			
			if (score > max) {
				max = score;
				maxId = id;
				maxBases = bases;
			}
			
			writer.write(id + "\tscore:\t" + score);
			writer.write(bases);
			
			System.out.println(id + "\tscore:\t" + score);
			System.out.println(bases);
			
			id = reader.readLine();
			bases = null;
			if (id != null) {
				bases = reader.readLine();
			}
			count += 1;
			if ((count % 1000) == 0) {
				System.out.println("Processed: " + count + " haplotypes.");
			}
			long e = System.currentTimeMillis();
			System.out.println("count: " + count + (e-s));
		}
		
		reader.close();
		writer.close();
		
		System.out.println(maxId + "\tscore:\t" + max + "\n" + maxBases);
	}
	
	private double getReadLikelihood(SAMRecord read, byte[] haplotype) {
		hmm.setValues(read.getReadLength() + 1, haplotype.length + 1);
		
		double score1 = hmm.subComputeReadLikelihoodGivenHaplotypeLog10(haplotype, read.getReadBases(), read.getBaseQualities(), insertionGOP, deletionGOP, overallGCP, 0, true, 0);
		double score2 = hmm.subComputeReadLikelihoodGivenHaplotypeLog10(haplotype, rc.reverseComplement(read.getReadBases()), rc.reverse(read.getBaseQualities()), insertionGOP, deletionGOP, overallGCP, 0, true, 0);
		
		return Math.max(score1, score2);
//		return score1;
	}

	public void run() {
		VdjPairHmm hmm = new VdjPairHmm();
		
		hmm.initialize(50, 5000);
		
//		String hap = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
		String hap = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATC";
		String read  = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
//		String read = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";
		
		long s = System.currentTimeMillis();
		double likelihood = 0;
		
		byte[] hapBases = hap.getBytes();
		byte[] readBases = read.getBytes();
//		byte[] readQuals = new byte[] { 40, 40, 30, 25, 39 };
		byte[] readQuals = new byte[readBases.length];
		Arrays.fill(readQuals, (byte) 40);
		
		hmm.setValues(readBases.length + 1, hapBases.length + 1);
		
		likelihood = hmm.subComputeReadLikelihoodGivenHaplotypeLog10(hapBases, readBases, readQuals, insertionGOP, deletionGOP, overallGCP, 0, true, 0);
		
		long e = System.currentTimeMillis();
		
		System.out.println("likelihood: " + likelihood);
		
//		System.out.println("Elapsed msecs:" + (e-s));
		
//		System.out.println("likelihood: " + likelihood);
	}
	
	static class VdjPairHmm extends LoglessPairHMM {
		void setValues(int paddedReadLen, int paddedHapLen) {
			this.paddedReadLength = paddedReadLen;
			this.paddedHaplotypeLength = paddedHapLen;
		}
	}
	
	public static void main(String[] args) throws Exception {
		System.out.println("Starting HMM");
		
//		String readsBam = args[0];
//		String haplotypeFasta = args[1];
//		String output = args[2];
		
		ReadLikelihoodCalc c = new ReadLikelihoodCalc(48);
		c.run();

		
		/*
		String readsBam = "/home/lmose/dev/vdj/scoring/assembly.reads.bam";
		String haplotypeFasta = "/home/lmose/dev/vdj/scoring/map_to_v.fa";
		String output = "/home/lmose/dev/vdj/scoring/haplotype_scores.fa";
		
		c.scoreHaplotypes(readsBam, haplotypeFasta, output);
		
		System.out.println("Done");
		*/
	}
}
