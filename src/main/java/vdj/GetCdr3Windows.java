package vdj;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class GetCdr3Windows {
	private ReverseComplementor rc = new ReverseComplementor();
	private Pattern regex = Pattern.compile("(TT[TC]|TA[CT])(TT[CT]|TA[TC]|CA[TC]|GT[AGCT]|TGG)(TG[TC])(([GA][AGCT])|TC)[AGCT]([ACGT]{3}){5,32}TGGG[GCT][GCT]");

	public void run(String input, String output, int windowSize) throws IOException {
		
		BufferedReader reader = new BufferedReader(new FileReader(input));
		BufferedWriter writer = new BufferedWriter(new FileWriter(output, false));
		
		
		String id = reader.readLine();
		String bases = null;
		if (id != null) {
			bases = reader.readLine();
		}
		
		while (id != null && bases != null) {

			String seq = null;
			if (windowSize <= 0) {
				seq = getWindow3(bases);
			} else {
				seq = getWindow(bases, windowSize);
			}
			
			if (seq != null) {
				writer.write(id + "_1\n" + seq + "\n");
			}
			
			if (windowSize <= 0) {
				seq = getWindow3(rc.reverseComplement(bases));	
			} else {
				seq = getWindow(rc.reverseComplement(bases), windowSize);
			}
			
			if (seq != null) {
				writer.write(id + "_2\n" + seq + "\n");
			}
			
			id = reader.readLine();
			bases = null;
			if (id != null) {
				bases = reader.readLine();
			}
		}

		reader.close();
		writer.close();
	}
	
	private String getWindow(String bases, int windowSize) {
		String window = null;
		Matcher m = regex.matcher(bases);
		if (m.find()) {
			int start = m.start();
			int end = m.end();
			
			String cdr3 = bases.substring(start, end);
			
			int cdr3Len = end-start;
			if (cdr3Len > windowSize) {
				System.out.println("cdr3 len > window size: " + bases);
				System.exit(-1);
			}
			
			int padding = windowSize - cdr3Len;
			
			int newEnd = end + padding/2;
			int newLength = newEnd-start;
			int newStart = start - (windowSize - (newLength));
			
			if (newStart < 0) {
				System.out.println("Insufficient padding to left of CDR3: " + bases);
			} else if (newEnd >= bases.length()) {
				System.out.println("Insuffient padding to right of CDR3: " + bases);
			} else {
				window = bases.substring(newStart, newEnd);
			}
		}
		
		return window;
	}
	
	private String getWindow2(String bases) {
		String window = null;
		Matcher m = regex.matcher(bases);
		if (m.find()) {
			int start = m.start();
			int end = m.end();
			
//			String cdr3 = bases.substring(start, end);
			
			int cdr3Len = end-start;
			
			int MAX_LEFT_WINDOW_PAD = 500;
			int MIN_LEFT_WINDOW_PAD = 50;
			int RIGHT_WINDOW_PAD = 40;
			
			if (start < MIN_LEFT_WINDOW_PAD) {
				System.out.println("Insufficient padding to left of CDR3: " + bases);
			} else if (start > MAX_LEFT_WINDOW_PAD) {
				System.out.println("Too much padding to left of CDR3: " + bases);
			} else if (bases.length() < end +  RIGHT_WINDOW_PAD) {
				System.out.println("Insuffient padding to right of CDR3: " + bases);
			} else {
				window = bases.substring(0, end +  RIGHT_WINDOW_PAD);
			}		
		}
		
		return window;
	}

	private String getWindow3(String bases) {
		
		String window = null;
		
		if (bases.length() >= 500) {
			Matcher m = regex.matcher(bases);
			if (m.find()) {
				int start = m.start();
				int end = m.end();
				
	//			String cdr3 = bases.substring(start, end);
				
				int cdr3Len = end-start;
				
				int MAX_LEFT_WINDOW_PAD = 500;
				int MIN_LEFT_WINDOW_PAD = 50;
				int MAX_CDR3_END = 400;
				int RIGHT_WINDOW_PAD = 80;
				
				if (start < MIN_LEFT_WINDOW_PAD) {
					System.out.println("Insufficient padding to left of CDR3: " + bases);
				} else if (start > MAX_LEFT_WINDOW_PAD) {
					System.out.println("Too much padding to left of CDR3: " + bases);
				} else if (end > MAX_CDR3_END) {
					System.out.println("CDR3 too far away from contig start: " + bases);
				} else if (bases.length() < end +  RIGHT_WINDOW_PAD) {
					System.out.println("Insuffient padding to right of CDR3: " + bases);
				} else {
					//window = bases.substring(0, end +  RIGHT_WINDOW_PAD);
					window = bases.substring(0, 500);
				}
			}
		}
		
		return window;
	}
	
	public static void main(String[] args) throws Exception {
		
		String input = args[0];
		String output = args[1];
		int windowSize = Integer.parseInt(args[2]);
		
		GetCdr3Windows cdr3 = new GetCdr3Windows();
		cdr3.run(input,  output, windowSize);
		
		System.out.println("Done.");
		
		/*
		String bases = "TTAAATCAGAAGTTCAAGGACAAGGCCACATTGACTGTAGACAAATCCTCCAACACAGCCTACATGCAACTCAGCAGCCCGACATCTGAGGACTCTGCGGTCTATTACTGTGCAAGAGGGGGGACCTACTATTGGTACGACTATGCTATGGACTACTGGGGTCAAGGAACCTCAGTCACCGTCTCCTCAGCCAAAACAACACCCCCATCTGTCTATCCAC";
		
		
		GetCdr3Windows cdr3 = new GetCdr3Windows();
		String window = cdr3.getWindow2(bases, 80);
		System.out.println(window);
		*/
	}
}
