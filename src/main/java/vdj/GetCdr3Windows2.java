package vdj;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class GetCdr3Windows2 {
	private ReverseComplementor rc = new ReverseComplementor();
	private Pattern regex = Pattern.compile("(TT[TC]|TA[CT])(TT[CT]|TA[TC]|CA[TC]|GT[AGCT]|TGG)(TG[TC])(([GA][AGCT])|TC)[AGCT]([ACGT]{3}){5,32}TGGG[GCT][GCT]");
	
	private Set<String> windows = new HashSet<String>();

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
//				writer.write(id + "_1\n" + seq + "\n");
				windows.add(seq);
			}
			
			if (windowSize <= 0) {
				seq = getWindow3(rc.reverseComplement(bases));	
			} else {
				seq = getWindow(rc.reverseComplement(bases), windowSize);
			}
			
			if (seq != null) {
//				writer.write(id + "_2\n" + seq + "\n");
				windows.add(seq);
			}
			
			id = reader.readLine();
			bases = null;
			if (id != null) {
				bases = reader.readLine();
			}
		}
		
		Map<String,Set<String>> aaGroups = new HashMap<String,Set<String>>();
		
		for (String window : windows) {
			String aa = AminoAcids.convert(window);
			Set<String> group = aaGroups.get(aa);
			if (group == null) {
				group = new HashSet<String>();
				group.add(window);
				aaGroups.put(aa, group);
			} else {
				group.add(window);
			}
		}
		
		for (String aa : aaGroups.keySet()) {
			StringBuffer buf = new StringBuffer();
			buf.append(aa);
			for (String window : aaGroups.get(aa)) {
				buf.append('\t');
				buf.append(window);
			}
			
			writer.write(buf.toString() + '\n');
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
//				System.exit(-1);
				return null;
			}
			
			int padding = windowSize - cdr3Len;
			
			int startPadding = padding / 2;
			
			// Keep CDR3 in frame
			startPadding = startPadding - (startPadding % 3);
			int endPadding = padding - startPadding;
			
			start = start - startPadding;
			end = end + endPadding;
						
			if (start < 0) {
				System.out.println("Insufficient padding to left of CDR3: " + bases);
			} else if (end >= bases.length()) {
				System.out.println("Insuffient padding to right of CDR3: " + bases);
			} else {
				window = bases.substring(start, end);
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
		
		GetCdr3Windows2 cdr3 = new GetCdr3Windows2();
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
