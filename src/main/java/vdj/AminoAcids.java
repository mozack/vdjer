package vdj;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

public class AminoAcids {

	private Map<String, String> entries = new HashMap<String, String>();
	
	private static AminoAcids aminos = new AminoAcids();
	
	static {
		aminos.init();
	}
	
	private void add(String aminoAcid, String codonsStr) {
		String[] codons = codonsStr.split(",");
		
		for (String codon : codons) {
			entries.put(codon, aminoAcid);
		}
	}
	
	public static String convert(String bases) {
		StringBuffer seq = new StringBuffer();
		for (int i=0; i<bases.length()-2; i+=3) {
			String amino = aminos.entries.get(bases.substring(i,i+3));
			if (amino != null) {
				seq.append(amino);
			} else {
				seq.append('X');
			}
		}
		
		return seq.toString();
	}
	
	public static void convertFasta(String input, String output) throws IOException  {
		BufferedReader reader = new BufferedReader(new FileReader(input));
		BufferedWriter writer = new BufferedWriter(new FileWriter(output, false));
		
		String id = reader.readLine();
		String bases = null;
		if (id != null) {
			bases = reader.readLine();
		}
		
		long s = System.currentTimeMillis();
		while (id != null && bases != null) {
						
			String aminoSeq = convert(bases);
			
			writer.write(id + '\n' + aminoSeq + '\n');
			
			id = reader.readLine();
			bases = null;
			if (id != null) {
				bases = reader.readLine();
			}

		}
		
		reader.close();
		writer.close();
	}
	
	private void init() {
		aminos.add("A","GCA,GCC,GCG,GCT");
		aminos.add("C","TGC,TGT");
		aminos.add("D","GAC,GAT");
		aminos.add("E","GAA,GAG");
		aminos.add("F","TTC,TTT");
		aminos.add("G","GGA,GGC,GGG,GGT");
		aminos.add("H","CAC,CAT");
		aminos.add("I","ATA,ATC,ATT");
		aminos.add("K","AAA,AAG");
		aminos.add("L","TTA,TTG,CTA,CTC,CTG,CTT");
		aminos.add("M","ATG");
		aminos.add("N","AAC,AAT");
		aminos.add("P","CCA,CCC,CCG,CCT");
		aminos.add("Q","CAA,CAG");
		aminos.add("R","CGA,CGC,CGG,CGT,AGA,AGG");
		aminos.add("S","TCA,TCC,TCG,TCT,AGC,AGT");
		aminos.add("T","ACA,ACC,ACG,ACT");
		aminos.add("V","GTA,GTC,GTG,GTT");
		aminos.add("W","TGG");
		aminos.add("Y","TAC,TAT");
		aminos.add(".","TAA,TAG,TGA");
	}
	
	public static void main(String[] args) throws Exception {
		String mode = args[0];
		
		if (mode.equals("fasta")) {
			String input = args[1];
			String output = args[2];
			
			convertFasta(input, output);
		} else if (mode.equals("bases")) {
			String bases = args[1];
			System.out.println(convert(bases));
		}
	}
}
