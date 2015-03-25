package vdj;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Reads a fasta file, outputs in frame sequence (i.e. sequence lacking a stop codon)
 * Shifts 0,1 or 2 bases as appropriate.  Sequences containing a stop codon are filtered 
 * @author lmose
 */
public class FilterStopCodons {
	
	ReverseComplementor rc = new ReverseComplementor(); 

	public void run(String input, String output, int dir) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(input));
		BufferedWriter writer = new BufferedWriter(new FileWriter(output, false));
		
		
		
		String id = reader.readLine();
		String bases = null;
		if (id != null) {
			bases = reader.readLine();
		}
		
		while (id != null && bases != null) {
		
			String seq = getInFrameSequence(bases, dir);
			
			if (seq != null) {
				writer.write(id + "\n" + seq + "\n");
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
	
	private String getInFrameSequence(String bases, int dir) {
		
		if (dir == 1 || dir == 3) {
			if (isInFrame(bases)) {
				return bases;
			}
			
			if (isInFrame(bases.substring(1))) {
				return bases.substring(1);
			}
			
			if (isInFrame(bases.substring(2))) {
				return bases.substring(2);
			}
		}
		
		if (dir == 2 || dir == 3) {
			bases = rc.reverseComplement(bases);
			
			if (isInFrame(bases)) {
				return bases;
			}
			
			if (isInFrame(bases.substring(1))) {
				return bases.substring(1);
			}
			
			if (isInFrame(bases.substring(2))) {
				return bases.substring(2);
			}
		}
		
		return null;
	}
	
	private boolean isInFrame(String bases) {
		for (int i=0; i<bases.length()-6; i+=3) {
			String codon = bases.substring(i, i+3);
			if (codon.equals("TAG") || codon.equals("TAA") || codon.equals("TGA")) {
				return false;
			}
		}
		
		return true;
	}
	
	public static void main(String[] args) throws Exception {
		String input = args[0];
		String output = args[1];
		int rc = Integer.parseInt(args[2]);  // 1 = forward.  2 = reverse.  3 = both
		
		FilterStopCodons f = new FilterStopCodons();
		f.run(input, output, rc);
	}
}


