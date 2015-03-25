package vdj;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class ReadExtractorFull {
	
	private ReverseComplementor rc = new ReverseComplementor();

	public void run(String inFile1, String inFile2, String outFile) throws IOException {
		
		BufferedWriter output = new BufferedWriter(new FileWriter(outFile, false));

		extract(inFile1, output);
		extract(inFile2, output);
		output.close();
	}
	
	private void extract(String inFile, BufferedWriter output) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(inFile));
		
		int index = 0;
		String bases = "";
		String quals = "";
		
		String line = reader.readLine();
		while (line != null) {
			if (index == 1) {
				bases = line;
			} else if (index == 3) {
				quals = line;
			}
			index += 1;
			
			if (index == 4) {
				index = 0;
				output(bases, quals, output);
				bases = "";
				quals = "";
			}
			
			line = reader.readLine();
		}
		
		reader.close();
	}
	
	private void output(String bases, String quals, BufferedWriter output) throws IOException {
		output.append(getOutputString(false, bases, quals));
		
		// We don't know the orientation of the read, so reverse complement it as well.
		output.append(getOutputString(false, rc.reverseComplement(bases), rc.reverse(quals)));
	}
	
	private String getOutputString(boolean negativeStrandFlag, String bases, String qualities) {
		StringBuffer readBuffer = new StringBuffer();
		readBuffer.append(negativeStrandFlag ? "1" : "0");
		
		readBuffer.append(bases);
		readBuffer.append(qualities);
		
		return readBuffer.toString();
	}
	
	public static void main(String[] args) throws IOException {
		String in1 = args[0];
		String in2 = args[1];
		String out = args[2];
	
		new ReadExtractorFull().run(in1, in2, out);
	}
}
