package vdj;

import java.io.BufferedReader;
import java.io.FileReader;

public class RC {

	public static void run(String input) throws Exception {
		ReverseComplementor rc = new ReverseComplementor();
		BufferedReader reader = new BufferedReader(new FileReader(input));
		
		String line = reader.readLine();
		while (line != null) {
			String[] fields = line.split("\\s");
			fields[0] = rc.reverseComplement(fields[0]);
			
			StringBuffer buf = new StringBuffer();
			for (int i=0; i<fields.length; i++) {
				buf.append(fields[i]);
				buf.append('\t');
			}
			
			System.out.println(buf.toString());
			line = reader.readLine();
		}
		
		reader.close();
	}
	
	public static void main(String[] args) throws Exception {
		run(args[0]);
	}
}
