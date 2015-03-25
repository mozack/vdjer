package vdj;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

/**
 * Extracts all reads from a BAM file into abra/vdj assembly format.
 * @author lmose
 *
 */
public class ReadExtractorFull2 {
	
	private ReverseComplementor rc = new ReverseComplementor();
	
	private static int READ_LENGTH = 50;

	public void run(String inFile1, String outFile) throws IOException {
		
		BufferedWriter output = new BufferedWriter(new FileWriter(outFile, false));

		extract(inFile1, output);
		output.close();
	}
	
	private void extract(String inFile, BufferedWriter output) throws IOException {
		
		SAMFileReader reader = new SAMFileReader(new File(inFile));
		reader.setValidationStringency(ValidationStringency.SILENT);
		for (SAMRecord read : reader) {
			String bases = read.getReadString();
			
			while (bases.length() < READ_LENGTH) {
				bases += 'N';
			}
			
			String quals = read.getBaseQualityString();
			while (quals.length() < READ_LENGTH) {
				quals += '#';
			}
			
			output(bases, quals, output);
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
		String out = args[1];
	
		new ReadExtractorFull2().run(in1, out);
	}
}
