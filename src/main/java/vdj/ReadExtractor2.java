package vdj;

import static vdj.VdjRegion.loadRegions;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * Extracts all reads mapping to specified regions or that are unmapped or have mate unmapped
 * along with their mates.
 * @author lmose
 *
 */
public class ReadExtractor2 {
	
	private List<VdjRegion> regions;
	private Map<String, ReadPair> readPairs = new HashMap<String, ReadPair>();
	private static final CompletedReadPair COMPLETED_READ_PAIR = new CompletedReadPair();
	private BufferedWriter output;
	private ReverseComplementor rc = new ReverseComplementor();
	private int numReads = 0;

	public void extract(String inFile, String outFile, String regionFile) throws IOException {
		regions = loadRegions(regionFile);
		SAMFileReader reader = new SAMFileReader(new File(inFile));
		reader.setValidationStringency(ValidationStringency.SILENT);
		
		output = new BufferedWriter(new FileWriter(outFile, false));
		
		// Identify reads that overlap the specified regions
		int readCount = 0;
		for (VdjRegion region : regions) {
			CloseableIterator<SAMRecord> iter = reader.queryOverlapping(region.chromosome, region.start, region.stop);
			while (iter.hasNext()) {
				SAMRecord read = iter.next();
				processRead(read);
				readCount += 1;
			}
			iter.close();
		}
		
		System.out.println("Processed: " + readCount + " aligned reads");
		
		// Identify reads where either end is unmapped or only one end is in a BCR region
		readCount = 0;
			
		for (SAMRecord read : reader) {
//			if ((read.getReadUnmappedFlag() || read.getMateUnmappedFlag()) || readPairs.containsKey(read.getReadName())) {
			if (readPairs.containsKey(read.getReadName())) {
				processRead(read);
			}
		}
		
		reader.close();
		output.close();
		
		System.out.println("Processed: " + numReads + " read pairs");
	}
	
	private void processRead(SAMRecord read) throws IOException {
		ReadPair readPair = readPairs.get(read.getReadName());
		if (readPair == null) {
			readPair = new InProgressReadPair();
			readPairs.put(read.getReadName(), readPair);
		}
		
		if (!readPair.isComplete() && !read.getCigarString().contains("H")) {
			
			if (read.getFirstOfPairFlag()) {
				readPair.setRead1(read);
			} else {
				readPair.setRead2(read);
			}
			
			if (readPair.isReadyForOutput()) {
				output(readPair);
				// Store marker indicating this fragment is complete
				readPairs.put(read.getReadName(), COMPLETED_READ_PAIR);
				numReads += 1;
			}
		}
	}
	
	private String getOutputString(boolean negativeStrandFlag, String bases, String qualities) {
		StringBuffer readBuffer = new StringBuffer();
		readBuffer.append(negativeStrandFlag ? "1" : "0");
		
		readBuffer.append(bases);
		readBuffer.append(qualities);
		
		return readBuffer.toString();
	}
	
	private void output(ReadPair readPair) throws IOException {
		output(readPair.getRead1());
		output(readPair.getRead2());
	}
	
	private void output(SAMRecord read) throws IOException {
		output.append(getOutputString(false, read.getReadString(), read.getBaseQualityString()));
		
		// We don't know the orientation of the read, so reverse complement it as well.
		output.append(getOutputString(false, rc.reverseComplement(read.getReadString()), rc.reverse(read.getBaseQualityString())));
	}
		
	interface ReadPair {
		SAMRecord getRead1();
		SAMRecord getRead2();
		public void setRead1(SAMRecord read1);
		public void setRead2(SAMRecord read2);
		boolean isReadyForOutput();
		boolean isComplete();
	}
	
	static class InProgressReadPair implements ReadPair {
		SAMRecord read1;
		SAMRecord read2;
		
		public SAMRecord getRead1() {
			return read1;
		}
		
		public SAMRecord getRead2() {
			return read2;
		}
		
		public void setRead1(SAMRecord read1) {
			this.read1 = read1;
		}
		
		public void setRead2(SAMRecord read2) {
			this.read2 = read2;
		}
		
		public boolean isReadyForOutput() {
			return read1 != null && read2 != null;
		}
		
		public boolean isComplete() {
			return false;
		}
	}
	
	static class CompletedReadPair implements ReadPair {
		
		public SAMRecord getRead1() {
			throw new UnsupportedOperationException();
		}
		
		public SAMRecord getRead2() {
			throw new UnsupportedOperationException();
		}
		
		public void setRead1(SAMRecord read1) {
			throw new UnsupportedOperationException();
		}
		
		public void setRead2(SAMRecord read2) {
			throw new UnsupportedOperationException();
		}
		
		public boolean isComplete() {
			return true;
		}
		
		public boolean isReadyForOutput() {
			return false;
		}
	}
	
	public static void main(String[] args) throws Exception {
		String input = args[0];
		String output = args[1];
		String regions = args[2];
		
		new ReadExtractor2().extract(input, output, regions);
	}
}
