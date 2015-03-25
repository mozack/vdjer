package vdj;


import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
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
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;



/**
 * Extracts all reads mapping to specified regions or that are unmapped or have mate unmapped
 * along with their mates.
 * @author lmose
 *
 */
public class CopyOfVdjReadExtractor {
	
	private Set<String> readIds = new HashSet<String>();
	private Set<String> reads = new HashSet<String>();
	private int numReads = 0;
	
	private static final int CONSTANT_START = 113256204;
	private static final int CONSTANT_END = 113422730;
	private static final int VDJ_START = CONSTANT_END+1;
	private static final int VDJ_END = 116009954;

	public void extract(String inFile, String outBam) throws IOException {
		SAMFileReader reader = new SAMFileReader(new File(inFile));
		reader.setValidationStringency(ValidationStringency.SILENT);
		
		SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
		
		SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(
				reader.getFileHeader(), false, new File(outBam)); 
		
		CloseableIterator<SAMRecord> iter = reader.queryOverlapping("chr12", VDJ_START, VDJ_END);
		while (iter.hasNext()) {
			SAMRecord read = iter.next();
			String id = read.getReadName() + "___" + (read.getFirstOfPairFlag() ? "1" : "2");
			if (!readIds.contains(id)) {
				reads.add(read.getReadName());
				readIds.add(id);
				writer.addAlignment(read);
				numReads += 1;
			}
		}
		iter.close();
		
		iter = reader.queryOverlapping("chr12", CONSTANT_START, CONSTANT_END);
		while (iter.hasNext()) {
			SAMRecord read = iter.next();
			String id = read.getReadName() + "___" + (read.getFirstOfPairFlag() ? "1" : "2");
			if (!readIds.contains(id)) {
				reads.add(read.getReadName());
				readIds.add(id);
				writer.addAlignment(read);
				numReads += 1;
			}
		}
		iter.close();
		
		iter = reader.queryUnmapped();
		while (iter.hasNext()) {
			SAMRecord read = iter.next();
			if (reads.contains(read.getReadName())) {
				String id = read.getReadName() + "___" + (read.getFirstOfPairFlag() ? "1" : "2");
				if (!readIds.contains(id)) {
					readIds.add(id);
					writer.addAlignment(read);
					numReads += 1;
				}
			}
		}
		iter.close();
		
		reader.close();
		writer.close();
		
		System.out.println("Processed: " + numReads + " read pairs");
	}
	
	
	public static void main(String[] args) throws Exception {
		String input = args[0];
		String output = args[1];
		
		new CopyOfVdjReadExtractor().extract(input, output);
	}
}
