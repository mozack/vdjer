package vdj;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.ValidationStringency;

import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;


/** 
 * Gathers read names from readIdBam and extracts records from readBam matching the read ids
 * into outFile
 * 
 * @author lmose
 *
 */
public class CustomReadExtractor {
	
	private Set<String> readIds = new HashSet<String>();
	
	private Set<String> readIdWithSequences = new HashSet<String>();

	public void run(String readIdBam, String readBam, String outFile) throws IOException {
		SAMFileReader idReader = new SAMFileReader(new File(readIdBam));
		idReader.setValidationStringency(ValidationStringency.SILENT);

		for (SAMRecord read : idReader) {
			readIds.add(read.getReadName());
		}
		
		idReader.close();
		
		System.out.println("Loaded: " + readIds.size() + " read ids.");
		
		SAMFileReader reader = new SAMFileReader(new File(readBam));
		reader.setValidationStringency(ValidationStringency.SILENT);
		
		SAMFileWriterFactory writerFactory = new SAMFileWriterFactory();
		
		SAMFileWriter writer = writerFactory.makeSAMOrBAMWriter(
				reader.getFileHeader(), false, new File(outFile));

		for (SAMRecord read : reader) {
			String uniqueId = read.getReadName() + "_" + read.getReadString();
			if (readIds.contains(read.getReadName()) && !(readIdWithSequences.contains(uniqueId))) {
				writer.addAlignment(read);
				readIdWithSequences.add(uniqueId);
			}
		}
		
		reader.close();
		writer.close();
	}
	
	public static void main(String[] args) throws IOException {
		String readIdBam = args[0];
		String readBam = args[1];
		String outFile = args[2];
		
//		String readIdBam = "/home/lmose/dev/vdj/extraction/selected_subset_pass2.k15.sort.bam";
//		String readBam = "/home/lmose/dev/vdj/extraction/igh_mapsplice_pairs.bam";
//		String outFile = "/home/lmose/dev/vdj/extraction/igh_mapsplice_check.bam";
		
		new CustomReadExtractor().run(readIdBam, readBam, outFile);
	}
}
