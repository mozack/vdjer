package vdj;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;

public class ReferenceExtractor {

	ReverseComplementor rc = new ReverseComplementor(); 
	
	public void extract(String reference, String outFile, String regionFile, boolean rcInBothDirections) throws IOException {
		
		BufferedWriter output = new BufferedWriter(new FileWriter(outFile, false));
		CompareToReference2 c2r = new CompareToReference2();
		c2r.init(reference);
		List<VdjRegion> regions = VdjRegion.loadRegions(regionFile, true);
		
		for (VdjRegion region : regions) {
			String seq = c2r.getSequence(region.chromosome, region.start+1, region.stop-region.start-1);
		
			if (rcInBothDirections) {
				output.write(">" + region.name + "_+\n");
				output.write(seq + "\n");
				
				seq = rc.reverseComplement(seq);
				output.write(">" + region.name + "_-\n");
				output.write(seq + "\n");
				
			} else {
				if (region.strand == '-') {
					seq = rc.reverseComplement(seq);
				}
				
				output.write(">" + region.name + "\n");
				output.write(seq + "\n");
			}
		}
		
		output.close();
		
		System.out.println("Done.");
	}
	
	public static void main(String[] args) throws Exception {
		String reference = args[0];
		String outFile = args[1];
		String regionFile = args[2];
		boolean rc = false;
		
		if (args.length > 3) {
			rc = Boolean.parseBoolean(args[3]);
		}
		
		ReferenceExtractor re = new ReferenceExtractor();
		re.extract(reference, outFile, regionFile, rc);
	}
}
