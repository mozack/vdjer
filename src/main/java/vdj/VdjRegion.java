package vdj;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class VdjRegion {
	String chromosome;
	int start;
	int stop;
	char strand = '0';
	String name;
	
	VdjRegion(String chromosome, int start, int stop) {
		this.chromosome = chromosome;
		this.start = start;
		this.stop = stop;
	}
	
	VdjRegion(String chromosome, int start, int stop, char strand, String name) {
		this(chromosome, start, stop);
		this.name = name;
		this.strand = strand;
	}
	
	public static List<VdjRegion> loadRegions(String regionFile) throws IOException {
		return loadRegions(regionFile, false);
	}
	
	public static List<VdjRegion> loadRegions(String regionFile, boolean isStranded) throws IOException {
		List<VdjRegion> regions = new ArrayList<VdjRegion>();
		BufferedReader reader = new BufferedReader(new FileReader(regionFile));
		
		String line = reader.readLine();
		while (line != null) {
			String[] fields = line.split("\\s");
			
			VdjRegion region = null;
			
			if (isStranded) {
				region = new VdjRegion(fields[0], Integer.parseInt(fields[1]), Integer.parseInt(fields[2]), fields[5].charAt(0), fields[3]);
			} else {
				region = new VdjRegion(fields[0], Integer.parseInt(fields[1]), Integer.parseInt(fields[2]));
			}
			
			regions.add(region);
			line = reader.readLine();
		}
		
		reader.close();
		
		return regions;
	}
}