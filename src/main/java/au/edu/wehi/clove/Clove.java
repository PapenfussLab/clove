package au.edu.wehi.clove;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.StringTokenizer;
import java.util.Map.Entry;
import java.util.TreeSet;


import htsjdk.samtools.util.Tuple;

import htsjdk.samtools.*;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalList;
import htsjdk.samtools.util.SamLocusIterator;


class EventIterator implements Iterator<Event> {
	
	private Iterator<GenomicNode> myNodeIterator;
	private GenomicNode currentNode;
	private int nextEventIndex;
	private HashSet<Event> skipEvents;
	private Event currentEvent;

	public EventIterator(Iterator<GenomicNode> nodeIterator, HashSet<Event> skip){
		myNodeIterator = nodeIterator;
		currentNode = nodeIterator.next();
		nextEventIndex = 0;
		skipEvents = skip;
	}
	
	@Override
	public boolean hasNext() {
		System.err.println("Method in EventIterator should not be used!");
		return false;
	}

	@Override
	public Event next() {
		if(currentNode == null){
			return null;
		}
		if(nextEventIndex < currentNode.getEvents().size()) {
			currentEvent = currentNode.getEvents().get(nextEventIndex);
			nextEventIndex ++;
			if(skipEvents.contains(currentEvent))
				return this.next();
			else 
				return currentEvent; 
		}
		nextEventIndex = 0;
		currentNode = (myNodeIterator.hasNext()? myNodeIterator.next() : null);
		while(currentNode!=null && currentNode.getEvents().size() == 0){
			if (myNodeIterator.hasNext())
				currentNode = myNodeIterator.next();
			else
				currentNode = null;
		}	
		if(currentNode != null)
			return this.next();
		else 
			return null;
	}
	
	public GenomicCoordinate getInsertionCoordinate() {
		if(currentEvent.getNode(true) == currentNode){
			return currentEvent.getC1();
		} else {
			return currentEvent.getC2();
		}
	}

	@Override
	public void remove() {
		System.err.println("Method in EventIterator should not be used!");
	}
	
}

public class Clove {
	private static int softclipLength5prime(SAMRecord s){
		int sc_start = s.getAlignmentStart() - s.getUnclippedStart();
		return sc_start;
	}
	private static int softclipLength3prime(SAMRecord s){
		int sc_end = s.getUnclippedEnd() - s.getAlignmentEnd();
		return sc_end;
	}
	private static boolean isInteresingSoftclip(SAMRecord s, int start, int end, String orientation){
		if(orientation.equals("+") && s.getAlignmentEnd() <= end + 2 && softclipLength5prime(s) > 4){
			return true;
		} 
		if(orientation.equals("-") && s.getAlignmentStart() >= start -2 && softclipLength3prime(s) > 4) {
			return true;
		}
//		if(s.getAlignmentStart()-2 <= breakpoint && s.getAlignmentEnd()+2>=breakpoint && (softclipLength3prime(s) > 4 || softclipLength5prime(s) > 4))
//			return true;
		return false;
	}
	private static boolean isAlignedAcrossBreakpoint(SAMRecord s, int breakpointPosition){
		if(s.getAlignmentStart() < breakpointPosition-5 && s.getAlignmentEnd() > breakpointPosition+5 && s.getAlignmentStart() - s.getUnclippedStart() < 10 && s.getUnclippedEnd() - s.getAlignmentEnd() < 10)
			return true;
		return false;
	}
	private static double fetchStuff(
			String chr, int start, int end, String orientation,  SAMFileReader bamFile) {
		
		SAMRecordIterator iter = bamFile.queryOverlapping(chr, start, end);
		
		int properAlignment = 0,  softclipped = 0;
		
		int breakpointPosition = (start + end)/2;
		
		for(SAMRecord s; iter.hasNext();){
			s = iter.next();
			if(isInteresingSoftclip(s, start, end, orientation)){
				softclipped ++;
			} else if (isAlignedAcrossBreakpoint(s, breakpointPosition)){
				properAlignment++;
			}
		}
		
		return (double)softclipped/(softclipped+properAlignment);
	}
	
	private static void addEventToNodeList(Event e, Hashtable<String, TreeSet<GenomicNode>> genomicNodes, boolean useFirstCoordinate){
		GenomicCoordinate coordinate = (useFirstCoordinate? e.getC1() : e.getC2());
		TreeSet<GenomicNode> nodeSet;
		if(! genomicNodes.containsKey(coordinate.getChr())){
			nodeSet = new TreeSet<GenomicNode>();
			genomicNodes.put(coordinate.getChr(), nodeSet);
		} else {
			nodeSet = genomicNodes.get(coordinate.getChr()) ;
		}
		GenomicNode newNode = new GenomicNode(coordinate, e);
		e.setNode(newNode, useFirstCoordinate);
		nodeSet.add(newNode);

	}
	
	private static String generateNodeLabel(GenomicNode n){
		return n.getStart().getChr()+"_"+n.getStart().getPos()+"_"+n.getEnd().getPos();
	}
	private static void graphVisualisation(String outputFilename, Hashtable<String, TreeSet<GenomicNode>> genomicNodes) throws IOException{
		HashSet<Event> eventsWritten = new HashSet<Event>();
		FileWriter output = new FileWriter(outputFilename);
		output.write("digraph g {\n");
		for(Entry<String, TreeSet<GenomicNode>> tableEntry: genomicNodes.entrySet()) {
			if(!tableEntry.getKey().equals("ecoli"))
				continue;
			output.write("{rank=same; ");
			for(GenomicNode n : tableEntry.getValue()){
				output.write(generateNodeLabel(n)+"; ");
			}
			output.write("}\n");
			for(GenomicNode n : tableEntry.getValue()){
				for(Event e: n.getEvents()) {
					if(eventsWritten.contains(e))
						continue;
					else
						eventsWritten.add(e);
					String l1 = generateNodeLabel(e.getNode(true)), l2 = generateNodeLabel(e.getNode(false));
					switch(e.getType()){
					case DEL: output.write(l1+"->"+l2+"[label=\"DEL\"];\n"); break;
					case TAN: output.write(l1+"->"+l2+"[label=\"TAN\" arrowtail=normal arrowhead=none dir=both];\n"); break;
					case COMPLEX_INVERSION: 
						Event ee = ((ComplexEvent)e).getEventsInvolvedInComplexEvent()[0];
						l1 = generateNodeLabel(ee.getNode(true));
						l2 = generateNodeLabel(ee.getNode(false));
						//if (l1.equals("ecoli_28204143_28204160")){
						//if (l1.equals("ecoli_24118329_24118583")){
							System.out.println("COMPLEX_INV:");
							System.out.println("l1:\t"+generateNodeLabel(ee.getNode(true))+"\t");
							System.out.println("l2:\t"+generateNodeLabel(ee.getNode(false))+"\t");
						//}
						output.write(l1+"->"+l2+"[label=\"COMPLEX_INV\"  dir=both];\n"); break;//should we differentiate COMPLEX_INV and INV?
					case COMPLEX_TRANSLOCATION: 
					case COMPLEX_DUPLICATION:
						GenomicNode insNode = e.getNode(true);
						String label = (e.getType() == EVENT_TYPE.COMPLEX_DUPLICATION? "DUP" : "TRANS");
						l1 = generateNodeLabel(insNode);
						Event[] allEvents = ((ComplexEvent)e).getEventsInvolvedInComplexEvent();
						for(int i=0;i<allEvents.length;i++){
							if(allEvents[i].getNode(true) == insNode){
								l2 = generateNodeLabel(allEvents[i].getNode(false));
								output.write(l1+"->"+l2+"[label=\""+label+"\" arrowtail=odiamond arrowhead=normal dir=both];\n");
							} else if(allEvents[i].getNode(false) == insNode){
								l2 = generateNodeLabel(allEvents[i].getNode(true));
								output.write(l1+"->"+l2+"[label=\""+label+"\" arrowtail=odiamond arrowhead=normal dir=both];\n");
							}
							
						}
						break;
					default: output.write(l1+"->"+l2+"[label=\""+e.getType()+"\"];\n");
					}
				}
			}
		}
		output.write("}\n");
		output.flush();
		output.close();
	}
	
	
	private static void compareToGoldStandard(String goldFileName, Hashtable<String, TreeSet<GenomicNode>> genomicNodes, int margin, boolean compareStrictly) throws IOException {
		boolean checkAgain = true;
		if(oldFns.size() == 0){
			checkAgain = false;
		}
		
		BufferedReader gold = new BufferedReader(new FileReader(goldFileName));
		String goldLine = gold.readLine();
		String currentChromosome = goldLine.replace(":","\t").split( "\t")[1];
		System.out.println("Working on first chromosome: "+currentChromosome);
		//Iterator<GenomicNode> iter = genomicNodes.get("gi|260447279|gb|CP001637.1|").iterator();
		
		HashSet<Event> skip = new HashSet<Event>();
		HashSet<Event> tryAgain = new HashSet<Event>();
		HashSet<String> recalledOnce = new HashSet<String>();
		Iterator<GenomicNode> iter = genomicNodes.get(currentChromosome).iterator();
		EventIterator events = new EventIterator(iter, skip);
		Event e = events.next();
		
		Hashtable<EVENT_TYPE, int[]> statsByType = new Hashtable<EVENT_TYPE, int[]>();
		for(EVENT_TYPE t: EVENT_TYPE.values()){
			//the convention used below is: TP index 0, FP 1, and FN 2
			statsByType.put(t, new int[4]);
		}
			//static conversion table
			Hashtable<String, EVENT_TYPE> typeConversion = new Hashtable<String, EVENT_TYPE>();
		{
			typeConversion.put("INVERSION", EVENT_TYPE.COMPLEX_INVERSION);
			typeConversion.put("DELETION", EVENT_TYPE.DEL);
			typeConversion.put("TANDEM", EVENT_TYPE.TAN);
			//typeConversion.put("INSERTION", EVENT_TYPE.INS);
			typeConversion.put("INSERTION", EVENT_TYPE.COMPLEX_DUPLICATION);
			typeConversion.put("TRANSLOCATION", EVENT_TYPE.COMPLEX_TRANSLOCATION);
			typeConversion.put("INVERTED_TRANSLOCATION", EVENT_TYPE.COMPLEX_INVERTED_TRANSLOCATION);
			typeConversion.put("INVERTED_INSERTION", EVENT_TYPE.COMPLEX_INVERTED_DUPLICATION);
			typeConversion.put("DUPLICATION", EVENT_TYPE.COMPLEX_DUPLICATION);
			typeConversion.put("INVERTED_DUPLICATION", EVENT_TYPE.COMPLEX_INVERTED_DUPLICATION);
			typeConversion.put("INTERCHROMOSOMAL_TRANSLOCATION", EVENT_TYPE.COMPLEX_INTERCHROMOSOMAL_TRANSLOCATION);
			typeConversion.put("INTERCHROMOSOMAL_INVERTED_TRANSLOCATION", EVENT_TYPE.COMPLEX_INTERCHROMOSOMAL_INVERTED_TRANSLOCATION);
			typeConversion.put("INTERCHROMOSOMAL_DUPLICATION", EVENT_TYPE.COMPLEX_INTERCHROMOSOMAL_DUPLICATION);
			typeConversion.put("INTERCHROMOSOMAL_INVERTED_DUPLICATION", EVENT_TYPE.COMPLEX_INTERCHROMOSOMAL_INVERTED_DUPLICATION);
		}
		
		
		while(goldLine != null ){
			StringTokenizer st = new StringTokenizer(goldLine, ":-\t ");
			String type = st.nextToken();
			String chr = st.nextToken();
			//TODO: what if the SVs are on a different chromosome?
			if(!currentChromosome.equals(chr)) {
				while(e!=null){
					if(tryAgain.contains(e)){
						System.out.println("FPFPFPFP: "+e);
						statsByType.get(e.getType())[1]++;
					} else {
						tryAgain.add(e);
					}
					e = events.next();
				}
				
				currentChromosome = chr;
				System.out.println("Working on chromosome: "+currentChromosome);
				iter = genomicNodes.get(currentChromosome).iterator();
				events = new EventIterator(iter, skip);
				e = events.next();
			}

			int start = Integer.parseInt(st.nextToken());
			//int end = Integer.parseInt(st.nextToken());
			if(type.equals("SNP") || type.equals("JOIN") || type.equals("TRANSLOCATION_DELETION")){
				goldLine = gold.readLine();
				continue;
			}
			if(e==null){
				//System.out.println("DEFINITE FN: "+goldLine);
				goldLine = gold.readLine();
				continue;
			} 
			
			GenomicCoordinate compare;
			if(e.getType() == EVENT_TYPE.COMPLEX_INVERTED_TRANSLOCATION || e.getType() == EVENT_TYPE.COMPLEX_INVERTED_DUPLICATION 
					|| e.getType() == EVENT_TYPE.COMPLEX_DUPLICATION || e.getType() == EVENT_TYPE.COMPLEX_TRANSLOCATION
					|| e.getType() == EVENT_TYPE.COMPLEX_INTERCHROMOSOMAL_INVERTED_TRANSLOCATION || e.getType() == EVENT_TYPE.COMPLEX_INTERCHROMOSOMAL_INVERTED_DUPLICATION
					|| e.getType() == EVENT_TYPE.COMPLEX_INTERCHROMOSOMAL_DUPLICATION || e.getType() == EVENT_TYPE.COMPLEX_INTERCHROMOSOMAL_TRANSLOCATION) {
				compare = ((ComplexEvent)e).getInsertionPoint();
				tryAgain.add(e); // do not attempt single coordinate events twice.
			} else if(tryAgain.contains(e)){
				compare = events.getInsertionCoordinate();
			} else {
				compare = events.getInsertionCoordinate();
			}
			if(!compare.getChr().equals(chr)) {
				//Fusion on different chr than goldLine
				if(compare.getChr().compareTo(chr) < 0){
					System.out.println("DEFAULT FP? "+e);
					e = events.next();
					continue;
				} else {
					System.out.println("DEFAULT FN?");
					goldLine= gold.readLine();
					continue;
				}
			}

			//GenomicCoordinate compare = (e.getType() == EVENT_TYPE.COMPLEX_INVERTED_TRANSLOCATION || e.getType() == EVENT_TYPE.COMPLEX_INVERTED_DUPLICATION || e.getType() == EVENT_TYPE.COMPLEX_DUPLICATION || e.getType() == EVENT_TYPE.COMPLEX_TRANSLOCATION? ((ComplexEvent)e).getInsertionPoint() : e.getC1());
			if(compare.distanceTo(new GenomicCoordinate(chr, start)) > margin) {
				if(compare.compareTo(new GenomicCoordinate(chr, start)) < 0) {
					//half TP?
					if(tryAgain.contains(e) || compareStrictly){
						//System.out.println("FP: "+e);
						statsByType.get(e.getType())[1]++;
					} else {
						tryAgain.add(e);
					}
					e = events.next();
				} else {
					if(!recalledOnce.contains(goldLine) || compareStrictly){
						System.out.println("FN: "+goldLine);
						if(checkAgain && !oldFns.contains(goldLine)){
							System. out.println("New FN: "+goldLine);
						} else if (!checkAgain){
							oldFns.add(goldLine);
						}
						statsByType.get(typeConversion.get(type))[2]++;
					}
					goldLine = gold.readLine();
					continue;
				}
			} else {
				if(typeConversion.get(type) == e.getType()){
//				if(type.equals("INVERSION") && e.getType()==EVENT_TYPE.COMPLEX_INVERSION || type.equals("DELETION") && e.getType()==EVENT_TYPE.DEL
//						|| type.equals("TANDEM") && e.getType()==EVENT_TYPE.TAN || type.equals("INSERTION") && e.getType()==EVENT_TYPE.COMPLEX_DUPLICATION
//						|| type.equals("TRANSLOCATION") && e.getType()==EVENT_TYPE.COMPLEX_TRANSLOCATION || type.equals("INSERTION") && e.getType()==EVENT_TYPE.INS) {
					//System.out.println("TP: "+e+" "+goldLine);
					if(recalledOnce.contains(goldLine)){
						//redundant TP?
						System.out.println("Redundant TP!");
					} else {
						statsByType.get(e.getType())[0]++;
					}
					//goldLine = gold.readLine();
					
					
				} else {
					//System.out.println("Half TP: Type mismatch: "+e+" "+goldLine);
					if(recalledOnce.contains(goldLine)){
						//redundant HTP?
						System.out.println("Redundant HTP!");
					} else {
						statsByType.get(e.getType())[3]++;
					}	
					
				}
				recalledOnce.add(goldLine);
				skip.add(e);
				e = events.next();
			}
			
		}
		while(goldLine!=null){
			if(! goldLine.contains("SNP") && ! goldLine.contains("TRANSLOCATION_DELETION") && (!recalledOnce.contains(goldLine) || compareStrictly)){
				String type = (new StringTokenizer(goldLine)).nextToken();
				System.out.println("FN: "+goldLine);
				statsByType.get(typeConversion.get(type))[2]++;
			}
			goldLine = gold.readLine();
		}
		e = events.next();
		
		int tps=0, fps=0, fns=0, htps=0;
		System.out.println("Stats:\tTP\tHalf TP\tFP\tFN\tSen\tSpe");
		for(EVENT_TYPE t: EVENT_TYPE.values()){
			int[] stats = statsByType.get(t);
			double sen = (stats[0]+stats[2]==0? 0: (double)stats[0]/(stats[0]+stats[2]));
			double spe = (stats[0]+stats[1]==0? 0: (double)stats[0]/(stats[0]+stats[1]));
			System.out.println(t+"\t"+stats[0]+"\t"+stats[3]+"\t"+stats[1]+"\t"+stats[2]+"\t"+sen+"\t"+spe);
			tps+=stats[0]; fps+=stats[1]; fns+=stats[2]; htps +=stats[3];
		}
		System.out.println("Total\t"+tps+"\t"+htps+"\t"+fps+"\t"+fns+"\t"+((double)(tps+htps)/(tps+htps+fns))+"\t"+((double)(tps+htps)/(tps+htps+fps)));
		System.out.println("Accuracy:\t"+((double)tps/(tps+fps+fns))+"\t"+((double)(tps+htps)/(tps+fps+fns+htps)));
		gold.close();
	}
	
	private static void reportEventComposition(Hashtable<String, TreeSet<GenomicNode>> genomicNodes) {
		Hashtable<EVENT_TYPE, Integer> eventCounts = new Hashtable<EVENT_TYPE, Integer>();
		int selfRef = 0;
		for(EVENT_TYPE t: EVENT_TYPE.values()){
			eventCounts.put(t, 0);
		}
		HashSet<Event> skipEvents = new HashSet<Event>();
		for(Entry<String, TreeSet<GenomicNode>> tableEntry: genomicNodes.entrySet()) {
			for(GenomicNode n: tableEntry.getValue()){
				for(Event e: n.getEvents()){
					if(skipEvents.contains(e))
						continue;
					Integer i = eventCounts.get(e.getType()) + 1;
					eventCounts.put(e.getType(), i);
					if(e.otherNode(n) == n){
						selfRef++;
					} else
						skipEvents.add(e);
				}
			}
		}
		for(EVENT_TYPE t: EVENT_TYPE.values()){
			System.out.println(t+": "+eventCounts.get(t));
		}
		System.out.println("Self refs: "+selfRef);
	}
	
	
	
	//private static double getReadDepth(String str, String chr, int start, int end){
	private static double getReadDepth(SAMFileReader samReader, String chr, int start, int end){
		
		if(start >= end){
			return -1;
		}
		
		//SAMFileReader  samReader=new  SAMFileReader(new  File(str));
        String chromId=chr;
        int chromStart=start;
        int chromEnd=end;
        int pos=0;
        int depth=0;
        int total=0;
        int count=0;
        Interval  interval=new  Interval(chromId,chromStart,chromEnd);
        IntervalList  iL=new  IntervalList(samReader.getFileHeader());
        iL.add(interval);

        SamLocusIterator  sli=new  SamLocusIterator(samReader,iL,true);

        for(Iterator<SamLocusIterator.LocusInfo> iter=sli.iterator(); iter.hasNext();){
            SamLocusIterator.LocusInfo  locusInfo=iter.next();
            //pos = locusInfo.getPosition();
            depth = locusInfo.getRecordAndPositions().size();
            total+=depth;
            count++;
            //System.out.println("POS="+pos+" depth:"+depth);
            }
        //System.out.println("total: "+total+"\tcount: "+count);
        sli.close();
        //samReader.close();
        
        return (double)total/count;
    }
	
	public static void createVCFHeader(PrintWriter output) throws IOException{
		  //file format
        output.write("##fileformat=VCFv4.2\n");
        
        //fileDate
        DateFormat dateFormat = new SimpleDateFormat ("yyyyMMdd");
        Date date = new Date();
        output.write("##fileDate="+dateFormat.format(date)+"\n");
        
        //INFO
        output.write("##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"PE confidence interval around END\">\n");
        output.write("##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"PE confidence interval around POS\">\n");
        output.write("##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for END coordinate in case of a translocation\">\n");
        output.write("##INFO=<ID=START>,Number=1,Type=Integer,Description=\"Start position of the interval (for certain types only)\">\n");
        output.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">\n");
        output.write("##INFO=<ID=PE,Number=1,Type=Integer,Description=\"Paired-end support of the structural variant\">\n");
        output.write("##INFO=<ID=MAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of paired-ends\">\n");
        output.write("##INFO=<ID=SR,Number=1,Type=Integer,Description=\"Split-read support\">\n");
        output.write("##INFO=<ID=SRQ,Number=1,Type=Float,Description=\"Split-read consensus alignment quality\">\n");
        output.write("##INFO=<ID=CONSENSUS,Number=1,Type=String,Description=\"Split-read consensus sequence\">\n");
        output.write("##INFO=<ID=CT,Number=1,Type=String,Description=\"Paired-end signature induced connection type\">\n");
        output.write("##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n");
        output.write("##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">\n");
        output.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n");
        output.write("##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Type of approach used to detect SV\">\n");
        output.write("##INFO=<ID=ADP,Number=1,Type=Integer,Description=\"Average Read Depth\">\n");
        
        //FILTER
        output.write("##FILTER=<ID=LowQual,Description=\"PE support below 3 or mapping quality below 20.\">\n");
        
        //FORMAT
        output.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        output.write("##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Log10-scaled genotype likelihoods for RR,RA,AA genotypes\">\n");
        output.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n");
        output.write("##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Per-sample genotype filter\">\n");
        output.write("##FORMAT=<ID=RC,Number=1,Type=Integer,Description=\"Raw high-quality read counts for the SV\">\n");
        output.write("##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# high-quality reference pairs\">\n");
        output.write("##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# high-quality variant pairs\">\n");
        output.write("##FORMAT=<ID=RR,Number=1,Type=Integer,Description=\"# high-quality reference junction reads\">\n");
        output.write("##FORMAT=<ID=RV,Number=1,Type=Integer,Description=\"# high-quality variant junction reads\">\n");

        //ALT
        output.write("##ALT=<ID=DEL,Description=\"Deletion\">\n");
        output.write("##ALT=<ID=TAN,Description=\"Tandem Duplication\">\n");
        output.write("##ALT=<ID=INV,Description=\"Inversion\">\n");
        output.write("##ALT=<ID=INS,Description=\"Insertion\">\n");
        output.write("##ALT=<ID=DUP,Description=\"Complex Duplication\">\n");
        output.write("##ALT=<ID=TRA,Description=\"Complex Translocation\">\n");
        output.write("##ALT=<ID=CIV,Description=\"Complex Inversion\">\n");
        output.write("##ALT=<ID=CVT,Description=\"Complex Inverted Translocation\">\n");
        output.write("##ALT=<ID=CVD,Description=\"Complex Inverted Duplication\">\n");
        output.write("##ALT=<ID=CIT,Description=\"Complex Interchromosomal Translocation\">\n");
        output.write("##ALT=<ID=CID,Description=\"Complex Interchromosomal Duplication\">\n");
        output.write("##ALT=<ID=IVT,Description=\"Complex Inverted Interchromosomal Translocation\">\n");
        output.write("##ALT=<ID=IVD,Description=\"Complex Inverted Interchromosomal Duplication\">\n");

        //Header Line
        output.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n");
	}
	
	enum SV_ALGORITHM {SOCRATES, DELLY, CREST, GUSTAF, BEDPE};
	
	
	static ArrayList<String> oldFns = new ArrayList<String>();
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		//Start Time
		long startTime = System.nanoTime();
	
		if(args.length < 8){
			System.err.println("Options (all mandatory -- input can be specified more than once):" +
					"\n\t-i <list of breakpoints> <algorithm (Socrates/Delly/Crest/Gustaf/BEDPE)>" +
					"\n\t-b <BAM file> \n\t-c <mean coverage> <coverage>" +
					"\n\t-o <output filename> [default: CLOVE.vcf]");
			System.exit(0);
		}
		
		/*parse the options from the command line */
		int argindex = 0;
		ArrayList<Tuple<BufferedReader, SV_ALGORITHM>> inputs = new ArrayList<Tuple<BufferedReader,SV_ALGORITHM>>();
		SAMFileReader  samReader = null;
		double mean = 0;
		double interval= 0;
		String goldStandard = null;
		String outputVCF = "CLOVE.vcf";
		while (argindex < args.length){
			if (args[argindex].equals("-i")){
				try{
					BufferedReader input = new BufferedReader(new FileReader(args[argindex + 1]));
					SV_ALGORITHM algorithm = SV_ALGORITHM.valueOf(args[argindex + 2].toUpperCase());
					inputs.add(new Tuple<BufferedReader, Clove.SV_ALGORITHM>(input, algorithm));
					argindex += 3;
				} catch (IllegalArgumentException e){
					System.err.println("Unable to parse input breakpoints.");
					System.exit(1);
				}
			} else if (args[argindex].equals("-b")){
				try {
					samReader=new  SAMFileReader(new  File(args[argindex + 1]));
					argindex += 2;
				} catch (IllegalArgumentException e){
					System.err.println("Unable to load bam file.");
					System.exit(1);
				}
			} else if(args[argindex].equals("-c")){
				try{
					mean = Double.parseDouble(args[argindex + 1]);
					interval = 2*Double.parseDouble(args[argindex + 2]);
					argindex += 3;
				} catch (IllegalArgumentException e){
					System.err.println("Unable to parse coverage and std.");
					System.exit(1);
				}
			} else if (args[argindex].equals("-d")){
				goldStandard = args[argindex + 1];
				argindex += 2;
			} else if(args[argindex].equals("-o")){
				outputVCF = args[argindex + 1];
				argindex += 2;
			}
			
			else {
				System.err.println("Unknown option: "+args[argindex]);
				throw new IllegalArgumentException();
			}
		}

		
		/*
		 * parse the entire input file and collect all events in list
		 */
		ArrayList<Event> allEvents = new ArrayList<Event>();

		String line;
		int count = 0;
		
		for(Tuple<BufferedReader, SV_ALGORITHM> input_tuple : inputs){
			System.out.println("Reading input...");
			BufferedReader input = input_tuple.a;
			SV_ALGORITHM algorithm = input_tuple.b;
			while ((line = input.readLine()) != null){
				//TODO: make # check algorithm specific?
				if(line.startsWith("#"))
					continue;
				Event e;
				switch(algorithm){
				case SOCRATES: 	e = Event.createNewEventFromSocratesOutputLatest(line, count++); 	break;
				case DELLY: 	e = Event.createNewEventFromDellyOutputLatest(line);break;
				case CREST:		e = Event.createNewEventFromCrestOutputLatest(line, count++); 		break;
				case GUSTAF: e = Event.createNewEventFromGustafOutput(line);	  if(e.size()<50) continue; break;
				case BEDPE: 	e = Event.createNewEventFromBEDPE(line); break;
				default:		e = null;
				}
				allEvents.add(e);
			}
			input.close();
		}
		System.out.println("Total events: "+allEvents.size());
		
		/*VCF Header*/
		//PrintWriter writer = new PrintWriter("/Users/schroeder/Downloads/VCF.txt", "UTF-8");
		//PrintWriter writer = new PrintWriter("/home/adrianto/Downloads/VCF.txt", "UTF-8");
		//PrintWriter writer = new PrintWriter("/home/users/allstaff/schroeder/tools/GenotypeBreakpoints/VCF.txt", "UTF-8");
		
		PrintWriter writer = new PrintWriter(outputVCF, "UTF-8");
		createVCFHeader(writer);
		
		/*
		 * Create nodes data structure that combines close events into the same 
		 * genomic node
		 */
		Hashtable<String, TreeSet<GenomicNode>> genomicNodes = new Hashtable<String, TreeSet<GenomicNode>>();
		
		//parse all events and create new nodes 
		for (Event e: allEvents){
			addEventToNodeList(e, genomicNodes, true);
			addEventToNodeList(e, genomicNodes, false);
		}
		
		//establish distance for "close" events according to algorithm
		int maxDistanceForNodeMerge = 15;
//		switch(algorithm){
//		case SOCRATES: 	maxDistanceForNodeMerge = 15; break;
//		case DELLY:		maxDistanceForNodeMerge = 100; break;
//		case CREST:		maxDistanceForNodeMerge = 15; break;
//		case GUSTAF:	maxDistanceForNodeMerge = 15; break;
//		default:		System.err.println("Node merge distance set to 0!");
//		}
		//static parameter to classify single inversions as FP or TP
		final boolean classifySimpleInversion = false;
		
		//iterate through node sets and merge nodes where necessary
		//also checks each node for redundant members
		//TODO: handle redundant members
		for(Entry<String, TreeSet<GenomicNode>> tableEntry: genomicNodes.entrySet()) {
			int nodesMerged = 0;
			GenomicNode[] staticList = new GenomicNode[tableEntry.getValue().size()];
			tableEntry.getValue().toArray(staticList);
			if(staticList.length == 0)
				break;
			GenomicNode lastNode = staticList[0], currentNode = null;
			for(int i = 1; i < staticList.length; i++){
				currentNode = staticList[i];
				if(currentNode.getStart().distanceTo(lastNode.getEnd()) < maxDistanceForNodeMerge){
//					System.out.println("Merging Node at "+lastNode.getStart()+" with node at "+currentNode.getStart());
//					for(Event e:lastNode.getEvents()){
//						System.out.println(e);
//					}
//					System.out.println("----");
//					for(Event e:currentNode.getEvents()){
//						System.out.println(e);
//					}
					lastNode.mergeWithNode(currentNode);	
					if(!tableEntry.getValue().remove(currentNode))
						System.out.println("NO CAN DO");
					nodesMerged++;
				} else {
					lastNode.checkForRedundantEvents(maxDistanceForNodeMerge);
					lastNode = currentNode;
				}
			}
			lastNode.checkForRedundantEvents(maxDistanceForNodeMerge);
			System.out.println("Nodes merged: "+nodesMerged);
		}
		System.out.println("Events merged: "+GenomicNode.global_event_merge_counter);
		
		//String goldStandard = args[1].substring(0, 22)+"_2.fa";
		//String goldStandard = "/home/users/allstaff/schroeder/GenotypeBreakpoints/data/ecoli/SV_list_2.txt";
		

		//compareToGoldStandard(goldStandard, genomicNodes, 150, true);
		if(goldStandard != null)
			compareToGoldStandard(goldStandard, genomicNodes, 150, false);
		
		String tempInfo = null;
		//iterate through node sets again, and genotype events
		for(Entry<String, TreeSet<GenomicNode>> tableEntry: genomicNodes.entrySet()) {
			System.out.println("Nodes on chr:"+tableEntry.getValue().size());
			for(GenomicNode currentNode: tableEntry.getValue()){
				//iterate through all event-event pairing in this node and assess for complex events
				Event e1, e2;
				HashSet<Event> removeEvents = new HashSet<Event>();
				HashSet<ComplexEvent> newComplexEvents = new HashSet<ComplexEvent>();
				ComplexEvent newComplexEvent = null;
				for(int i=0; i<currentNode.getEvents().size(); i++){
					e1 = currentNode.getEvents().get(i);
					for(int j=0; j<currentNode.getEvents().size(); j++){
						e2 = currentNode.getEvents().get(j);
						if(e1 == e2 || removeEvents.contains(e2) || removeEvents.contains(e1) 
								|| e1.otherNode(currentNode) == currentNode || e2.otherNode(currentNode) == currentNode)
							continue;
						switch(e1.getType()){
							//inversions
							case INV1: {
								if(e2.getType() == EVENT_TYPE.INV2 && Event.sameNodeSets(e1, e2)){
									//System.out.println("Complex inversion between "+e1+" and "+e2);
									GenomicCoordinate invstart = (e1.getC1().compareTo(e1.getC2()) < 0? e1.getC1() : e1.getC2());
									GenomicCoordinate invend   = (e2.getC2().compareTo(e2.getC1()) < 0? e2.getC1() : e2.getC2());
									//System.out.println(e1.getC1()+"\t"+e1.getC2()+"\t"+e2.getC1()+"\t"+e2.getC2()+"\t"+invstart+"\t"+invend);
									newComplexEvent = new ComplexEvent(invstart, invend, EVENT_TYPE.COMPLEX_INVERSION, (new Event[] {e1, e2}), currentNode);
									//newComplexEvent.setId(e1.getId());
									newComplexEvent.setCoord(invstart);
									newComplexEvent.setId(e1.getId()+"+"+e2.getId());
									newComplexEvent.setRef(e1.getRef());
									newComplexEvent.setAlt("<CIV>");
									newComplexEvent.setFilter(e1.getFilter());
									newComplexEvent.setQual(e1.getQual());
//									tempInfo = e1.getInfo();
//									//System.out.println(tempInfo+"\n");
//									if (tempInfo.contains("SVTYPE")){
//										tmpOld = tempInfo.substring(tempInfo.indexOf("SVTYPE=")+7, tempInfo.indexOf("SVTYPE=")+10);
//										tmpNew = newComplexEvent.getAlt();
//										tempInfo=tempInfo.replace(tmpOld, tmpNew);
//										tmpOld = tempInfo.substring(tempInfo.indexOf("CHR2=")+5, tempInfo.indexOf(";",tempInfo.indexOf("CHR2=")));
//										tmpNew = invend.getChr();
//										tempInfo=tempInfo.replace(tmpOld, tmpNew);
//										tmpOld = tempInfo.substring(tempInfo.indexOf(";END=")+5, tempInfo.indexOf(";CT",tempInfo.indexOf(";END=")));
//										tmpNew = Integer.toString(invend.getPos());
//										tempInfo.replace(tmpOld, tmpNew);
//										newComplexEvent.setInfo(tempInfo);
//									}
									double readDepth = getReadDepth(samReader, invstart.getChr(), invstart.getPos(), invend.getPos());
									//tempInfo="SVTYPE="+newComplexEvent.getAlt().substring(1, 4)+"; CHR2="+invend.getChr()+"; END="+Integer.toString(invend.getPos());
									tempInfo="SVTYPE="+newComplexEvent.getAlt().substring(1, 4)+"; CHR2="+invend.getChr()+"; END="+Integer.toString(invend.getPos())+"; ADP="+readDepth;
									newComplexEvent.setInfo(tempInfo);
																									
									//writer.write(newComplexEvent.getC1().getChr()+"\t"+invstart+"\t"+newComplexEvent.getId()+"\t"+newComplexEvent.getRef()+"\t"+newComplexEvent.getAlt()+"\t"+newComplexEvent.getQual()+"\t"+newComplexEvent.getFilter()+"\t"+newComplexEvent.getInfo()+"\n");
									//currentNode?
									//System.out.println(currentNode.getStart().toString());
									//	System.out.println(currentNode.getEnd().toString());
								}
								else if(e2.getType() == EVENT_TYPE.INV2 ){
									GenomicNode other1 = e1.otherNode(currentNode), other2 = e2.otherNode(currentNode);
									if(other1.compareTo(other2) > 0){
										Event e3 = other1.existsDeletionEventTo(other2);
										if(e3 != null){
											 GenomicCoordinate invstart = (e2.getNode(true) == currentNode ? e2.getC2() : e2.getC1() ),
											                	 invend   = (e1.getNode(true) == currentNode ? e1.getC2() : e1.getC1() ),
																				 insert   = (e1.getNode(true) == currentNode ? e1.getC1() : e1.getC2() );
											 newComplexEvent = new ComplexEvent(invstart, invend, EVENT_TYPE.COMPLEX_INVERTED_TRANSLOCATION, (new Event[] {e1, e2, e3}), currentNode, insert);
											 //newComplexEvent.setId(e1.getId());
											 newComplexEvent.setCoord(invstart);
											 newComplexEvent.setId(e1.getId()+"+"+e2.getId());
											 newComplexEvent.setRef(e1.getRef());
											 newComplexEvent.setAlt("<CVT>");
											 newComplexEvent.setFilter(e1.getFilter());
											 newComplexEvent.setQual(e1.getQual());
//											 tempInfo = e1.getInfo();
//											 //System.out.println(tempInfo+"\n");
//											 if (tempInfo.contains("SVTYPE")){
//												 tmpOld = tempInfo.substring(tempInfo.indexOf("SVTYPE=")+7, tempInfo.indexOf("SVTYPE=")+10);
//												 tmpNew = newComplexEvent.getAlt();
//												 tempInfo=tempInfo.replace(tmpOld, tmpNew);
//												 tmpOld = tempInfo.substring(tempInfo.indexOf("CHR2=")+5, tempInfo.indexOf(";",tempInfo.indexOf("CHR2=")));
//												 tmpNew = invend.getChr();
//												 tempInfo=tempInfo.replace(tmpOld, tmpNew);
//												 tmpOld = tempInfo.substring(tempInfo.indexOf(";END=")+5, tempInfo.indexOf(";CT",tempInfo.indexOf(";END=")));
//												 tmpNew = Integer.toString(invend.getPos());
//												 tempInfo.replace(tmpOld, tmpNew);
//												 newComplexEvent.setInfo(tempInfo);
//											 }
											 double readDepth = getReadDepth(samReader, invstart.getChr(), invstart.getPos(), invend.getPos());
											 //tempInfo="SVTYPE="+newComplexEvent.getAlt().substring(1, 4)+"; CHR2="+invend.getChr()+"; END="+Integer.toString(invend.getPos());
											 tempInfo="SVTYPE="+newComplexEvent.getAlt().substring(1, 4)+"; CHR2="+invend.getChr()+"; END="+Integer.toString(invend.getPos())+"; ADP="+readDepth;
											 newComplexEvent.setInfo(tempInfo);
											 //writer.write(newComplexEvent.getC1().getChr()+"\t"+invstart+"\t"+newComplexEvent.getId()+"\t"+newComplexEvent.getRef()+"\t"+newComplexEvent.getAlt()+"\t"+newComplexEvent.getQual()+"\t"+newComplexEvent.getFilter()+"\t"+newComplexEvent.getInfo()+"\n");
										} else {
											//System.out.println("INVDUP!"+e1+e2);
											GenomicCoordinate invstart = (e2.getNode(true) == currentNode ? e2.getC2() : e2.getC1() ),
															  invend   = (e1.getNode(true) == currentNode ? e1.getC2() : e1.getC1() ),
															  insert   = (e1.getNode(true) == currentNode ? e1.getC1() : e1.getC2() );
											newComplexEvent = new ComplexEvent(invstart, invend, EVENT_TYPE.COMPLEX_INVERTED_DUPLICATION, (new Event[] {e1, e2}), currentNode, insert);
											//newComplexEvent.setId(e1.getId());
											newComplexEvent.setCoord(invstart);
											newComplexEvent.setId(e1.getId()+"+"+e2.getId());
											newComplexEvent.setRef(e1.getRef());
											newComplexEvent.setAlt("<CVD>");
											newComplexEvent.setFilter(e1.getFilter());
											newComplexEvent.setQual(e1.getQual());
//											tempInfo = e1.getInfo();
//											//System.out.println(tempInfo+"\n");
//											if (tempInfo.contains("SVTYPE")){
//												tmpOld = tempInfo.substring(tempInfo.indexOf("SVTYPE=")+7, tempInfo.indexOf("SVTYPE=")+10);
//												tmpNew = newComplexEvent.getAlt();
//												tempInfo=tempInfo.replace(tmpOld, tmpNew);
//												tmpOld = tempInfo.substring(tempInfo.indexOf("CHR2=")+5, tempInfo.indexOf(";",tempInfo.indexOf("CHR2=")));
//												tmpNew = invend.getChr();
//												tempInfo=tempInfo.replace(tmpOld, tmpNew);
//												tmpOld = tempInfo.substring(tempInfo.indexOf(";END=")+5, tempInfo.indexOf(";CT",tempInfo.indexOf(";END=")));
//												tmpNew = Integer.toString(invend.getPos());
//												tempInfo.replace(tmpOld, tmpNew);
//												newComplexEvent.setInfo(tempInfo);
//											}
											double readDepth = getReadDepth(samReader, invstart.getChr(), invstart.getPos(), invend.getPos());
											//tempInfo="SVTYPE="+newComplexEvent.getAlt().substring(1, 4)+"; CHR2="+invend.getChr()+"; END="+Integer.toString(invend.getPos());
//											tempInfo="SVTYPE="+newComplexEvent.getAlt().substring(1, 4)+"; CHR2="+invend.getChr()+"; END="+Integer.toString(invend.getPos())+"; ADP="+readDepth;
											tempInfo="SVTYPE="+newComplexEvent.getAlt().substring(1, 4)+"; CHR2="+invstart.getChr()+"; START="+Integer.toString(invstart.getPos())+"; END="+Integer.toString(invend.getPos())+"; ADP="+readDepth;
											newComplexEvent.setInfo(tempInfo);
											newComplexEvent.setCoord(insert);
											//writer.write(newComplexEvent.getC1().getChr()+"\t"+invstart+"\t"+newComplexEvent.getId()+"\t"+newComplexEvent.getRef()+"\t"+newComplexEvent.getAlt()+"\t"+newComplexEvent.getQual()+"\t"+newComplexEvent.getFilter()+"\t"+newComplexEvent.getInfo()+"\n");
										}
									}
								}
								else {
									//unknown pairing
								}
								break;
							}
							//duplications and translocations
							case DEL: {
								if(e2.getType() == EVENT_TYPE.TAN){
									GenomicNode other1 = e1.otherNode(currentNode), other2 = e2.otherNode(currentNode);
									if(other1.compareTo(other2) < 0 && currentNode.compareTo(other1) < 0
											|| other2.compareTo(other1) < 0 && other1.compareTo(currentNode) < 0){
										Event e3 = other1.existsDeletionEventTo(other2);
										if(e3 != null){
											//System.out.println("Translocation between "+e1+ " and "+ e2);
											//intrachromosomal translocations are actually ambiguous as to where they come from and got to
											//as a convention, we call the smaller bit the translocated one inserted into the larger bit.
											if(e1.size() < e3.size()) {
												//area under e1 is translocated
												GenomicCoordinate transtart, tranend;
												if(e1.getC1().compareTo(e1.getC2()) < 0) { // don't assume ordered coordinates
													transtart = e1.getC1();
													tranend   = e1.getC2();
												} else {
													transtart = e1.getC2();
													tranend   = e1.getC1();
												}
												GenomicCoordinate traninsert = (e2.getNode(true) == currentNode? e2.getC2() : e2.getC1());
												GenomicNode hostingNode = (e2.getNode(true) == currentNode? e2.getNode(false) : e2.getNode(true));
												newComplexEvent = new ComplexEvent(transtart, tranend, EVENT_TYPE.COMPLEX_TRANSLOCATION, (new Event[] {e1, e2, e3}), hostingNode, traninsert);
												//newComplexEvent.setId(e1.getId());
												newComplexEvent.setCoord(transtart);
												newComplexEvent.setId(e1.getId()+"+"+e2.getId());
												newComplexEvent.setRef(e1.getRef());
												newComplexEvent.setAlt("<TRA>");
												newComplexEvent.setFilter(e1.getFilter());
												newComplexEvent.setQual(e1.getQual());
//												tempInfo = e1.getInfo();
//												//System.out.println(tempInfo+"\n");
//												if (tempInfo.contains("SVTYPE")){
//													tmpOld = tempInfo.substring(tempInfo.indexOf("SVTYPE=")+7, tempInfo.indexOf("SVTYPE=")+10);
//													tmpNew = newComplexEvent.getAlt();
//													tempInfo=tempInfo.replace(tmpOld, tmpNew);
//													tmpOld = tempInfo.substring(tempInfo.indexOf("CHR2=")+5, tempInfo.indexOf(";",tempInfo.indexOf("CHR2=")));
//													tmpNew = tranend.getChr();
//													tempInfo=tempInfo.replace(tmpOld, tmpNew);
//													tmpOld = tempInfo.substring(tempInfo.indexOf(";END=")+5, tempInfo.indexOf(";CT",tempInfo.indexOf(";END=")));
//													tmpNew = Integer.toString(tranend.getPos());
//													tempInfo.replace(tmpOld, tmpNew);
//													newComplexEvent.setInfo(tempInfo);
//												}
												double readDepth = getReadDepth(samReader, transtart.getChr(), transtart.getPos(), tranend.getPos());
												//tempInfo="SVTYPE="+newComplexEvent.getAlt().substring(1, 4)+"; CHR2="+tranend.getChr()+"; END="+Integer.toString(tranend.getPos());
												tempInfo="SVTYPE="+newComplexEvent.getAlt().substring(1, 4)+"; CHR2="+transtart.getChr()+"; START="+Integer.toString(transtart.getPos())+"; END="+Integer.toString(tranend.getPos())+"; ADP="+readDepth;
												newComplexEvent.setInfo(tempInfo);
												newComplexEvent.setCoord(traninsert);
												//writer.write(newComplexEvent.getC1().getChr()+"\t"+newComplexEvent.getC1().getPos()+"\t"+newComplexEvent.getId()+"\t"+newComplexEvent.getRef()+"\t"+newComplexEvent.getAlt()+"\t"+newComplexEvent.getQual()+"\t"+newComplexEvent.getFilter()+"\t"+newComplexEvent.getInfo()+"\n");
											} else {
												//area under e3 is translocated
												GenomicCoordinate transtart, tranend;
												if(e3.getC1().compareTo(e3.getC2()) < 0) {
													transtart = e3.getC1();
													tranend   = e3.getC2();
												} else {
													transtart = e3.getC2();
													tranend   = e3.getC1();
												}
												GenomicCoordinate traninsert = (e2.getNode(true) == currentNode? e2.getC1() : e2.getC2());
												newComplexEvent = new ComplexEvent(transtart, tranend, EVENT_TYPE.COMPLEX_TRANSLOCATION, (new Event[] {e1, e2, e3}), currentNode, traninsert);
												newComplexEvent.setCoord(transtart);
												//newComplexEvent.setId(e1.getId());
												newComplexEvent.setId(e1.getId()+"+"+e2.getId());
												newComplexEvent.setRef(e1.getRef());
												newComplexEvent.setAlt("<TRA>");
												newComplexEvent.setFilter(e1.getFilter());
												newComplexEvent.setQual(e1.getQual());
//												tempInfo = e1.getInfo();
//												//System.out.println(tempInfo+"\n");
//												if (tempInfo.contains("SVTYPE")){
//													tmpOld = tempInfo.substring(tempInfo.indexOf("SVTYPE=")+7, tempInfo.indexOf("SVTYPE=")+10);
//													tmpNew = newComplexEvent.getAlt();
//													tempInfo=tempInfo.replace(tmpOld, tmpNew);
//													tmpOld = tempInfo.substring(tempInfo.indexOf("CHR2=")+5, tempInfo.indexOf(";",tempInfo.indexOf("CHR2=")));
//													tmpNew = tranend.getChr();
//													tempInfo=tempInfo.replace(tmpOld, tmpNew);
//													tmpOld = tempInfo.substring(tempInfo.indexOf(";END=")+5, tempInfo.indexOf(";CT",tempInfo.indexOf(";END=")));
//													tmpNew = Integer.toString(tranend.getPos());
//													tempInfo.replace(tmpOld, tmpNew);
//													newComplexEvent.setInfo(tempInfo);
//												}
												double readDepth = getReadDepth(samReader, transtart.getChr(), transtart.getPos(), tranend.getPos());
												//tempInfo="SVTYPE="+newComplexEvent.getAlt().substring(1, 4)+"; CHR2="+tranend.getChr()+"; END="+Integer.toString(tranend.getPos());
												tempInfo="SVTYPE="+newComplexEvent.getAlt().substring(1, 4)+"; CHR2="+transtart.getChr()+"; START="+Integer.toString(transtart.getPos())+"; END="+Integer.toString(tranend.getPos())+"; ADP="+readDepth;
												newComplexEvent.setInfo(tempInfo);
												newComplexEvent.setCoord(traninsert);
												//writer.write(newComplexEvent.getC1().getChr()+"\t"+newComplexEvent.getC1().getPos()+"\t"+newComplexEvent.getId()+"\t"+newComplexEvent.getRef()+"\t"+newComplexEvent.getAlt()+"\t"+newComplexEvent.getQual()+"\t"+newComplexEvent.getFilter()+"\t"+newComplexEvent.getInfo()+"\n");
											}
										}
										else {
											//System.out.println("Duplication between "+e1+ " and "+ e2);
											GenomicCoordinate dupstart, dupend, insert;
											if(other1.compareTo(other2) < 0){
												//duplicated bit is downstream of currentNode
												dupstart = (e1.getNode(true) == currentNode? e1.getC2() : e1.getC1());
												dupend   = (e2.getNode(true) == currentNode? e2.getC2() : e2.getC1());
												insert   = (e1.getNode(true) == currentNode? e1.getC1() : e1.getC2());
												newComplexEvent = new ComplexEvent(dupstart, dupend, EVENT_TYPE.COMPLEX_DUPLICATION, (new Event[] {e1, e2}), currentNode, insert);
												newComplexEvent.setCoord(dupstart);
												//newComplexEvent.setId(e1.getId());
												newComplexEvent.setId(e1.getId()+"+"+e2.getId());
												newComplexEvent.setRef(e1.getRef());
												newComplexEvent.setAlt("<DUP>");
												newComplexEvent.setFilter(e1.getFilter());
												newComplexEvent.setQual(e1.getQual());
//												tempInfo = e1.getInfo();
//												//System.out.println(tempInfo+"\n");
//												if (tempInfo.contains("SVTYPE")){
//													tmpOld = tempInfo.substring(tempInfo.indexOf("SVTYPE=")+7, tempInfo.indexOf("SVTYPE=")+10);
//													tmpNew = newComplexEvent.getAlt();
//													tempInfo=tempInfo.replace(tmpOld, tmpNew);
//													tmpOld = tempInfo.substring(tempInfo.indexOf("CHR2=")+5, tempInfo.indexOf(";",tempInfo.indexOf("CHR2=")));
//													tmpNew = dupend.getChr();
//													tempInfo=tempInfo.replace(tmpOld, tmpNew);
//													tmpOld = tempInfo.substring(tempInfo.indexOf(";END=")+5, tempInfo.indexOf(";CT",tempInfo.indexOf(";END=")));
//													tmpNew = Integer.toString(dupend.getPos());
//													tempInfo.replace(tmpOld, tmpNew);
//													newComplexEvent.setInfo(tempInfo);
//												}
												double readDepth = getReadDepth(samReader, dupstart.getChr(), dupstart.getPos(), dupend.getPos());
												//tempInfo="SVTYPE="+newComplexEvent.getAlt().substring(1, 4)+"; CHR2="+dupend.getChr()+"; END="+Integer.toString(dupend.getPos());
												tempInfo="SVTYPE="+newComplexEvent.getAlt().substring(1, 4)+"; CHR2="+dupstart.getChr()+"; START="+Integer.toString(dupstart.getPos())+"; END="+Integer.toString(dupend.getPos())+"; ADP="+readDepth;
												newComplexEvent.setInfo(tempInfo);
												newComplexEvent.setCoord(insert);
												//writer.write(newComplexEvent.getC1().getChr()+"\t"+newComplexEvent.getC1().getPos()+"\t"+newComplexEvent.getId()+"\t"+newComplexEvent.getRef()+"\t"+newComplexEvent.getAlt()+"\t"+newComplexEvent.getQual()+"\t"+newComplexEvent.getFilter()+"\t"+newComplexEvent.getInfo()+"\n");
											} else {
												//duplication upstream of currentNode
												dupstart = (e2.getNode(true) == currentNode? e2.getC2() : e2.getC1());
												dupend   = (e1.getNode(true) == currentNode? e1.getC2() : e1.getC1());
												insert   = (e2.getNode(true) == currentNode? e2.getC1() : e2.getC2());
												newComplexEvent = new ComplexEvent(dupstart, dupend, EVENT_TYPE.COMPLEX_DUPLICATION, (new Event[] {e1, e2}), currentNode, insert);
												newComplexEvent.setCoord(dupstart);
												//newComplexEvent.setId(e1.getId());
												newComplexEvent.setId(e1.getId()+"+"+e2.getId());
												newComplexEvent.setRef(e1.getRef());
												newComplexEvent.setAlt("<DUP>");
												newComplexEvent.setFilter(e1.getFilter());
												newComplexEvent.setQual(e1.getQual());
//												tempInfo = e1.getInfo();
//												//System.out.println(tempInfo+"\n");
//												if (tempInfo.contains("SVTYPE")){
//													tmpOld = tempInfo.substring(tempInfo.indexOf("SVTYPE=")+7, tempInfo.indexOf("SVTYPE=")+10);
//													tmpNew = newComplexEvent.getAlt();
//													tempInfo=tempInfo.replace(tmpOld, tmpNew);
//													tmpOld = tempInfo.substring(tempInfo.indexOf("CHR2=")+5, tempInfo.indexOf(";",tempInfo.indexOf("CHR2=")));
//													tmpNew = dupend.getChr();
//													tempInfo=tempInfo.replace(tmpOld, tmpNew);
//													tmpOld = tempInfo.substring(tempInfo.indexOf(";END=")+5, tempInfo.indexOf(";CT",tempInfo.indexOf(";END=")));
//													tmpNew = Integer.toString(dupend.getPos());
//													tempInfo.replace(tmpOld, tmpNew);
//													newComplexEvent.setInfo(tempInfo);
//												}
												double readDepth = getReadDepth(samReader, dupstart.getChr(), dupstart.getPos(), dupend.getPos());
												//tempInfo="SVTYPE="+newComplexEvent.getAlt().substring(1, 4)+"; CHR2="+dupend.getChr()+"; END="+Integer.toString(dupend.getPos());
												tempInfo="SVTYPE="+newComplexEvent.getAlt().substring(1, 4)+"; CHR2="+dupstart.getChr()+"; START="+Integer.toString(dupstart.getPos())+"; END="+Integer.toString(dupend.getPos())+"; ADP="+readDepth;
												newComplexEvent.setInfo(tempInfo);
												newComplexEvent.setCoord(insert);
												//writer.write(newComplexEvent.getC1().getChr()+"\t"+newComplexEvent.getC1().getPos()+"\t"+newComplexEvent.getId()+"\t"+newComplexEvent.getRef()+"\t"+newComplexEvent.getAlt()+"\t"+newComplexEvent.getQual()+"\t"+newComplexEvent.getFilter()+"\t"+newComplexEvent.getInfo()+"\n");
											}
											
										}
									}
								}
								break;
							}
							//interchromosomal events
							case ITX1: {
								if(e2.getType() == EVENT_TYPE.ITX2) {
									GenomicNode other1 = e1.otherNode(currentNode), other2 = e2.otherNode(currentNode);
									if(! other1.getStart().onSameChromosome(other2.getStart()) )
										break;
									GenomicCoordinate eventStart, eventEnd, eventInsert;
									if(currentNode.compareTo(other1) < 0 && other1.getEnd().compareTo(other2.getStart()) < 0 ) {
										eventStart = (e1.getNode(true)==currentNode? e1.getC2() : e1.getC1()); 
										eventEnd = (e2.getNode(true)==currentNode? e2.getC2() : e2.getC1());
									} else if (	currentNode.compareTo(other1) > 0 && other1.getStart().compareTo(other2.getEnd()) > 0) {
										eventStart = (e2.getNode(true)==currentNode? e2.getC2() : e2.getC1());
										eventEnd = (e1.getNode(true)==currentNode? e1.getC2() : e1.getC1());  
									} else {
										break;
									}
									eventInsert = (e1.getNode(true)==currentNode? e1.getC1() : e1.getC2());
									if(eventStart.compareTo(eventEnd) >= 0){
										System.out.println("Fishes!");
									}
									Event e3 = other1.existsDeletionEventTo(other2);
									if(e3 != null){
										newComplexEvent = new ComplexEvent(eventStart, eventEnd, EVENT_TYPE.COMPLEX_INTERCHROMOSOMAL_TRANSLOCATION, new Event[] {e1, e2, e3}, currentNode, eventInsert);
										newComplexEvent.setCoord(eventStart);
										//newComplexEvent.setId(e1.getId());
										newComplexEvent.setId(e1.getId()+"+"+e2.getId());
										newComplexEvent.setRef(e1.getRef());
										newComplexEvent.setAlt("<CIT>");
										newComplexEvent.setFilter(e1.getFilter());
										newComplexEvent.setQual(e1.getQual());
//										tempInfo = e1.getInfo();
//										//System.out.println(tempInfo+"\n");
//										if (tempInfo.contains("SVTYPE")){
//											tmpOld = tempInfo.substring(tempInfo.indexOf("SVTYPE=")+7, tempInfo.indexOf("SVTYPE=")+10);
//											tmpNew = newComplexEvent.getAlt();
//											tempInfo=tempInfo.replace(tmpOld, tmpNew);
//											tmpOld = tempInfo.substring(tempInfo.indexOf("CHR2=")+5, tempInfo.indexOf(";",tempInfo.indexOf("CHR2=")));
//											tmpNew = eventEnd.getChr();
//											tempInfo=tempInfo.replace(tmpOld, tmpNew);
//											tmpOld = tempInfo.substring(tempInfo.indexOf(";END=")+5, tempInfo.indexOf(";CT",tempInfo.indexOf(";END=")));
//											tmpNew = Integer.toString(eventEnd.getPos());
//											tempInfo.replace(tmpOld, tmpNew);
//											newComplexEvent.setInfo(tempInfo);
//										}
										double readDepth = getReadDepth(samReader, eventStart.getChr(), eventStart.getPos(), eventEnd.getPos());
										tempInfo="SVTYPE="+newComplexEvent.getAlt().substring(1, 4)+"; CHR2="+eventStart.getChr()+"; START="+Integer.toString(eventStart.getPos())+"; END="+Integer.toString(eventEnd.getPos())+"; ADP="+readDepth;
										newComplexEvent.setInfo(tempInfo);
										newComplexEvent.setCoord(eventInsert);
										//writer.write(newComplexEvent.getC1().getChr()+"\t"+newComplexEvent.getC1().getPos()+"\t"+newComplexEvent.getId()+"\t"+newComplexEvent.getRef()+"\t"+newComplexEvent.getAlt()+"\t"+newComplexEvent.getQual()+"\t"+newComplexEvent.getFilter()+"\t"+newComplexEvent.getInfo()+"\n");
									} else {
										newComplexEvent = new ComplexEvent(eventStart, eventEnd, EVENT_TYPE.COMPLEX_INTERCHROMOSOMAL_DUPLICATION, new Event[] {e1, e2}, currentNode, eventInsert);
										newComplexEvent.setCoord(eventStart);
										//newComplexEvent.setId(e1.getId());
										newComplexEvent.setId(e1.getId()+"+"+e2.getId());
										newComplexEvent.setRef(e1.getRef());
										newComplexEvent.setAlt("<CID>");
										newComplexEvent.setFilter(e1.getFilter());
										newComplexEvent.setQual(e1.getQual());
//										tempInfo = e1.getInfo();
//										//System.out.println(tempInfo+"\n");
//										if (tempInfo.contains("SVTYPE")){
//											tmpOld = tempInfo.substring(tempInfo.indexOf("SVTYPE=")+7, tempInfo.indexOf("SVTYPE=")+10);
//											tmpNew = newComplexEvent.getAlt();
//											tempInfo=tempInfo.replace(tmpOld, tmpNew);
//											tmpOld = tempInfo.substring(tempInfo.indexOf("CHR2=")+5, tempInfo.indexOf(";",tempInfo.indexOf("CHR2=")));
//											tmpNew = eventEnd.getChr();
//											tempInfo=tempInfo.replace(tmpOld, tmpNew);
//											tmpOld = tempInfo.substring(tempInfo.indexOf(";END=")+5, tempInfo.indexOf(";CT",tempInfo.indexOf(";END=")));
//											tmpNew = Integer.toString(eventEnd.getPos());
//											tempInfo.replace(tmpOld, tmpNew);
//											newComplexEvent.setInfo(tempInfo);
//										}
										double readDepth = getReadDepth(samReader, eventStart.getChr(), eventStart.getPos(), eventEnd.getPos());
										tempInfo="SVTYPE="+newComplexEvent.getAlt().substring(1, 4)+"; CHR2="+eventStart.getChr()+"; START="+Integer.toString(eventStart.getPos())+"; END="+Integer.toString(eventEnd.getPos())+"; ADP="+readDepth;
										newComplexEvent.setInfo(tempInfo);
										newComplexEvent.setCoord(eventInsert);
										//writer.write(newComplexEvent.getC1().getChr()+"\t"+newComplexEvent.getC1().getPos()+"\t"+newComplexEvent.getId()+"\t"+newComplexEvent.getRef()+"\t"+newComplexEvent.getAlt()+"\t"+newComplexEvent.getQual()+"\t"+newComplexEvent.getFilter()+"\t"+newComplexEvent.getInfo()+"\n");
									}
								}
								break;
							}
							case INVTX1: {
								if(e2.getType() == EVENT_TYPE.INVTX2) {
									GenomicNode other1 = e1.otherNode(currentNode), other2 = e2.otherNode(currentNode);
									if(other1.getStart().onSameChromosome(other2.getStart()) && other1.getEnd().compareTo(other2.getStart()) > 0){
										GenomicCoordinate eventStart = (e2.getNode(true)==currentNode? e2.getC2() : e2.getC1()),
										 	eventEnd = (e1.getNode(true)==currentNode? e1.getC2() : e1.getC1()),
											eventInsert = (e1.getNode(true)==currentNode? e1.getC1() : e1.getC2());
									if(eventStart.compareTo(eventEnd) >= 0){
										System.out.println("Fishes!");
									}
									Event e3 = other1.existsDeletionEventTo(other2);
									if(e3 != null){
										newComplexEvent = new ComplexEvent(eventStart, eventEnd, EVENT_TYPE.COMPLEX_INTERCHROMOSOMAL_INVERTED_TRANSLOCATION, new Event[] {e1, e2, e3}, currentNode, eventInsert);
										newComplexEvent.setCoord(eventStart);
										//newComplexEvent.setId(e1.getId());
										newComplexEvent.setId(e1.getId()+"+"+e2.getId()+"+"+e3.getId());
										newComplexEvent.setRef(e1.getRef());
										newComplexEvent.setAlt("<IVT>");
										newComplexEvent.setFilter(e1.getFilter());
										newComplexEvent.setQual(e1.getQual());
//										tempInfo = e1.getInfo();
//										//System.out.println(tempInfo+"\n");
//										if (tempInfo.contains("SVTYPE")){
//											tmpOld = tempInfo.substring(tempInfo.indexOf("SVTYPE=")+7, tempInfo.indexOf("SVTYPE=")+10);
//											tmpNew = newComplexEvent.getAlt();
//											tempInfo=tempInfo.replace(tmpOld, tmpNew);
//											tmpOld = tempInfo.substring(tempInfo.indexOf("CHR2=")+5, tempInfo.indexOf(";",tempInfo.indexOf("CHR2=")));
//											tmpNew = eventEnd.getChr();
//											tempInfo=tempInfo.replace(tmpOld, tmpNew);
//											tmpOld = tempInfo.substring(tempInfo.indexOf(";END=")+5, tempInfo.indexOf(";CT",tempInfo.indexOf(";END=")));
//											tmpNew = Integer.toString(eventEnd.getPos());
//											tempInfo.replace(tmpOld, tmpNew);
//											newComplexEvent.setInfo(tempInfo);
//										}
										double readDepth = getReadDepth(samReader, eventStart.getChr(), eventStart.getPos(), eventEnd.getPos());
										tempInfo="SVTYPE="+newComplexEvent.getAlt().substring(1, 4)+"; CHR2="+eventStart.getChr()+"; START="+Integer.toString(eventStart.getPos())+"; END="+Integer.toString(eventEnd.getPos())+"; ADP="+readDepth;
										newComplexEvent.setInfo(tempInfo);
										newComplexEvent.setCoord(eventInsert);
										//writer.write(newComplexEvent.getC1().getChr()+"\t"+newComplexEvent.getC1().getPos()+"\t"+newComplexEvent.getId()+"\t"+newComplexEvent.getRef()+"\t"+newComplexEvent.getAlt()+"\t"+newComplexEvent.getQual()+"\t"+newComplexEvent.getFilter()+"\t"+newComplexEvent.getInfo()+"\n");
									} else {
										newComplexEvent = new ComplexEvent(eventStart, eventEnd, EVENT_TYPE.COMPLEX_INTERCHROMOSOMAL_INVERTED_DUPLICATION, new Event[] {e1, e2}, currentNode, eventInsert);
										newComplexEvent.setCoord(eventStart);//
										//newComplexEvent.setId(e1.getId());
										newComplexEvent.setId(e1.getId()+"+"+e2.getId());
										newComplexEvent.setRef(e1.getRef());
										newComplexEvent.setAlt("<IVD>");
										newComplexEvent.setFilter(e1.getFilter());
										newComplexEvent.setQual(e1.getQual());
										tempInfo = e1.getInfo();
										//System.out.println(tempInfo+"\n");
//										if (tempInfo.contains("SVTYPE")){
//											tmpOld = tempInfo.substring(tempInfo.indexOf("SVTYPE=")+7, tempInfo.indexOf("SVTYPE=")+10);
//											tmpNew = newComplexEvent.getAlt();
//											tempInfo=tempInfo.replace(tmpOld, tmpNew);
//											tmpOld = tempInfo.substring(tempInfo.indexOf("CHR2=")+5, tempInfo.indexOf(";",tempInfo.indexOf("CHR2=")));
//											tmpNew = eventEnd.getChr();
//											tempInfo=tempInfo.replace(tmpOld, tmpNew);
//											tmpOld = tempInfo.substring(tempInfo.indexOf(";END=")+5, tempInfo.indexOf(";CT",tempInfo.indexOf(";END=")));
//											tmpNew = Integer.toString(eventEnd.getPos());
//											tempInfo.replace(tmpOld, tmpNew);
//											newComplexEvent.setInfo(tempInfo);
//										}
										double readDepth = getReadDepth(samReader, eventStart.getChr(), eventStart.getPos(), eventEnd.getPos());
										tempInfo="SVTYPE="+newComplexEvent.getAlt().substring(1, 4)+"; CHR2="+eventStart.getChr()+"; START="+Integer.toString(eventStart.getPos())+"; END="+Integer.toString(eventEnd.getPos())+"; ADP="+readDepth;
										newComplexEvent.setInfo(tempInfo);
										newComplexEvent.setCoord(eventInsert);
										//writer.write(newComplexEvent.getC1().getChr()+"\t"+newComplexEvent.getC1().getPos()+"\t"+newComplexEvent.getId()+"\t"+newComplexEvent.getRef()+"\t"+newComplexEvent.getAlt()+"\t"+newComplexEvent.getQual()+"\t"+newComplexEvent.getFilter()+"\t"+newComplexEvent.getInfo()+"\n");
									}
								}
							}
							break;
						}
							
							default: //don't even attempt other types
						}
						//check if a new complex event has been generated
						if(newComplexEvent != null){
							//-> add events to cleanup and break current loop
							newComplexEvents.add(newComplexEvent);
							for(Event e: newComplexEvent.getEventsInvolvedInComplexEvent()){
								removeEvents.add(e);
							}
							newComplexEvent = null;
							break; //break the for-j loop, as this guy is already paired
						}
					}
				}
				//all event pairings have been investigated 
				//-> clean up some stuff by removing events and adding the new complex ones.
				for(Event e: removeEvents){
					e.getNode(true).getEvents().remove(e);
					e.getNode(false).getEvents().remove(e);
				}
				for(Event e: newComplexEvents){
					e.getNode(true).getEvents().add(e);
				}
			}
		}
		
		//while we're at it: let's run through the nodes again!
		//this time for output
		int totalEvents = 0;
		for(Entry<String, TreeSet<GenomicNode>> tableEntry: genomicNodes.entrySet()) {
			//System.out.println("Working on Entry: "+tableEntry.toString());
			for(GenomicNode currentNode: tableEntry.getValue()){
				if(currentNode.getEvents().size() > 1){
//					System.out.println("Node might be shifty: "+currentNode.getEvents().size()+" members!");
//					System.out.println(currentNode.getEvents().get(0)+"  "+currentNode.getEvents().get(1));
				}
				totalEvents += currentNode.getEvents().size();
				HashSet<Event> skipEvents = new HashSet<Event>(), deleteEvents = new HashSet<Event>(), newEvents = new HashSet<Event>();
				for(Event e: currentNode.getEvents()){
					if(skipEvents.contains(e))
						continue;
					//if(currentNode.getEvents().size() < 2 && e instanceof ComplexEvent ){//&& e.otherNode(currentNode) != currentNode){// (e.getType() == EVENT_TYPE.COMPLEX_INTERCHROMOSOMAL_DUPLICATION || e.getType()==EVENT_TYPE.COMPLEX_INTERCHROMOSOMAL_TRANSLOCATION)){
						e.processAdditionalInformation(); //TODO: this is a bit of a sly hack to classify insertions in Socrates... not sure how to do it more transparently. 	
						switch(e.getType()) {
						case INV1: 
						case INV2:
							skipEvents.add(e);
							//deleteEvents.add(e);
							if(classifySimpleInversion) {
								ComplexEvent e2 = new ComplexEvent(e.getC1(), e.getC2(), EVENT_TYPE.COMPLEX_INVERSION, new Event[] {e}, currentNode);
								e = e2;
								newEvents.add(e2);
								e2.setId(e2.getId());
								e2.setRef(e2.getRef());
								e2.setAlt("<INV>");
								e2.setFilter(e2.getFilter());
								e2.setQual(e2.getQual());
								e2.setInfo(e2.getInfo());
							}
							else {
								e.setAlt("<INV>");
								e.setFailFilter();
							}
							break;
						case DEL:
							//check for deletion
							//double readDepth = meanReadDepth(reader, e.getC1().getPos()+1, e.getC2().getPos()-1);
							double readDepth = getReadDepth(samReader, e.getC1().getChr(), e.getC1().getPos()+1, e.getC2().getPos()-1);
							skipEvents.add(e);
							if(readDepth < 0 || readDepth > mean-interval){
								//deleteEvents.add(e);
								e.setFailFilter(); 
							} else {
								//System.out.print("read depth for event: "+readDepth+"\t");
							}
							e.setAlt("<DEL>");
							e.setInfo(e.getInfo()+"; ADP="+readDepth );
							break;
						case TAN:
							//double readDepth = meanReadDepth(reader, e.getC1().getPos()+1, e.getC2().getPos()-1);
							readDepth = getReadDepth(samReader, e.getC1().getChr(), e.getC1().getPos(), e.getC2().getPos());
							skipEvents.add(e);
//							//double flank = (meanReadDepth(reader, e.getC1().getPos()-200, e.getC1().getPos()) + meanReadDepth(reader, e.getC2().getPos(), e.getC2().getPos()+200))/2;
							if( readDepth < 0 || readDepth < mean+interval){
								//System.out.println("\t\t\t\t\t\tNot proper duplication!!");
								//deleteEvents.add(e);
								e.setFailFilter();
							} else {
								//System.out.print("read depth for event: "+readDepth+"\t");
							}
							e.setAlt("<DUP>");
							e.setInfo(e.getInfo()+"; ADP="+readDepth );
							break;
						case COMPLEX_DUPLICATION:
						case COMPLEX_INVERTED_DUPLICATION:
						case COMPLEX_INTERCHROMOSOMAL_DUPLICATION:
						case COMPLEX_INTERCHROMOSOMAL_INVERTED_DUPLICATION:
							if(e.getC2().getPos() - e.getC1().getPos() < 50){
								//too small for RD check
								break;
							}
//							readDepth = getReadDepth(samReader, e.getC1().getChr(), e.getC1().getPos(), e.getC2().getPos());
//							if(readDepth < mean+interval){
//								deleteEvents.add(e);
//								skipEvents.add(e);
//								continue;
//							} else {
//								System.out.print("read depth for event: "+readDepth+"\t");
//							}
							break;
						case ITX1:
						case ITX2:
						case INVTX1:
						case INVTX2:
							e.setFailFilter();
							e.setAlt(Event.altVCF(e.getType()));
							skipEvents.add(e);
							break;
						}
						
						//System.out.println(e);
					//}
					if(e.otherNode(currentNode) == currentNode){
						skipEvents.add(e);
						//System.out.println("Self reference: "+e);
					} else {
						e.otherNode(currentNode).getEvents().remove(e);
					}
				}
				currentNode.getEvents().addAll(newEvents);
				for(Event e: deleteEvents){
					e.getNode(true).getEvents().remove(e);
					e.getNode(false).getEvents().remove(e);
				}
			}
		}
		//System.out.println("Total events: "+totalEvents);
		
		//compareToGoldStandard(goldStandard, genomicNodes, 150, true);
		if(goldStandard != null)
			compareToGoldStandard(goldStandard, genomicNodes, 150, false);
	
		//graphVisualisation("data/simul_ecoli_graph.gv", genomicNodes);
		
		count=0;
		//reportEventComposition(genomicNodes);
		/*VCF Output*/
		HashSet<Event> skipEvents = new HashSet<Event>();
		for(Entry<String, TreeSet<GenomicNode>> tableEntry: genomicNodes.entrySet()) {
			for(GenomicNode currentNode: tableEntry.getValue()){
				for(Event e: currentNode.getEvents()){
					if(skipEvents.contains(e)){
						continue;
					} else {
						skipEvents.add(e);
					}
					try{
						writer.write(e.toVcf()+"\n");
					} catch (NullPointerException npe) {
						System.out.println("VCF fail");
					}
				}
			}
		}
		
		writer.close();		
		samReader.close();	
		//End Time
		long endTime = System.nanoTime();
		System.out.println("Took "+(endTime - startTime) + " ns"); 
	}

	

}
