package au.edu.wehi.clove;

import java.util.HashSet;
import java.util.StringTokenizer;


enum EVENT_TYPE {INS, INV1, INV2, DEL, TAN, INVTX1, INVTX2, ITX1, ITX2, XXX, COMPLEX_INVERSION, COMPLEX_INVERTED_DUPLICATION, COMPLEX_DUPLICATION, COMPLEX_TRANSLOCATION, COMPLEX_INVERTED_TRANSLOCATION, COMPLEX_INTERCHROMOSOMAL_TRANSLOCATION, COMPLEX_INTERCHROMOSOMAL_DUPLICATION, COMPLEX_INTERCHROMOSOMAL_INVERTED_TRANSLOCATION, COMPLEX_INTERCHROMOSOMAL_INVERTED_DUPLICATION};

public class Event {

	private GenomicCoordinate c1, c2;
	private EVENT_TYPE type;
	private GenomicNode[] myNodes;
	private String additionalInformation;
	private GenomicCoordinate coord;
	private String id;
	private String ref;
	private String alt;
	private String qual;
	private String filter;
	private String info;
	private HashSet<Clove.SV_ALGORITHM> calledBy;
	private int calledTimes;
		
	public Event(GenomicCoordinate c1, GenomicCoordinate c2, EVENT_TYPE type){
		if(c1.compareTo(c2) < 0){
			this.c1 = c1;
			this.c2 = c2;
		} else {
			this.c1 = c2;
			this.c2 = c1;
		}
		this.type = type;
		myNodes = new GenomicNode[2];
		this.info="";
		this.calledBy = new HashSet<Clove.SV_ALGORITHM>();
		this.calledTimes = 0;
	}
	
	public Event(GenomicCoordinate c1, GenomicCoordinate c2, EVENT_TYPE type, String additionalInformation){
		this(c1,c2,type);
		this.additionalInformation = additionalInformation;
	}

	/*Create event with VCF Info*/
	public Event(GenomicCoordinate c1, GenomicCoordinate c2, EVENT_TYPE type, String id, String ref, String alt, String qual, String filter, String info, HashSet<Clove.SV_ALGORITHM> calledBy, int calledTimes){
		if(c1.compareTo(c2) < 0){
			this.c1 = c1;
			this.c2 = c2;
		} else {
			this.c1 = c2;
			this.c2 = c1;
		}
		this.coord = c1;
		this.type = type;
		myNodes = new GenomicNode[2];
		this.id=id;
		this.ref=ref;
		this.alt=alt;
		this.qual=qual;
		this.filter=filter;
		this.info=info;
		this.calledBy = calledBy;
		this.calledTimes = calledTimes;
	}
	
	/*
	 * Static function to handle the particularities of Socrates output, and convert it into a general
	 * purpose Event.
	 */
	public static Event createNewEventFromSocratesOutput(String output){
		String line = output.replace("\t\t", "\tX\t");
		StringTokenizer t = new StringTokenizer(line);
		String chr1 = t.nextToken(":");
		int p1 = Integer.parseInt(t.nextToken(":\t"));
		String o1 = t.nextToken("\t");
		t.nextToken("\t");
		String chr2 = t.nextToken("\t:");
		int p2 = Integer.parseInt(t.nextToken("\t:"));
		String o2 = t.nextToken("\t");
		
		GenomicCoordinate c1 = new GenomicCoordinate(chr1, p1);
		GenomicCoordinate c2 = new GenomicCoordinate(chr2, p2);
		EVENT_TYPE type = classifySocratesBreakpoint(c1, o1, c2, o2);
		
		//look for additional information at the end of the call
		int i = 0;
		while(i<19 && t.hasMoreTokens()){
			i++;
			t.nextToken();
		}
		String additionalComments = (t.hasMoreTokens()? t.nextToken() : "");
		if(additionalComments.startsWith("Inserted sequence")){
			String insert = additionalComments.substring("Inserted sequence: ".length());
			return new Event(c1, c2, type, insert);
		}
		
		return new Event(c1, c2, type);
	}
	/*
	 * Static function to handle the particularities of Socrates output, and convert it into a general
	 * purpose Event.
	 */
	public static Event createNewEventFromSocratesOutputLatest(String output, int count){
		String line = output.replace("\t\t", "\tX\t");
		StringTokenizer t = new StringTokenizer(line);
		String chr1 = t.nextToken(":");
		int p1 = Integer.parseInt(t.nextToken(":\t"));
		String o1 = t.nextToken("\t");
		t.nextToken("\t");
		String chr2 = t.nextToken("\t:");
		int p2 = Integer.parseInt(t.nextToken("\t:"));
		String o2 = t.nextToken("\t");
		
		String id="SOC"+Integer.toString(count);
		String ref=".";
		String qual=".";
		String filter="PASS";
		
		GenomicCoordinate c1 = new GenomicCoordinate(chr1, p1);
		GenomicCoordinate c2 = new GenomicCoordinate(chr2, p2);
		EVENT_TYPE type = classifySocratesBreakpoint(c1, o1, c2, o2);
		
		//look for additional information at the end of the call
		int i = 0;
		while(i<19 && t.hasMoreTokens()){
			i++;
			t.nextToken();
		}
		String additionalComments = (t.hasMoreTokens()? t.nextToken() : "");
		if(additionalComments.startsWith("Inserted sequence")){
			String insert = additionalComments.substring("Inserted sequence: ".length());
			return new Event(c1, c2, type, insert);
		}
		
		String alt=altVCF(type);
		String info="SVTYPE="+alt.substring(1, 4)+"; CHR2="+chr2+"; END="+p2;
				
		return new Event(c1, c2, type, id, ref, alt, qual, filter, info, new HashSet<Clove.SV_ALGORITHM>() {{add(Clove.SV_ALGORITHM.SOCRATES);}}, 1);
	}
	/*
	 * Function to classify a line of Socrates output into a genomic event type.
	 * The distinctions between INV1/2 etc are arbitrary, and have to be consistent across all the inputs.
	 */
	private static EVENT_TYPE classifySocratesBreakpoint(GenomicCoordinate c1, String o1, GenomicCoordinate c2, String o2){
		if(c1.onSameChromosome(c2)){
			if(o1.equals(o2)){
				if(o1.equals("+"))
					return EVENT_TYPE.INV1;
				else
					return EVENT_TYPE.INV2;
			} else if (o1.equals("+") && c1.compareTo(c2) < 0 || o1.equals("-") && c1.compareTo(c2) >=0 ){
				return EVENT_TYPE.DEL;
			} else if (o1.equals("-") && c1.compareTo(c2) < 0 || o1.equals("+") && c1.compareTo(c2) >=0 ){
				return EVENT_TYPE.TAN;
			} else {
				return EVENT_TYPE.XXX;
			}
		} else if(o1.equals(o2)) {
			if(o1.equals("+"))
				return EVENT_TYPE.INVTX1;
			else
				return EVENT_TYPE.INVTX2;
		} else if(o1.equals("+") &&  c1.compareTo(c2) < 0 || o1.equals("-") && c1.compareTo(c2) >= 0){
			return EVENT_TYPE.ITX1;
		} else {
			return EVENT_TYPE.ITX2;
		}
	}
	/*
	 * Static function to handle the particularities of Delly output, and convert it into a general
	 * purpose Event.
	 */
	public static Event createNewEventFromDellyOutput(String output){
		StringTokenizer t = new StringTokenizer(output, "\t:");
		String chr1 = t.nextToken();
		String chr2 = chr1;
		int p1 = Integer.parseInt(t.nextToken());
		int p2 = Integer.parseInt(t.nextToken());
		t.nextToken();
		t.nextToken();
		t.nextToken();
		String typeT;
		String tempT = t.nextToken(); 
		typeT = tempT.substring(1,tempT.indexOf("_"));
		if (typeT.equals("Inversion")){
			typeT = tempT.substring(1,(tempT.indexOf("_")+2));
		}
		
		GenomicCoordinate c1 = new GenomicCoordinate(chr1, p1);
		GenomicCoordinate c2 = new GenomicCoordinate(chr2, p2);
		EVENT_TYPE type = classifyDellyBreakpoint(c1, c2, typeT);
		
		//System.out.println(chr1 +"\t"+ p1 +"\t"+ p2 +"\t" + type +"\t"+ typeT);
		
		return new Event(c1, c2, type);
	}

	/*
	 * Function to classify a line of Delly output into a genomic event type.
	 * The distinctions between INV1/2 etc are arbitrary, and have to be consistent across all the inputs.
	 * c1 and c2 are always the same chromosome
	 */
	private static EVENT_TYPE classifyDellyBreakpoint(GenomicCoordinate c1, GenomicCoordinate c2, String t){
		if(t.equals("Inversion_0")){
			return EVENT_TYPE.INV1;
		} else if (t.equals("Inversion_1")){
			return EVENT_TYPE.INV2;
		} else if (t.equals("Deletion")){
			return EVENT_TYPE.DEL;
		} else if (t.equals("Duplication")){
			return EVENT_TYPE.TAN;
		} else {
			return EVENT_TYPE.XXX;
		}
	}
	
	
	/*
	 * Static function to handle the particularities of Delly output, and convert it into a general
	 * purpose Event.
	 */
	public static Event createNewEventFromDellyOutputLatest(String output){
		String[] bits = output.split("\t");
		String chr1 = bits[0];
		int p1 = Integer.parseInt(bits[1]);
		String[] moreBits = bits[7].split(";");
		String chr2 = moreBits[5].replace("CHR2=", "");
		int p2 = Integer.parseInt(moreBits[6].replace("END=", ""));
		String o = moreBits[7].replace("CT=", "");
		String o1 = (Integer.parseInt(o.split("to")[0]) == 3? "+" : "-");
		String o2 = (Integer.parseInt(o.split("to")[1]) == 3? "+" : "-");
		
		String id=bits[2];
		String ref=bits[3];
		String alt=bits[4];
		String qual=bits[5];
		String filter=bits[6];
		String info=bits[7];
		
		GenomicCoordinate c1 = new GenomicCoordinate(chr1, p1);
		GenomicCoordinate c2 = new GenomicCoordinate(chr2, p2);
		EVENT_TYPE type = classifySocratesBreakpoint(c1, o1, c2, o2);
		
		//System.out.println(chr1 +"\t"+ p1 +"\t"+ p2 +"\t" + type +"\t"+ typeT);
		
		return new Event(c1, c2, type, id, ref, alt, qual, filter, info, new HashSet<Clove.SV_ALGORITHM>() {{add(Clove.SV_ALGORITHM.DELLY);}}, 1);
		//return new Event(c1, c2, type);
		
	}
	/*
	 * Function to classify a line of Delly output into a genomic event type.
	 * The distinctions between INV1/2 etc are arbitrary, and have to be consistent across all the inputs.
	 * c1 and c2 are always the same chromosome
	 */
	private static EVENT_TYPE classifyDellyBreakpointLatest(GenomicCoordinate c1, GenomicCoordinate c2){
		return null;
	}
	
	public static Event createNewEventFromBEDPE (String output){
		String[] bits = output.split("\t");
		String chr1 = bits[0];
		int p1 = Integer.parseInt(bits[1]);
		String chr2 = bits[3];
		int p2 = Integer.parseInt(bits[4]);
		String o1 = bits[8];
		String o2 = bits[9];
		String qual=bits[7];
		String id=bits[6];
		String ref="";
		String alt="";
		String info= (bits.length>10? bits[10]: "");
		String filter = "PASS";
		GenomicCoordinate c1 = new GenomicCoordinate(chr1, p1);
		GenomicCoordinate c2 = new GenomicCoordinate(chr2, p2);
		EVENT_TYPE type = classifySocratesBreakpoint(c1, o1, c2, o2);
		return new Event(c1, c2, type, id, ref, alt, qual, filter, info, new HashSet<Clove.SV_ALGORITHM>() {{add(Clove.SV_ALGORITHM.BEDPE);}}, 1);
	}
	
	
	public static Event createNewEventFromCrestOutput(String output) {
		StringTokenizer t = new StringTokenizer(output, "\t");
		
		String chr1 = t.nextToken();
		int p1 = Integer.parseInt(t.nextToken());
		String o1 = t.nextToken();
		t.nextToken();
		String chr2 = t.nextToken();
		int p2 = Integer.parseInt(t.nextToken());
		String o2 = t.nextToken();
		
		GenomicCoordinate c1 = new GenomicCoordinate(chr1, p1);
		GenomicCoordinate c2 = new GenomicCoordinate(chr2, p2);
		
		t.nextToken();
		EVENT_TYPE type = classifyCrestBreakpoint(t.nextToken(), chr1, chr2, o1, o2);
		
		if(type == EVENT_TYPE.COMPLEX_INVERSION){
			return new ComplexEvent(c1, c2, type, null, null);
		}
		return new Event(c1, c2, type);
	}
	
	public static Event createNewEventFromCrestOutputLatest(String output, int count) {
		StringTokenizer t = new StringTokenizer(output, "\t");
		
		String chr1 = t.nextToken();
		int p1 = Integer.parseInt(t.nextToken());
		String o1 = t.nextToken();
		t.nextToken();
		String chr2 = t.nextToken();
		int p2 = Integer.parseInt(t.nextToken());
		String o2 = t.nextToken();
		
		String id="CRT"+Integer.toString(count);
		String ref=".";
		String qual=".";
		String filter="PASS";
		
		GenomicCoordinate c1 = new GenomicCoordinate(chr1, p1);
		GenomicCoordinate c2 = new GenomicCoordinate(chr2, p2);
		
		t.nextToken();
		EVENT_TYPE type = classifyCrestBreakpoint(t.nextToken(), chr1, chr2, o1, o2);
		
		if(type == EVENT_TYPE.COMPLEX_INVERSION){
			return new ComplexEvent(c1, c2, type, null, null);
		}
		
		String alt=altVCF(type);
		String info="SVTYPE="+alt.substring(1, 4)+"; CHR2="+chr2+"; END="+p2;
		
		return new Event(c1, c2, type, id, ref, alt, qual, filter, info, new HashSet<Clove.SV_ALGORITHM>() {{add(Clove.SV_ALGORITHM.CREST);}}, 1);
	}
	
	private static EVENT_TYPE classifyCrestBreakpoint(String t, String c1, String c2, String o1, String o2){
		if(t.equals("DEL")){
			return EVENT_TYPE.DEL;
		} else if (t.equals("INS")){
			return EVENT_TYPE.TAN;
		} else if (t.equals("INV")){
			return EVENT_TYPE.COMPLEX_INVERSION;
		} else if(t.equals("ITX")){
			if(o1.equals("+"))
				return EVENT_TYPE.INV1;
			else 
				return EVENT_TYPE.INV2;
		} else if (t.equals("CTX")) {
//			if(o1.equals(o2)) {
//				if(o1.equals("+"))
//					return EVENT_TYPE.INVTX1;
//				else
//					return EVENT_TYPE.INVTX2;
//			} else if(o1.equals("+") &&  c1.compareTo(c2) < 0 || o1.equals("-") && c1.compareTo(c2) >= 0){
//				return EVENT_TYPE.ITX1;
//			} else {
//				return EVENT_TYPE.ITX2;
//			}
			if(o1.equals(o2)) {
				if(o1.equals("+")) {
					if(c1.compareTo(c2) < 0)
						return EVENT_TYPE.ITX1;
					else
						return EVENT_TYPE.ITX2;
				}
				else
					return EVENT_TYPE.XXX;
//			} else if(o1.equals("+") &&  c1.compareTo(c2) < 0 || o1.equals("-") && c1.compareTo(c2) >= 0){
			} else if(c1.compareTo(c2) < 0 && o1.equals("-")){
				return EVENT_TYPE.INVTX2;
			} else if(c1.compareTo(c2) < 0 && o1.equals("+")){
				return EVENT_TYPE.INVTX1;
			}
			else if (c1.compareTo(c2) >= 0 && o1.equals("+")){
				return EVENT_TYPE.INVTX1;
			} else {
				return EVENT_TYPE.INVTX2;
			}
		} else {
			return EVENT_TYPE.XXX;
		}
	}


	public static Event createNewEventFromGustafOutput(String output) {
		StringTokenizer t = new StringTokenizer(output, "\t");
		
		String chr1 = t.nextToken();
		t.nextToken();
		EVENT_TYPE type = classifyGustafBreakpoint(t.nextToken());
		int p1 = Integer.parseInt(t.nextToken());
		int p2 = Integer.parseInt(t.nextToken());
		t.nextToken();
		String o = t.nextToken();
		if(type==EVENT_TYPE.INV1 && o.equals("-"))
			type = EVENT_TYPE.INV2;
		t.nextToken();
		String info = t.nextToken();
		StringTokenizer i = new StringTokenizer(info, "=;");
		i.nextToken();
		i.nextToken();
		String d = i.nextToken();
		String chr2;
		if(d.equals("endChr")){
			chr2 = i.nextToken();
			i.nextToken();
			p2 = Integer.parseInt(i.nextToken());
		} else if (d.equals("size")) {
			chr2 = chr1;
		} else {
			System.err.println("Confusion in the Gustaf camp!");
			chr2=null;
		}
		
		GenomicCoordinate c1 = new GenomicCoordinate(chr1, p1);
		GenomicCoordinate c2 = new GenomicCoordinate(chr2, p2);
		
		return new Event(c1, c2, type);
	}
	private static EVENT_TYPE classifyGustafBreakpoint(String t){
		if(t.equals("deletion")){
			return EVENT_TYPE.DEL;
		} else if (t.equals("duplication")){
			return EVENT_TYPE.TAN;
		} else if (t.equals("inversion")){
			return EVENT_TYPE.INV1;
		} else if(t.equals("ITX")){
			return EVENT_TYPE.INV1;
		} else if (t.equals("insertion")) {
			return EVENT_TYPE.INS;
		} else {
			return EVENT_TYPE.XXX;
		}
	}
	
	public GenomicCoordinate getC1() {
		return c1;
	}

	public GenomicCoordinate getC2() {
		return c2;
	}

	public EVENT_TYPE getType() {
		return type;
	}
	
	public void setNode(GenomicNode n, boolean firstCoordinate){
		if(firstCoordinate)
			myNodes[0] = n;
		else
			myNodes[1] = n;
	}
	
	public GenomicNode getNode(boolean firstCoordinate){
		if(firstCoordinate)
			return myNodes[0];
		else
			return myNodes[1];
	}
	
	
	public static boolean sameNodeSets(Event e1, Event e2){
		if(e1.myNodes[0] == e2.myNodes[0] && e1.myNodes[1] == e2.myNodes[1] 
				|| e1.myNodes[0] == e2.myNodes[1] && e1.myNodes[1] == e2.myNodes[0])
			return true;
		return false;		
	}
	
	@Override
	public String toString() {
		if(c1.onSameChromosome(c2)){
			return c1.getChr()+":"+c1.getPos()+"-"+c2.getPos()+" "+type;
		} else {
			return c1+"<->"+c2+" "+type;
		}
	}
	
	public GenomicNode otherNode(GenomicNode node){
		if(myNodes[0] == node) {
			return myNodes[1];
		}
		if(myNodes[1] == node) {
			return myNodes[0];
		}
		System.err.println("otherNode: query node is not assiciated with Event!");
		return null;
	}
	
	public int size() {
		return c1.distanceTo(c2);
	}
	
	public void processAdditionalInformation(){
		if(this.additionalInformation!= null && this.additionalInformation.matches("[ACGT]+") && myNodes[0] == myNodes[1]){
			this.type = EVENT_TYPE.INS;
		}
	}
	
	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public String getQual() {
		return qual;
	}

	public void setQual(String qual) {
		this.qual = qual;
	}

	public String getAlt() {
		return alt;
	}

	public void setAlt(String alt) {
		this.alt = alt;
	}

	public String getFilter() {
		return filter;
	}

	public void setFilter(String filter) {
		this.filter = filter;
	}

	public String getRef() {
		return ref;
	}

	public void setRef(String ref) {
		this.ref = ref;
	}

	public String getInfo() {
		return info;
	}

	public void setInfo(String info) {
		this.info = info;
	}

	public GenomicCoordinate getCoord() {
		return coord;
	}

	public void setCoord(GenomicCoordinate newCoord) {
		coord = newCoord;
	}

	public static String altVCF(EVENT_TYPE type){
		if(type.equals(EVENT_TYPE.DEL)){
			return "<DEL>";
		} else if(type.equals(EVENT_TYPE.INS)){
			return "<INS>";
		} else if(type.equals(EVENT_TYPE.TAN)){
			return "<TAN>";
		} else if(type.equals(EVENT_TYPE.INV1)){
			return "<INV>";
		} else if(type.equals(EVENT_TYPE.INV2)){
			return "<INV>";
		} else if(type.equals(EVENT_TYPE.INVTX1)){
			return "<INV>";
		} else if(type.equals(EVENT_TYPE.INVTX2)){
			return "<INV>";
		} else if(type.equals(EVENT_TYPE.ITX1)){
			return "<TRA>";
		} else if(type.equals(EVENT_TYPE.ITX2)){
			return "<TRA>";
		} else if(type.equals(EVENT_TYPE.COMPLEX_DUPLICATION)){
			return "<DUP>";
		} else if(type.equals(EVENT_TYPE.COMPLEX_INVERTED_TRANSLOCATION)){
			return "<CVT>";
		} else if(type.equals(EVENT_TYPE.COMPLEX_INVERTED_DUPLICATION)){
			return "<CVD>";
		} else if(type.equals(EVENT_TYPE.COMPLEX_TRANSLOCATION)){
			return "<TRA>";
		} else if(type.equals(EVENT_TYPE.COMPLEX_INVERSION)){
			return "<CIV>";
		} else if(type.equals(EVENT_TYPE.COMPLEX_INTERCHROMOSOMAL_INVERTED_TRANSLOCATION)){
			return "<IVT>";
		} else if(type.equals(EVENT_TYPE.COMPLEX_INTERCHROMOSOMAL_INVERTED_DUPLICATION)){
			return "<IVD>";
		} else {
			return "<XXX>";
		} 
	}
	
	public HashSet<Clove.SV_ALGORITHM> getCalledBy() {
		return calledBy;
	}
	public void addCaller(HashSet<Clove.SV_ALGORITHM> caller){
		this.calledBy.addAll(caller);
	}

	public int getCalledTimes() {
		return calledTimes;
	}
	public void increaseCalls(int inc){
		this.calledTimes += inc;
	}

	public String toVcf() {
		return this.getCoord().getChr()+"\t"+this.getCoord().getPos()+"\t"+this.getId()+"\t"
				+this.getRef()+"\t"+this.getAlt()+"\t"+this.getQual()+"\t"+this.getFilter()
				+"\t"+this.getInfo()+";SUPPORT="+this.calledBy.size()+","+this.calledTimes;
	}
	
	public void setFailFilter(){
		this.filter = "FAIL";
	}
}
