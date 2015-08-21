package au.edu.wehi.clove;


class GenomicCoordinate implements Comparable<GenomicCoordinate>{
	private String chr;
	private int pos;
	public GenomicCoordinate(String chr, int pos){
		this.chr = chr;
		this.pos = pos;
	}
	@Override
	public int compareTo(GenomicCoordinate other) {
		if(this.onSameChromosome(other)){
			if(this.pos < other.pos)
				return -1;
			if(this.pos > other.pos)
				return 1;
			return this.chr.compareTo(other.chr);
		}
		return compareChromosomes(other);
	}
	public String getChr() {
		return chr;
	}
	public int getPos() {
		return pos;
	}
	public boolean onSameChromosome(GenomicCoordinate other){
		if (this.chr.equals(other.chr))
			return true;
		return false;
	}
	public int distanceTo(GenomicCoordinate other){
		if(!onSameChromosome(other))
			return Integer.MAX_VALUE;
		return Math.abs(this.pos - other.pos);
	}
	@Override
	public String toString() {
		return chr+":"+pos;
	}
	
	private int compareChromosomes(GenomicCoordinate other){
		return this.chr.compareTo(other.chr); 
	}
}