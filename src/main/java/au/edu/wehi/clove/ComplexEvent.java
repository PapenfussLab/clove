package au.edu.wehi.clove;


public class ComplexEvent extends Event{

	private Event[] eventsInvolvedInComplexEvent;
	private GenomicCoordinate insertionPoint;
	public ComplexEvent(GenomicCoordinate c1, GenomicCoordinate c2,
			EVENT_TYPE type, Event[] involvedEvents, GenomicNode hostingNode) {
		super(c1, c2, type);
		this.eventsInvolvedInComplexEvent = involvedEvents;
		super.setNode(hostingNode, true);
		super.setNode(hostingNode, false);
		
	}
	public ComplexEvent(GenomicCoordinate c1, GenomicCoordinate c2,
			EVENT_TYPE type, Event[] involvedEvents, GenomicNode hostingNode, GenomicCoordinate insertionPoint) {
		this(c1,c2,type,involvedEvents,hostingNode);
		this.insertionPoint = insertionPoint;
	}

	public Event[] getEventsInvolvedInComplexEvent(){
		return this.eventsInvolvedInComplexEvent;
	}
	
	public GenomicCoordinate getInsertionPoint() {
		return this.insertionPoint;
	}
	
	@Override
	public String toString(){
		if(this.getType() == EVENT_TYPE.COMPLEX_TRANSLOCATION || this.getType() == EVENT_TYPE.COMPLEX_DUPLICATION 
				|| this.getType() == EVENT_TYPE.COMPLEX_INTERCHROMOSOMAL_DUPLICATION || this.getType() == EVENT_TYPE.COMPLEX_INTERCHROMOSOMAL_TRANSLOCATION
				|| this.getType() == EVENT_TYPE.COMPLEX_INVERTED_DUPLICATION || this.getType() == EVENT_TYPE.COMPLEX_INVERTED_TRANSLOCATION
				|| this.getType() == EVENT_TYPE.COMPLEX_INTERCHROMOSOMAL_INVERTED_DUPLICATION || this.getType() == EVENT_TYPE.COMPLEX_INTERCHROMOSOMAL_INVERTED_TRANSLOCATION){
			return this.getInsertionPoint()+" "+super.toString();
		} else {
			return super.toString();
		}
	}
	
}
