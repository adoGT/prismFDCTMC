//==============================================================================
//	
//	Copyright (c) 2002-
//	Authors:
//	* Dave Parker <david.parker@comlab.ox.ac.uk> (University of Oxford)
//	
//------------------------------------------------------------------------------
//	
//	This file is part of PRISM.
//	
//	PRISM is free software; you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation; either version 2 of the License, or
//	(at your option) any later version.
//	
//	PRISM is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//	
//	You should have received a copy of the GNU General Public License
//	along with PRISM; if not, write to the Free Software Foundation,
//	Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//	
//==============================================================================

package explicit;

import java.lang.UnsupportedOperationException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import prism.PrismLog;
import prism.PrismException;
import prism.ModelType;

/**
 * Simple explicit-state representation of a FDCTMC.
 */
public class FDCTMCSimple extends CTMCSimple implements FDCTMC {
	protected List<FDEvent> fdEvents;
	protected Map<Integer,Map<Integer,Map<Integer,String>>> transToSynchLabel;

	// Constructors

	/**
	 * Constructor: empty FDCTMC.
	 */
	public FDCTMCSimple() throws PrismException {
		super();
		initialise();
	}

	/**
	 * Constructor: new FDCTMC with fixed number of states.
	 */
	public FDCTMCSimple(int numStates) throws PrismException {
		super(numStates);
		initialise();
	}

	/**
	 * Copy constructor.
	 */
	public FDCTMCSimple(FDCTMCSimple fdctmc) {
		super(fdctmc);
		this.fdEvents = new ArrayList<FDEvent>(fdctmc.getNumFDEvents());
		for (int i = 0; i < fdctmc.getNumFDEvents(); ++i) {
			this.fdEvents.add(new FDEvent(fdctmc.getFDEvent(i)));
		}
		this.transToSynchLabel = fdctmc.transToSynchLabel; //TODO perform deep copy if needed
	}

	/**
	 * TODO: implement properly Construct a FDCTMC from an existing one and a
	 * state in permutation, i.e. in which state index i becomes index
	 * permut[i]. Note: have to build new Distributions from scratch anyway to
	 * do this, so may as well provide this functionality as a constructor.
	 */
	public FDCTMCSimple(FDCTMCSimple fdctmc, int permut[]) 
	{
		super(fdctmc, permut);
		this.fdEvents = new ArrayList<FDEvent>(fdctmc.getNumFDEvents());
		for (int i = 0; i < fdctmc.getNumFDEvents(); ++i) {
			this.fdEvents.add(new FDEvent(fdctmc.getFDEvent(i), permut));
		}
		
		transToSynchLabel = new HashMap<Integer,Map<Integer,Map<Integer,String>>> (fdctmc.getNumFDEvents()+1);
		if(fdctmc.transToSynchLabel == null) return; 
		
		for(int i : fdctmc.transToSynchLabel.keySet()) {
			Map<Integer,Map<Integer,String>> sources = new HashMap<Integer,Map<Integer,String>>();
			for(int src: fdctmc.transToSynchLabel.get(i).keySet()) {
				Map<Integer,String> destinations = new HashMap<Integer,String>(fdctmc.transToSynchLabel.get(i).get(src).size());
				for(Entry<Integer,String> dest: fdctmc.transToSynchLabel.get(i).get(src).entrySet())
					destinations.put(permut[dest.getKey()],dest.getValue());
				sources.put(permut[src], destinations);
			}
			transToSynchLabel.put(i, sources);
		}
	}

	private void initialise() throws PrismException {
		this.fdEvents = new ArrayList<FDEvent>();
		transToSynchLabel = new HashMap<Integer,Map<Integer,Map<Integer,String>>> (getNumFDEvents()+1);
		transToSynchLabel.put(-1, new HashMap<Integer, Map<Integer,String>>());
		for (int i= 0; i <getNumFDEvents();++i) {
			transToSynchLabel.put(i, new HashMap<Integer, Map<Integer,String>>());
		}
	}

	@Override
	public int addState() {
		return super.addState();
	}

	@Override
	public void addStates(int numToAdd) {
		super.addStates(numToAdd);

		for (FDEvent fdEvent : fdEvents)
			fdEvent.addStates(numToAdd);
	}

	/**
	 * Add to the probability for a transition.
	 */
	public void addToProbability(int i, int j, double prob, double delay,
			int module, int fdIndex, String label) {
		if (fdIndex < 0)
			addToProbability(i, j, prob);
		else {
			FDEvent fdEvent = getFDEvent(module, fdIndex);
			if (fdEvent == null) {
				fdEvent = new FDEvent(label, numStates, delay, module, fdIndex);
				addFDEvent(fdEvent);
			}
			fdEvent.addToProbability(i, j, prob);
		}
	}

	// Accessors (for ModelSimple, overrides CTMCSimple)

	@Override
	public ModelType getModelType() {
		return ModelType.FDCTMC;
	}

	@Override
	public List<FDEvent> getAllFDEvents() {
		return fdEvents;
	}
        
	public List<Integer> getAllFDEventsIndexes() {
		List<Integer> result = new ArrayList<>(fdEvents.size());
		for(int i = 0; i < fdEvents.size(); ++i) {
			result.add(i);
		}
		return result;
	}

	public int getNumFDEvents() {
		if (fdEvents == null)
			return 0;
		return fdEvents.size();
	}

	public void addFDEvent(FDEvent fdEvent) {
		fdEvents.add(fdEvent);
	}

	public FDEvent getFDEvent(int i) {
		return fdEvents.get(i);
	}

	public FDEvent getFDEvent(int module, int index) {
		for (FDEvent fdEvent : fdEvents) {
			if (fdEvent.getModule() == module && fdEvent.getIndex() == index)
				return fdEvent;
		}
		return null;
	}

	@Override
	public boolean isFDEventActive(FDEvent fdEvent, int state) {
		throw new UnsupportedOperationException("Not implemented");
	}

	@Override
	public void exportToPrismExplicitTra(PrismLog out) throws PrismException {
		super.exportToPrismExplicitTra(out);


		out.println(this);
		out.println("FDEvents: " + getNumFDEvents());
	}



	@Override
	public String toString() {
		return "FDCTMCSimple [fdEvents=" + fdEvents + ", transToSynchLabel="
				+ transToSynchLabel + ", trans=" + trans + ", initialStates="
				+ initialStates + "]";
	}

	@Override
	public void addSynchLabel(int fixedD, int src, int dest, String label) throws PrismException {
		if(transToSynchLabel.get(fixedD) == null)
			transToSynchLabel.put(fixedD, new HashMap<Integer, Map<Integer,String>>());
		
		Map<Integer, String> map;
		if(transToSynchLabel.get(fixedD).get(src) == null) {
			map = new HashMap<Integer, String>(1);
			map.put(dest, label);
			transToSynchLabel.get(fixedD).put(src, map);
			return;
		}
		
		map = transToSynchLabel.get(fixedD).get(src);
		
		if(map.get(dest) == null){
			map.put(dest,label);
			return;
		}
		
		if(map.get(dest) != label)
			throw new PrismException("Multiple synchronization labels from state " + src + " to state " + dest + "!");

	}

	@Override
	public void clearSynchLabels() {
		transToSynchLabel.clear();
	}

	@Override
	public Map<Integer, String> getSychLabelsForState(int fixedD, int state) {
		if(transToSynchLabel.get(fixedD) == null)
			return null;
		return transToSynchLabel.get(fixedD).get(state);
	}

}
