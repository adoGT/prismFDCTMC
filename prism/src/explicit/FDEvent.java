package explicit;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.*;
import java.util.Map.Entry;
import prism.PrismException;

public class FDEvent extends DTMCSimple {
	private String label;
	private double delay;
	private double weight;
	private BitSet active;
	private int module;
	private int index;

	// Constructors

	/**
	 * Constructor: empty FDEvent.
	 */
	public FDEvent() {
		super();
		label = "";
		delay = 0;
		module = 0;
		index = 0;
		weight = 1;
		clearActive();
	}

	/**
	 * Constructor: new FDEvent with fixed number of states.
	 */
	public FDEvent(int numStates) {
		super(numStates);
		label = "";
		delay = 0;
		module = 0;
		index = 0;
		weight = 1;
		clearActive();
	}

	public FDEvent(String label, int numStates, double delay, int module, int index) {
		super(numStates);
		this.label = label;
		this.delay = delay;
		this.module = module;
		this.index = index;
		weight = 1;
		clearActive();
	}

	/**
	 * Copy constructor.
	 */
	public FDEvent(FDEvent fdEvent) {
		super(fdEvent);
		this.label = fdEvent.label;
		this.delay = fdEvent.delay;
		this.weight = fdEvent.weight;
		this.module = fdEvent.module;
		this.index = fdEvent.index;
		clearActive();
		this.active.or(fdEvent.active);
	}

	/**
	 * Copy constructor.
	 */
	public FDEvent(FDEvent fdEvent, int permut[]) {
		super(fdEvent, permut);
		this.label = fdEvent.label;
		this.delay = fdEvent.delay;
		this.weight = fdEvent.weight;
		this.module = fdEvent.module;
		this.index = fdEvent.index;
		clearActive();
		int min = (numStates < permut.length ? numStates : permut.length);
		for (int i = 0; i < min; i++) {
			if (fdEvent.isActive(i))
				active.set(permut[i]);
		}
		// this.active.or(fdEvent.active);
	}

	public int getNumberOfSteps(double interval) throws PrismException {
		System.out.println("Delay: " + delay + " interval: " + interval
				+ " res: "
				+ (new Double(Math.floor(delay / interval))).intValue());
        if(delay < interval)
        	throw new PrismException("Delay(" + delay + " is smaller than discretization step(" + interval + " ).");
		return (int) Math.round(delay / interval);
	}

	/**
	 * Add to the probability for a transition.
	 */
	public void addToProbability(int i, int j, double prob) {
		super.addToProbability(i, j, prob);
		setActive(i);
	}

	private void clearActive() {
		active = new BitSet(numStates);
	}

	public void setActive(int state) {
		active.set(state);
	}

	public void setWeight(double weight) {
		this.weight = weight;
	}

	public void setDelay(double delay) {
		this.delay = delay;
	}

	public boolean isActive(int state) {
		return active.get(state);
	}

	public String getLabel() {
		return label;
	}
	
	public BitSet getActive() {
		return active;
	}

	public double getWeight() {
		return weight;
	}

	public int getModule() {
		return module;
	}

	public int getIndex() {
		return index;
	}

	public double getDelayTime() {
		return delay;
	}

	@Override
	public String toString() {
		return "FDEvent{" + "label=" + label + ", delay=" + delay + ", weight=" + weight
				+ ", active=" + active + ", module=" + module + ", index="
				+ index + ", super=" + super.toString() + '}';
	}

}
