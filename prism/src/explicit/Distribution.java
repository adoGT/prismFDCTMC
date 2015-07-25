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

import java.util.*;
import java.util.Map.Entry;

import prism.PrismUtils;

/**
 * Explicit representation of a probability distribution.
 * Basically, a mapping from (integer-valued) state indices to (non-zero, double-valued) probabilities. 
 */
public class Distribution implements Iterable<Entry<Integer,Double>>
{
	private HashMap<Integer,Double> map;
	
	/**
	 * Create an empty distribution.
	 */
	public Distribution()
	{
		map = new HashMap<>();
	}
        
        public Distribution(double[] values)
	{
		map = new HashMap<>();
                for(int i=0; i<values.length; ++i)
                    set(i, values[i]);
	}
	
	/**
	 * Copy constructor.
	 */
	public Distribution(Distribution distr)
	{
		this();
		Iterator<Entry<Integer,Double>> i = distr.iterator();
		while (i.hasNext()) {
			Map.Entry<Integer,Double> e = i.next();
			add(e.getKey(), e.getValue());
		}
	}
	
	/**
	 * Construct a distribution from an existing one and a state index permutation,
	 * i.e. in which state index i becomes index permut[i].
	 * Note: have to build the new distributions from scratch anyway to do this,
	 * so may as well provide this functionality as a constructor.
	 */
	public Distribution(Distribution distr, int permut[])
	{
		this();
		Iterator<Entry<Integer,Double>> i = distr.iterator();
		while (i.hasNext()) {
			Map.Entry<Integer,Double> e = i.next();
			add(permut[e.getKey()], e.getValue());
		}
	}
	
	/**
	 * Clear all entries of the distribution.
	 */
	public void clear()
	{
		map.clear();
	}
	
	/**
	 * Add 'prob' to the probability for state 'j'.
	 * Return boolean indicating whether or not there was already
	 * non-zero probability for this state (i.e. false denotes new transition).
	 */
	public boolean add(int j, double prob)
	{
		Double d = (Double) map.get(j);
		if (d == null) {
			map.put(j, prob);
			return false;
		} else {
			set(j, d + prob);
			return true;
		}
	}
        
        public Distribution add(Distribution distr) {
            for(Map.Entry<Integer,Double> item : distr.map.entrySet()) {
                this.add(item.getKey(), item.getValue());
            }
            return this;
        }

	/**
	 * Set the probability for state 'j' to 'prob'.
	 */
	public void set(int j, double prob)
	{
		if (prob == 0.0)
			map.remove(j);
		else
			map.put(j, prob);
	}

	/**
	 * Get the probability for state j. 
	 */
	public double get(int j)
	{
		Double d;
		d = (Double) map.get(j);
		return d==null ? 0.0 : d.doubleValue();
	}
	
	/**
	 * Returns true if index j is in the support of the distribution. 
	 */
	public boolean contains(int j)
	{
		return map.get(j) != null;
	}
	
	/**
	 * Returns true if all indices in the support of the distribution are in the set. 
	 */
	public boolean isSubsetOf(BitSet set)
	{
		Iterator<Entry<Integer,Double>> i = iterator();
		while (i.hasNext()) {
			Map.Entry<Integer,Double> e = i.next();
			if (!set.get((Integer) e.getKey()))
				return false;
		}
		return true;
	}
	
	/**
	 * Returns true if all at least one index in the support of the distribution is in the set. 
	 */
	public boolean containsOneOf(BitSet set)
	{
		Iterator<Entry<Integer,Double>> i = iterator();
		while (i.hasNext()) {
			Map.Entry<Integer,Double> e = i.next();
			if (set.get((Integer) e.getKey()))
				return true;
		}
		return false;
	}
	
	/**
	 * Get an iterator over the entries of the map defining the distribution.
	 */
	public Iterator<Entry<Integer,Double>> iterator()
	{
		return map.entrySet().iterator();
	}

	/**
	 * Returns true if the distribution is empty.
	 */
	public boolean isEmpty()
	{
		return map.isEmpty();
	}
	
	/**
	 * Get the size of the support of the distribution.
	 */
	public int size()
	{
		return map.size();
	}

	public Set<Double> getSortedValues()
	{
		Set<Double> toReturn = new TreeSet<Double>();
		for(Double value : map.values()) 
		{
			toReturn.add(value);
		}
		return toReturn;
	}
	
	public String sortedValues() 
	{
		String output = "";
		for(Double value : getSortedValues())
		{
			output = value.toString().concat(" " + output);
		}
		return output;
	}
	
	/**
	 * Get the sum of the probabilities in the distribution.
	 */
	public double sum()
	{
		double d = 0.0;
                for(Map.Entry<Integer,Double> item : map.entrySet()) {
                    d += item.getValue();
                }
		return d;
	}
	
	/**
	 * Get the sum of all the probabilities in the distribution except for state j.
	 */
	public double sumAllBut(int j)
	{
		double d = 0.0;
		Iterator<Entry<Integer,Double>> i = iterator();
		while (i.hasNext()) {
			Map.Entry<Integer,Double> e = i.next();
			if (e.getKey() != j)
				d += e.getValue();
		}
		return d;
	}
	
	/**
	 * Create a new distribution, based on a mapping from the state indices
	 * used in this distribution to a different set  of state indices.
	 */
	public Distribution map(int map[])
	{
		Distribution distrNew = new Distribution();
		Iterator<Entry<Integer,Double>> i = iterator();
		while (i.hasNext()) {
			Map.Entry<Integer,Double> e = i.next();
			distrNew.add(map[e.getKey()], e.getValue());
		}
		return distrNew;
	}
        
        public void addDistr(Distribution distr) {
            for(Map.Entry<Integer, Double> entry : distr) {
                add(entry.getKey(), entry.getValue());
            }
        }
        
        /**
         * Create a new distribution by multiply every probability by weight.
         * @param weight 
         */
        public Distribution multiply(double weight)
        {
            Distribution distrNew = new Distribution();
            for(Map.Entry<Integer, Double> entry : map.entrySet()) {
                distrNew.add(entry.getKey(), weight*entry.getValue());
            }
            return distrNew;
        }
        
        public void multiplyThis(double weight)
        {
            for(Map.Entry<Integer, Double> entry : map.entrySet()) {
                entry.setValue(weight*entry.getValue());
            }
        }
        
        
        public void minMax(Distribution distr, boolean min) {
        	if(min){
                    for( Map.Entry<Integer,Double> prob: this) {
                            if(prob.getValue() > distr.get(prob.getKey()))
                            {
                                    prob.setValue(distr.get(prob.getKey()));
                            }
                    }        		
        	}else{
                    for(Map.Entry<Integer,Double> prob: distr) {
                            if(this.get(prob.getKey()) < prob.getValue()){
                                    this.set(prob.getKey(),prob.getValue());
                            }
                    }
        	}
        }

        public boolean diff(Distribution dist, double epsilon) {
        	Set<Integer> keys = new HashSet<>();
        	keys.addAll(map.keySet());
        	keys.addAll(dist.map.keySet());
        	for( Integer key : keys) {
        		if(Math.abs(map.get(key)-dist.get(key))>epsilon)
        			return true;
        	}
        	return false;
        }
        
	@Override
	public boolean equals(Object o)
	{
                if (this == o) return true;
                if (o == null) return false;
		Double d1, d2;
		Distribution d = (Distribution) o;
                if (d.size() != size())
			return false;
		Iterator<Entry<Integer,Double>> i = iterator();
		while (i.hasNext()) {
			Map.Entry<Integer,Double> e = i.next();
			d1 = e.getValue();
			d2 = d.map.get(e.getKey());
			if (d2 == null || !PrismUtils.doublesAreClose(d1, d2, 1e-12, false))
				return false;
		}
		return true;
	}
	
	@Override
	public int hashCode()
	{
		// Simple hash code
		return map.hashCode() + 15;
	}

	@Override
	public String toString()
	{
		return map.toString();
	}
        
        public String toJSON()
        {
            String result = "[ ";
            
            boolean first = true;
            for(Map.Entry<Integer, Double> entry : this) 
            {
                if(first) 
                    first = false;
                else
                    result += ",";
                
                result += "{ \"state\": " + entry.getKey() + ", \"prob\": " + entry.getValue() + " }";
            }
            result += " ]";
            return result;
        }
}
