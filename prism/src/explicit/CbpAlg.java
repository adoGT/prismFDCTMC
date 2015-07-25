//==============================================================================
//	
//	Copyright (c) 2015-
//	Authors:
//	* Lubos Korenciak <lkorenciak@gmail.com> (Masaryk University)
//	* Adrian Farmadin <a.farmadin@gmail.com> (Masaryk University)
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

import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import prism.PrismException;
import prism.PrismLog;

/**
 * 
 * 
 * @author Adrian
 */
public class CbpAlg {
    
    
	protected PrismLog mainLog;
        
	protected FDCTMCSimple fdctmc;
        
    // stores main steps of computation and final result
	protected FDMatrix<PositionKey> transProb;
    // temporarily stores transProb, for it's update between steps
	protected FDMatrix<PositionKey> transProbCopy;
	protected FDMatrix<PositionKey> transProbResult;
        
    //list of all possible subsets of fix-delays and active states in subset  
	protected List<FDSubset> MTEgroups;
        
	protected double intervalLength;
	protected DTMCSimple uniformisedDTMC;
	protected double uniformizationRate;
	protected Map<FDEvent, Integer> numberOfSteps;

	public CbpAlg(PrismLog mainLog) throws PrismException {

		// TODO we are now computing with DTMCSimple and multiplying
		// distributions with
		// quite complicated objects, but also taking adnvantage of sparse
		// distributions and matrices
		// for some reason CTMCSimple is computing matrix-vector multiplication
		// by
		// by method vmMult from very basic object DTMCExplicit, resulting in
		// multiplication
		// of large arrays containing possibly a lot of zeroes. But maybe this
		// could be faster using
		// some heuristics. It should be tested!
        this.mainLog = mainLog;
		fdctmc = null;
		transProb = new FDMatrix<>();
		transProbCopy = new FDMatrix<>();
		transProbResult = new FDMatrix<>();
		MTEgroups = new ArrayList<>();
		numberOfSteps = new HashMap<FDEvent, Integer>();
		intervalLength = 0.0;
		uniformizationRate = 0.0;
	}

	/**
	 * Copy constructor.
	 */
	public CbpAlg(CbpAlg cbp) {
		this.uniformizationRate = cbp.getUniformizationRate();
		this.uniformisedDTMC = cbp.getUniformisedDTMC();
		this.fdctmc = cbp.fdctmc;
		this.transProb = FDMatrix.copyFDMatrixPositionKey(cbp.transProb);
		this.transProbCopy = FDMatrix.copyFDMatrixPositionKey(cbp.transProbCopy);
		this.transProbResult = FDMatrix.copyFDMatrixPositionKey(cbp.transProbResult);
		this.MTEgroups = new ArrayList<>(cbp.MTEgroups.size());
		for (FDSubset fdSubset : cbp.MTEgroups) {
			this.MTEgroups.add(new FDSubset(fdSubset));
		}
	}
	
	public Distribution runCBP(double intervalLength, int numberOfIntervals,
			double error) throws PrismException // lengthOfInterval,
												// numberOfIntervals, error
    {
		System.out.println("Num of poisson steps: " + getNumberOfPoissonStepsFromError(intervalLength,
                             numberOfIntervals, error, uniformizationRate));
		
        compute(
            numberOfIntervals,
            intervalLength,
            getNumberOfPoissonStepsFromError(intervalLength,
                             numberOfIntervals, error, uniformizationRate));
        return getResult();
    }

	public void setFdctmc(FDCTMCSimple fdctmc) throws PrismException {
		if (fdctmc == null)
			throw new IllegalArgumentException("CBPalg: fdctmc can't be null.");
		this.fdctmc = fdctmc;
		// Build (implicit) uniformised DTMC
		uniformizationRate = fdctmc.getMaxExitRate();// getDefaultUniformisationRate();
		uniformisedDTMC = fdctmc.buildUniformisedDTMC(uniformizationRate);
		computeMTEgroups();
        computeNumberOfSteps();
	}

	public void setInterval(double interval) {
		if (interval <= 0)
			throw new IllegalArgumentException(
					"Interval length must be greater than zero, but was "
							+ interval);
		intervalLength = interval;
	}

	protected void computeNumberOfSteps() throws PrismException {
		if (fdctmc == null)
			throw new NullPointerException("fdctmc was not initialized");
		if(intervalLength <= 0)
			throw new PrismException("Interval length must be larger than 0");
		for (FDEvent event : fdctmc.getAllFDEvents()) {
			numberOfSteps.put(event, event.getNumberOfSteps(intervalLength));
		}
	}

	public double getUniformizationRate() {
		return uniformizationRate;
	}

	public DTMCSimple getUniformisedDTMC() {
		return uniformisedDTMC;
	}

	public void compute(int numOfSteps, double intervalLength, int numOfPoissSteps)
			throws PrismException {
		this.intervalLength = intervalLength;
		initialize();

		System.out.println(fdctmc);
		
		System.out.println("GROUPS: " + MTEgroups);

		
		for (int i = 1; i <= numOfSteps; ++i) {
			Map<Integer , Set<PositionKey>> sameKeys = new HashMap<>();
			executeExponentialTransitions(intervalLength, numOfPoissSteps);

//			if(i > numOfSteps/2) {
//				HashSet<PositionKey> keys = new HashSet<>();
//				
//				for(PositionKey key2: transProb.keySet()) {
//					keys.add(key2);
//				}
//				
//				for(PositionKey key: keys) {
//					for(PositionKey key2: transProb.keySet()) {
//						if(!key.equals(key2) && key.hashCode()==key2.hashCode()) {
//							/*System.out.println("############################" + key.hashCode() + "==" + key2.hashCode());
//							System.out.println("key1" + key );
//							System.out.println("key2" + key2 );*/
//							Set setKey = sameKeys.get(key.hashCode());
//							if(setKey==null) {
//								setKey = new HashSet<PositionKey>();
//								sameKeys.put(key.hashCode(), setKey);
//							}
//							setKey.add(new PositionKey(key));
//							setKey.add(new PositionKey(key2));
//						}
//					} 
//				}
//	
//				for(int a : sameKeys.keySet())
//					if(sameKeys.get(a).size()>1)
//						System.out.println("HASH: " + a + " KEYS: " + sameKeys.get(a));
//			}
			
			executeFD();
			updateWaitingFD();
			
		}
		

		

	}

	public Distribution getResult() {
		return transProb.getDistr();
	}

	public List<FDSubset> getMTEGroups() {
		return MTEgroups;
	}

	public FDMatrix getTransProb() {
		return transProb;
	}

	void computeMTEgroups() {
		if (fdctmc.getNumFDEvents() > 25)
			throw new IllegalArgumentException(
					"Too many fix-delays. We support max 25 fix-delays.");

		long numSubGroups = Math.round((Math.pow(2, fdctmc.getNumFDEvents())));
		long longKey = 0;
		BitSet key, activeStates = new BitSet(fdctmc.getNumStates()), andStates = new BitSet(
				fdctmc.getNumStates());

		do {
			activeStates.set(0, fdctmc.getNumStates(), true);
			key = BitSet.valueOf(new long[] { longKey });

			for (int i = 0; i < fdctmc.getNumFDEvents(); ++i) {
				andStates.set(0, fdctmc.getNumStates(), false);

				andStates.or(fdctmc.getFDEvent(i).getActive());

				if (!key.get(i))
					andStates.flip(0, fdctmc.getNumStates());

				activeStates.and(andStates);

			}
			List<Integer> states = new ArrayList<>();

			int i = activeStates.nextSetBit(0);
			while (i != -1) {
				states.add(i);
				i = activeStates.nextSetBit(i + 1);
			}

			if (!states.isEmpty())
				MTEgroups.add(new FDSubset(key, states));

			longKey++; // get next subset of fdelays
		} while (longKey < numSubGroups);
	}

	protected void initialize() throws PrismException {
		Distribution pi = new Distribution(); // TODO maybe incorporate to later
												// code / no need to create
												// distribution
		if (fdctmc.getNumInitialStates() == 1) {
			pi.add(fdctmc.getFirstInitialState(), 1);
		} else {
			// pi.add(0, 1);
			throw new PrismException("Not one initial state: "
					+ fdctmc.getNumInitialStates() + ".");
		}

		for (FDSubset mTEGroup : MTEgroups) {
			BitSet fdSubset = mTEGroup.getActiveFDs();
			// key in matrix, remaining time on fdelay
			List<Integer> position = new ArrayList<>(fdctmc.getNumFDEvents());
			Distribution distr = new Distribution(); // value to store

			for (int i = 0; i < fdctmc.getNumFDEvents(); ++i) {
				if (fdSubset.get(i)) {
					position.add(numberOfSteps.get(fdctmc.getFDEvent(i)));
				} else
					position.add(0);
			}

			for (int i : mTEGroup.getActiveStates()) {
				distr.add(i, pi.get(i));
			}
			transProb.put(new PositionKey(position), distr);
		}
	}

	public void executeExponentialTransitions(double intervalLength,
			int numOfPoissSteps) {
		transProbResult.clear();
		FDMatrix<PositionKey> transProbTmp = null;
		double weight = Math.exp(-uniformizationRate * intervalLength);

		for (Map.Entry<PositionKey, Distribution> entry : transProb)
			transProbResult.put(new PositionKey(entry.getKey()), new Distribution(entry.getValue())
					.multiply(weight));

		
		for (int j = 1; j <= numOfPoissSteps; ++j) {
			weight = weight * uniformizationRate * intervalLength / j;
			// copy
			transProbTmp = transProbCopy;
			transProbCopy = transProb;
			transProb= transProbTmp; 
			transProb.clear();
			
			for (Map.Entry<PositionKey, Distribution> entry : transProbCopy) {
				Distribution distr2 = uniformisedDTMC
						.multiplyDistribution(entry.getValue());
				if (this instanceof CbpAlgLow)
					((CbpAlgLow)this).update(entry.getKey(), distr2, false, weight);
				else 
					update(entry.getKey(), distr2, false, weight);
			}
		}
		transProbTmp = transProb;
		transProb = transProbResult;
		transProbResult = transProbTmp;
	}

	private void update(PositionKey position, Distribution distr,
			boolean restart, double weight) {
                PositionKey position2 = new PositionKey(fdctmc.getNumFDEvents());

		for (FDSubset mTEGroup : MTEgroups) {
			position2.clear();

			for (int i = 0; i < mTEGroup.getActiveFDs().size(); ++i)
				if (mTEGroup.getActiveFDs().get(i))
					if (position.get(i) > 0) {
						position2.set(i, position.get(i));
					} else {
						position2.set(i,
								numberOfSteps.get(fdctmc.getFDEvent(i)) + 1);
					}
                                
			Distribution distr1 = transProbResult.getAndCreate(new PositionKey(position2));
			Distribution distr2 = transProb.getAndCreate(new PositionKey(position2));
			for (int j : mTEGroup.getActiveStates()) {
				if (distr.get(j) == 0)
					continue;
				distr1.add(j, distr.get(j) * weight);
				distr2.add(j, distr.get(j));
			}
		}
	}

	public void executeFD() throws PrismException {
		if (fdctmc.getNumFDEvents() == 0)
			return;

        Set<PositionKey> keys = new HashSet<>();
        for( PositionKey key : transProb.keySet()) {
            if(key.numOfOnes()!=0) {
                keys.add(new PositionKey(key));
            }
        }
		
        Set<PositionKey> newKeys = new HashSet<>();
        
        for (int i = fdctmc.getNumFDEvents(); i >= 1; --i) {
            for(Iterator<PositionKey> it = keys.iterator(); it.hasNext(); ) {
                PositionKey key = it.next();

		        if(key.numOfOnes()==i) {
		            List<Integer> activeFds = key.getActiveFd();
		            double overalWeight = getFdWeigth(activeFds);
		            for(int active : activeFds) {
		                double weight = fdctmc.getFDEvent(active).getWeight() / overalWeight;
		                Distribution distr = fdctmc.getFDEvent(active).multiplyDistribution(transProb.get(key));
						newKeys.addAll(updateByFD(distr, key, active, weight));
		            }
		            
		            if(!activeFds.isEmpty()) {
		                transProb.remove(key);
		            }
		        }
            }
            keys.addAll(newKeys);
        	newKeys.clear();
		}
	}


	private Set<PositionKey> updateByFD(Distribution distr, PositionKey position, int fDIndex, double weight) {
		Set<PositionKey> newKeysToProcess = new HashSet<PositionKey>();
		
                PositionKey position2 = new PositionKey(fdctmc.getNumFDEvents());
                Distribution distr2 = new Distribution();
			
                for (FDSubset mTEGroup : MTEgroups) {
                        position2.clear();
			for (int j = 0; j < fdctmc.getNumFDEvents(); ++j) {
				if (mTEGroup.getActiveFDs().get(j)) {
					if (j != fDIndex && position.get(j) != 0)
						position2.set(j, position.get(j));
					else
						position2.set(j, numberOfSteps.get(fdctmc.getFDEvent(j)) + 1);
				}
			}
			
			
			distr2.clear();
			
			for (int j : mTEGroup.getActiveStates()) {
				if (distr.get(j) == 0)
					continue;
				distr2.add(j, distr.get(j) * weight);
			}
			
			if(distr2.sum() == 0) continue;
			
			PositionKey position3 = new PositionKey(position2);
			
            transProb.getAndCreate(position3).add(distr2); 
            
            if(position3.numOfOnes()>0) newKeysToProcess.add(position3);
		}
		
		return newKeysToProcess;
	}


    protected void updateWaitingFD() {
    	FDMatrix<PositionKey> result = new FDMatrix<>(transProb.size());

		for (Map.Entry<PositionKey, Distribution> entry : transProb) {
                        //TODO it is possible to do smartly- only change the key/ possible extension with different structure of transProb
			PositionKey position2 = new PositionKey(fdctmc.getNumFDEvents()); 

			for (int i = 0; i < fdctmc.getNumFDEvents(); ++i) {
				if (entry.getKey().get(i) > 1) {
					position2.set(i, entry.getKey().get(i) - 1);
				}
			}
			result.put(position2, entry.getValue());
		}
		transProb = result;
	}
        /**
	 * Calculates the number of Poisson steps to avoid accumulation of too large
	 * error.
	 * 
	 * @param lambdaT
	 *            - uniformization rate, multiplied by length of interval
	 * @param numOfIntervals
	 *            - number of same intervals to which the time bound is divided
	 * @param epsilon
	 *            - error
	 * @return
	 * @throws PrismException 
	 */
	public static int getNumberOfPoissonStepsFromError(double delta,
			int numOfIntervals, double epsilon, double rate) throws PrismException {
		// computattion of normal poisson
		Double poissonCDF = 0.0;
		double lambdaTWorkingCopy = 1;
		long factorial = 1;
		Double error = (1-epsilon) * Math.exp(rate * delta * numOfIntervals); 

		double lambdaT = rate * delta;
		
		if(error.isInfinite() || error.isNaN() || error == 0)
			throw new PrismException("Too big numbers: (1-epsilon) * Math.exp(lambdaT) is infinity, NaN, or 0.");

		int j = 0;
		while (Math.pow(poissonCDF, numOfIntervals) < error) {
			poissonCDF += lambdaTWorkingCopy / factorial;
			lambdaTWorkingCopy = lambdaTWorkingCopy * lambdaT;
			j++;
			factorial = factorial * j;
			if(poissonCDF.isInfinite() || poissonCDF.isNaN() || poissonCDF == 0)
				throw new IllegalArgumentException("Too big numbers, computation reached infinity, NaN, or 0.");
			}
		return (j-1);
	}
	
	public static double getPoissonError(double delta,
			int numOfIntervals, double numOfPoissonSteps, double rate) throws PrismException {
		// computattion of normal poisson
		Double poissonCDF = 0.0;
		double lambdaTWorkingCopy = 1;
		long factorial = 1;
		Double weight = Math.exp(-rate * delta * numOfIntervals); 

		double lambdaT = rate * delta;
		
		if(weight.isInfinite() || weight.isNaN() || weight == 0)
			throw new PrismException("Too big numbers: (1-epsilon) * Math.exp(lambdaT) is infinity, NaN, or 0.");

		for(int j=0;j<= numOfPoissonSteps;){
			poissonCDF += lambdaTWorkingCopy / factorial;
			lambdaTWorkingCopy = lambdaTWorkingCopy * lambdaT;
			j++;
			factorial = factorial * j;
			if(poissonCDF.isInfinite() || poissonCDF.isNaN() || poissonCDF == 0)
				throw new IllegalArgumentException("Too big numbers, computation reached infinity, NaN, or 0.");
		}

		return 1-Math.pow(poissonCDF,numOfIntervals) * weight;
	}

    protected double getFdWeigth(List<Integer> fds) {
        double weight = 0;
        for(int fd : fds) {
            weight += fdctmc.getFDEvent(fd).getWeight();
        }
        return weight;
    }
	
}
