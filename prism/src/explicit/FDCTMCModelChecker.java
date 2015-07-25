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

import java.io.File;
import java.math.BigInteger;
import java.util.*;

import explicit.rewards.MCRewards;
import explicit.rewards.StateRewardsArray;
import parser.State;
import parser.ast.*;
import parser.type.*;
import prism.*;

/**
 * Explicit-state model checker for continuous-time Markov chains (CTMCs).
 * 
 * This uses various bits of functionality of DTMCModelChecker, so we inherit from that.
 * (This way DTMCModelChecker picks up any options set on this one.) 
 */
public class FDCTMCModelChecker extends CTMCModelChecker
{
	// Steady-state/transient probability computation

	/**
	 * Compute transient probability distribution (forwards).
	 * Start from initial state (or uniform distribution over multiple initial states).
	 */
	public StateValues doTransient(FDCTMCSimple fdctmc, double time) throws PrismException
	{
		return doTransient(fdctmc, time, (StateValues) null);
	}

	/**
	 * Compute transient probability distribution (forwards).
	 * Optionally, use the passed in file initDistFile to give the initial probability distribution (time 0).
	 * If null, start from initial state (or uniform distribution over multiple initial states).
	 * @param fdctmc The FDCTMC
	 * @param t Time point
	 * @param initDistFile File containing initial distribution
	 */
	public StateValues doTransient(FDCTMCSimple fdctmc, double t, File initDistFile) throws PrismException
	{
		StateValues initDist = readDistributionFromFile(initDistFile, fdctmc);
		return doTransient(fdctmc, t, initDist);
	}

	/**
	 * Compute transient probability distribution (forwards).
	 * Optionally, use the passed in vector initDist as the initial probability distribution (time 0).
	 * If null, start from initial state (or uniform distribution over multiple initial states).
	 * For reasons of efficiency, when a vector is passed in, it will be trampled over,
	 * so if you wanted it, take a copy. 
	 * @param fdctmc The FDCTMC
	 * @param t Time point
	 * @param initDist Initial distribution (will be overwritten)
	 */
	public StateValues doTransient(FDCTMCSimple fdctmc, double t, StateValues initDist) throws PrismException
	{
		StateValues initDistNew = null, probs = null;

		// Build initial distribution (if not specified)
		if (initDist == null) {
			initDistNew = new StateValues(TypeDouble.getInstance(), new Double(0.0), fdctmc);
			double initVal = 1.0 / fdctmc.getNumInitialStates();
			for (int in : fdctmc.getInitialStates()) {
				initDistNew.setDoubleValue(in, initVal);
			}
		} else {
			initDistNew = initDist;
		}
                double discretizationStep = 0.01, error = 0.01;
                int numOfIntervals = (int) Math.round( t/discretizationStep);
                
		// Compute transient probabilities
                CbpAlg cbpAlg = new CbpAlg(mainLog);
                cbpAlg.setInterval(discretizationStep);
                cbpAlg.setFdctmc(fdctmc);
                Distribution distr = cbpAlg.runCBP(discretizationStep, numOfIntervals, error);
                double[] res = new double[fdctmc.getNumStates()];
                for(int i = 0; i<fdctmc.getNumStates(); ++i) {
                    res[i] = distr.get(i);
                }
                
		probs = StateValues.createFromDoubleArray(res, fdctmc);

		return probs;
	}

	
}
