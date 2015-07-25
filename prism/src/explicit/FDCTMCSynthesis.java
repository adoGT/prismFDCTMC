package explicit;

import java.util.BitSet;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;
import java.util.Vector;
import java.util.Map.Entry;

import prism.PrismException;
import explicit.rewards.FDCTMCRewardsSimple;
import explicit.rewards.MCRewards;
import explicit.rewards.MDPRewardsSimple;
import explicit.rewards.StateRewardsArray;

public class FDCTMCSynthesis {

	FDCTMCSimple model;
	FDCTMCRewardsSimple fDCTMCRewards; 
	List<Integer> targetStatesFDCTMC;
	
	private MDPSimple mdp;
	private MDPRewardsSimple mdpRewards;
	
	private Map<Integer,Integer> fDCTMCStateToFDIndex;
	private Map<Integer,Set<Integer>> fDIndexToPotato; //FD to potato
	private Map<Integer,Set<Integer>> fDIndexToOutStates;
	private Map<Integer,Integer> fDIndexToEntry; //FD to entry
	private Map<Integer,Integer> fDCTMCStateToMDPState;
	private Map<Integer,Integer> mDPStateToFDCTMCState;	
	
	double uniformizationRate;
	DTMCSimple uniformisedDTMC;
	
	public void computeSynthesis(FDCTMCSimple model, FDCTMCRewardsSimple fDCTMCRewards, 
			List<Integer> targetStatesFDCTMC, double error) throws PrismException{
		this.model = model;
		this.fDCTMCRewards = fDCTMCRewards;
		this.targetStatesFDCTMC = targetStatesFDCTMC;

		initialise();
		
		
		/////////////////////add choices to MDP and fill also transition rewards
	
		uniformizationRate = model.getMaxExitRate();// getDefaultUniformisationRate();
		uniformisedDTMC = model.buildUniformisedDTMC(uniformizationRate);

		
		
		//computation of relevant constants
		double boundExpectedCost;
		double boundNumOfSteps;
		double alpha;
		double d1;

		
		//////////// boundExpectedCost

		//we use values of FD from models and solve model
		processNonPotatoStates();
		
		//compute choices and transition costs corresponding to entry to potatoes
		
		//first we have to evaluate bound on expected cost from every state
		//we will use initial values of FDs fill corresponding MDP and solve it
		
		for(Entry<Integer,Integer> entry : fDIndexToEntry.entrySet()) {

			processPotatoState(entry.getKey(), entry.getValue(), model.getFDEvent(entry.getKey()).getDelayTime(),
					model.getFDEvent(entry.getKey()).getDelayTime(), 1e-11); //TODO error should be computed accordingly, i.e. try some error, then test perturbation bound, if it is OK, let it be, otherwise make another iteration with smaller error
		}

		//model check created MDP
		MDPModelChecker mDPModelChecker = new MDPModelChecker(); 
		BitSet target = new BitSet(mdp.getNumStates());
		for (int i :targetStatesFDCTMC) {
    		target.set(fDCTMCStateToMDPState.get(i), true);
		}
		
		ModelCheckerResult result;
		result = mDPModelChecker.computeReachRewards(mdp, mdpRewards, target, true, null, null);

		//take maximum from all MDP states
		boundExpectedCost = 0.0;
		for(int i=0; i<mdp.getNumStates();i++){
			if(boundExpectedCost < result.soln[i]) 
				boundExpectedCost = result.soln[i];
		}
		
		System.out.println("BOUND ON EXPECTED COST FROM EACH STATE: " + boundExpectedCost);
		

		///////// boundNumOfSteps
		
		//we check whether each entry of potato has FD transition with positive impulse cost 
		double minimalOneStepCost = Double.MAX_VALUE;
		for(int i=0;i<model.getNumStates();i++){
			if(targetStatesFDCTMC.contains(i)) continue;
			if(minimalOneStepCost > fDCTMCRewards.getStateReward(i)/uniformizationRate)
				minimalOneStepCost = fDCTMCRewards.getStateReward(i)/uniformizationRate;
		}
		
		for(Entry<Integer,Integer> entry : fDIndexToEntry.entrySet()) {
			double currentImpulseCost = 0.0;
			Iterator<Entry<Integer,Double>> it;
			it = model.getFDEvent(entry.getKey()).getTransitionsIterator(entry.getValue());
			while(it.hasNext()) {
				Entry<Integer,Double> entry2 =it.next();
				currentImpulseCost = currentImpulseCost 
									+ entry2.getValue() * fDCTMCRewards.getTransitionReward(entry.getKey(), 
											entry.getValue() , entry2.getKey());
			}
			if(currentImpulseCost<=0)
				throw new PrismException("State " + entry.getValue() 
						+ " is entry of potato but does not have positive cost (" + currentImpulseCost +") on FD transition!");
			if(minimalOneStepCost > currentImpulseCost )
				minimalOneStepCost = currentImpulseCost;
		}
		
//		System.out.println("MINIMAL ONE STEP COST IS: " + minimalOneStepCost);
		
		boundNumOfSteps = boundExpectedCost/minimalOneStepCost;
	
		System.out.println("BOUND ON EXPECTED NUMBER OF STEPS FROM EACH STATE: "+ boundNumOfSteps);
		
		
		/////////////// alpha
		
		alpha = error /(boundNumOfSteps * (1+boundExpectedCost*mdp.getNumStates()));
		double alpha2 = 1/(2*boundNumOfSteps*mdp.getNumStates());
		if(alpha2 < alpha ) alpha = alpha2;
		
		System.out.println("ALPHA: " + alpha);
		
		////////////// d1
		
		//computation of minimal and maximal nonzero costs/rewards
		double maxReward = 0.0; // = 2;
		double minReward = Double.MAX_VALUE; // = 1;
		double minBranchingPst = Double.MAX_VALUE; // = 0.5;
		
		//exponential transitions
		for(int i =0; i<model.getNumStates();i++) {
			Iterator<Entry<Integer,Double>> it;
			it = model.getTransitionsIterator(i);
			while(it.hasNext()) {
				Entry<Integer,Double> entry2 =it.next();
				if(minBranchingPst > entry2.getValue() && entry2.getValue() > 0)
					minBranchingPst = entry2.getValue();
				if(minReward > fDCTMCRewards.getTransitionReward(-1,i, entry2.getKey()) && fDCTMCRewards.getTransitionReward(-1,i, entry2.getKey()) > 0)
					minReward=fDCTMCRewards.getTransitionReward(-1,i, entry2.getKey());
				if(maxReward < fDCTMCRewards.getTransitionReward(-1,i, entry2.getKey()))
					maxReward = fDCTMCRewards.getTransitionReward(-1,i, entry2.getKey());
			}
		}
		
		//FD transitions
		for(int j : model.getAllFDEventsIndexes()){
			for(int i =0; i<model.getNumStates();i++) {
				Iterator<Entry<Integer,Double>> it;
				it = model.getFDEvent(j).getTransitionsIterator(i);
				while(it.hasNext()) {
					Entry<Integer,Double> entry2 =it.next();
					if(minBranchingPst > entry2.getValue() && entry2.getValue() > 0)
						minBranchingPst = entry2.getValue();
					if(minReward > fDCTMCRewards.getTransitionReward(j,i, entry2.getKey()) && fDCTMCRewards.getTransitionReward(j,i, entry2.getKey()) > 0)
						minReward=fDCTMCRewards.getTransitionReward(j,i, entry2.getKey());
					if(maxReward < fDCTMCRewards.getTransitionReward(j,i, entry2.getKey()))
						maxReward = fDCTMCRewards.getTransitionReward(j,i, entry2.getKey());
				}
			}
		}
		
		//state rewards
		for(int i =0; i<model.getNumStates();i++) {
			if(minReward > fDCTMCRewards.getStateReward(i)&& fDCTMCRewards.getStateReward(i) > 0)
				minReward=fDCTMCRewards.getStateReward(i);
			if(maxReward < fDCTMCRewards.getStateReward(i))
				maxReward = fDCTMCRewards.getStateReward(i);
		}
		
		System.out.println("\n\tmaxReward: " + maxReward + "\n\tminReward: " + minReward + "\n\tminBranchingPst: "+ minBranchingPst);
		
		//end of computation of minimal and maximal nonzero costs/rewards
		
		
		d1 = 2*(uniformizationRate + 1) * maxReward; 
		if (d1< 2* uniformizationRate) d1 = 2* uniformizationRate;
		
		System.out.println("D1: " + d1);
		
		//compute delta and dmax for every potato    		
		double delta = alpha/d1;
		double kappa = (error * delta * minReward)/(2* mdp.getNumStates()*(1+Math.pow(boundExpectedCost,2)));
		double dMax  = boundExpectedCost/(Math.pow(minBranchingPst, model.getNumStates()-mdp.getNumStates()+model.getNumFDEvents())*minReward);
		double dMax2 = (Math.E * Math.abs(Math.log(alpha/2)))/(uniformizationRate*minBranchingPst); 
		if(dMax2 > dMax) dMax = dMax2;
		
		long timer;
		

		for(Entry<Integer,Integer> entry : fDIndexToEntry.entrySet()) {

			System.out.println("Computation of potato for FD " + entry.getKey());
			
			   			
//    		delta =  0.009;//0.000009;    		
//    		kappa = 1e-9; //TODO pri 14 zacina pocitat extremne pomaly!!! mozno nestaci presnost double!!!
//    		dMax  = 24.8;
//    		System.out.println("\n\tDELTA: " + delta + "\n\tKAPPA: " + kappa + "\n\tD_MAX: "+ dMax + "\n");

			processPotatoState(entry.getKey(), entry.getValue(), delta, dMax, kappa);
		} //done computation of this potato
		

//		System.out.println("\nBUILT MDP "+ mdp);
//		System.out.println("BUILT MDPRewards: "+ mdpRewards + "\n");

		System.out.println("\nComputing optimal strategy for created MDP.\n");
		timer = System.currentTimeMillis();
		
		//model check created MDP
		result = mDPModelChecker.computeReachRewards(mdp, mdpRewards, target, true, null, null);

		timer =System.currentTimeMillis() - timer;
		System.out.println("\nOptimal strategy computed in " + (timer/1000.0) + " s.\n");
		
		System.out.println("RESULT in initial state: " + result.soln[mdp.getFirstInitialState()] );

//		mDPModelChecker.expReachStrategy(mdp, mdpRewards, mdp.getFirstInitialState(), target, true, result.lastSoln);
	//TODO doesn't work I don't know wht is wrong	
	}
	

	//copied from newer version of prism
    //temporary solution- when we learn how to compute accumulated rate rewards more efficiently, we will not need this anymore
	/**
	 * Perform cumulative reward computation.
	 * Compute, for each state of {@ctmc}, the expected rewards accumulated until {@code t}
	 * when starting in this state and using reward structure {@code mcRewards}.
	 * @param ctmc The CTMC
	 * @param mcRewards The rewards
	 * @param t Time bound //TODO
	 */
	private double computeCumulativeRewards(DTMC dtmc, double uniformisationRate, 
			MCRewards mcRewards, double t, double error, int initState, Set<Integer> targetStates) throws PrismException
	{
//		ModelCheckerResult res = null;
		int i, n, iters;
		double soln[], soln2[], tmpsoln[], sum[];
//		long timer;
		// Fox-Glynn stuff
		FoxGlynn fg;
		int left, right;
		double q, qt, acc, weights[], totalWeight;

		// Optimisation: If t = 0, this is easy.
		if (t == 0) {
//			res = new ModelCheckerResult();
//			res.soln = new double[ctmc.getNumStates()];
//			return res;
			return 0;
		}

		// Start backwards transient computation
//		timer = System.currentTimeMillis();
//		mainLog.println("\nStarting backwards cumulative rewards computation...");

		// Store num states
		n = dtmc.getNumStates();

		// Get uniformisation rate; do Fox-Glynn
//		q = ctmc.getDefaultUniformisationRate();
		q = uniformisationRate; //TODO was like above- maybe default uniformisation rate is needed instead of max exit rate!
		qt = q * t;
//		mainLog.println("\nUniformisation: q.t = " + q + " x " + t + " = " + qt);
//		acc = termCritParam / 8.0;
		acc = error / 8.0;
		fg = new FoxGlynn(qt, 1e-300, 1e+300, acc);
		left = fg.getLeftTruncationPoint();
		right = fg.getRightTruncationPoint();
		if (right < 0) {
			throw new PrismException("Overflow in Fox-Glynn computation (time bound too big?)");
		}
		weights = fg.getWeights();
		totalWeight = fg.getTotalWeight();
		for (i = left; i <= right; i++) {
			weights[i - left] /= totalWeight;
		}

		// modify the poisson probabilities to what we need for this computation
		// first make the kth value equal to the sum of the values for 0...k
		for (i = left+1; i <= right; i++) {
			weights[i - left] += weights[i - 1 - left];
		}
		// then subtract from 1 and divide by uniformisation constant (q) to give mixed poisson probabilities
		for (i = left; i <= right; i++) {
			weights[i - left] = (1 - weights[i - left]) / q;
		}
//		mainLog.println("Fox-Glynn (" + acc + "): left = " + left + ", right = " + right);

		// Build (implicit) uniformised DTMC
		DTMC dtmcUnif = dtmc;

		// Create solution vector(s)
		soln = new double[n];
		soln2 = new double[n];

		// Initialise solution vectors.
		for (i = 0; i < n; i++)
			soln[i] = mcRewards.getStateReward(i);

		// do 0th element of summation (doesn't require any matrix powers)
		sum = new double[n];
		if (left == 0) {
			for (i = 0; i < n; i++)
				sum[i] += weights[0] * soln[i];
		} else {
			for (i = 0; i < n; i++)
				sum[i] += soln[i] / q;
		}

		// Start iterations
		iters = 1;
		while (iters <= right) {
			// Matrix-vector multiply
			dtmcUnif.mvMult(soln, soln2, null, false);
			// Swap vectors for next iter
			tmpsoln = soln;
			soln = soln2;
			soln2 = tmpsoln;
			// Add to sum
			if (iters >= left) {
				for (i = 0; i < n; i++)
					sum[i] += weights[iters - left] * soln[i];
			} else {
				for (i = 0; i < n; i++)
					sum[i] += soln[i] / q;
			}
			iters++;
		}

		// Finished backwards transient computation
//		timer = System.currentTimeMillis() - timer;
//		mainLog.print("Backwards transient cumulative rewards computation");
//		mainLog.println(" took " + iters + " iters and " + timer / 1000.0 + " seconds.");

		// Return results
//		res = new ModelCheckerResult();
//		res.soln = sum;
//		res.lastSoln = soln2;
//		res.numIters = iters;
//		res.timeTaken = timer / 1000.0;
//		res.timePre = 0.0;
		return sum[initState];
	}
	
	/**
	 * Creates mdp and mdpRewards. Fills states to mdp.
	 * Creates and fully fills all remaining global variables.
	 * @throws PrismException 
	 */
	private void initialise() throws PrismException{
		////////////////////creation of MDP
		mdp = new MDPSimple();
        
        //find potatoes
        
        fDCTMCStateToFDIndex = new HashMap<Integer,Integer>(model.getNumStates());
		fDIndexToPotato = new HashMap<Integer,Set<Integer>>(); //FD to potato
		fDIndexToOutStates = new HashMap<Integer,Set<Integer>>(model.getNumFDEvents());
		
		Set<Integer> potato;
		FDEvent fD;
        for(int j=0; j<model.getNumFDEvents();++j) {
        	fD = model.getFDEvent(j); 
        	potato = new HashSet<Integer>();
        	for(int i=0;i<model.getNumStates();++i){
        		if (fD.isActive(i)){
        			if(fDCTMCStateToFDIndex.get(i) != null && fDCTMCStateToFDIndex.get(i) != fD.getIndex() )
        				throw new PrismException("Two active FDs in state " + i + "! Not supported.");
        			fDCTMCStateToFDIndex.put(i, j);
        			potato.add(i);
        		}
        	}
        	if(!potato.isEmpty()) fDIndexToPotato.put(j,potato);
        	fDIndexToOutStates.put(j, new HashSet<Integer>());
        }
        
        System.out.println("fDCTMCStateToFDIndex: " + fDCTMCStateToFDIndex);
        System.out.println("fDIndexToPotato: " + fDIndexToPotato);
        
        
        //fill necessary structures
		fDIndexToEntry = new HashMap<Integer,Integer>(model.getNumFDEvents()); //FD to entry
		fDCTMCStateToMDPState = new HashMap<Integer,Integer>();
		mDPStateToFDCTMCState = new HashMap<Integer,Integer>();

		
		
		if(model.getNumInitialStates() != 1)
			System.out.println("WARNING: detected multiple initial states. Using the firs one.");


		//add newly found state
		int initialState = model.getFirstInitialState();
		int newMDPState = mdp.addState();
		mdp.addInitialState(newMDPState);
		fDCTMCStateToMDPState.put(initialState, newMDPState);
		mDPStateToFDCTMCState.put(newMDPState, initialState);
		
		//process initial state
		if(fDCTMCStateToFDIndex.get(initialState)!=null){
			//set entry of potato to corresponding FD
			fDIndexToEntry.put(fDCTMCStateToFDIndex.get(initialState), initialState);
		}
			
		//structures for traversing
		Stack<Integer> toDoStates = new Stack<Integer>(); //grey states
		Set<Integer> doneStates = new HashSet<Integer>(model.getNumStates()); //black states
		Integer currentState, nextState;
		
		toDoStates.push(model.getFirstInitialState());
		
		while (!toDoStates.isEmpty()) {
			currentState = toDoStates.pop();
			doneStates.add(currentState);
			
			//explore exponential successors 
			Iterator<Entry<Integer,Double>> it = model.getTransitions(currentState).iterator(); //iterate over distribution put to stack, check conditions!!!
			while (it.hasNext()) {
				nextState = ((Entry<Integer,Double>)it.next()).getKey();
				
				if(!(doneStates.contains(nextState) || toDoStates.contains(nextState))){ //not found yet
					toDoStates.push(nextState); //process traversing information
				
					//check whether there is an active FD event
					if(fDCTMCStateToFDIndex.get(nextState) == null){ //there is none, it is an MDP state
			    		newMDPState = mdp.addState();
						fDCTMCStateToMDPState.put(nextState, newMDPState);
						mDPStateToFDCTMCState.put(newMDPState, nextState);    						
					}
				}
				
				//check whether this transition goes to entry of potato
				if(fDCTMCStateToFDIndex.get(nextState) != null &&  				//some FD must be active in nextState
						(fDCTMCStateToFDIndex.get(currentState) == null 		//and was not active before
																				//or before was active different FD
						  || (fDCTMCStateToFDIndex.get(currentState) != null 
						  		&& fDCTMCStateToFDIndex.get(currentState)!=fDCTMCStateToFDIndex.get(nextState))))		
				{
					if(!(doneStates.contains(nextState) || toDoStates.contains(nextState))){ //not found yet
			    		newMDPState = mdp.addState();
						fDCTMCStateToMDPState.put(nextState, newMDPState);
						mDPStateToFDCTMCState.put(newMDPState, nextState);
					}
					
					//check whether there is different entry to potato
					if(fDIndexToEntry.get(fDCTMCStateToFDIndex.get(nextState)) != null 
							&& fDIndexToEntry.get(fDCTMCStateToFDIndex.get(nextState)) != nextState)
						throw new PrismException("There are ,multiple entries to potato!");
					fDIndexToEntry.put(fDCTMCStateToFDIndex.get(nextState), nextState);
				}
				
				//check whether transition leaves potato
				if(fDCTMCStateToFDIndex.get(currentState) != null &&  			//some FD must be active in currentState
						(fDCTMCStateToFDIndex.get(nextState) == null 			//and was not active in nextState
																				//or different FD is active in nextState
						  || (fDCTMCStateToFDIndex.get(nextState) != null 
						  		&& fDCTMCStateToFDIndex.get(currentState)!=fDCTMCStateToFDIndex.get(nextState))))
					fDIndexToOutStates.get(fDCTMCStateToFDIndex.get(currentState)).add(nextState);
			}
			
			//explore FD successors
			if(fDCTMCStateToFDIndex.get(currentState) == null) continue; //no FD
			
			it = model.getFDEvent(fDCTMCStateToFDIndex.get(currentState)).getTransitionsIterator(currentState);
			while (it.hasNext()) {
				nextState = ((Entry<Integer,Double>)it.next()).getKey();
				
				if(!(doneStates.contains(nextState) || toDoStates.contains(nextState))){ //not found yet
					toDoStates.push(nextState); //process traversing information

					//check whether there is an active FD event
					if(fDCTMCStateToFDIndex.get(nextState) == null){ //there is none, it is an MDP state
			    		newMDPState = mdp.addState();
						fDCTMCStateToMDPState.put(nextState, newMDPState);
						mDPStateToFDCTMCState.put(newMDPState, nextState);    						
					}
				}
				
				//check whether this transition goes to entry of potato
				if(fDCTMCStateToFDIndex.get(nextState) != null) {
					if(!(doneStates.contains(nextState) || toDoStates.contains(nextState))){ //not found yet
			    		newMDPState = mdp.addState();
						fDCTMCStateToMDPState.put(nextState, newMDPState);
						mDPStateToFDCTMCState.put(newMDPState, nextState);
					}
					
					//check whether there is already different entry to potato
					if(fDIndexToEntry.get(fDCTMCStateToFDIndex.get(nextState)) != null 
							&& fDIndexToEntry.get(fDCTMCStateToFDIndex.get(nextState)) != nextState)
						throw new PrismException("There are ,multiple entries to potato!");
					fDIndexToEntry.put(fDCTMCStateToFDIndex.get(nextState), nextState);

				}
				
				//this transition leaves potato because it is FD
				fDIndexToOutStates.get(fDCTMCStateToFDIndex.get(currentState)).add(nextState);
			}
		}
		//done traversing
		
		
		System.out.println("fDIndexToEntry: "+ fDIndexToEntry);
		System.out.println("fDCTMCStateToMDPState: " + fDCTMCStateToMDPState);
		System.out.println("mDPStateToFDCTMCState: " + mDPStateToFDCTMCState);
		System.out.println("fDIndexToOutStates: "+ fDIndexToOutStates);

		System.out.println("MDP "+ mdp);
		
		mdpRewards = new MDPRewardsSimple(mdp.getNumStates());
	}
	
	/**
	 * Adds one choice to every non-potato state to mdp
	 * and corresponding rewards to mdpRewards.
	 */
    private void processNonPotatoStates() {
		//choices corresponding to non-entry to potato estates
		
		Iterator<Entry<Integer,Double>> it;
		double overallCost;
		int src;
		for(int i=0;i<mdp.getNumStates();++i) {
			src = mDPStateToFDCTMCState.get(i);
			if(fDIndexToPotato.containsKey(src)) continue; //entry-of potato -> continue
				
			overallCost = (1/uniformizationRate) * fDCTMCRewards.getStateReward(src); //add rate cost
			
			it = uniformisedDTMC.getTransitionsIterator(src);
			Distribution distr = new Distribution();
			while (it.hasNext()) {
				Entry<Integer,Double> entry = ((Entry<Integer,Double>)it.next());

				distr.add(fDCTMCStateToMDPState.get(entry.getKey()), entry.getValue());
				// add transition cost, i.e. trans_probability * impulse_cost
				overallCost += entry.getValue() * fDCTMCRewards.getTransitionReward(-1, src, entry.getKey());
			}
			mdp.addChoice(i, distr);
			mdpRewards.addToTransitionReward(i, 0, overallCost);
		}
	}
    
	public void processPotatoState(int fdIndex, int entryToPotato, double delta, double dMax, double kappa) throws PrismException {
	
		//clear choices/rewards for this state in mdp/mdpRewards (new values will be computed)
		mdp.clearChoices(fDCTMCStateToMDPState.get(entryToPotato));
		mdpRewards.clearRewards(fDCTMCStateToMDPState.get(entryToPotato));
		
		DTMCSimple expDTMC, fDDTMC, expTransRewDTMC, fDTransRewDTMC;
		StateRewardsArray expDTMCRewards;
		Set<Integer> dTMCOutStates; 
		int numOfPoissonSteps = CbpAlg.getNumberOfPoissonStepsFromError(delta, ((Double)(dMax/delta)).intValue(), kappa, uniformizationRate);
	
		System.out.println("NUM of poiss steps: " + numOfPoissonSteps);
		
		//initialization
		expDTMC = new DTMCSimple();
		fDDTMC = new DTMCSimple();
		expTransRewDTMC = new DTMCSimple();
		fDTransRewDTMC = new DTMCSimple();
		expDTMCRewards = new StateRewardsArray(fDIndexToPotato.get(fdIndex).size() 
												+ fDIndexToOutStates.get(fdIndex).size());
		dTMCOutStates = new HashSet<Integer>();
	
		Vector<Integer> dTMCStatesToFDCTMCState = new Vector<Integer>(fDIndexToOutStates.get(fdIndex).size());
		Map<Integer,Integer> fDCTMCStateToDTMCState = new HashMap<Integer,Integer>(fDIndexToOutStates.get(fdIndex).size());
		int lastIndex, indexOfFirstNonPotatoState;
		
		//add states
		for(int state: fDIndexToPotato.get(fdIndex)){
			lastIndex = expDTMC.addState();
			fDDTMC.addState();
			expTransRewDTMC.addState();
			fDTransRewDTMC.addState();
			dTMCStatesToFDCTMCState.add(state);
			fDCTMCStateToDTMCState.put(state, lastIndex);
			
			expDTMCRewards.setStateReward(lastIndex, fDCTMCRewards.getStateReward(state));
		}
		indexOfFirstNonPotatoState = expDTMC.getNumStates();
		
		for(int state: fDIndexToOutStates.get(fdIndex)){
			lastIndex = expDTMC.addState();
			fDDTMC.addState();
			expTransRewDTMC.addState();
			fDTransRewDTMC.addState();        			
			dTMCStatesToFDCTMCState.add(state);
			fDCTMCStateToDTMCState.put(state, lastIndex);
	
			dTMCOutStates.add(lastIndex);
	//		expDTMCRewards.setStateReward(lastIndex, 0); //TODO if exception may be uncommented
		}
		
		
		//add exp transitions
		Iterator<Entry<Integer,Double>> it;
		for(int state: fDIndexToPotato.get(fdIndex)){
			it = uniformisedDTMC.getTransitionsIterator(state);
			while(it.hasNext()) {
				Entry<Integer,Double> entry2 = ((Entry<Integer,Double>)it.next());
				expDTMC.addToProbability(fDCTMCStateToDTMCState.get(state),
						fDCTMCStateToDTMCState.get(entry2.getKey()),
						entry2.getValue());
				//if trans reward is positive add it to reward kernel 
				if(fDCTMCRewards.getTransitionReward(-1,state,entry2.getKey()) > 0)
					expTransRewDTMC.addToProbability(fDCTMCStateToDTMCState.get(state),
							fDCTMCStateToDTMCState.get(entry2.getKey()),
							entry2.getValue() * fDCTMCRewards.getTransitionReward(-1,state,entry2.getKey()));
			}
		}
	
		//adding self loops
		for(int i = indexOfFirstNonPotatoState;i < expDTMC.getNumStates();i++){
			expDTMC.addToProbability(i,i,1.0);
		}
		
		//add FD transitions
		for(int state: fDIndexToPotato.get(fdIndex)){
			it = model.getFDEvent(fdIndex).getTransitionsIterator(state);
			while(it.hasNext()) {
				Entry<Integer,Double> entry2 = ((Entry<Integer,Double>)it.next());
				fDDTMC.addToProbability(fDCTMCStateToDTMCState.get(state),
						fDCTMCStateToDTMCState.get(entry2.getKey()),
						entry2.getValue());
				//if trans reward is positive add it to reward kernel 
				if(fDCTMCRewards.getTransitionReward(fdIndex,state,entry2.getKey()) > 0)
					fDTransRewDTMC.addToProbability(fDCTMCStateToDTMCState.get(state),
							fDCTMCStateToDTMCState.get(entry2.getKey()),
							entry2.getValue() * fDCTMCRewards.getTransitionReward(fdIndex,state,entry2.getKey()));
	
			}
		}
		
		System.out.println("\n expDTMC: "+ expDTMC);
		System.out.println("\n fDDTMC: "+ fDDTMC);
		
		System.out.println("\n expTransRewDTMC: "+ expTransRewDTMC);
		System.out.println("\n fDTransRewDTMC: "+ fDTransRewDTMC);
	
		System.out.print("\n expDTMCRewards[");
		for(int i=0;i<(expDTMC.getNumStates()-1);i++)
			System.out.print(expDTMCRewards.getStateReward(i) + ", ");
		System.out.println(expDTMCRewards.getStateReward(expDTMC.getNumStates()-1) + "]");
	
		
		double rateReward;
		
		
		//////////////////computing and filling rate rewards
		
		System.out.println("\nComputing rate costs.\n");
		long timer = System.currentTimeMillis();
		
		for(int i=1;i * delta <= dMax;++i){
			rateReward = computeCumulativeRewards(expDTMC, uniformizationRate, 
					expDTMCRewards, i * delta, kappa, fDCTMCStateToDTMCState.get(fdIndex),
					dTMCOutStates);
			
			mdpRewards.addToTransitionReward(fDCTMCStateToMDPState.get(fdIndex), 
												i-1, rateReward);
		}
		
		timer = System.currentTimeMillis()-timer;
		System.out.println("\nRate costs computed in " + (timer/1000.0) + " s.\n");
		
		
	//	System.out.println("COSTS: " + mdpRewards);
		
		
		
		/////////////////computing transition rewards
		
		System.out.println("\nComputation of transition probabilities and transition costs.\n");
		timer = System.currentTimeMillis();
		
		double[] temp;
		
		double[] transProb = new double[expDTMC.getNumStates()];
		double[] transProb2 = new double[expDTMC.getNumStates()];
		double[] resultProb = new double[expDTMC.getNumStates()];
		
		//setting initial distribution
		transProb[fDCTMCStateToDTMCState.get(fDIndexToEntry.get(fdIndex))]= 1.0;
		double exponential = Math.exp(-uniformizationRate * delta);
	
		
	
	//	
	//	System.out.println("transProb:" + transProb[0] + ", " + transProb[1] + ", " + transProb[2] + ", " +transProb[3]);
	//	System.out.println("resultProb:" + resultProb[0] + ", " + resultProb[1] + ", " + resultProb[2] + ", " +resultProb[3]);
		
		double overallReward = 0.0;
	
	
		
		//for each interval do
		for(int i=1;i * delta <= dMax;++i){
			double lambdaT = 1;
	
			
			//zero step
			//add to resultProb
			for(int k = 0; k< resultProb.length;k++)
				resultProb[k] = transProb[k] * exponential * lambdaT;        			
	
			for(int j=1;j<=numOfPoissonSteps;++j) {
				double oneStepReward = 0.0; 
	//			System.out.println("===transProb:" + transProb[0] + ", " + transProb[1] + ", " + transProb[2] + ", " +transProb[3]);
				lambdaT = lambdaT * (uniformizationRate * delta /j); 
				
				//first compute exp. transition rewards
				expTransRewDTMC.vmMult(transProb, transProb2);
				for(int k = 0; k< transProb.length;k++)
					oneStepReward = transProb2[k];
				overallReward = overallReward + oneStepReward * exponential * lambdaT;
				
				expDTMC.vmMult(transProb, transProb2);
				//add to resultProb
				for(int k = 0; k< resultProb.length;k++)
					resultProb[k] = resultProb[k] + (transProb2[k] * exponential * lambdaT);
	
				//swap vectors
				temp = transProb2;
				transProb2 = transProb;
				transProb = temp;
	
			}
	
			//filling rewards to MDP
			
			//compute fixed delay rewards
			double fixedDelayReward = 0.0;
			fDTransRewDTMC.vmMult(resultProb, transProb2);
			for(int k = 0; k< resultProb.length;k++)
				fixedDelayReward = fixedDelayReward + transProb2[k];
	
			mdpRewards.addToTransitionReward(fDCTMCStateToMDPState.get(entryToPotato), 
					i-1, fixedDelayReward + overallReward);
	
					
	//		System.out.println("transProb:" + transProb[0] + ", " + transProb[1] + ", " + transProb[2] + ", " +transProb[3]);
	//		System.out.println("resultProb:" + resultProb[0] + ", " + resultProb[1] + ", " + resultProb[2] + ", " +resultProb[3]);
			
			
			//fire fixed delays 
			fDDTMC.vmMult(resultProb, transProb2);
			
	//		System.out.println("transProb2:" + transProb2[0] + ", " + transProb2[1] + ", " + transProb2[2] + ", " +transProb2[3]);
			

			//filling new choice for MDP
			Distribution distr = new Distribution();
			
			//copy probability from non-potato states to mdp
			//We should ensure that distribution sums to exactly one. For the last state
			//that has nonzero probability we set its probability to 1 - sum of the rest.
			double sum =0.0;
			double notAddedValue = 0;
			int lastNonZeroIndex=-1;
			double prob = 0.0;
//			System.out.println("indexOfFirstNonPotatoState: "+ indexOfFirstNonPotatoState + " expDTMC.getNumStates(): " + expDTMC.getNumStates());
			for(int j=indexOfFirstNonPotatoState;j<expDTMC.getNumStates(); ++j) {
				prob = resultProb[j]+transProb2[j];
				if(prob > 0){
					if(lastNonZeroIndex != -1){
						distr.set(fDCTMCStateToMDPState.get(dTMCStatesToFDCTMCState.get(lastNonZeroIndex)), notAddedValue);
						sum = sum + notAddedValue;
					}
					lastNonZeroIndex = j;
					notAddedValue = prob;
				}
	
			}
			distr.set(fDCTMCStateToMDPState.get(dTMCStatesToFDCTMCState.get(lastNonZeroIndex)), 1-sum);
			
			

//			prob2 = resultProb[2]+transProb2[2];
////			prob2 = prob2 + resultProb[2]+transProb2[2];
//			distr.add(fDCTMCStateToMDPState.get(dTMCStatesToFDCTMCState.get(2)), 0);			
			
			
	//		System.out.println("CREATING DISTRIBUTION: " + distr);
			
			//add built distribution as a choice to mdp
			mdp.addChoice(fDCTMCStateToMDPState.get(fDIndexToEntry.get(fdIndex)), distr);
			
			//set resultProb as initial distribution for next interval
			temp = resultProb;
			resultProb = transProb;
			transProb = temp;
			
		}
		
		
		
		timer = System.currentTimeMillis()-timer;
		System.out.println("\nTransition probabilities and costs computed in " + (timer/1000.0) + " s.\n");
		
	} //done computation of this potato
}
