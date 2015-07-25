package explicit;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import cern.colt.matrix.DoubleFactory1D;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import prism.PrismException;
import prism.PrismLog;

/**
 * Class for partial kernel computation. Partial kernel holds values for steps from 1.
 * Values for step 0 are not included (Kernel in 0 is identity and is added in different
 * computation part).
 * 
 */
public class ComputeKernelForRegStateColt {
	protected int regState, numOfPoissSteps, numberOfIntervals, inputNumberOfIntervals, numberOfStepsToComputeKernel;
	protected List<DoubleMatrix1D> kernelLocalMin, kernelLocalMax, kernelGlobalMin, kernelGlobalMax;
	protected CbpAlgExt cbpAlg;
	protected double intervalLength, error, poissError;
	protected BitSet regStates;
	protected int newState;
	protected boolean continueComputeMin, continueComputeMax;
	protected double deletedProbMin, sumTotalMin, sumTotalNoRegMin, deletedProbabilityMax;
	protected int step;
	

	
	public ComputeKernelForRegStateColt(PrismLog mainLog, int regState, FDCTMCSimple fdctmc, double t,
			int numOfPoissSteps, int numberOfIntervals, int inputNumberOfIntervals, double error, BitSet regStates) throws PrismException {
		this.regState = regState;
		this.regStates = regStates;
		this.intervalLength = t;
		this.numOfPoissSteps = numOfPoissSteps;
		this.numberOfIntervals = numberOfIntervals;
		this.inputNumberOfIntervals = inputNumberOfIntervals;
		this.error = error;
		kernelLocalMin = new ArrayList<>();
		kernelLocalMax = new ArrayList<>();
		kernelGlobalMin = new ArrayList<>();
		kernelGlobalMax = new ArrayList<>();
		cbpAlg = new CbpAlgExt(mainLog);
		cbpAlg.setInterval(t);
		System.out.println("FDCTMC original: " + fdctmc);
		cbpAlg.setFdctmc(prepareFDCTMC(fdctmc));
		//fixUniformisedDTMC();
		System.out.println("FDCTMC: " + cbpAlg.fdctmc);
		cbpAlg.initialize();
		continueComputeMin = true;
		continueComputeMax = true;
	}

	protected FDCTMCSimple prepareFDCTMC(FDCTMCSimple fdctmc) {
		FDCTMCSimple result = new FDCTMCSimple(fdctmc);

		result.removeAllInitialStates();
		result.addInitialState(regState);
		return result;
	}
	
	protected void fixUniformisedDTMC() {
		cbpAlg.uniformisedDTMC.getTransitions(newState).set(regState, cbpAlg.uniformisedDTMC.getTransitions(newState).get(newState));
		cbpAlg.uniformisedDTMC.getTransitions(newState).set(newState, 0);
	}
	
	protected void distrToMatrix(Distribution distr, boolean min) {
//		System.out.println("STEP: " + step + " DISTR " + distr);
		DoubleMatrix1D local = DoubleFactory1D.sparse.make(cbpAlg.fdctmc.numStates, 0);
		DoubleMatrix1D global = DoubleFactory1D.sparse.make(cbpAlg.fdctmc.numStates, 0);
		for(Map.Entry<Integer,Double> entry : distr) {
				if(regStates.get(entry.getKey())) {
					global.set(entry.getKey(), entry.getValue());
				} else {
					local.set(entry.getKey(), entry.getValue());
				}
		}
		
		if(min) {
			kernelGlobalMin.add(global);
			kernelLocalMin.add(local);
		} else {
//			kernelGlobalMax.add(global);
//			kernelLocalMax.add(local);
		}
	}
	
	protected boolean checkConditions(){
		// min
		double probSumMin = cbpAlg.resultMin.sum();

//		System.out.println("XXXXXXXX Step = " + step + " Interval length: "+ intervalLength 
//				+ " Num of Intervals: " + numberOfIntervals + " One step error: " + ((1-deletedProbMin-probSumMin)/step) 
//				+ " Error: " + (1-deletedProbMin - probSumMin) + " Poisson error: " + poissError);
		
		// if we need to try smaller discretizationStep?
//		System.out.println("Min nutna: (1-(" + deletedProbabilityMin + "-" + probSumMin + "))/" + step +")" + "<=" + "(" + error + "/" + numberOfIntervals + ")");
//		System.out.println("Min nutna:" + ((1-deletedProbabilityMin-probSumMin)/step) + "<=" + (error/numberOfIntervals));

//		if((1-deletedProbabilityMin - probSumMin)/step>(error/numberOfIntervals))
		if((1-deletedProbMin - probSumMin)>error)
				return false;		

		if((1-deletedProbMin - probSumMin) < 0 )
			throw new IllegalArgumentException("Sum of probability is greater thatn 1: " + (deletedProbMin + probSumMin));
			
		// if we can stop
//		System.out.println("Min zastavujuca: ((1-" + sumMin + ")/" + step + ")" + ">?" + "(" + error + "/" + numberOfIntervals + ")");
//		System.out.println("Min zastavujuca: " + ((1-sumMin)/step) + ">?" + (error/numberOfIntervals));
//		if(((1-sumMin)/step)<=(error/numberOfIntervals)){
//		if(((1-deletedProbMin)*(inputNumberOfIntervals- (step / numberOfStepsToComputeKernel)))<=error){
//		if(((1-deletedProbMin))<=error){
		if(step >= numberOfIntervals){
			//System.out.println("STOP");
			continueComputeMin = false;	
			continueComputeMax = false;
		}
		
		//max
		
		return true;
	}
	
	/**
	 * Process step. If addKernelStep == true, add result to partial kernels.
	 */
	protected void processStep(boolean addKernelStep) { 
		if(addKernelStep){
			distrToMatrix(cbpAlg.resultMin, true);
			deletedProbMin += cbpAlg.transProb.get(new PositionKeyMinMax(cbpAlg.fdctmc.getNumFDEvents(), true)).sum();
	        cbpAlg.transProb.put(new PositionKeyMinMax(cbpAlg.fdctmc.getNumFDEvents(), true), new Distribution());

//	        System.out.println("DeletedProb: " + deletedProbMin + "\tsum of transProb: "+ cbpAlg.resultMin.sum());

	        
//			distrToMatrix(cbpAlg.resultMax, false);
//	        cbpAlg.transProb.put(new PositionKeyMinMax(cbpAlg.fdctmc.getNumFDEvents(), false), new Distribution());
		}

	}
	
	public boolean computeKernel() throws PrismException {
		poissError = CbpAlg.getPoissonError(intervalLength, numberOfIntervals, numOfPoissSteps, cbpAlg.fdctmc.getMaxExitRate());
		
		System.out.println("########## PoissErro: " + poissError + " = " 
		 + "\nIntervalLength: " +  intervalLength + "\nNumberOfIntervals: " + numberOfIntervals 
		 + "\nNumOfPoissSteps: " + numOfPoissSteps + "\nRate: " + cbpAlg.fdctmc.getMaxExitRate());
		deletedProbMin = 0; sumTotalMin=0; sumTotalNoRegMin=0; deletedProbabilityMax = 0; step = 1;
		numberOfStepsToComputeKernel = numberOfIntervals/inputNumberOfIntervals; 
		do {
			cbpAlg.executeStep(numOfPoissSteps);
			cbpAlg.executeFixedD(cbpAlg.transProb, cbpAlg.fdctmc.getAllFDEventsIndexes());
			cbpAlg.updateWaitingFD();
			cbpAlg.computeResults();
			if(!checkConditions()) {
				return false;
			}
			processStep((step % numberOfStepsToComputeKernel) ==0); 
			++step;
        } while(continueComputeMin || continueComputeMax);
		return true;
	}
	
	public void computeRest(int toMin, int toMax) {
		
	}
	
}
