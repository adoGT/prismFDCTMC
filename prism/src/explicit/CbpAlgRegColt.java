package explicit;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

import prism.PrismException;
import prism.PrismLog;
import cern.colt.function.DoubleDoubleFunction;
import cern.colt.matrix.DoubleFactory2D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.doublealgo.Formatter;

public class CbpAlgRegColt {
	protected PrismLog mainLog;
	protected FDCTMCSimple fdctmc;
	protected BitSet regStates;
	protected double discretizationStep, error, time;
	protected double inputDiscretizationStep;
	protected int inputNumberOfIntervals;
	protected int numOfPoissSteps, numberOfIntervals;
	protected List<DoubleMatrix2D> kernelLocalMin, kernelLocalMax, kernelGlobalMin, kernelGlobalMax;
	protected List<ComputeKernelForRegStateColt> ckfrsList;
	
	public CbpAlgRegColt(PrismLog mainLog) throws PrismException {
		this.mainLog = mainLog;
		fdctmc = null;
		regStates = null;
		kernelLocalMin = new ArrayList<>();
		kernelLocalMax = new ArrayList<>();
		kernelGlobalMin = new ArrayList<>();
		kernelGlobalMax = new ArrayList<>();
		ckfrsList = new ArrayList<>();
	}
	
	public CbpAlgRegColt(CbpAlgRegColt cbp) {
		this.mainLog = cbp.mainLog;
		this.fdctmc = cbp.fdctmc;
		this.regStates = cbp.regStates;
		this.kernelLocalMin = cbp.kernelLocalMin;
		this.kernelLocalMax = cbp.kernelLocalMax;
		this.kernelGlobalMin = cbp.kernelGlobalMin;
		this.kernelGlobalMax = cbp.kernelGlobalMax;
	}
	
	protected void findRegStates() {
		regStates = new BitSet(fdctmc.getNumStates());
		regStates.set(0, fdctmc.getNumStates(), true);
		for (int i = 0; i < fdctmc.getNumFDEvents(); ++i) {
			regStates.andNot(fdctmc.getFDEvent(i).getActive());
		}
		System.out.println("RegStates: " + regStates);
	}
	
	public void setFdctmc(FDCTMCSimple fdctmc) {
		if (fdctmc == null)
			throw new IllegalArgumentException("CBPalg: fdctmc can't be null.");
		this.fdctmc = fdctmc;
		findRegStates();
	}
	
	protected boolean computeKernels() throws PrismException {
		int i = regStates.nextSetBit(0);
		while (i != -1) {
			if(!computeKernel(i))
				return false;
			i = regStates.nextSetBit(i + 1);
		}
		return true;
	}
	
	protected boolean computeKernel(int i) throws PrismException {
		ComputeKernelForRegStateColt ckfrs = new ComputeKernelForRegStateColt(mainLog, i, fdctmc, discretizationStep, numOfPoissSteps,
				numberOfIntervals, inputNumberOfIntervals, error, regStates);
		if(!ckfrs.computeKernel())
			return false;
		ckfrsList.add(i, ckfrs);
		return true;
	}
	
	protected void computeRest() {
		int maxStepsMax = 0, maxStepsMin = 0;
		for(ComputeKernelForRegStateColt ckfrs : ckfrsList) {
			if(ckfrs.kernelGlobalMax.size() > maxStepsMax)
				maxStepsMax = ckfrs.kernelGlobalMax.size();
			if(ckfrs.kernelGlobalMin.size() > maxStepsMin)
				maxStepsMin = ckfrs.kernelGlobalMin.size();
		}

		for(ComputeKernelForRegStateColt ckfrs : ckfrsList) {
			ckfrs.computeRest(maxStepsMin, maxStepsMax);
		}
		
		kernelGlobalMax.add(DoubleFactory2D.sparse.identity(fdctmc.numStates));
		kernelLocalMax.add(DoubleFactory2D.sparse.identity(fdctmc.numStates));
		kernelGlobalMin.add(DoubleFactory2D.sparse.identity(fdctmc.numStates));
		kernelLocalMin.add(DoubleFactory2D.sparse.identity(fdctmc.numStates));
		
		for(int i=1; i<=maxStepsMax; ++i){
			kernelLocalMax.add(DoubleFactory2D.sparse.make(fdctmc.numStates, fdctmc.numStates, 0));
			kernelGlobalMax.add(DoubleFactory2D.sparse.make(fdctmc.numStates, fdctmc.numStates, 0));
			for(ComputeKernelForRegStateColt ckfrs : ckfrsList) {
				kernelLocalMax.get(i).viewRow(ckfrs.regState).assign(ckfrs.kernelLocalMax.get(i-1));
				kernelGlobalMax.get(i).viewRow(ckfrs.regState).assign(ckfrs.kernelGlobalMax.get(i-1));
			}
		}
		
		for(int i=1; i<=maxStepsMin; ++i){
			kernelLocalMin.add(DoubleFactory2D.sparse.make(fdctmc.numStates, fdctmc.numStates, 0));
			kernelGlobalMin.add(DoubleFactory2D.sparse.make(fdctmc.numStates, fdctmc.numStates, 0));
			for(ComputeKernelForRegStateColt ckfrs : ckfrsList) {
				kernelLocalMin.get(i).viewRow(ckfrs.regState).assign(ckfrs.kernelLocalMin.get(i-1));
				kernelGlobalMin.get(i).viewRow(ckfrs.regState).assign(ckfrs.kernelGlobalMin.get(i-1));
			}
		}
			
	}
	
	protected void computeParams() throws PrismException {
		this.numberOfIntervals = (int) Math.round( time/this.discretizationStep);
		this.numOfPoissSteps = CbpAlg.getNumberOfPoissonStepsFromError(discretizationStep,
                numberOfIntervals, (error * 0.1), //this is to devote some portion of error to poisson error
                fdctmc.getMaxExitRate());		
	}
	
	public void runCBP(double time, double error, double discretizationStep) throws PrismException {
		this.time = time;
		this.error = error;
		this.discretizationStep = discretizationStep;
		this.inputDiscretizationStep = discretizationStep;
		this.inputNumberOfIntervals = (int) Math.round( time/this.discretizationStep);
		clear();
		computeParams();
		while(!computeKernels()) {
			clear();
			this.discretizationStep /= 2;
			computeParams();
		}
		computeRest();
	}
	
	protected void clear() {
		kernelGlobalMax.clear();
		kernelGlobalMin.clear();
		kernelLocalMax.clear();
		kernelLocalMin.clear();
		ckfrsList.clear();
	}
	
	public void resultInStep() {
		/*System.out.println("Kernel sizes: \nkernelGlobalMin:" + kernelGlobalMin.size()
				+ "\nkernelLocalMin:" + kernelLocalMin.size() + "\nkernelGlovalMax:"
				+ kernelGlobalMax.size() + "\nkernelLocalMax: " + kernelLocalMax.size());*/
		/*System.out.println("Global: ");
		for(DoubleMatrix2D m : kernelGlobalMax) {
			System.out.println(m);
		}System.out.println("Local: ");
		for(DoubleMatrix2D m : kernelLocalMax) {
			System.out.println(m);
		}*/
//		List<DoubleMatrix2D> stepsMax = new ArrayList<>();
		List<DoubleMatrix2D> stepsMin = new ArrayList<>();
//		stepsMax.add(DoubleFactory2D.sparse.identity(fdctmc.numStates));
		stepsMin.add(DoubleFactory2D.sparse.identity(fdctmc.numStates));
//		double oldError =0;
		for(int i=1; i<=inputNumberOfIntervals; ++i) {

			//System.out.println("step " + i + " : ");
//			DoubleMatrix2D tmpMax = DoubleFactory2D.sparse.make(fdctmc.numStates, fdctmc.numStates,0);
			DoubleMatrix2D tmpMin = DoubleFactory2D.sparse.make(fdctmc.numStates, fdctmc.numStates,0);
//			for(int j=1; j<=i && j<kernelGlobalMax.size() ; ++j) {
				//System.out.println("Pos 0x0: " + kernelGlobalMax.get(j).get(0, 0));
//				tmpMax.assign(kernelGlobalMax.get(j).zMult(stepsMax.get(i-j), null), plus);
				//System.out.println(stepsMax.get(i-j) + "\n * \n " + kernelGlobalMax.get(j) + " = " + tmpMax);
//			}
			for(int j=1; j<=i && j<kernelGlobalMin.size() ; ++j) {
				tmpMin.assign(kernelGlobalMin.get(j).zMult(stepsMin.get(i-j), null), plus);
				//System.out.println(stepsMin.get(i-j) + "\n * \n " + kernelGlobalMin.get(j) + " = " + tmpMin);
			}
//			if(i<kernelLocalMax.size()) {
//				tmpMax.assign(kernelLocalMax.get(i), plus);
//			}
			if(i<kernelLocalMin.size()) {
				tmpMin.assign(kernelLocalMin.get(i), plus);
			}
//			stepsMax.add(tmpMax);
			stepsMin.add(tmpMin);
			
			
//			System.out.println("\n" + stepsMin.get(i));
			
//			double sum =0;
//			for(int k=0; k< stepsMin.get(i).columns(); k++){
//				 sum += stepsMin.get(i).get(0,k);
//			}
//			System.out.println(i +": \terror: "+ (1-sum) + " \tdiff: " + ((1-sum) -oldError));
//			oldError = (1-sum);
			
			
		}
		System.out.println("Min " + inputNumberOfIntervals + " (sum= " + stepsMin.get(inputNumberOfIntervals).zSum() + "): " + stepsMin.get(inputNumberOfIntervals));
//		System.out.println("Max " + numberOfIntervals + " : " + stepsMax.get(numberOfIntervals));
	}
	
	DoubleDoubleFunction plus = new DoubleDoubleFunction() {
	    public double apply(double a, double b) { return a+b; }
	};
	
}
