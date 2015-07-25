package explicit;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.HashSet;
import java.util.Iterator;

import prism.PrismException;
import prism.PrismLog;

public class CbpAlgQuickFix extends CbpAlg {
    
    protected FDMatrix<PositionKeyMinMax> transProb;
    protected FDMatrix<PositionKeyMinMax> transProbCopy;
    protected FDMatrix<PositionKeyMinMax> transProbResult;
    protected Distribution resultMin, resultMax;
    
    
    private static final int EMPTY_SET = -72;
    
    public CbpAlgQuickFix(PrismLog mainLog) throws PrismException {
        super(mainLog);
        transProb = new FDMatrix<>();
        transProbCopy = new FDMatrix<>();
        transProbResult = new FDMatrix<>();
    }
    
    protected void clear() {
    	transProb.clear();
    	transProbCopy.clear();
    	transProbResult.clear();
    }
    
    @Override
    public void compute(int numOfSteps, double t, int numOfPoissSteps)
			throws PrismException {
        intervalLength = t;
        initialize();

        
        for (int i = 1; i <= numOfSteps; ++i) 
        {
            executeStep(numOfPoissSteps);
            executeFixedD(transProb, fdctmc.getAllFDEventsIndexes()); 
            updateWaitingFD();
         //   System.out.println(",{ \"step\" : " + i + ", \"transProbMin\": " + getResult(true).toJSON() 
           //         + ", \"transProbMax\": " + getResult(false).toJSON() + " }");
        }
       // System.out.println("]}");

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
    
    public Distribution getResult() {
        return transProb.getDistr();
    }
    
    @Override
    protected void updateWaitingFD() {
    	FDMatrix<PositionKeyMinMax> result = new FDMatrix<>(transProb.size());

		for (Map.Entry<PositionKeyMinMax, Distribution> entry : transProb) {
                        //TODO it is possible to do smartly- only change the key/ possible extension with different structure of transProb
			PositionKeyMinMax position2 = new PositionKeyMinMax(entry.getKey()); 

			for (int i = 0; i < fdctmc.getNumFDEvents(); ++i) {
				if (entry.getKey().get(i) > 1) {
					position2.set(i, entry.getKey().get(i) - 1);
				}
			}
			result.put(position2, entry.getValue());
		}
		transProb = result;
	}
    
   
    
    @Override
    public void initialize() throws PrismException {
    	clear();
        Distribution pi = new Distribution(); 
        
        if (fdctmc.getNumInitialStates() == 1) {
                pi.add(fdctmc.getFirstInitialState(), 1);
        } else {
                throw new PrismException("Not one initial state: "
                                + fdctmc.getNumInitialStates() + ".");
        }
        
        Distribution distr = new Distribution(); 
        List<Integer> position = new ArrayList<>(fdctmc.getNumFDEvents());
        
        for (FDSubset mTEGroup : MTEgroups) {
                BitSet fdSubset = mTEGroup.getActiveFDs();
                // key in matrix, remaining time on fdelay
                position.clear();
                distr.clear();

                for (int i = 0; i < fdctmc.getNumFDEvents(); ++i) {
                        if (fdSubset.get(i)) {
                                position.add(numberOfSteps.get(fdctmc.getFDEvent(i)));
                        } else
                                position.add(0);
                }

                for (int i : mTEGroup.getActiveStates()) {
                    distr.set(i, pi.get(i));
                }
                
                transProb.put(new PositionKeyMinMax(position, true, PositionKeyMinMax.IN_TICK), new Distribution(distr));
                transProb.put(new PositionKeyMinMax(position, false, PositionKeyMinMax.IN_TICK), new Distribution(distr));
                
            }
    }
    
    private PositionKeyMinMax getPostition(FDSubset fdSubset, PositionKeyMinMax oldKey, int fd) {
        PositionKeyMinMax resultPosition = new PositionKeyMinMax(oldKey.size(), oldKey.isMin());
        resultPosition.nextGroup();
        for (int i = 0; i < fdctmc.getNumFDEvents(); ++i) {
            if (fdSubset.getActiveFDs().get(i)) {
                if( oldKey.inGroup(i) == PositionKeyMinMax.NOT_IN_GROUP ) {
                    resultPosition.setInLastGroup(i, numberOfSteps.get(fdctmc.getFDEvent(i))+1);
                } else {
                    resultPosition.set(i, 
                        oldKey.get(i), 
                        oldKey.inGroup(i));
                }
                          
                // TODO: should be in new group??
                if(fd==i && fd!=EMPTY_SET) {
                    resultPosition.setInLastGroup(i, numberOfSteps.get(fdctmc.getFDEvent(i))+1);
                }
            }
        }
        return resultPosition;
    }
    
    
    protected void executeStep(int numberOfPoissonSteps) {
        //mainLog.println("CbpAlgExt.executeStep(" + numberOfPoissonSteps +")", PrismLog.VL_ALL);
    	
        FDMatrix<PositionKeyMinMax> transProbRedundant = new FDMatrix<>(transProb.size());
        for (Iterator<Map.Entry<PositionKeyMinMax, Distribution>> it =  transProb.entrySet().iterator(); it.hasNext(); ) {
            Map.Entry<PositionKeyMinMax, Distribution> entry = it.next();
            
            double weight = Math.exp(-uniformizationRate * intervalLength);
            //mainLog.println("Position: " + entry.getKey() + "\nDistribution: " + entry.getValue() + "\nWeight: " + weight + "\ntoExplode: " + toExplodeBetweenTicks, PrismLog.VL_ALL);

        	transProbCopy.clear();
        	transProbCopy.put(new PositionKeyMinMax(entry.getKey()), new Distribution(entry.getValue()));
            
        	for(int i=0; true; ) {
            	//mainLog.println("Step(" + i + ")", PrismLog.VL_ALL);
            	transProbResult = (FDMatrix<PositionKeyMinMax>)FDMatrix.copyFDMatrixPositionKeyMinMax(transProbCopy);
            	
                for (Map.Entry<PositionKeyMinMax, Distribution> resEntry : transProbResult) {
                    transProbRedundant.add(resEntry.getKey(), resEntry.getValue().multiply(weight));
                }
        		i++;
                if(i>numberOfPoissonSteps) break;
        		weight = weight * uniformizationRate * intervalLength / i;
        		transProbCopy = executeExp(transProbCopy);
                //mainLog.println("Weight: " + weight + "\nTransProbRedundant: " + transProbRedundant + "\nEND Step(" + i +")", PrismLog.VL_ALL);
            }

            it.remove(); //TODO I believe that this saves some memory, if not, or it is too time consuming can be erased
            
          }

        transProb = transProbRedundant;
        //mainLog.println("END - CbpAlgExt.executeStep(" + numberOfPoissonSteps +")", PrismLog.VL_ALL);
    }
    

    private FDMatrix<PositionKeyMinMax> executeExp(FDMatrix<PositionKeyMinMax> map) {
        FDMatrix<PositionKeyMinMax> result = new FDMatrix<>(map.size());
        for (Map.Entry<PositionKeyMinMax, Distribution> entry : map) {
            Distribution distr = uniformisedDTMC.multiplyDistribution(entry.getValue());
            update(distr, result, entry.getKey(), EMPTY_SET, 1);
        }
        return result;
    }
    

    private Set<PositionKeyMinMax> update(Distribution distr, FDMatrix<PositionKeyMinMax> result, PositionKeyMinMax position, int fdSet, double weight) {
 		Set<PositionKeyMinMax> newKeysToProcess = new HashSet<PositionKeyMinMax>();
     	for(FDSubset fdSubset : MTEgroups) {
             PositionKeyMinMax newPosition = getPostition(fdSubset, position, fdSet);
             
 			Distribution distr2 = new Distribution();
 			
             for (int j : fdSubset.getActiveStates()) {
                 if (distr.get(j) == 0)
                     continue;
                 distr2.add(j, distr.get(j)*weight);
             }

 			if(distr2.isEmpty()) continue;
 			
 			result.getAndCreate(newPosition).add(distr2);
             
             if(newPosition.numOfOnes()>0) newKeysToProcess.add(newPosition);

         }
         return newKeysToProcess;
     }
     
     protected void executeFixedD(FDMatrix<PositionKeyMinMax> map, List<Integer> fd) {
         if (fdctmc.getNumFDEvents() == 0)
 			return;
         
         Set<PositionKeyMinMax> keys = new HashSet<>();
         for( PositionKeyMinMax key : map.keySet()) {
             if(key.numOfOnes()!=0) {
                 keys.add(new PositionKeyMinMax(key));
             }
         }

         Set<PositionKeyMinMax> newKeys = new HashSet<>();
         
         for (int i = fdctmc.getNumFDEvents(); i >= 1; --i) {
             for(Iterator<PositionKeyMinMax> it = keys.iterator(); it.hasNext(); ) {
                 PositionKeyMinMax key = it.next();
                 
                 if(key.numOfOnes()==i) {
                     List<Integer> activeFds = key.getActiveFd();
                     activeFds.retainAll(fd);
                     double overalWeight = getFdWeigth(activeFds);
                     for(int active : activeFds) {
                         double weight = fdctmc.getFDEvent(active).getWeight() / overalWeight;
                         Distribution distr = fdctmc.getFDEvent(active).multiplyDistribution(map.get(key));
                         newKeys.addAll(update(distr, map, key, active, weight));
                     }
                     
                     if(!activeFds.isEmpty()) {
                         map.remove(key);
                         it.remove();
                     }
                 }
             }
             keys.addAll(newKeys);
             newKeys.clear();
         }
     }
    
    
    public Distribution getResult(boolean min) {
    	if (min)
    		return resultMin;
		else
			return resultMax;
    }
    	
    public void computeResults() {
    	resultMin = new Distribution();
    	resultMax = new Distribution();
        for (Map.Entry<PositionKeyMinMax, Distribution> entry : transProb) {
            if(entry.getKey().isMin()==true)
                resultMin.add(entry.getValue());
            else
            	resultMax.add(entry.getValue());
        }
    }
    
    
}
