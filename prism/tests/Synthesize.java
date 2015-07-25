import java.io.*; 
import java.util.Map.*;
import java.util.*;

import parser.*;
import parser.ast.*;
import prism.*;
import simulator.*;
import explicit.*;
import explicit.rewards.*;


public class Synthesize {
  

    public static void main(String[] args)
	{
    	
		try {
                    
            if (args.length != 2){
            	System.out.println("usage: java Synchronize <FDCTMC.pm> <error>");
            	return;
            }


            double error=0;
            error = Double.parseDouble(args[1]);
            if(error>1.0 || error <= 0)
            	throw new PrismException("Error " + args[1] + " is too high or too low.");
            
            
            // Simple example: parse a PRISM file from a file, construct the model and export to a .tra file
			PrismLog mainLog = new PrismPrintStreamLog(System.out);
            //mainLog.setVerbosityLevel(PrismLog.VL_ALL);
            Prism prism = new Prism(mainLog, mainLog);
			
            ModulesFile modulesFile = prism.parseModelFile(new File(args[0]));
			
            UndefinedConstants undefinedConstants = new UndefinedConstants(modulesFile, null);
              
            modulesFile.setUndefinedConstants(undefinedConstants.getMFConstantValues());
			    
//            System.out.println("MODULES FILE: " + modulesFile);
            
            
            ConstructModel constructModel = new ConstructModel(prism.getSimulator(), mainLog);

//            undefinedConstants.defineUsingConstSwitch(args[2]);

            FDCTMCSimple model = (FDCTMCSimple) constructModel.constructModel(modulesFile);
            
           
            boolean foundTarget= false;
            Expression targetExpression = null;
            for(int i =0; i< modulesFile.getLabelList().size();i++){
            	if(modulesFile.getLabelList().getLabelName(i).equals("target")){
            		foundTarget= true;
            		targetExpression = modulesFile.getLabelList().getLabel(i);
            	}	
            }

            if(!foundTarget)
            	throw new PrismException("Model does not contain \"target\" label specigying target states!");
            
            
    		List<Integer> targetStatesFDCTMC = new ArrayList<Integer>();
    		
    		List<State> statesList = model.getStatesList();
			for (int i = 0; i < model.getNumStates(); i++) {
				if (targetExpression.evaluateBoolean(modulesFile.getConstantValues(), statesList.get(i))) {
					targetStatesFDCTMC.add(i);
				}
			}
            
		
            System.out.println("PARSED MODEL: " + model); //TODO preco ma FD udalost ako label waiting cost
			System.out.println("TARGET STATES: " + targetStatesFDCTMC);

            
            
    		// Get reward info
            ConstructRewards constructRewards = new ConstructRewards(mainLog);
            
            FDCTMCRewardsSimple fDCTMCRewards = null;
            RewardStruct rewStruct = null; // Reward struct object
    		if (modulesFile == null)
    			throw new PrismException("No model file to obtain reward structures");
    		if (modulesFile.getNumRewardStructs() == 0)
    			throw new PrismException("Model has no rewards specified");
    		if (modulesFile.getNumRewardStructs() > 1) {
    			throw new PrismException("Model has multiple reward structures specified, please specify only one.");
    		}
    			rewStruct = modulesFile.getRewardStruct(0);
    		if (rewStruct == null)
    			throw new PrismException("Invalid reward structure index \"" + 0 + "\"");
            
    		fDCTMCRewards = constructRewards.buildFDCTMCRewardStructure((FDCTMC) model, rewStruct, modulesFile.getConstantValues());
    		
    		model.clearSynchLabels();
    		
    		//check whether every non-target state has positive state reward
    		for(int i=0;i<model.getNumStates();i++){
    			if(targetStatesFDCTMC.contains(i)) continue;
    			if(fDCTMCRewards.getStateReward(i)<= 0)
    				throw new PrismException("State " + i + " das not have positive state reward: " + fDCTMCRewards.getStateReward(i));
    		}
    		
    		System.out.println(fDCTMCRewards.toString());
            
            FDCTMCSynthesis synthesis = new FDCTMCSynthesis();
            
            synthesis.computeSynthesis(model, fDCTMCRewards, targetStatesFDCTMC, error);
            
            
        } catch (FileNotFoundException e) {
			System.out.println("Error: " + e.getMessage());
			System.exit(1);
		} catch (PrismException e) {
			System.out.println("Error: " + e.getMessage());
			System.exit(1);
		}
	}
    

}
