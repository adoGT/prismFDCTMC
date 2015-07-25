package explicit.rewards;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import prism.PrismException;

public class FDCTMCRewardsSimple extends StateRewardsArray implements MCRewards {

	protected Map<Integer,Map<Integer,Map<Integer,Double>>> transReward;
	
	public FDCTMCRewardsSimple(int numStates) {
		super(numStates);
		transReward = new HashMap<Integer,Map<Integer,Map<Integer,Double>>> ();
	}
	

	/**
	 * adds transition reward
	 * @param fixedDelay
	 * @param src
	 * @param dest
	 * @param rew
	 */
	public void addTransitionReward(int fixedDelay, int src, int dest, double rew){
		if(transReward.get(fixedDelay) == null){
			transReward.put(fixedDelay,new HashMap<Integer,Map<Integer,Double>>());
		}
			
		
		Map<Integer, Double> map;
		if(transReward.get(fixedDelay).get(src) == null) {
			map = new HashMap<Integer, Double>(1);
			map.put(dest, rew);
			transReward.get(fixedDelay).put(src, map);
			return;
		}
		
		map = transReward.get(fixedDelay).get(src);
		if(map.get(dest) == null){
			map.put(dest,rew);
			return;
		}

		map.put(dest,rew + map.get(dest));

	}
	

	/**
	 * @param fixedDelay
	 * @param src
	 * @param dest
	 * @return transition reward
	 */
	public double getTransitionReward(int fixedDelay, int src, int dest)
	{
		if(transReward.get(fixedDelay) != null 
				&& transReward.get(fixedDelay).get(src) != null
				&& transReward.get(fixedDelay).get(src).get(dest) !=null)
			return transReward.get(fixedDelay).get(src).get(dest);
		else
			return 0.0;
	}


	@Override
	public String toString() {
		return "FDCTMCRewardsSimple [transReward=" + transReward
				+ ", stateRewards=" + Arrays.toString(stateRewards) + "]";
	}
	
	
	
}
