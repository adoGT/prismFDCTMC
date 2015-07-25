package explicit;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;


public class DTMCDeadlineSimple extends DTMCSimple{

	
	// deadline expiration transition matrix (distribution list) 
	protected List<Distribution> deadlineExpirationToActive;
	protected List<Distribution> deadlineExpirationToNonActive;
	protected int numActiveDeadlineStates = 0;
	protected int permut[];
	
	
	public DTMCDeadlineSimple(int numStates)
	{
		super(numStates);
		
		numActiveDeadlineStates = 0;
		deadlineExpirationToActive = null;
		deadlineExpirationToNonActive = null;
	}
	
	public DTMCDeadlineSimple(DTMCDeadlineSimple dtmc, List<Integer> deadlineTrans)
	{
		this(dtmc.numStates);
		
		//compute permutation
		permut = new int[numStates];
		int drift = 0;
		for (int i= 0; i< deadlineTrans.size();i++)
		{
			if(deadlineTrans.get(i) < 0)
			{
				drift++;
				permut[i] = numStates - drift;
			}
			else
			{
				permut[i] = i - drift;
			}
		}

		for(int i = deadlineTrans.size(); i<numStates;i++)
		{
			permut[i] = i - drift;
		}

		numActiveDeadlineStates = deadlineTrans.size()-drift;
		
		//use permutation to rearrange states so that the active deadline states will be in front
		copyFrom(dtmc, permut);
		for (int i = 0; i < numStates; i++) {
			trans.set(permut[i], new Distribution(dtmc.trans.get(i), permut));
		}
		numTransitions = dtmc.numTransitions;
		

		//initialization of deadlineExpiration transition matrix
		//actually we can imagine deadline expiration as doing probabilistic step when we reached 
		//deadline and then instantaneously changing state according to deadlineTrans		
		deadlineExpirationToActive = new ArrayList<Distribution>(); // for active states		
		deadlineExpirationToNonActive = new ArrayList<Distribution>(); // for active states		

		
		//economical storage of transposed deadlineTranse which are already rearranged by permut
		//this is needed to perform computation of deadlineExpiration quickly
		List<Set<Integer>> deadlineTransBackward = new ArrayList< Set<Integer>>();		


		//Initialization of the list
		for(int i=0;i<numStates;i++)
		{
			deadlineTransBackward.add(null);
		}
		
		int j =0;		
		for(int i=0;i<deadlineTrans.size();i++)
		{
			j = deadlineTrans.get(i);
			if(j >= 0)
			{
				if(deadlineTransBackward.get(permut[i]) == null)
				{
					deadlineTransBackward.add(permut[i], new HashSet<Integer>());
				}
				deadlineTransBackward.get(permut[i]).add(permut[j]);
			}
		}
	
		Distribution dist, newDist;
		double temp;
		for(int i =0; i<numActiveDeadlineStates; i++)
		{
			 dist = trans.get(i);
			 newDist = new Distribution();
			 for(j =0;j<numActiveDeadlineStates;j++)
			 {
				 temp=0;
				 for(Integer index : deadlineTransBackward.get(j))
				 {
					 temp += dist.get(index);
				 }
				 newDist.add(j, temp);
			 }
			 deadlineExpirationToActive.add(i, newDist);
		}


		
		for(int i= numActiveDeadlineStates; i< numStates; i++)
		{
			 dist = trans.get(i);
			 newDist = new Distribution();
			 for(j =0;j<numActiveDeadlineStates;j++)
			 {
				 //adding P2, resp. probability of changing from active to non-active deadline state
				 temp= trans.get(i-numActiveDeadlineStates).get(j); 
				 
				 for(Integer index : deadlineTransBackward.get(j))
				 {
					 temp += dist.get(index);
				 }
				 newDist.add(j, temp);
			 }
			 deadlineExpirationToNonActive.add(i-numActiveDeadlineStates, newDist);
		}

	}
	
	public int[] getPermut()
	{
		return permut;
	}
	
	public int getNumActiveDeadlineStates()
	{
		return numActiveDeadlineStates;
	}
}
