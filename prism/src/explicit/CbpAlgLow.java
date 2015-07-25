package explicit;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import prism.PrismException;
import prism.PrismLog;

public class CbpAlgLow extends CbpAlg {

//	public CbpAlgLow() {
//		throw new NotImplementedException();
//	}
	
	public CbpAlgLow(PrismLog mainLog) throws PrismException {
        super(mainLog);
	}

	
	public void update(PositionKey position, Distribution distr,
			boolean restart, double weight) {
		for (FDSubset mTEGroup : MTEgroups) {
			List<Integer> position2 = new ArrayList<>(fdctmc.getNumFDEvents());
			for (int i = 0; i < fdctmc.getNumFDEvents(); ++i)
				position2.add(0);

			for (int i = 0; i < mTEGroup.getActiveFDs().size(); ++i)
				if (mTEGroup.getActiveFDs().get(i))
					if (position.get(i) > 0) {
						position2.set(i, position.get(i));
					} else {
						position2.set(i,
								numberOfSteps.get(fdctmc.getFDEvent(i)));
					}
			Distribution distr1 = transProbResult.get(new PositionKey(position2));
			Distribution distr2 = transProb.get(new PositionKey(position2));
			for (int j : mTEGroup.getActiveStates()) {
				if (distr.get(j) == 0)
					continue;
				distr1.add(j, distr.get(j) * weight);
				distr2.add(j, distr.get(j));
			}
            transProbResult.put(new PositionKey(position2), distr1);
			transProb.put(new PositionKey(position2), distr2);
		}
	}
	
}
