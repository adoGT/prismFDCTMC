package explicit;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

/**
 * TODO vymysliet efektivne prechadzanie aktivnych fd
 * @author ado
 */
public class FDSubset {
    private BitSet activeFDs;
    private List<Integer> activeStates;
    
    public FDSubset() {
        activeFDs = new BitSet();
        activeStates = new ArrayList<Integer>();
    }
    
    public FDSubset(BitSet FDs, List<Integer> states) {
        activeFDs = FDs;
        activeStates = states;
    }
    
    /**
     * Copy constructor
     */
    public FDSubset(FDSubset fdSubst) {
        this.activeFDs = (BitSet) fdSubst.activeFDs.clone();
        for(int i : fdSubst.getActiveStates())
            this.activeStates.add(i);
    }
    
    public BitSet getActiveFDs() {
        return activeFDs;
    }
    
    public List<Integer> getActiveStates() {
        return activeStates;
    }

    @Override
    public String toString() {
        return "FDSubset{" + "activeFDs=" + activeFDs + ", activeStates=" + activeStates + '}';
    }
    
    
}