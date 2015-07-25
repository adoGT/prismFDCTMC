package explicit;

import java.util.ArrayList;
import java.util.List;

import prism.PrismException;

/**
 *   TODO: zrecyklovat index pre group
 * @author ado
 */
public class PositionKeyMinMax extends PositionKey {

    public static final int IN_TICK = -1;
    public static final int NOT_IN_GROUP = -45;
    public static final int NEXT_GROUP = -72;
    
    private boolean min;
    
    private List<Integer> groups;
    private int lastGroup;
    
    public PositionKeyMinMax(Integer size, boolean min) {
        super(size);
        this.min=min;
        lastGroup = NEXT_GROUP;
        groups = new ArrayList<>(size);
        for(int i = 0; i<size; i++)
            groups.add(NOT_IN_GROUP);
    }
    
    public PositionKeyMinMax(Integer size) {
        super(size);
        groups = new ArrayList<>(size);
        for(int i = 0; i<size; i++)
            groups.add(NOT_IN_GROUP);
        min=false;
        lastGroup = NEXT_GROUP;
    }
    
    public PositionKeyMinMax(List<Integer> tokens, boolean min, int group) {
        super(tokens);
        groups = new ArrayList<>(tokens.size());
        for(int i = 0; i<tokens.size(); i++) {
            if(tokens.get(i)==0) {
                groups.add(NOT_IN_GROUP);
            } else {
                groups.add(group);
            }
        }
        this.min = min;
        lastGroup = NEXT_GROUP;
    }
   
    
    public PositionKeyMinMax(PositionKeyMinMax marking) {
        super(marking);
        min = marking.min;
        groups = new ArrayList<>(marking.size());
        for(int i = 0;i < marking.groups.size();i++) {
            groups.add(marking.groups.get(i));
        }
        lastGroup = marking.lastGroup;        
    }
    
    public void nextGroup() {
        lastGroup = NEXT_GROUP;
    }
    
    public int inGroup(int index) {
        if(tokens.size()<=index)
            return NOT_IN_GROUP;
        if(tokens.get(index)==0)
            return NOT_IN_GROUP;
        return groups.get(index);
    }
    
    public List<Integer> getFDInGroup(int group){
        List<Integer> result = new ArrayList<>();
        for(int i = 0; i<groups.size(); i++) {
            if(groups.get(i)==group) 
                result.add(i);
        }
        return result;
    }
    
    /**
     * Return groups of delays, in which at least one delay is going to explode
     * in next tick
     * @param inTick
     * @return 
     */
    public List<Integer> toExplode(boolean inTick) {
        List<Integer> result = new ArrayList<>();
        for(int i=0; i<groups.size(); i++) {
            if(!inTick && groups.get(i)==IN_TICK)
                continue;
            if(!result.contains(groups.get(i)) && tokens.get(i)==1) {
                result.add(groups.get(i));
            }
            
        }
        return result;
    }
    
    public boolean isMin() {
        return min;
    }
    
    @Override
    public void add(int value) {
        super.add(value);
        if(lastGroup==NEXT_GROUP) {
            lastGroup = tokens.size()-1;
        }
        groups.add((value==0) ? NOT_IN_GROUP : lastGroup);
    }
    
    @Override
    public void set(int index, int value) {
        super.set(index, value);
        if(value==0) {
            changeGroup(index, NOT_IN_GROUP);
        }
    }
    
    public int set(int index, int value, int group) {
        super.set(index, value);
        if(value==0) {
            return changeGroup(index, NOT_IN_GROUP);
        }
        return changeGroup(index, group);
    }
    
    public void setInLastGroup(int index, int value) {
        super.set(index, value);
        if(value==0) {
            changeGroup(index, NOT_IN_GROUP);
            return;
        }
        if(lastGroup==NEXT_GROUP) {
            lastGroup = index;
            removeIndex(index);
        }
        
        lastGroup = changeGroup(index, lastGroup);
    }
    
    public int changeGroup(int index, int group) {
        if(groups.get(index)==group) {
            return group;
        }
        
        // the group of index will change
        removeIndex(index);
        
        if(group==NOT_IN_GROUP) 
            return NOT_IN_GROUP;
        
        if(groups.contains(group)) {
            if(index<group) {
                groups.set(index, index);
                for(int i=group; i < groups.size(); ++i) {
                    if(groups.get(i)==group)
                        groups.set(i, index);
                }
                return index;
            } else {
                groups.set(index, group);
                return group;
            }
        } else {
            groups.set(index, index);
            return index;
        }
        
    }
    
    private void removeIndex(int index) {
        groups.set(index, NOT_IN_GROUP);
        if(groups.contains(index)) {
            int newIndex = NEXT_GROUP;
            for(int i = index+1; i<groups.size(); ++i) {
                if(groups.get(i)==index) {
                    if(newIndex==NEXT_GROUP)
                        newIndex=i;
                    groups.set(i, newIndex);
                }
            }
        }
    }
    
    @Override
    public int compareTo(PositionKey mark) {
        if(mark==null)
            throw  new NullPointerException("PositionKeyMinMax.compareTo: mark is null!");
        
        PositionKeyMinMax marking = (PositionKeyMinMax) mark;
        if(marking.min != min)
            if(min)
                return -3;
            else
                return 3;
        
        if(super.compareTo(mark)!=0)
            return super.compareTo(mark);
        
        if(marking.groups.size()!=groups.size())
            return 5;
        
        for(int i=0; i<groups.size(); ++i) {
            if(groups.get(i).compareTo(marking.groups.get(i))!=0)
                return groups.get(i).compareTo(marking.groups.get(i));
        }
        
        return 0;
    }

	@Override
	public int hashCode() {
		final int prime = 29;
		int result = super.hashCode();
		result = prime * result + ((groups == null) ? 0 : groups.hashCode());
		result = prime * result + (min ? 1231 : 1237);
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (!super.equals(obj))
			return false;
		if (getClass() != obj.getClass())
			return false;
		PositionKeyMinMax other = (PositionKeyMinMax) obj;
		if (groups == null) {
			if (other.groups != null) {
				return false;
                        }
		} else if (!groups.equals(other.groups)) {
			return false;
                }
		if (min != other.min) {
			return false;
                }
		return true;
	}

    @Override
    public String toString() {
        return "M[" + min +"] T[" + tokens + "] G[" + groups + "]";
    }
    
}
