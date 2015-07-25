package explicit;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 *   
 * @author ado
 */
public class PositionKey implements Comparable<PositionKey>, Iterable<Integer> {

    protected List<Integer> tokens;
   
    public PositionKey(Integer size) {
        tokens = new ArrayList<Integer>(size);
        for(int i = 0; i<size; ++i) {
            tokens.add(0);
        }
    }
   
    public PositionKey(List<Integer> marking) {
        tokens = new ArrayList<Integer>(marking.size());
        for(int i = 0;i < marking.size();i++) {
            tokens.add(marking.get(i));
        }
    }
    
    public PositionKey(PositionKey marking) {
        tokens = new ArrayList<Integer>(marking.size());
        for(int i = 0;i < marking.size();i++) {
            tokens.add(marking.tokens.get(i));
        }
    }
    
    public void clear() {
        for(int i =0; i < tokens.size(); ++i) {
            tokens.set(i, 0);
        }
    }
    
    public void updateWaitingFD() {
        for(int i =0; i < tokens.size(); ++i) {
            if (tokens.get(i) > 1) {
                tokens.set(i, tokens.get(i) - 1);
            }
        }
    }
   
   
    public void fill(Integer index, Integer value) {
        tokens.set(index,value);
    }

    public Integer get(int index) {
        if(tokens.size()<=index)
            return 0;
        return tokens.get(index);
    }
    
    public void add(int value) {
        tokens.add(value);
    }
    
    public void set(int index, int value) {
        tokens.set(index, value);
    }
   
    public int size() {
        return tokens.size();
    }
    
    public List<Integer> getActiveFd() {
        List<Integer> result = new ArrayList<>();
        for(int i = 0; i < tokens.size(); i++) {
            if(tokens.get(i)==1)
                result.add(i);
        }
        return result;
    }
    
    public boolean onlyZero() {
        for(int pos : tokens) 
            if(pos!=0)
                return false;
        
        return true;
    }
    
    public int numOfOnes() {
        int result = 0;

        for (int i : tokens)
                if (i == 1)
                        ++result;

        return result;
    }
    
    public void updateTicks() {
        for(int i = 0; i < tokens.size(); i++)
            if(tokens.get(i)>1)
                tokens.set(i, tokens.get(i)-1);
                
    }
   
    @Override
    public int compareTo(PositionKey marking) {
        if(marking==null)
            throw new NullPointerException("PositionKey.compareTo: marking was null!");
        
        int keyLen = tokens.size();
        
        if(tokens.size() != marking.size()) {
            if(tokens.size() < marking.size()) {
                for(int i = tokens.size(); i < marking.size(); ++i) {
                    if(marking.get(i)!=0)
                        return 8;
                }
            } else {
                keyLen = marking.size();
                for(int i = marking.size(); i < tokens.size(); ++i) {
                    if(tokens.get(i)!=0)
                        return 8;
                }
            }
        }        
        
        int numOfOnesThis=0, numOfOnesMarking=0, sumThis=0, sumMarking=0,
                tmpThis=0, tmpMarking=0;
        for(int i=0; i<keyLen; ++i) {
            tmpThis = tokens.get(i);
            tmpMarking = marking.get(i);
            if(tmpThis==1)
                numOfOnesThis++;
            if(tmpMarking==1)
                numOfOnesMarking++;
            sumThis += tmpThis;
            sumMarking += tmpMarking;
        }
       
        if(sumThis==0 || sumMarking==0) {
            if(sumThis==sumMarking)
                return 0;
            if(sumThis==0)
                return -1;
            else
                return 1;
        }
        
        if(numOfOnesThis>numOfOnesMarking){
            return -2;
        }
        
        if(numOfOnesThis<numOfOnesMarking){
            return 2;
        }
        
        for(int i=0; i<keyLen; ++i) {
            if(tokens.get(i).compareTo(marking.get(i))!=0)
                return tokens.get(i).compareTo(marking.get(i));
        }
        
        return 0;
    }

    @Override
    public int hashCode() {
        return ((tokens == null) ? 0 : tokens.hashCode());
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj)
            return true;
        if (obj == null)
            return false;
        if (getClass() != obj.getClass())
            return false;
        PositionKey other = (PositionKey) obj;
        if (tokens == null) {
            if (other.tokens != null)
                return false;
        } else if (!tokens.equals(other.tokens))
            return false;
        return true;
    }

    @Override
    public String toString() {
        return "Key [" + tokens + "]";
    }

    @Override
    public Iterator<Integer> iterator() {
        return tokens.iterator();
    }
   
   
   
}
