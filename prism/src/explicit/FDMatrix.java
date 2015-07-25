package explicit;

import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;
import java.util.Map;
import java.util.HashMap;

public class FDMatrix<T> implements Iterable<Map.Entry< T, Distribution>> {
     private HashMap< T, Distribution> matrix;
     
     public FDMatrix() {
         matrix = new HashMap<>();
     }
  
     public FDMatrix(int size) {
         matrix = new HashMap<>();
     }
     
     /**
      * Copy constructor
      */
    /* public FDMatrix(FDMatrix<T> fdMatrix) {
         matrix = new HashMap<>();
         for(Map.Entry<T, Distribution> entry : fdMatrix)
             this.put(new PositionKey(entry.getKey()), new Distribution(entry.getValue()));
     }*/
     
     public static FDMatrix copyFDMatrixPositionKey(FDMatrix<PositionKey> fdMatrix) {
         FDMatrix<PositionKey> copy = new FDMatrix<>(fdMatrix.size());
         for(Map.Entry<PositionKey, Distribution> entry : fdMatrix)
             copy.put(new PositionKey(entry.getKey()), new Distribution(entry.getValue()));
         
         return copy;
     }
     
     public static FDMatrix copyFDMatrixPositionKeyMinMax(FDMatrix<PositionKeyMinMax> fdMatrix) {
         FDMatrix<PositionKeyMinMax> copy = new FDMatrix<>(fdMatrix.size());
         for(Map.Entry<PositionKeyMinMax, Distribution> entry : fdMatrix)
             copy.put(new PositionKeyMinMax(entry.getKey()), new Distribution(entry.getValue()));
         
         return copy;
     }
     
     public static void minMax(FDMatrix<PositionKeyMinMax> first, FDMatrix<PositionKeyMinMax> second) {
         for (Iterator<Map.Entry<PositionKeyMinMax, Distribution>> it =  first.entrySet().iterator(); it.hasNext(); ) {
            Map.Entry<PositionKeyMinMax, Distribution> entry = it.next();
            if(entry.getKey().isMin() && !second.existsDistr(entry.getKey()))
                it.remove();
         }
         
         for(Map.Entry<PositionKeyMinMax, Distribution> entry : second) {
             if(first.existsDistr(entry.getKey()))
            	 first.get(entry.getKey()).minMax(entry.getValue(), entry.getKey().isMin());
             else {
            	 if(!entry.getKey().isMin()){
                     first.put(entry.getKey(), entry.getValue()); //if this causes exception make new distribution
            	 }
             }
         }
     }
     
     public static boolean diffReg(double epsilon, FDMatrix<PositionKeyMinMax> first, FDMatrix<PositionKeyMinMax> second) {
         Set<PositionKeyMinMax> keys = new HashSet<>();
         for(PositionKeyMinMax key : first.keySet())
        	 keys.add(new PositionKeyMinMax(key));
         for(PositionKeyMinMax key : second.keySet())
        	 keys.add(new PositionKeyMinMax(key));
         for(PositionKeyMinMax key : keys) {
        	if( first.get(key).diff(second.get(key), epsilon))
        		return true;
         }
    	 return false;
     }
     
     public void multiply(double weight) {
         for (Iterator<Map.Entry<T, Distribution>> it =  matrix.entrySet().iterator(); it.hasNext(); ) {
             Map.Entry<T, Distribution> entry = it.next();
             	entry.getValue().multiplyThis(weight);
         }
     }
     
     public void add(FDMatrix<T> toAdd) {
         for(Map.Entry<T, Distribution> entry : toAdd) {
             add(entry.getKey(), entry.getValue());
         }
     }
     
     public void clear() {
         matrix.clear();
     }
     
     public int size() {
    	 return matrix.size();
     }
     
     public boolean isEmpty(){
    	 return matrix.isEmpty();
     }
     
     public void put(T position, Distribution distr) {
         if(!distr.isEmpty() && distr.sum() != 0)
            matrix.put(position, distr);
         else
             matrix.remove(position);
     }
     
     public Set<Map.Entry<T, Distribution>> entrySet() {
         return matrix.entrySet();
     }
     
     public Set<T> keySet() {
         return matrix.keySet();
     }
     
     public void remove(T position) {
         matrix.remove(position);
     }
     
     
     /**
      * See if we need to create always new Distribution
      * @param position
      * @return 
      */
     public Distribution get(T position) {  
         Distribution distr = matrix.get(position);
         if(distr==null){
             distr = new Distribution();
         }
         return distr;   
     }
     
     
     /**
      * @param position
      * @return 
      */
     public Distribution getAndCreate(T position) {  
         Distribution distr = matrix.get(position);
         if(distr==null){
             distr = new Distribution();
             matrix.put(position,distr);
         }
         return distr;   
     }
     
     public boolean existsDistr(T position) {
         Distribution distr = matrix.get(position);
         if(distr==null) {
             return false;
         }
         return true;
     }
     
     public Distribution getDistr() {
         Distribution result = new Distribution();
         for(Map.Entry<T, Distribution> entry : matrix.entrySet()) {
             result.addDistr(entry.getValue());
         }
         return result;
     }
     
     public void add(T position, Distribution distr) {
         Distribution inDistr = matrix.get(position);
         if(inDistr==null)
             inDistr = new Distribution();
         inDistr.add(distr);
         matrix.put(position, inDistr);
     }
     
     public double getDistrSum() {
         
         return getDistr().sum();
     }
     
    @Override
    public Iterator<Map.Entry<T, Distribution>> iterator() {
        return matrix.entrySet().iterator();
    }

    @Override
    public String toString() {
        String result = "FDMatrix{";
        for(Map.Entry<T, Distribution> entry : matrix.entrySet()) {
            result += "\n" + entry.getKey() + " : " + entry.getValue();
         }
        result += "}";
        return result;
    }
    
    @Override
    public boolean equals(Object o) {
        if ( this == o ) return true;
        if ( !(o instanceof FDMatrix) ) return false;
        FDMatrix<T> cmp = (FDMatrix<T>) o;
        if(size() != cmp.size())
            return false;
        Distribution d1, d2;
        Iterator<Map.Entry<T, Distribution>> it = iterator();
        while(it.hasNext()) {
            Map.Entry<T, Distribution> entry = it.next();
            d1 = entry.getValue();
            d2 = cmp.matrix.get(entry.getKey());;
            if(!d1.equals(d2))
                return false;
        }
        return true;
        
    }
    
    
     
}