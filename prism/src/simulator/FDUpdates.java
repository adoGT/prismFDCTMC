package simulator;


import parser.ast.Updates;

public class FDUpdates {
    
    private Updates updates;
    private String fDelay;
    
    
    public FDUpdates(Updates updates, String delay) {
        this.updates = updates;
        this.fDelay = delay;
    }
    
    public void setFDelay(String delay) {
        this.fDelay = delay;
    }
    
    public String getFDelay() {
        return fDelay;
    }
    
    public void setUpdates(Updates updates) {
        this.updates = updates;
    }
    
    public Updates getUpdates() {
        return updates;
    }
}