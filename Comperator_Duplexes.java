package Targeting;

import java.util.Comparator;

import BasicComponents.Gene;

public class Comperator_Duplexes implements Comparator<Object>{
    public int compare(Object obj1, Object obj2)
    {
    	Gene gene1 = (Gene)obj1;
    	Gene gene2 = (Gene)obj2;

    	if (gene1._duplexes.size() < gene2._duplexes.size())
    		return 1;
    	else if(gene1._duplexes.size() > gene2._duplexes.size())
    		return -1;
    	else
    		return 0;
    	
    	//sorted from high to low
    }

}
