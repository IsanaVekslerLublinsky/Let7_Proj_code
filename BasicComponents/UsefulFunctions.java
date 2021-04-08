package BasicComponents;

import java.math.BigDecimal;
import java.util.ArrayList;

public class UsefulFunctions {
	
	public static int findMirWithName(ArrayList<miRNA> mirs,String mirName, int startingIndex)
	{//1-human,2-viral
		int i=startingIndex;
		while(i<mirs.size())
		{
			if (mirs.get(i).getName().equals(mirName))
				return i;
			i++;
		}

		return -1;
	}
	
	//given nameLine
	public static int findGeneWithName(ArrayList<Gene> genes,String geneName, int startingIndex)
	{//1-human,2-viral
		int i=startingIndex;
		while(i<genes.size())
		{
			if (genes.get(i)._nameLine.equals(geneName))
				return i;
			i++;
		}

		return -1;
	}

	public static int compareInts(int num1, int num2)
	{
		if (num1<num2) return -1;
		else if (num1==num2) return 0;
		else return 1;
	}
	
	public static String reverse(String toFlip)
	{
		String ans = "";
		for (int i=0;i<toFlip.length();i++){
			ans = toFlip.charAt(i)+ans;
		}
		return ans;
	}
	
	public static String getComplement(String sequence)
	{

		String ans = sequence.toLowerCase();
		ans = ans.replace('a','T');
		ans = ans.replace('t','A');
		ans = ans.replace('u','A');
		ans = ans.replace('g','C');
		ans = ans.replace('c','G');
		return ans;
	}
	
	//copied from StringUtils *************************************
	 public static int countMatches(String str, String sub) 
	 {
        if (isEmpty(str) || isEmpty(sub)) {
            return 0;
        }
        int count = 0;
        int idx = 0;
        while ((idx = str.indexOf(sub, idx)) != -1) {
            count++;
            idx += sub.length();
        }
        return count;
	 }
	 
	 
	 public static boolean isEmpty(String str) {
	        return str == null || str.length() == 0;
	    }
	 
	 public static double round(double d, int decimalPlace)
	 {
			BigDecimal bd = new BigDecimal(Double.toString(d));
		    bd = bd.setScale(decimalPlace,BigDecimal.ROUND_HALF_UP);
		    return bd.doubleValue();
	 }
}
