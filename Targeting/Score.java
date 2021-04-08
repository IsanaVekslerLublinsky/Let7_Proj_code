package Targeting;

public class Score
{
	int _matches=0;//matches
	int _mismatches=0;//mismatches (starting after the first matching point)
	int _GU=0;//GU
	int _bulges=0;//bulges on the mRNA
	int _bulges_miR=0;
	int _indexStart=-1;//index on the miRNA from the first match
	
	int _totalErrors=0;//mismatches + GU + bulges 
	
	int _maxConsecutiveMatches = 0; //not interrupted by anything
	int _maxConsecutiveMatchesAndGUs = 0;//matches+GUs
	int _currConsecutiveMatches = 0; //the last block of matches
	int _currConsecutiveMatchesAndGUs = 0;
	
	int _bestIndex_i,_bestIndex_j; //for easier trace back
	public String[] _alignment;
	public String _1D;
	
	public int _currIndex_i,_currIndex_j;
	
	public int _startOfmiRNA; //for the non seed region, the actual start of the paired alignment
	
	public Score(int matches)
	{
		_matches = matches;
		_alignment = new String[4];
		_1D = "";
		
		for(int i=0;i<4;i++)
		{
			_alignment[i]="";
		}
	}
	
	public void changeLines(String line1, String line2,String line3,String line4)
	{
	   _alignment[0] = line1;
	   _alignment[1] = line2;
	   _alignment[2] = line3;
	   _alignment[3] = line4;
	}
	
	public String toString()
	{
		String ans = "";
		ans = _matches + ":" + _GU + ":" + _mismatches + ":" + _bulges + "\n";
		for(int i=0;i<4;i++)
			ans = ans + _alignment[i] + "\n"; 
		
		return ans;
	}
}
