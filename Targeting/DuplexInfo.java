package Targeting;

import java.util.StringTokenizer;

//import let7_targeting.DuplexesFileReader.DuplexesPair;

import BasicComponents.Gene;
import BasicComponents.miRNA;

public class DuplexInfo
{
	public Gene _gene;
	public miRNA _miRNA;
	public String _geneName;
	public String[] _alignment,_alignemntStraight;
	public String _targetSeq;
	public int _locOnTarget;
	public double _energy;
	
	public String _targetPart;//that participates in the alignment
	public int _targetPart_start;
	public int _targetPart_end;
	
	public int _startOnTarget,_endOnTarget; //the actual coordinated 5' to 3'
	public int _startOnMiRNA,_endOnMiRNA; //the actual coordinated 5' to 3'
	
	public String _duplex_1D,_duplex_1D_trimmed;
	public int _matches;
	public int _maxConsecutiveMatches,_maxConsecutiveMatchesAndGus;
	public String[] _seedType;//type, seed-like,stand_alone_bulge_G,stand_alone_GU,stand_alone_perfect
	
	public String[] _threeSuppl;//type,isSuppl
	public int _central;//0-no,1-perfect,2-imperfect
	
	public String _typeSeed,_typeCentral;
	public Score _bestScore;
	
	public boolean _overlapsWithSeed;
	public DuplexInfo _overlappingSeed;
	
	public String _seedContained; //if the see 1d contained in the full duplex 1d
	
	public DuplexInfo _fullDuplex; //for the seed duplex to keep everything
	
	public String _mRNA; //the actual part that matched
	
	public String getCentral() 
	{
		if(_central==0)
			return "none";
		else if(_central==1)
			return "perfect";
		else
			return "imperfect";
	}

	public void parseSeedType(String token) 
	{
		//[1m6mer8m;1GU;1b;, T, F, F, F]
		String temp = token.substring(1,token.length()-1);
		StringTokenizer st = new StringTokenizer(temp,",");
		System.out.println(token + " " + st.countTokens());
		_seedType = new String[5];
		if(st.countTokens()==5)
			_seedType[0] =  st.nextToken().trim();
		
		for(int i=1;i<5;i++)
			_seedType[i] = st.nextToken().trim();
	}

	/*public void parseThreeSuppl(String value) {
		//TODO
		
		if(value.equals("true"))
			_threeSuppl = true;
		else
			_threeSuppl = false;
			
	}

	public void parseCentral(String value) 
	{
		if(value.equals("none"))
			_central = 0;
		else if(value.equals("perfect"))
			_central = 1;
		else
			_central = 2;
	}*/

	public String getType(int option) 
	{
		
		if(option==1)
			return 
					(_seedType[0] + "\t" + _seedType[1] + "\t" + _seedType[2] + "\t" + _seedType[3] + "\t" + _seedType[4] + "\t" + _threeSuppl[0] + "\t" + _threeSuppl[1] + "\t" + _central);
		
		String ans = "";
		if(_seedType[2].equals("T")||_seedType[3].equals("T")||_seedType[4].equals("T"))
		{
			ans = "seed = " + _seedType[0];
			if(_threeSuppl[1].equals("T"))
				ans = ans + " + 3suppl";
		}
		else if(_seedType[1].equals("T") && _threeSuppl[1].equals("T"))
		{
			ans = "seed = " + _seedType[0] + " + 3suppl";
		}
		else if(_central==1)
			ans = "central perfect";
		else if(_central==2)
			ans = "central imperfect";
		else
			ans = "other";
		
		return ans;
	}
	
	////////////////////////////////////////////////////////////////
	boolean _seed2_8 = false;
	boolean _seed2_7 = false;
	boolean _seed2_8_GU = false;
	boolean _seed2_8_b = false;
	boolean _seed2_8_b_G = false; //G bulge between pos 5 and 6
	boolean _seed2_7_b_G = false;
	boolean _inPair = false;
	//DuplexesPair _pair = null;
	
	int _type = 3;//1-perfect, 2-standAlone, not perfect, 3 - not stand alone, 4- doesn't pass
	
	//checks the seed based on the _duplex_1D and bestScore info
	public void checkSeedType()
	{
		//seed2-8
		if(_bestScore._matches==7 && _bestScore._totalErrors==0)
		{
			_seed2_8 = true;
			_type = 1;
		}
		//seed 2-7
		else if(_bestScore._matches==6 && _bestScore._totalErrors==0)
		{
			_seed2_7 = true;
			_type = 1;
		}
		else if(_duplex_1D.startsWith("111111"))
		{
			_seed2_7 = true;
			_type = 1;
		}
		else if(_bestScore._matches==6 && _bestScore._GU==1 && _bestScore._totalErrors==1)
		{
			//seed 2-8, 1 GU
			_seed2_8_GU = true;
			_type = 2;
		}
		
		else if(_bestScore._matches==7 && _bestScore._bulges==1 && _bestScore._totalErrors==1)
		{
			//seed 2-8, 1 b
			_seed2_8_b = true;
			_type = 2;
		}
		//seed 2-8, 1 b, G pos 5-6
		else if(_bestScore._matches==7 && _bestScore._bulges==1 && _bestScore._totalErrors==1 && _duplex_1D.charAt(4)=='4')
		{
			_seed2_8_b_G = true;
			_type = 2;
		}
		else if(_bestScore._matches>=6 && _bestScore._totalErrors==1)
		{
			//seed 2-7, 1 b, 1 mismatch
			//seed 2-8, 1 mismatch
			_type = 3;
		}
		else if(_bestScore._matches==5 && _bestScore._GU==1)
		{
			//seed 2-7, 1 GU
			_type = 3;
		}
		else if(_bestScore._matches==4 && _bestScore._GU==1)
		{
			//2-6, 1 GU
			_type=3;
		}
		else if(_bestScore._matches==4 && _bestScore._totalErrors==0)
		{
			//2-5
			_type=3;
		}
		else
			_type=4;
	}

	public String getSeedType() 
	{
		if(_seed2_8)
			return"seed_2_8";
		if(_seed2_7)
			return "seed_2_7";
		if(_seed2_8_GU)
			return "seed_2_8_GU";
		if(_seed2_8_b_G)
			return "seed_2_8_b_G";
		if(_seed2_8_b)
			return "seed_2_8_b";
		if(_seed2_7_b_G)
			return "seed_2_7_b_G";
		
		return "other";
	}

	
	public void checkNonSeedType()
	{
		// 1 - perfect 6mer and more
		// 2 - perfect 6mer with up to two GUs and more
		// 3 - other
		if(_bestScore._maxConsecutiveMatches>=6)
			_type = 1;
		else if(_bestScore._maxConsecutiveMatchesAndGUs>=7)
			_type = 2;
		else 
			_type = 3;
	}

	public void checkNonSeedType_13_16()
	{
		if(_bestScore._maxConsecutiveMatches>=4)
			_type = 1;
		else 
			_type = 3;
	}
	
	public void checkCentralMatchType()
	{
		if(_bestScore._matches>=11)
			_type=1;
		else if(_bestScore._matches>=10 && _bestScore._totalErrors<=1)
			_type=2;
		else if(_bestScore._matches>=9 && _bestScore._totalErrors<=2)
			_type=3;
		else
			_type=4;
	}
	
	public String getNonSeedType() 
	{
		if(_type==1)
			return "perfect";
		else if(_type==2)
			return "GUs";
		else
			return "other";	
	}
	
	
	public String getAlignment() 
	{
		String ans = "";
		
		for(int i=3;i>=0;i--)
		{
			ans = ans + _alignment[i] + "\n";
		}
		return ans;
	}

	//new May 2020
	int _seedTypeNew = -1;
	int _loopSize = 0;
	int _seedStart = -1;
	int _seedEnd = -1;
	int nonSeedStart = -1;
	int nonSeedEnd = -1;
	int _nonSeedType = -1;
	
	public void classifyDuplex() {
		
		//based on duplex_1D
		//get to the first position that is not 4 or 5 -- the beginning of miR
		int curPos = 0;
		int curPosOnTheMir = 0;
		
		while(_duplex_1D.charAt(curPos)=='4' || _duplex_1D.charAt(curPos)=='5')
			curPos++;
		
		curPos=curPos+1; //get to the second pos on the miRNA
		_seedStart = curPos;
		//curPosOnTheMir = 1;
		
		if(_duplex_1D.substring(curPos,curPos+6).equals("111111")) { //2-7 perfect
			//System.out.println(_duplex_1D + " " + curPos);
			
			//System.out.println(_duplex_1D.substring(curPos,curPos+6));
			_seedEnd = _seedStart+6-1;
			_seedTypeNew = 1;
			curPosOnTheMir = 6;
		}
		else if(_duplex_1D.substring(curPos,curPos+3).equals("111")){ //2-4 perfect and 5-7 1b or 5-8GU/mismatch
			curPos=curPos+3;
			if(_duplex_1D.substring(curPos,curPos+3).contains("4") || _duplex_1D.substring(curPos,curPos+3).contains("5")) 
			{
				int temp = count(_duplex_1D.substring(curPos,curPos+4),'1');
				if(temp==3)
					_seedTypeNew = 2;
				
				_seedEnd = _seedStart + 7 -1;
				curPosOnTheMir = 6;
			}
			else {
				int temp_1 = count(_duplex_1D.substring(curPos,curPos+4),'1');
				int temp_2 = count(_duplex_1D.substring(curPos,curPos+4),'2');
				int temp_3 = count(_duplex_1D.substring(curPos,curPos+4),'3');
				if(temp_1==3 && (temp_2==1 || temp_3==1))
					_seedTypeNew = 2;
				
				_seedEnd = _seedStart + 7 -1;
				curPosOnTheMir = 7;
			}
			
		}
		else { //2-4 1b + 5-7 perfect or 2-4 mismatch/GU + 5-8 perfect
			//TODO
			if(_duplex_1D.substring(curPos,curPos+3).contains("4") || _duplex_1D.substring(curPos,curPos+3).contains("5"))
			{
				int temp = count(_duplex_1D.substring(curPos,curPos+7),'1');
				if(temp==6)
					_seedTypeNew = 3;
				
				_seedEnd = _seedStart + 7 -1;
				curPosOnTheMir = 6;
			}
			else {
				int temp_1 = count(_duplex_1D.substring(curPos,curPos+3),'1');
				int temp_2 = count(_duplex_1D.substring(curPos,curPos+3),'2');
				int temp_3 = count(_duplex_1D.substring(curPos,curPos+3),'3');
				int temp_4 = count(_duplex_1D.substring(curPos+3,curPos+7),'1');
				if(temp_1==2 && (temp_2==1 || temp_3==1) && temp_4==4) 
					_seedTypeNew = 3;
				
				_seedEnd = _seedStart + 7 -1;
				curPosOnTheMir = 7;
			}
		}
		
		//11-16 region
		curPos = _seedEnd+1;
		while(curPosOnTheMir<10) {
			char temp = _duplex_1D.charAt(curPos);
			if(temp=='1' || temp=='2' || temp=='3' || temp=='6')
				curPosOnTheMir++;
		
			if(curPosOnTheMir<10)
				curPos++;
		}
		
		nonSeedStart = curPos;
		
		//now curPos is on the 11nt of the miR
		//11-13
		int temp_1_1 = count(_duplex_1D.substring(curPos,curPos+3),'1');
		int temp_1_2 = count(_duplex_1D.substring(curPos,curPos+3),'2');
		
		curPos++;
		while(curPosOnTheMir<11) {
			char temp = _duplex_1D.charAt(curPos);
			if(temp=='1' || temp=='2' || temp=='3' || temp=='6')
				curPosOnTheMir++;
		
			if(curPosOnTheMir<11)
				curPos++;
		}
		//12-14
		int temp_2_1 = count(_duplex_1D.substring(curPos,curPos+3),'1');
		
		curPos++;
		while(curPosOnTheMir<12) {
			char temp = _duplex_1D.charAt(curPos);
			if(temp=='1' || temp=='2' || temp=='3' || temp=='6')
				curPosOnTheMir++;
		
			if(curPosOnTheMir<12)
				curPos++;
		}
		//12-14
		int temp_3_1 = count(_duplex_1D.substring(curPos,curPos+3),'1');
		int temp_3_2 = count(_duplex_1D.substring(curPos,curPos+3),'1');
		
		
		if(temp_1_1==3 || temp_2_1==3 || temp_3_1==3)
			_nonSeedType=1;
		else if(temp_1_1+temp_1_2==3 || temp_3_1+temp_3_2==3 || temp_3_1+temp_3_2==3)
			_nonSeedType=2;
		
		//Loop
		//count everything except 6
		
		_loopSize = countOtherThan(_duplex_1D.substring(_seedEnd+1,nonSeedStart),'6');
		
	}

	private int count(String stringLong, char ch) {
		int ans=0;
		for(int i=0;i<stringLong.length();i++)
			if(stringLong.charAt(i)==ch)
				ans++;
		return ans;
	}
	
	private int countOtherThan(String stringLong, char ch) {
		int ans=0;
		for(int i=0;i<stringLong.length();i++)
			if(stringLong.charAt(i)!=ch)
				ans++;
		return ans;
	}
}
