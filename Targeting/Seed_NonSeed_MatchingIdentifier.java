package Targeting;

public class Seed_NonSeed_MatchingIdentifier 
{
	
	//the targetSeq is already reversedComplemented
	Score[][] scoresTable;
	char[] miRNA_arr;
	char[] target_arr;
	
	int actual_miR_Length;
	int actual_target_Length;
	
	//for the seed region
	int mismatches_max=1; //TODO change to 0
	int GU_max=1;
	int bulges_max=1;
	int maxTotalErrors = 1;
	int minPosForBulge_Mismatch = 0;//after 3 matches on the seed (0,1,2)
	
	//for the non seed / central
	int mismatches_max_NS=1;
	int GU_max_NS=2;
	int bulges_max_NS=1;
	int bulges_miR_max_NS=1;
	int maxTotalErrors_NS = 2;
	
	int Inf=100;
	
	Score bestScore;
	//no overhangs allowed on either side
	public Seed_NonSeed_MatchingIdentifier(int miRNA_length, int target_length) 
	{
		scoresTable = new Score[miRNA_length+1][target_length+1];
	}
	
	//option 1- seed, 2-nonSeed
	public void processSeed(String miRNA, String targetSeq) 
	{
		//System.out.println(miRNA + " " + targetSeq);
		
		//TODO - make it more efficient
		miRNA_arr = (" " + miRNA).toLowerCase().toCharArray();
		target_arr = (" " + targetSeq).toLowerCase().toCharArray();
		
		actual_miR_Length = miRNA.length();
		actual_target_Length = targetSeq.length();
		
		//init
		for(int i=0;i<actual_miR_Length+1;i++)
			for(int j=0;j<actual_target_Length+1;j++)
			{
				//System.out.println(i + " " + j);
				scoresTable[i][j] = new Score(-Inf);
			}
		
		scoresTable[0][0]._matches = 0;
		
		///
		bestScore = scoresTable[0][0];
		
		for(int i=1;i<actual_miR_Length+1;i++)
		{
			//Since bulges are not allowed in the miRNA it starts with j=1, otherwise j=1
			for(int j=i;j<actual_target_Length+1;j++)
			{
				Score score = compareNodesSeed(i,j);
				
				score._currIndex_i = i;
				score._currIndex_j = j;
				
				scoresTable[i][j] = score;
				
			}
		}
		
		
	}

	
	//NEW -  if a full seed 2-7 is achieved, and position 8 cannot form match, do not try to extend it further
	public Score compareNodesSeed(int i, int j)
	{
		//System.out.println(i + " " + j + " " + miRNA_arr[i] + " " + target_arr[j]);
		Score ans = new Score(-Inf);
		Score temp1,temp2;
		
		int d; //0-match,1-wobble,2-mismatch
		if (match(miRNA_arr[i],target_arr[j])) 
			d=0;
		else if(wobble(miRNA_arr[i],target_arr[j]))
			d=1;
		else d=2;
		
		temp1 = scoresTable[i-1][j-1];
		temp2 = scoresTable[i][j-1];
		//Score temp3 = scoresTable[i-1][j]; //no bulges in miRNA are allowed

		if(temp1._matches==-Inf && temp2._matches==-Inf)
		{
			//System.out.println("returning");
			ans._matches = -Inf;
			return ans;
		}
		
		int exp1,exp2;
		
		if(d==0) //d=0
		{
			if(temp1._matches!=-Inf)
				exp1 = temp1._matches+1;
			else
				exp1 = temp1._matches;
		}
		else //wobble/mismatches
		{
			if(temp1._totalErrors<maxTotalErrors)
			{
				if(d==1)
				{
					if(temp1._GU<GU_max)
						exp1 = temp1._matches;
					else
						exp1 = -Inf;
				}
				else 
				{
					if(temp1._mismatches<mismatches_max && i>=minPosForBulge_Mismatch)
						exp1 = temp1._matches;
					else
						exp1 = -Inf;
				}
			}
			else
				exp1 = -Inf;
		}
				
		if(temp2._totalErrors<maxTotalErrors && temp2._bulges<bulges_max && i>=minPosForBulge_Mismatch)
		{
			exp2 = temp2._matches;
		}
		else
			exp2 = -Inf;

		
		//NEW
		//if both are equal, but one of them is a full seed and the other is not take the seed
		if(exp1==exp2)
		{
			if(temp1._1D.contains("111111"))
				exp2=-Inf;
			else if(temp2._1D.contains("111111"))
				exp1=-Inf;
		}
		
		int tempy = Math.max(exp1, exp2);
			
		ans._matches = tempy;

		if(tempy==-Inf)
			return ans;
		
		if (ans._matches == exp1)
		{
			ans._bestIndex_i = i-1;
			ans._bestIndex_j = j-1;
			
			if(d==0 || d==1) //d=0
			{	
				if(d==0)
				{
					ans._GU = temp1._GU;
					ans.changeLines(temp1._alignment[0]+" ", temp1._alignment[1]+miRNA_arr[i], temp1._alignment[2]+target_arr[j], temp1._alignment[3]+" ");
					ans._totalErrors = temp1._totalErrors;
					ans._1D = temp1._1D + "1";
					
					ans._currConsecutiveMatches = temp1._currConsecutiveMatches + 1;
					ans._maxConsecutiveMatches = Math.max(temp1._currConsecutiveMatches + 1, temp1._maxConsecutiveMatches);
					
					ans._currConsecutiveMatchesAndGUs = temp1._currConsecutiveMatchesAndGUs + 1;
					ans._maxConsecutiveMatchesAndGUs = Math.max(temp1._currConsecutiveMatchesAndGUs + 1, temp1._maxConsecutiveMatchesAndGUs);
				}
				if(d==1)
				{
					ans._GU = temp1._GU+1;
					ans._totalErrors = temp1._totalErrors+1;
					ans.changeLines(temp1._alignment[0]+" ", temp1._alignment[1]+miRNA_arr[i], temp1._alignment[2]+target_arr[j], temp1._alignment[3]+" ");
					ans._1D = temp1._1D + "2";
					
					ans._maxConsecutiveMatches = temp1._maxConsecutiveMatches;
					ans._currConsecutiveMatches = 0;
					
					ans._currConsecutiveMatchesAndGUs = temp1._currConsecutiveMatchesAndGUs + 1;
					ans._maxConsecutiveMatchesAndGUs = Math.max(temp1._currConsecutiveMatchesAndGUs + 1, temp1._maxConsecutiveMatchesAndGUs);
			
				}
				
				ans._mismatches = temp1._mismatches;
				if(temp1._indexStart==-1)
					ans._indexStart = i;
				else
					ans._indexStart = temp1._indexStart;
				
				if(ans._matches>bestScore._matches)
					bestScore = ans;
				else if(ans._matches==bestScore._matches && ans._GU>bestScore._GU)
					bestScore = ans;
			}
			else 
			{
				ans._mismatches = temp1._mismatches + 1;
				ans._totalErrors = temp1._totalErrors + 1;
			
				
				ans._GU = temp1._GU;
				ans._indexStart = temp1._indexStart;
				ans.changeLines(temp1._alignment[0]+miRNA_arr[i], temp1._alignment[1]+" ", temp1._alignment[2]+" ", temp1._alignment[3]+target_arr[j]);
				ans._1D = temp1._1D + "3";
				
				ans._maxConsecutiveMatches = temp1._maxConsecutiveMatches;
				ans._currConsecutiveMatches = 0;
				
				ans._maxConsecutiveMatchesAndGUs = temp1._maxConsecutiveMatchesAndGUs;
				ans._currConsecutiveMatchesAndGUs = 0;
			}
			
			ans._bulges = temp1._bulges;
		}
		
		else if (ans._matches == exp2)
		{
			ans._bestIndex_i = i;
			ans._bestIndex_j = j-1;
			
			ans._bulges = temp2._bulges + 1;
			ans._totalErrors = temp2._totalErrors + 1;
			
			ans._mismatches = temp2._mismatches;
			ans._GU = temp2._GU;
			ans._indexStart = temp2._indexStart;

			ans.changeLines(temp2._alignment[0]+" ", temp2._alignment[1]+" ", temp2._alignment[2]+" ", temp2._alignment[3]+target_arr[j]);

			if(target_arr[j]=='g' || target_arr[j]=='G')
				ans._1D = temp2._1D + "4";
			else
				ans._1D = temp2._1D + "5";
			
			ans._maxConsecutiveMatches = temp2._maxConsecutiveMatches;
			ans._currConsecutiveMatches = 0;
			
			ans._maxConsecutiveMatchesAndGUs = temp2._maxConsecutiveMatchesAndGUs;
			ans._currConsecutiveMatchesAndGUs = 0;

		}
		/*
		System.out.println(ans._matches);
		System.out.println(ans._line1);
		System.out.println(ans._line2);
		System.out.println(ans._line3);
		System.out.println(ans._line4);
		*/
		return ans;
		
	}
	
	///////////////////// Non seed region
	
	public void processNonSeed(String miRNA, String targetSeq) 
	{
		//ystem.out.println(miRNA + " " + targetSeq);
		
		miRNA_arr = (" " + miRNA).toLowerCase().toCharArray();
		target_arr = (" " + targetSeq).toLowerCase().toCharArray();
		
		actual_miR_Length = miRNA.length();
		actual_target_Length = targetSeq.length();
		
		//init
		scoresTable[0][0] = new Score(0);
		
		for(int i=1;i<actual_miR_Length+1;i++)
		{
			scoresTable[i][0] = new Score(0);
			//Score temp3 = scoresTable[i-1][0];
			//scoresTable[i][0].changeLines(temp3._alignment[0]+miRNA_arr[i], temp3._alignment[1]+" ", temp3._alignment[2]+" ", temp3._alignment[3]+" ");
			scoresTable[i][0]._startOfmiRNA = i;
		}
		
		//no overhangs are allowed in the target, since it is sliding window
		for(int j=1;j<actual_target_Length+1;j++)
			scoresTable[0][j] = new Score(-Inf);
		
		for(int i=1;i<actual_miR_Length+1;i++)
			for(int j=1;j<actual_target_Length+1;j++)
			{
				//System.out.println(i + " " + j);
				scoresTable[i][j] = new Score(-Inf);
			}
		
		///
		bestScore = scoresTable[0][0];
		
		for(int i=1;i<actual_miR_Length+1;i++)
		{
			//Since bulges are not allowed in the miRNA it starts with j=1, otherwise j=1
			for(int j=1;j<actual_target_Length+1;j++)
			{
				Score score = compareNodesNonSeed(i,j);
				
				score._currIndex_i = i;
				score._currIndex_j = j;
				
				scoresTable[i][j] = score;
			}
		}
	}

	
	
	
	//bulges are allowed on both sides
	public Score compareNodesNonSeed(int i, int j)
	{
		//System.out.println(i + " " + j + " " + miRNA_arr[i] + " " + target_arr[j]);
		Score ans = new Score(-Inf);
		
		Score temp1,temp2,temp3;
		
		int d; //0-match,1-wobble,2-mismatch
		if (match(miRNA_arr[i],target_arr[j])) 
			d=0;
		else if(wobble(miRNA_arr[i],target_arr[j]))
			d=1;
		else d=2;
		
		temp1 = scoresTable[i-1][j-1];
		temp2 = scoresTable[i][j-1];
		temp3 = scoresTable[i-1][j]; 

		//System.out.println(i + " " + j);
		
		if(temp1._matches==-Inf && temp2._matches==-Inf && temp3._matches==-Inf)
		{
			//System.out.println("returning");
			ans._matches = -Inf;
			return ans;
		}
		
		int exp1,exp2,exp3;
		
		if(d==0) //d=0
		{
			if(temp1._matches!=-Inf)
				exp1 = temp1._matches+1;
			else
				exp1 = temp1._matches;
		}
		else //wobble/mismatches
		{
			if(temp1._totalErrors<maxTotalErrors_NS)
			{
				if(d==1)
				{
					if(temp1._GU<GU_max_NS)
						exp1 = temp1._matches;
					else
						exp1 = -Inf;
				}
				else 
				{
					//no mismatches are allowed in the beginning
					if(temp1._indexStart>-1 && temp1._mismatches<mismatches_max_NS)
						exp1 = temp1._matches;
					else
						exp1 = -Inf;
				}
			}
			else
				exp1 = -Inf;
		}
		
		//bulges in the target
		if(temp2._totalErrors<maxTotalErrors_NS && temp2._bulges<bulges_max_NS && temp2._indexStart>-1)
		{
			exp2 = temp2._matches;
		}
		else
			exp2 = -Inf;

		if(temp3._totalErrors<maxTotalErrors_NS && temp3._bulges_miR<bulges_miR_max_NS)
		{
			exp3 = temp3._matches;
		}
		else
			exp3 = -Inf;

		
		int tempy = Math.max(exp1, exp2);
		tempy = Math.max(tempy, exp3);
		
		ans._matches = tempy;

		if(tempy==-Inf)
			return ans;
		
		if (ans._matches == exp1)
		{
			ans._bestIndex_i = i-1;
			ans._bestIndex_j = j-1;
			
			if(d==0 || d==1) //d=0
			{	
				if(d==0)
				{
					ans._GU = temp1._GU;
					ans.changeLines(temp1._alignment[0]+" ", temp1._alignment[1]+miRNA_arr[i], temp1._alignment[2]+target_arr[j], temp1._alignment[3]+" ");
					ans._totalErrors = temp1._totalErrors;
					ans._1D = temp1._1D + "1";
					
					ans._currConsecutiveMatches = temp1._currConsecutiveMatches + 1;
					ans._maxConsecutiveMatches = Math.max(temp1._currConsecutiveMatches + 1, temp1._maxConsecutiveMatches);
					
					ans._currConsecutiveMatchesAndGUs = temp1._currConsecutiveMatchesAndGUs + 1;
					ans._maxConsecutiveMatchesAndGUs = Math.max(temp1._currConsecutiveMatchesAndGUs + 1, temp1._maxConsecutiveMatchesAndGUs);
				}
				if(d==1)
				{
					ans._GU = temp1._GU+1;
					ans._totalErrors = temp1._totalErrors+1;
					ans.changeLines(temp1._alignment[0]+" ", temp1._alignment[1]+miRNA_arr[i], temp1._alignment[2]+target_arr[j], temp1._alignment[3]+" ");
					ans._1D = temp1._1D + "2";
					
					ans._maxConsecutiveMatches = temp1._maxConsecutiveMatches;
					ans._currConsecutiveMatches = 0;
					
					ans._currConsecutiveMatchesAndGUs = temp1._currConsecutiveMatchesAndGUs + 1;
					ans._maxConsecutiveMatchesAndGUs = Math.max(temp1._currConsecutiveMatchesAndGUs + 1, temp1._maxConsecutiveMatchesAndGUs);
				}
				
				ans._mismatches = temp1._mismatches;
				if(temp1._indexStart==-1)
					ans._indexStart = i;
				else
					ans._indexStart = temp1._indexStart;
				
				if(ans._matches>bestScore._matches)
					bestScore = ans;
				else if(ans._matches==bestScore._matches && ans._GU>bestScore._GU)
					bestScore = ans;
			}
			else 
			{
				ans._mismatches = temp1._mismatches + 1;
				ans._totalErrors = temp1._totalErrors + 1;
			
				
				ans._GU = temp1._GU;
				ans._indexStart = temp1._indexStart;
				ans.changeLines(temp1._alignment[0]+miRNA_arr[i], temp1._alignment[1]+" ", temp1._alignment[2]+" ", temp1._alignment[3]+target_arr[j]);
				ans._1D = temp1._1D + "3";
				
				ans._maxConsecutiveMatches = temp1._maxConsecutiveMatches;
				ans._currConsecutiveMatches = 0;
				
				ans._maxConsecutiveMatchesAndGUs = temp1._maxConsecutiveMatchesAndGUs;
				ans._currConsecutiveMatchesAndGUs = 0;
			}
			
			ans._bulges = temp1._bulges;
			ans._bulges_miR = temp1._bulges_miR;
			ans._startOfmiRNA = temp1._startOfmiRNA;
		}
		
		else if (ans._matches == exp2)
		{
			ans._bestIndex_i = i;
			ans._bestIndex_j = j-1;
			
			ans._bulges = temp2._bulges + 1;
			ans._totalErrors = temp2._totalErrors + 1;
			
			ans._bulges_miR = temp2._bulges_miR;
			ans._mismatches = temp2._mismatches;
			ans._GU = temp2._GU;
			ans._indexStart = temp2._indexStart;

			ans.changeLines(temp2._alignment[0]+" ", temp2._alignment[1]+" ", temp2._alignment[2]+" ", temp2._alignment[3]+target_arr[j]);

			if(target_arr[j]=='g' || target_arr[j]=='G')
				ans._1D = temp2._1D + "4";
			else
				ans._1D = temp2._1D + "5";
			
			ans._startOfmiRNA = temp2._startOfmiRNA;
			
			ans._maxConsecutiveMatches = temp2._maxConsecutiveMatches;
			ans._currConsecutiveMatches = 0;
			
			ans._maxConsecutiveMatchesAndGUs = temp2._maxConsecutiveMatchesAndGUs;
			ans._currConsecutiveMatchesAndGUs = 0;
		}
		else if (ans._matches == exp3)
		{
			ans._bestIndex_i = i-1;
			ans._bestIndex_j = j;
			
			ans._bulges_miR = temp3._bulges_miR + 1;
			ans._totalErrors = temp3._totalErrors + 1;
			
			ans._bulges = temp3._bulges;
			ans._mismatches = temp3._mismatches;
			ans._GU = temp3._GU;
			ans._indexStart = temp3._indexStart;

			ans.changeLines(temp3._alignment[0]+miRNA_arr[i], temp3._alignment[1]+" ", temp3._alignment[2]+" ", temp3._alignment[3] +" ");

			ans._1D = temp3._1D + "6";
			
			ans._startOfmiRNA = temp3._startOfmiRNA;
			
			ans._maxConsecutiveMatches = temp3._maxConsecutiveMatches;
			ans._currConsecutiveMatches = 0;
			
			ans._maxConsecutiveMatchesAndGUs = temp3._maxConsecutiveMatchesAndGUs;
			ans._currConsecutiveMatchesAndGUs = 0;
		}
		/*
		System.out.println(ans._matches);
		System.out.println(ans._line1);
		System.out.println(ans._line2);
		System.out.println(ans._line3);
		System.out.println(ans._line4);
		*/
		return ans;
		
	}
	
	////////////////////////////////////////////////////////////////////////
	private boolean match(char c1, char c2)//check complementarity
	{
		if ((c1=='a' && c2 =='t')||
			(c1=='t' && c2 =='a')||
			(c1=='a' && c2 =='u')||
			(c1=='u' && c2 =='a')||
			(c1=='g' && c2 =='c')||
			(c1=='c' && c2 =='g'))
			
			return true;
		else return false;
	}

	private boolean wobble(char c1, char c2)
	{
		if ((c1=='u' && c2 =='g')||
			(c1=='t' && c2 =='g')||
			(c1=='g' && c2 =='u')||
			(c1=='g' && c2 =='t'))
			return true;
		else return false;
	}

	
	
}
