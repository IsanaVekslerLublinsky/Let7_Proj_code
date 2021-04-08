package Targeting;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import java.util.LinkedList;
import java.util.StringTokenizer;


public class RNAduplex
{

	public String RNAduplexLoc = Parameters.RNAduplexLoc;
	
	private String _seq1, _seq2 ; //seq1 = miRNA

	public String _seq2_part; //the actual part that participates in the alignment
	public int _seq2_part_start;
	public int _seq2_part_end;
	
	public String _dotBracket;
	public double _energy_duplex = 0;

	String _seq;
	String _mfeLine;

	public String[] _alignment;
	public int _matches=0; //number of matches+GUs across the duplex

	public String[] _alignmentStraight;
	
	public RNAduplex()
	{
       //empty const
	}

	public void call_RNAcofold(String seq1 , String seq2) throws Exception
	{

		this._seq1 = seq1;
		this._seq2 = seq2;

		try{

		    String command = RNAduplexLoc + "RNAduplex";
	        Process child = Runtime.getRuntime().exec(command);

	        // Get output stream to write from it
	        OutputStream out = child.getOutputStream();

	        out.write((_seq1 + "\n" + _seq2+ "\n").getBytes());
	        out.close();

	        //System.out.println("TRYING " + _seq1 + " " + _seq2);

	        InputStream in = child.getInputStream();
	        int c;

	        String ans = "";
	        while ((c = in.read()) != -1 && c!='\n') {
	            ans = ans + ((char)c);
	        }
	        in.close();

	        //System.out.println(ans);
	        child.destroy();

	        _mfeLine = ans;

			int a = _mfeLine.lastIndexOf("(");
            int b = _mfeLine.lastIndexOf(")");

            if(ans.length()>0)
            {
	            _energy_duplex = Double.parseDouble(_mfeLine.substring(a + 1, b).trim());
	            //System.out.println(_mfeLine);
	            _dotBracket = _mfeLine.substring(0,a);
            }
            else
            	_energy_duplex = 0;
            //System.out.println("DOT " + _dotBracket);
            
		}
		catch(IOException ex){
			throw new Exception("Problem calling RNAduplex " + ex);

		}
	}//call_RNAduplex

		
	//cut the ends of the target to eliminate long overhang
	public void translateDotBracketToAlignment()
	{
		//((((((((..(.(((((.((((&.)))).)).).)).)..))))))))   1,22  :   6,30  (-18.70)

		//System.out.println("DOT " + _dotBracket);
		StringTokenizer st = new StringTokenizer(_dotBracket);
		String match = st.nextToken();
		String miRNA_indexes = st.nextToken();
		st.nextToken(); //:
		String mRNA_indexes = st.nextToken();

         int ampersand = match.indexOf('&');
         
         String miRNAmatch  = match.substring(0, ampersand);
         String mRNAmatch = reverse(match.substring(ampersand+1).trim());
         
          //instead of cutting the sequences
         //add ... on the dotBracket
         st = new StringTokenizer(miRNA_indexes,",");
         int begin1 = Integer.parseInt(st.nextToken());
         int end1 = Integer.parseInt(st.nextToken());

         //String miRNA  = _seq1.substring(begin-1,end);
         String miRNA = _seq1;
         
         
         
         st = new StringTokenizer(mRNA_indexes,",");
         int begin2 = Integer.parseInt(st.nextToken());
         int end2 = Integer.parseInt(st.nextToken());

         String mRNA = reverse(_seq2);
         
         //miRNA has to be completed
         for(int i=1;i<begin1;i++)
        	 miRNAmatch = "." + miRNAmatch;
         for(int i=end1;i<_seq1.length();i++)
        	 miRNAmatch = miRNAmatch + ".";
         
         //mRNA has to be completed only partially
         int toAdd = Math.min(_seq2.length()-end2, begin1-1); 
         for(int i=0;i<toAdd;i++)
         {
        	 mRNAmatch = "." + mRNAmatch;
         }
         end2 = end2+toAdd;
         
         toAdd = Math.min(_seq1.length()-end1, begin2-1); 
         for(int i=1;i<begin2;i++)
         {
        	 mRNAmatch = mRNAmatch + ".";
         }
         begin2 = begin2-toAdd;
         
         mRNA = reverse(_seq2.substring(begin2-1,end2));
         
         _seq2_part = _seq2.substring(begin2-1,end2);
         _seq2_part_start = begin2-1;
         _seq2_part_end = end2;
         
         LinkedList<Integer> miRNAstack = new LinkedList<Integer>();
         LinkedList<Integer> mRNAstack = new LinkedList<Integer>();
            //read "unclosed brackets" "(" in mRNA
         for (int i = 0; i < miRNA.length(); i++)
         {
            if (miRNAmatch.charAt(i)=='(')
                miRNAstack.add(new Integer(i));
            else if (miRNAmatch.charAt(i)==')')
                miRNAstack.removeLast();
         }
            //
        //System.out.println(mRNA)
        for (int i =0 ; i< mRNA.length(); i++)
        {
            if (mRNAmatch.charAt(i)==')')
                mRNAstack.add(new Integer(i));
            else if (mRNAmatch.charAt(i)=='(')
                mRNAstack.removeLast();
        }


        String[] lines = new String[4];
        for (int i = 0; i < 4; i++)
            lines[i] = "";


        int miRNA_start = -1;
        int mRNA_start = -1;
        int miRNA_end=0,mRNA_end=0;

        String miRNAsub,mRNAsub;
        
        while (miRNAstack.size() > 0)
        {
        	miRNA_end = miRNAstack.removeFirst().intValue();
        	mRNA_end = mRNAstack.removeFirst().intValue();
        	if (miRNA_end - miRNA_start==1)
        		miRNAsub="";
        	else
        		miRNAsub = miRNA.substring(miRNA_start+1,miRNA_end);

        	if (mRNA_start == mRNA_end)
        		mRNAsub="";
        	else
        		mRNAsub = mRNA.substring(mRNA_start+1,mRNA_end);


        	int diff1 = miRNAsub.length();
			int diff2 = mRNAsub.length();


        	if(diff1<diff2)
        	{
    			lines[0]= lines[0] +  miRNAsub + createString(' ',diff2-diff1) + ' ';
    			lines[1]= lines[1] + createString(' ' , diff2) + miRNA.charAt(miRNA_end);
    			lines[2]= lines[2] + createString(' ' , diff2) + mRNA.charAt(mRNA_end);
    			lines[3]= lines[3] + mRNAsub + ' ';

       		}
    		else
    		{
    			lines[0] = lines[0] + miRNAsub + ' ';
    			lines[1] = lines[1] + createString(' ' , diff1) + miRNA.charAt(miRNA_end);
    			lines[2] = lines[2] + createString(' ' , diff1) + mRNA.charAt(mRNA_end);
    			lines[3] = lines[3] + mRNAsub + createString(' ',diff1-diff2) + ' ';
    		}

        	miRNA_start = miRNA_end;
        	mRNA_start = mRNA_end;

        }
        if (miRNA_end<miRNA.length())
        	lines[0] = lines[0]+miRNA.substring(miRNA_end+1);
        if (mRNA_end<mRNA.length())
        	lines[3] = lines[3]+mRNA.substring(mRNA_end+1);

        _alignment = lines;
        straightLines();
	}
		
	//bring the alignment lines to the same length
	public void straightLines() {
		int maximalLength = Math.max(_alignment[0].length(),_alignment[1].length());
		maximalLength = Math.max(maximalLength, _alignment[2].length());
		maximalLength = Math.max(maximalLength, _alignment[3].length());
		
		char[] line1 = fillArray(maximalLength,_alignment[0]);
		char[] line2 = fillArray(maximalLength,_alignment[1]);
		char[] line3 = fillArray(maximalLength,_alignment[2]);
		char[] line4 = fillArray(maximalLength,_alignment[3]);
		
		_alignmentStraight = new String[4];
		_alignmentStraight[0] = new String (line1);
		_alignmentStraight[1] = new String (line2);
		_alignmentStraight[2] = new String (line3);
		_alignmentStraight[3] = new String (line4);					
	}
			
			
	private String reverse(String toFlip){
		String ans = "";
		for (int i=0;i<toFlip.length();i++){
			ans = toFlip.charAt(i)+ans;
		}
		return ans;
	}

	private String createString(char a , int diff)
	{
		String ans="";
		for(int i=0; i<diff;i++)
		{
			ans = ans + a;
		}
		return ans;
	}

	////***************************************************
	public String translateTo1D() 
	{
		//translates the alignment into a 1D string
		//1-match(AU/UA/GC/CG)
		//2-GU/UG
		//3-mismatch
		//4-bulge on mRNA G
		//5-bulge on mRNA other
		//6-bulge on miRNA
		
		//lines[0][1] - miRNA
		//lines[2][3] - mRNA
		String ans = "";
		_matches = 0;
		
		int maximalLength = Math.max(_alignment[0].length(),_alignment[1].length());
		maximalLength = Math.max(maximalLength, _alignment[2].length());
		maximalLength = Math.max(maximalLength, _alignment[3].length());
		
		char[] line1 = fillArray(maximalLength,_alignment[0]);
		char[] line2 = fillArray(maximalLength,_alignment[1]);
		char[] line3 = fillArray(maximalLength,_alignment[2]);
		char[] line4 = fillArray(maximalLength,_alignment[3]);
		
		for(int i=0;i<line1.length;i++)
		{
			if(line2[i]!=' ')//match/GU
			{
				if((line2[i]=='G' && line3[i]=='T') || (line2[i]=='T' && line3[i]=='G'))
					ans = ans + "2";
				else
				{
					ans = ans + "1";
					_matches++;
				}
			}
			else if(line1[i]!=' ' && line4[i]!=' ')//mismatch
				ans = ans + "3";
			else if(line4[i]=='G') //G bulge on mRNA
				ans = ans + "4";
			else if(line4[i]!=' ')
				ans = ans + "5";
			else
				ans = ans + "6";
		}
		
		return ans;
	}
	
	private char[] fillArray(int maximalLength, String string) 
	{
		char[] answer = new char[maximalLength];
		for(int i=0;i<maximalLength;i++)
		{
			answer[i] = ' ';
		}
		
		for(int i=0;i<string.length();i++)
		{
			answer[i] = string.charAt(i);
		}
		return answer;
	}
		
		
}

