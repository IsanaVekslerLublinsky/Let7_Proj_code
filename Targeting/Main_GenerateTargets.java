package Targeting;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;


import BasicComponents.Gene;
import BasicComponents.GeneFileParser;
import BasicComponents.UsefulFunctions;
import BasicComponents.miRNA;
import BasicComponents.miRNAFileParser;



public class Main_GenerateTargets {
	
	//java -Xmx20000m Targeting.Main_GenerateTargets
	//Slides across the target sequence, 
	//employs RNAduplex, prints the info into file
	
	public static void main(String[] args)throws Exception
	{
		GeneFileParser gene_parser = new GeneFileParser(Parameters.geneF);
		miRNAFileParser miRNAparser = new miRNAFileParser(Parameters.miRNAfile);
		
		Main_GenerateTargets generator = new Main_GenerateTargets(gene_parser);
		
		String riboSeqData = Parameters.riboSeqData;
		gene_parser.annotateGenesWithInfo(riboSeqData);
		
		String suffix = "";
		//generator._runMode = 1; //seed
		suffix = "_DP_seed.aln";
		
		generator.updatePositions();
		
		for(int i=0;i<miRNAparser._miRNAs.size();i++)
		{
			if(i==Parameters.miRNAindex)
			{
				miRNA mir = miRNAparser._miRNAs.get(i);
				//System.out.println(mir._bestEnergy);
				String outFile1 = Parameters.outFilePwd + "alignments_" + mir._name + suffix;
				String outFile2 = Parameters.outFilePwd + "alignments_" + mir._name + suffix + ".html";
				String outFile3 = Parameters.outFilePwd + "alignments_" + mir._name + suffix + ".csv";
				
				BufferedWriter out1 = new BufferedWriter(new FileWriter(outFile1));
				BufferedWriter out2 = new BufferedWriter(new FileWriter(outFile2));
				BufferedWriter out3 = new BufferedWriter(new FileWriter(outFile3));
				
				generator.writeToHTMLStart(out2, mir);
				generator.process(miRNAparser._miRNAs.get(i),out1,out2,out3);
				generator.writeToHTMLEnd(out2);
				out1.close();
				out2.close();
				out3.close();
			}
		}	
		
	}

	
	//int _runMode; //1-seed,2-nonSeed
	
	//for the seed parameters
	int _startPos_mir = 1; //2-8
	int _endPos_Mir = 8;
	int _mRNAfragmentSize = 10;
	int _mRNAfragmentSizeFull = 20; //additional 20
	
	GeneFileParser _gene_Parser;
	RNAduplex _duplexer;
	BufferedWriter _out1,_out2,_out3;//alignments,report
	
	Seed_NonSeed_MatchingIdentifier _identifier;
	
	public Main_GenerateTargets(GeneFileParser gene_parser) 
	{
		_gene_Parser =  gene_parser;
		_duplexer = new RNAduplex();	
	}
	
	public void updatePositions()
	{
		//seed 1,8,10
		//non-seed 10,18,11
		_startPos_mir = 1;
		_endPos_Mir = 7; //8
		_mRNAfragmentSize = 10;
			
		_identifier = new Seed_NonSeed_MatchingIdentifier(6,10); //7,10
		
	}
	
	//out2 - html
	public void process(miRNA miRNA,BufferedWriter out1,BufferedWriter out2,BufferedWriter out3) throws Exception
	{
		_out1 = out1;
		_out2 = out2;
		_out3 = out3;
		
		//String header = "miRname\tmiRIndex\tgeneName\tgeneIndex\tstartOnTar\tendOnTar\tstartOnMir\tendOnMir\tM\tMCM\tMCMGU\tGU\tMM\tB\tBM\tseed1D\taln3\taln2\taln1\taln0\ttype\t1D\taln3\taln2\taln1\taln0\tcontained\n";
		String header3 = "indexG\tindexT\tGene\tlog2FC\tpadj\t(1,1)\t(1,2)\t(1,3)\t(2,1)\t(2,2)\t(2,3)\t(3,1)\t(3,2)\t(3,3)\tlist";
		_out3.write(header3+"\n");
		
		//_out1.write(header+"\n");
		ArrayList<Gene> genes = _gene_Parser._fullGenes;
		for(int i=0;i<genes.size();i++)
		//for(int i=0;i<70;i++)	
		{
			//System.out.println("gene " + i);
			slideTheGene_DynamicPrograming(miRNA,genes.get(i),i);
		}
		
		
		//print the genes, collect UTRs of the same gene - sort, and filter 3UTRs with the same sites
		//for UTRs with no duplexes, keep one representative, for enrichment analysis
		Gene prevGene=null;
		ArrayList<Gene> group = new ArrayList<Gene>();
		for(int i=0;i<genes.size();i++)
		//for(int i=0;i<70;i++)	
		{
			Gene gene = genes.get(i);
			String geneName = gene._geneName;
			
			if(group.size()==0)
				group.add(gene);
			else {
				if(prevGene._geneName.equals(geneName)) {
					group.add(gene);
				}
				else
				{
					//System.out.println(indexGroup);
					sortAndPrintGroup(group);
					group.clear();
					group.add(gene);
				}
			}		
			prevGene = gene;
		}
	}

	
	public void sortAndPrintGroup(ArrayList<Gene> genes) throws Exception{
		Comperator_Duplexes comp = new Comperator_Duplexes(); //based on number of duplexes
		Object[] sortedGenes = genes.toArray();
		Arrays.sort(sortedGenes,comp);
		
		ArrayList<Gene> toPrint = new ArrayList<Gene>();
		toPrint.add((Gene)sortedGenes[0]);
		
		for(int i=1;i<sortedGenes.length;i++) {
			Gene gene_prev = (Gene)sortedGenes[0];
			Gene gene_curr = (Gene)sortedGenes[i];
			
			boolean contained = compareDuplexes(gene_prev._duplexes,gene_curr._duplexes);
			if(!contained)
			{
				System.out.println(gene_curr._geneName);
				toPrint.add(gene_curr);
			}
		}
		
		
		for(int i=0;i<toPrint.size();i++) {
			Gene gene = (Gene)sortedGenes[i];	
			printDuplexes(gene,gene._geneName,gene._duplexes);
			printDuplexesToHTML(gene, gene._geneName, gene._duplexes, _out2);
			printDuplexesToCSV(gene,_out3);
		}
	}
	
	private void printDuplexesToCSV(Gene gene, BufferedWriter out) throws Exception{
		
		//need to split the list based on sites combinations
		//Print to CSv all genes
		
		int index=0;
		if(gene._duplexes.size()>0)
			index=runningIndex;
		
			out.write(runningIndex_all + "\t" + index + "\t" + gene._geneName + "\t" + gene._log2FC + "\t" + gene._padj + "\t");
			String list = "";
			int[] groups = new int[9];
			for(int i=0;i<gene._duplexes.size();i++)
			{
				DuplexInfo cur = gene._duplexes.get(i);
				int seed = cur._fullDuplex._seedTypeNew;
				int nonseed = cur._fullDuplex._nonSeedType;
				if(nonseed==-1)
					nonseed=3;
				
				groups[(seed-1)*3+(nonseed-1)]++;
				
				list = list + cur._fullDuplex._seedTypeNew +  "," + cur._fullDuplex._nonSeedType + ";";
			}
			for(int i=0;i<groups.length;i++) {
				out.write(groups[i]+"\t");
			}
			out.write(list+"\n");
		//}
		
		
	}

	private boolean compareDuplexes(ArrayList<DuplexInfo> _duplexes, ArrayList<DuplexInfo> _duplexes2) {
		for(int i=0;i<_duplexes2.size();i++) {
			DuplexInfo dup_2 = _duplexes2.get(i);
			boolean found=false;
			for(int j=0;j<_duplexes.size()&& !found;j++) {
				DuplexInfo dup_1 = _duplexes.get(j);
				if(dup_1._fullDuplex._duplex_1D.equals(dup_2._fullDuplex._duplex_1D))
					found=true;
			}
			if(!found)
				return false;
		}
		return true;
		
	}

	
	private void slideTheGene_DynamicPrograming(miRNA miRNA,Gene gene,int index) throws Exception
	{
		String seq = gene._sequence;
		
		ArrayList<DuplexInfo> duplexes = new ArrayList<DuplexInfo>();
		String geneName = gene._nameLine;
		String miRNAseq = miRNA._sequence.substring(_startPos_mir,Math.min(_endPos_Mir,miRNA._sequence.length()));	
		
		String miRNAseqFull = miRNA._sequence;
		
		for(int i=_mRNAfragmentSizeFull;i<seq.length()-_mRNAfragmentSize;i++)
		{
			String targetSeq = seq.substring(i,i+_mRNAfragmentSize);
			String targetSeq_rev = UsefulFunctions.reverse(targetSeq);
			
			String targetSeqFull = seq.substring(i-_mRNAfragmentSizeFull,i+_mRNAfragmentSize+1);
			
			
			_identifier.processSeed(miRNAseq, targetSeq_rev);
			
			Score bestScore = _identifier.bestScore;
			String[] _alignment = bestScore._alignment;
			String duplex_1D = bestScore._1D;
			int matches = bestScore._matches;
			int maxConsMatches = bestScore._maxConsecutiveMatches;
			int maxConsMatchesAndGUs = bestScore._maxConsecutiveMatchesAndGUs;
			
			//it matched in reverse
			int targetStart = _mRNAfragmentSize - bestScore._currIndex_j;
			int targetEnd = _mRNAfragmentSize;
			String partSeq = targetSeq.substring(targetStart,targetEnd);
			
			DuplexInfo dupInfo = new DuplexInfo();
			dupInfo._miRNA = miRNA;
			dupInfo._gene =gene;
			dupInfo._geneName = geneName;
			dupInfo._alignment = _alignment;
			dupInfo._duplex_1D = duplex_1D;
			dupInfo._locOnTarget = i;
			dupInfo._matches = matches;
			dupInfo._maxConsecutiveMatches = maxConsMatches;
			dupInfo._maxConsecutiveMatchesAndGus = maxConsMatchesAndGUs;
			
			dupInfo._targetSeq = targetSeq;
			
			//TODO - check this
			dupInfo._targetPart = partSeq;
			dupInfo._targetPart_start = targetStart;
			dupInfo._targetPart_end = targetEnd - 1;
			
			dupInfo._startOnTarget = i + targetStart;
			dupInfo._endOnTarget = i + targetEnd -1;
			dupInfo._startOnMiRNA = 1;
			dupInfo._endOnMiRNA = bestScore._currIndex_i; //since the miRNA skips the first nt, and the dynamic programming adds " " in the beginning
			dupInfo._bestScore = bestScore;
			
			dupInfo.checkSeedType();
			
			if(dupInfo._type<4) {
				
				//Here - new Run RNAduplex to get the full alignment
				_duplexer.call_RNAcofold(miRNAseqFull, targetSeqFull);
				_duplexer.translateDotBracketToAlignment();
				
				DuplexInfo dupInfoFull = new DuplexInfo();
				dupInfoFull._targetSeq = targetSeqFull;
				dupInfoFull._alignment = _duplexer._alignment;
				dupInfoFull._alignemntStraight = _duplexer._alignmentStraight;
				dupInfoFull._energy = _duplexer._energy_duplex;
				dupInfoFull._duplex_1D = _duplexer.translateTo1D();
				dupInfo._fullDuplex = dupInfoFull;
				dupInfoFull._mRNA = _duplexer._seq2_part;
				dupInfoFull._targetPart = _duplexer._seq2_part;
				
				dupInfoFull.classifyDuplex();
				
				if(dupInfo._fullDuplex._seedTypeNew!=-1 && dupInfo._fullDuplex._loopSize<=14 && dupInfo._fullDuplex._targetPart.length()>=16)
				{
					//if(gene._geneName.equals("ins-35"))
					//	System.out.println(dupInfo._fullDuplex._targetPart);
					
					//Fix 17/4/2022
					if(duplexes.size()>0 && (duplexes.get(duplexes.size()-1)._fullDuplex._duplex_1D.equals(dupInfo._fullDuplex._duplex_1D) ||
							                duplexes.get(duplexes.size()-1)._fullDuplex._targetPart.contains(dupInfo._fullDuplex._targetPart) ||
											dupInfo._fullDuplex._targetPart.contains(duplexes.get(duplexes.size()-1)._fullDuplex._targetPart))) {
						//nothing
						//System.out.println()
					}
					else
						duplexes.add(dupInfo);
				}
			}
			
		}
		
		gene._duplexes = duplexes;
		
	}
	
	private void printDuplexes(Gene gene, String geneName,ArrayList<DuplexInfo> duplexes) throws Exception
	{
		DuplexInfo curr;
		
		//System.out.println(duplexes.size());
		
		for(int i=0;i<duplexes.size();i++)
		{
			curr = duplexes.get(i);
			//indexes of the miRNAs/genes
			_out1.write(curr._miRNA._name + "\t" + curr._miRNA._orderInList+ "\t");
			_out1.write(curr._geneName + "\t" + curr._gene._orderInList + "\t");
			
			//these are indexes that start from 0
			_out1.write(curr._startOnTarget + "\t" + curr._endOnTarget + "\t");
			
			_out1.write(curr._startOnMiRNA + "\t" + curr._endOnMiRNA + "\t");
			
			_out1.write(curr._matches + "\t");
			_out1.write(curr._bestScore._maxConsecutiveMatches + "\t");
			_out1.write(curr._bestScore._maxConsecutiveMatchesAndGUs + "\t");
			_out1.write(curr._bestScore._GU + "\t");
			_out1.write(curr._bestScore._mismatches + "\t");
			_out1.write(curr._bestScore._bulges + "\t");
			_out1.write(curr._bestScore._bulges_miR + "\t");
			
			String delim = "\t";//"\t"
			_out1.write(curr._duplex_1D + delim);
			_out1.write(curr._alignment[3]+delim);
			_out1.write(curr._alignment[2]+delim);
			_out1.write(curr._alignment[1]+delim);
			_out1.write(curr._alignment[0]+delim);
			
			/*
			if(_runMode==1)
				_out1.write(curr.getSeedType() + "\n");
			else
				_out1.write(curr.getNonSeedType() + "\n");
			*/
			_out1.write(curr._type + "\t");
			
			_out1.write(curr._fullDuplex._duplex_1D + delim);
			_out1.write(curr._fullDuplex._alignemntStraight[3] +delim);
			_out1.write(curr._fullDuplex._alignemntStraight[2] +delim);
			_out1.write(curr._fullDuplex._alignemntStraight[1] +delim);
			_out1.write(curr._fullDuplex._alignemntStraight[0] +delim);
			
			boolean contained = curr._fullDuplex._duplex_1D.substring(1).startsWith(curr._duplex_1D);
			_out1.write(contained+"\n");
			
			//_out1.write(curr._bestScore.toString()+"\n");
			//_out1.write(curr._bestScore.currIndex_i+"\n");
			
		}
		
		//_out2.write(gene._orderInList + "\t" + geneName + "\n");
	}
	
	public void writeToHTMLStart(BufferedWriter out, miRNA miRNA)throws Exception
	{
		String style= "<style type=\"text/css\"> .INFO1{color:blue;} .INFO2{color:green;} .INFO3{color:black;} .INFO4{color:red;}</style>";
		String html_first = "<HTML> <head>" + style+ "</head> <body>\n";
		
		html_first = html_first + "miRNA=" + miRNA._name + " (" + miRNA._sequence + ")" + "<br><br>";
        
		html_first = html_first + "1D_string_key: <br>  1-match <br> 2-GU <br> 3-mismatch <br>" +
		                          "4-mRNA_bulge_G <br>  5-mRNA_bulge_other <br> 6-miRNA_bulge <br><br>";
		
		
		html_first = html_first + "SeedType: <br>  1 : (2-7 full match) <br> 2 : (2-7 + 1mRNA bulge on 5-7) or (2-8 with 1GU/MM 5-8) <br> 3: (2-7 + 1mRNA bulge on 2-4) or (2-8 with 1GU/MM 2-4) <br>" +
		                          "NonSeedType: <br> 1 : match 11-13/12-14/13-16 <br> 2 : match+GU 11-13/12-14/13-16 <br> -1 : other <br>";	

		
		html_first = html_first + "<TABLE BORDER=1 cellpadding=10 align=left> \n";
		html_first = html_first + "<CAPTION>" + "<b> <span class=INFO1> let-7 Targeting Results </span> </b>" +" </CAPTION> <br>\n";
			
			
		html_first = html_first + "<tr> \n" +
		                          "<th>index</th>"+
						  		  "<th>gene </th> \n"+
		                          "<th>RiboSeqInfo</th>" +
						  			  "<th>Seed loc</th> \n" + 
						  			  //"<th>Seed Alignment</th> \n" +
						  			  "<th>Full alignment</th>" + 
						  			  "<th>Seed/NonSeed Type</th>" + "\n";
		
		out.write(html_first);
	}
	
	public void writeToHTMLEnd(BufferedWriter out)throws Exception{
		String html_last1 = "</TABLE> \n ";
		String html_last2 = "</body></HTML> ";
		out.write(html_last1 + html_last2);
	}
	
	static int runningIndex=0;
	static int runningIndex_all=0;
	
	public static void printDuplexesToHTML(Gene gene, String geneName,ArrayList<DuplexInfo> duplexes, BufferedWriter curOut) throws Exception
	{
		
		DuplexInfo curr;
		
		if(duplexes.size()>0)
			runningIndex++;
		
		runningIndex_all++;
		
		for(int i=0;i<duplexes.size();i++)
		{
			curr = duplexes.get(i);
			//for the first deplex 
			String ans = "";
			if(i==0) {
				ans = "<tr><th>" + runningIndex + "</th>" +  
				          "<th>" + curr._gene.getInfoForHTML() + "(" + curr._gene._orderInList +")" + "<br>" + "</th>" +
				          "<th>" + "<pre>" + getColoredPval(curr._gene._log2FC,curr._gene._pval,curr._gene._padj) + "</pre></th>"; 
			}
			else {
				ans = "<tr><th>" + "" + "</th>" +  
				           "<th>" + "" + "<br>"+ "</th>" + 
				           "<th>" + "" + "<br>"+ "</th>";
			}
			
			//indexes start with 0
			ans = ans  +     "<th>" + "[seedLoc = " + curr._startOnTarget + "," + curr._endOnTarget + "]"+ "<br>" + 
					                 "Seq=" + curr._fullDuplex._targetPart + "</th> \n"+
					                 "</th> \n" + 
                             
				             "<th ALIGN=left><pre>" +
				                            //getAlignmentNoColoring(curr._fullDuplex._alignment) + "<br><br>" + 
				                            getAlignmentColoring(curr._fullDuplex,curr._fullDuplex._alignemntStraight,curr._fullDuplex._duplex_1D) + "</pre></th> \n" +
                             "<th>" + "SeedType = " + curr._fullDuplex._seedTypeNew + "<br>" +
				                      "NonSeedType = " + curr._fullDuplex._nonSeedType + "</th> \n";
  				          
			
			//if(curr._fullDuplex._seedTypeNew!=-1)
			curOut.write("<td ALIGN=left>" +ans + "</td>\n");
  				             			
		}
		
	}

	

	
	private static String getAlignmentColoring(DuplexInfo fullDuplex, String[] alignment, String str1d) {//color int the alignment the miRNA seed region and the non seed region
		//System.out.println(alignment);
		//miRNA located in the line 0 and 1
		int start1=0,end1=0,start2=0,end2=0;
		int curr=-1;
		for(int i=0;i<Math.max(alignment[0].length(),alignment[1].length());i++) {
			if(alignment[0].charAt(i)!=' ' || alignment[1].charAt(i)!=' ')
			{
				curr++;
				if(curr==1)
					start1=i;
				else if(curr==6) //depends if it is perfect or not
					end1=i;
				else if (curr==10)
					start2=i;
				else if(curr==15)
					end2=i;
			}
			else
			{
				//System.out.println("in else");
			}
		}
		
		start1 = fullDuplex._seedStart;
		end1 = fullDuplex._seedEnd;
		
		String ans = "         " + alignment[3] + "<br>" + 
		             "mRNA  3' " + alignment[2] + "<br>" + 
		             "miRNA 5' " + alignment[1].substring(0, start1) + 
		            "<span class=INFO1>" + alignment[1].substring(start1, end1+1) + "</span>" + 
		                                   alignment[1].substring(end1+1, start2) + 
		            "<span class=INFO2>" + alignment[1].substring(start2, end2+1) + "</span>" +
		                                   alignment[1].substring(end2+1) + "<br>" +
		            
					"         " + alignment[0].substring(0, start1) + 
		            "<span class=INFO1>" + alignment[0].substring(start1, end1+1) + "</span>" + 
		                                   alignment[0].substring(end1+1, start2) + 
		            "<span class=INFO2>" + alignment[0].substring(start2, end2+1) + "</span>" +
		                                   alignment[0].substring(end2+1) + "<br><br>" + 
		            
                    "         " + str1d.substring(0, start1) + 
					"<span class=INFO1>" + str1d.substring(start1, end1+1) + "</span>" + 
					                       str1d.substring(end1+1, start2) + 
					"<span class=INFO2>" + str1d.substring(start2, end2+1) + "</span>" +
					                       str1d.substring(end2+1) + "<br>";
			
		return ans;
	}
	
	private static String getColoredPval(Double log2FC, Double pval, Double padj) {
		if(log2FC==null)
			return "";
		String ans = "";
		if (log2FC>=0.58 && padj<0.1)
			ans = "<span class=INFO4>" + "log2FC=" + round(log2FC,3) + "\n" + "padj=" + roundPadj(padj) + "\n" + "over-expressed" + "</span>";
		else
			ans = "log2FC=" + round(log2FC,3) + "\n" + "padj=" + roundPadj(padj);
		return ans;
	}
	
	
	 private static Double roundPadj(Double padj) {
		if(padj<0.001)
			return padj;
		else
			return round(padj,3);
	}

	public static double round(double d, int decimalPlace)
	 {
			BigDecimal bd = new BigDecimal(Double.toString(d));
		    bd = bd.setScale(decimalPlace,BigDecimal.ROUND_HALF_UP);
		    return bd.doubleValue();
	 }

}
