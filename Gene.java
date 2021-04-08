package BasicComponents;

import java.util.ArrayList;
import java.util.StringTokenizer;

import Targeting.DuplexInfo;

public class Gene 
{
	
	public String _nameLine,_sequence,_name;
	public int _orderInList;
	
	String _ensembl_gene_id="";
	String _ensembl_transcript_id="";
	String _UTR_start,_UTR_end ;//
	String _strand; //-1 or 1
	public String _exterName="";
	public String _function ="";
	public boolean _collected = false;

	public ArrayList<miRNA> _mirs = new ArrayList<miRNA>();
	public int totalNumOfHits =0;
	public int totalOfHumanTargetingMirs = 0;

	public boolean _contained = false; //for duplicates analysis
	public int _length;
	
	public ArrayList<String> _names;
	public boolean _probeFound = false;
	
	
	//////// New 2019 
	String _gene_id="";
	String _transcript_id="";
	public String _geneName = "";
	
	//New 2020 - Riboseq data
	public Double _padj,_log2FC, _pval;
	
	public ArrayList<DuplexInfo> _duplexes;
	
	public Gene(String nameLine, String sequence, int orderInList)
	{
		this._nameLine = nameLine;
		this._sequence = sequence;
		this._orderInList = orderInList;
	}

	public Gene(String name, int length, int orderInList)//for full gene
	{
		this._nameLine = name;
		_orderInList = orderInList;
		_sequence = "";
		_length = length;
	}
	
	//option 1 - WB_parasite (>caenorhabditis_elegans_prjna13758|WBGene00007080|AH6.2.1|sfxn-1.1)
	public Gene(String name, String sequence, int orderInList,int option)//for full gene
	{
		this._nameLine = name;
		this._orderInList = orderInList;
		this._sequence = sequence;
		
		parseNamingInfo(option);
	}
	
	private void parseNamingInfo(int option) {
		if(option==1) {
			StringTokenizer st = new StringTokenizer(_nameLine,"|");
			st.nextToken();
			_gene_id = st.nextToken();
			_transcript_id = st.nextToken();
			_geneName = st.nextToken();
		}
		
	}

	public String toString()
	{
		//String ans = _nameLine + "\n" + _sequence;
		String ans = _nameLine + "\n 3'UTR length:" + _sequence.length() + "\n";	
		return ans;
	}
	
	
	//
	public void parseHeader() 
	{
		//smed
		//>v31.010110|Name=mk4.010110.00;ID=2687650|-4526..13335|1|1|-4526..5142|
		StringTokenizer st = new StringTokenizer(_nameLine,"|");
		st.nextToken();
		String temp = st.nextToken();
		_geneName = temp.substring(temp.indexOf("=")+1,temp.indexOf(";"));
	}
	
	public void parseHeaderEnsembl() 
	{
		StringTokenizer st = new StringTokenizer(_nameLine,"|");
		_ensembl_gene_id = st.nextToken();
		_ensembl_transcript_id = st.nextToken();
		_exterName = st.nextToken();
	}
	
	
	//************************
	public void setGeneID(String id)
	{
		this._gene_id = id;
	}

	public String getGeneID()
	{
		return this._gene_id;
	}

	public void setEnsemblTranscriptID(String id)
	{
		this._ensembl_transcript_id = id;
	}

	public String getEnsemblTranscriptID()
	{
		return this._ensembl_transcript_id;
	}

	public void setStartEndUTR_and_strand(String start, String end, String strand)
	{
		_strand = strand;
		this._UTR_start = start;
		this._UTR_end = end;
	}

	public String getStartUTR()
	{
		return this._UTR_start;
	}

	public String getEndUTR()
	{
		return this._UTR_end;
	}

	public String getStrand()
	{
		return _strand;
	}

	public String getSequence()
	{
		return _sequence;
	}

	public void addTargetingMir(miRNA miRNA) 
	{
		if(!_mirs.contains(miRNA))
			_mirs.add(miRNA);
	}

	public String getInfoForHTML() {
		return _gene_id + ":" + _geneName;
	}
}
