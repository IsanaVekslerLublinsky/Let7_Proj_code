package BasicComponents;

import java.util.ArrayList;

public class miRNA
{
	
	public String _name, _accession;
	public String _sequence;
	public int _orderInList;
	
	//For the parsing alignments stage
	public int numOfTargetGenesBefore=0;
	public int numOfTargetGenesAfter=0;

	public ArrayList<Gene> _targetGenes;
	public ArrayList<String> _uniqueTargetGenes;

	public Double _bestEnergy;
	public int _type; //1-hairpin, 2-mature 
	public ArrayList<Sequence> _matchingHairpins;
	public ArrayList<Sequence> _matchingMatures;
	public int _matching_5p_3p; //1-5p,2-3p, 3-both
	
	public miRNA(String name, String accession, String sequence, int orderInList)
	{
		this._name = name;
		this._accession = accession;
		this._sequence = sequence;
		this._orderInList = orderInList;
	}

	public miRNA(Sequence sequence) 
	{
		this._name = sequence._header;
		this._sequence = sequence._sequence;
		this._orderInList = sequence._orderInList;
	}



	public String toString()
	{
		String ans = ">" + _name + "\n" + _sequence;	
		return ans;
	}
	
	public String getName()
	{
		return _name;
	}
	
	public void addTargetGene(Gene gene) 
	{

		if(_targetGenes==null)
			_targetGenes = new ArrayList<Gene>();

		if(_targetGenes.contains(gene))
			_targetGenes.add(gene);
	}
	
	public String getSequence()
	{
		return _sequence;
	}
}
