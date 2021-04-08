package BasicComponents;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;

import Targeting.RNAduplex;

public class miRNAFileParser
{

	public ArrayList<miRNA> _miRNAs;

	public miRNAFileParser(String file)throws IOException , Exception
	{
		parse_miRNA_File(new File(file));
	}
	
	public miRNAFileParser(File file)throws IOException , Exception
	{
		parse_miRNA_File(file);
	}

	public void parse_miRNA_File(File file)throws IOException , Exception
	{
		_miRNAs = new ArrayList<miRNA>();

		FileReader fr = new FileReader(file);
		BufferedReader reader = new BufferedReader(fr);
		String nameLine = reader.readLine();
		String sequenceLine;

		int numberOf_miRNA = 0;
		while(nameLine!=null)
		{
			sequenceLine = reader.readLine();
			nameLine = nameLine.substring(1);
			StringTokenizer st = new StringTokenizer(nameLine);
			String name = st.nextToken();
			String accession = "";
			if(st.hasMoreTokens())
				accession = st.nextToken();

			//System.out.println(name);

			miRNA mir = new miRNA(name, accession, sequenceLine, numberOf_miRNA);
			_miRNAs.add(mir);

			numberOf_miRNA++;
			nameLine = reader.readLine();
		}
		reader.close();
		fr.close();
		reader = null;
		fr = null;
		System.gc();

		
		System.out.println("Done parsing miRNAs file " + _miRNAs.size());
	}

	public void computeBestEnergies()throws Exception
	{
		RNAduplex cofold = new RNAduplex();
		for(int i =0; i<_miRNAs.size();i++)
		{
			miRNA mir = _miRNAs.get(i);
			cofold.call_RNAcofold(mir.getSequence(),UsefulFunctions.reverse(UsefulFunctions.getComplement(mir.getSequence())));
			mir._bestEnergy = cofold._energy_duplex;
		}
	}
	
	public ArrayList<miRNA> get_miRNAs()
	{
		return _miRNAs;
	}
	
	public miRNA getMiRNAWithName(String name)
	{
		for(int i=0;i<_miRNAs.size();i++)
		{
			if(_miRNAs.get(i)._name.equals(name) || _miRNAs.get(i)._name.contains(name+"("))
			{
				return _miRNAs.get(i);
			}
		}
		return null;
	}
}
