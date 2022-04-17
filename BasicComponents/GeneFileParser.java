package BasicComponents;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;

public class GeneFileParser {

	public ArrayList<Gene> _fullGenes = new ArrayList<Gene>();
	public ArrayList<Gene> _uniqueGenes = new ArrayList<Gene>();

	public GeneFileParser(String file)throws IOException , Exception
	{
		parse_genes_File(new File(file));
	}

	private void parse_genes_File(File file)throws IOException , Exception
	{
		FileReader fr = new FileReader(file);
		BufferedReader reader = new BufferedReader(fr);
		String nameLine = reader.readLine();
		String sequenceLine;

		int numberOfGene =0;
		Gene full_gene;

		while(nameLine!=null)
		{
			String fullSequence = "";
			sequenceLine = reader.readLine();
			while (sequenceLine!=null && !sequenceLine.startsWith(">"))
			{
				//System.out.println("in WHILE");
				fullSequence = fullSequence + sequenceLine;
				sequenceLine = reader.readLine();
			}

			String name = nameLine.substring(1);

			full_gene = new Gene(name, fullSequence,numberOfGene,1);
			_fullGenes.add(full_gene);
			numberOfGene++;

			nameLine = sequenceLine;
		}

		System.out.println("Done parsing genes file");
		System.out.println("All genes " + this._fullGenes.size());
		reader.close();
		fr.close();
		reader = null;
		fr = null;
		System.gc();
	}

	public ArrayList<Gene> geteFullGenes()
	{
		return this._fullGenes;
	}

	public ArrayList<Gene> geteUniqueGenes()
	{
		return this._uniqueGenes;
	}

	//binarySearch
	public ArrayList<Gene> getGenesWithID(String geneID)
	{
		//fullGenes are sorted by ENSG
		//name = ENSG00000006114(ENST00000339208)
		int low=0;
		int high = _fullGenes.size()-1;
		int mid;
		int loops = 0;
		ArrayList<Gene> answer = new ArrayList<Gene>();
        while( low <= high )
        {
        	
        	if (low+1==high)
 			 {
 				 if (_fullGenes.get(low).getGeneID().equals(geneID))
 					 return collectGenes(geneID,low);
 				 else if (_fullGenes.get(high).getGeneID().equals(geneID))
 					 return collectGenes(geneID,high);
 			 }
 			
        	
            mid = ( low + high ) / 2;
            if( _fullGenes.get(mid).getGeneID().compareTo(geneID) < 0 )
                low = mid + 1;
            else if( _fullGenes.get(mid).getGeneID().compareTo(geneID) > 0 )
                high = mid - 1;
            else
            {
            	return collectGenes(geneID, mid);
            }

            loops++;
        }

        return answer;
	}


	private ArrayList<Gene> collectGenes(String geneID, int mid) 
	{
		ArrayList<Gene> answer = new ArrayList<Gene>();
		
		//collect all sequences with the same geneID
		int startIndex = mid; 
		while(startIndex>=0 && _fullGenes.get(startIndex)._ensembl_gene_id.equals(geneID))
		{
			startIndex--;
		}
			
		startIndex++;
		
		int endIndex = mid;
		while(endIndex<_fullGenes.size() && _fullGenes.get(endIndex)._ensembl_gene_id.equals(geneID))
		{
			endIndex++;
		}
		endIndex--;
		
		for(int i=startIndex;i<=endIndex;i++)
		{
			if(_fullGenes.get(i)._collected == true)
				System.out.println("PROBLEM " + _fullGenes.get(i)._ensembl_gene_id);
			
			answer.add(_fullGenes.get(i));
			_fullGenes.get(i)._collected = true;
		}
		return answer;
		
	}

	

	//for briggsae
	public ArrayList<Gene> getGenesWithIDInHeader(String name)
	{
		ArrayList<Gene> ans = new ArrayList<Gene>();
		for(int i=0;i<_fullGenes.size();i++)
		{
			Gene gene = _fullGenes.get(i);
			//System.out.println(gene._name);
			
			if(gene._nameLine.startsWith(name) || gene._nameLine.contains("|" + name+"|") || gene._nameLine.endsWith("|" + name) || gene._nameLine.endsWith("|CELE-" + name))
			{
				//TODO - check this
				/*
				if(gene._probeFound)
					System.out.println("was found already " + name + " " + gene._name);
				*/
				ans.add(gene);
				gene._probeFound = true;
			}
			else if(ans.size()>0)
				return ans;
		}
		return ans;
	}

	//WB_id
	public ArrayList<Gene> getGenesWithIDInNames(String name)
	{
		ArrayList<Gene> ans = new ArrayList<Gene>();
		for(int i=0;i<_fullGenes.size();i++)
		{
			Gene gene = _fullGenes.get(i);
			//System.out.println(gene._names.size());
			boolean found = false;
			if(gene._nameLine.contains(name))
				found = true;
			if(!found && gene._names!=null)
			{
				for(int k=0;k<gene._names.size();k++)
				{	
					if(gene._names.get(k).equals(name))
					{
						found = true;
					}
				}
			}
			if(found)
			{
				gene._probeFound = true;
				ans.add(gene);
			}
			else if(!found && ans.size()>0) //to make sure all genes are returned
				return ans;
		}
		return ans;
	}
	
	
	public ArrayList<Gene> getGenesFromUTRsNames(ArrayList<Gene> utrNames) 
	{
		ArrayList<Gene> ans = new ArrayList<Gene>();
		for(int i=0;i<utrNames.size();i++)
		{
			ans.add(_fullGenes.get(utrNames.get(i)._orderInList));
		}
		return ans;
	}
	
	
	//this part is needed in DuplexFileReader in order to print let-7 targeting info
	public ArrayList<ArrayList<Integer>> combinedGenes;
	static ArrayList<String> combinedNames;
	
	public void combineIsoforms()
	{
		combinedGenes = new ArrayList<ArrayList<Integer>>();
		combinedNames = new ArrayList<String>();
		
		for(int i=0;i<_fullGenes.size();i++)
		{
			String geneName = _fullGenes.get(i)._nameLine;
			StringTokenizer stName = new StringTokenizer(geneName,"|");
			String name = stName.nextToken();
			int size = combinedNames.size();
				
			if(size>0 && name.equals(combinedNames.get(size-1)))
			{
				//need to add to the same gene
				combinedGenes.get(combinedGenes.size()-1).add(i);
			}
			else
			{
				ArrayList<Integer> newArr = new ArrayList<Integer>();
				newArr.add(i);
				combinedGenes.add(newArr);
				combinedNames.add(name);
			}		
		}
	}

	//new July 2020
	public void annotateGenesWithInfo(String riboSeqData) throws Exception{
		String outFile = Parameters.outFilePwd + "RiboSeqTable_statistics";
		BufferedWriter out = new BufferedWriter(new FileWriter(outFile));
		BufferedReader reader = new BufferedReader(new FileReader(riboSeqData));
		String header = reader.readLine(); //"gene.name","Geneid","ribo.log2FC","ribo.padj"
		String line = reader.readLine();
		while(line!=null) {
			
			StringTokenizer st = new StringTokenizer(line,",");
			String name = st.nextToken();
			String id = st.nextToken();
			String temp = st.nextToken();
			
			Double log2FC=0.1,padj=1.0;
			
			boolean proceed = true;
			if(temp.equals("NA"))
				proceed = false;
			else
				log2FC = Double.parseDouble(temp);
			
			temp = st.nextToken();
			if(temp.equals("NA"))
				proceed=false;
			else
				padj = Double.parseDouble(temp);
			
			//System.out.println(log2FC + " " + padj + " " + proceed);
			//if(proceed) {
				ArrayList<Gene> genes = getGenesWithName(name);
				out.write(name + " " + log2FC + " " + padj + " " + genes.size()+"\n");
				
				for(int i=0;i<genes.size();i++) {
					Gene cur = genes.get(i);
					cur._log2FC = log2FC;
					//cur._pval = pval; 
					cur._padj = padj;
				}
			//}
			line=reader.readLine();
		}
		reader.close();
		out.close();
		
	}

	private ArrayList<Gene> getGenesWithName(String name) {
		ArrayList<Gene> genes = new ArrayList<Gene>();
		for(int i=0;i<_fullGenes.size();i++) {
			Gene cur = _fullGenes.get(i); 
			if(cur._geneName.equals(name))
				genes.add(cur);
		}
		return genes;
	}
}
