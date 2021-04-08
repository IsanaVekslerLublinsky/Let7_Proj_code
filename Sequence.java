package BasicComponents;

public class Sequence
{
	public String _header,_sequence;
	public int _orderInList;
	public miRNA _miRNA;
	
	public Sequence(String header,String sequence, int counter)
	{
		_header = header;
		_sequence = sequence;
		_orderInList = counter;
	}
	
	public void addMiRNA(miRNA miR)
	{
		_miRNA = miR;
	}
}
