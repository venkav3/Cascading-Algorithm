// invoked with nodes[] array populated according to the graph tested, likewise appropriately set value of numNodes
// when populating nodes[] ignore all data elements, but populate degAll and nBors[], these are the only data required to determine components
public class Component {
	static public int Infty = 100000;
	static public int maxNodes = 500;
	
	
	public static int find_components(int numNodes, Node nodes[])
	{
		int label = 0;
		for(int i=0;i<numNodes;++i)
			{	// System.out.println(" i "+i);
			    nodes[i].component = -1; // initialize with dummy value
				nodes[i].depth = Infty;
			}
		for(int i = 0;i<numNodes;++i)
			if(nodes[i].component == -1) // begin depth-first-search from this node
			{
				nodes[i].component = label;
				nodes[i].depth = 0;
				DFS(numNodes, nodes, i, label);
				++label;
			}
		int numComponents = -1;
		for(int i=0;i<numNodes;++i)
			if(nodes[i].component > numComponents)
				numComponents = nodes[i].component;
		//for(int i =0;i<numNodes;++i)
			//System.out.println("node "+i+" comp "+nodes[i].component);
		return(numComponents+1);
		
		
	}// end_find_components()
	//-------------------------------------------------------------------------------------------------------------------------
	public static void DFS(int numNodes, Node nodes[], int s, int label) // depth-first-search from node s
	{
		for(int i = 0; i < nodes[s].degAll;++i)
		{
			int node2 = nodes[s].nBors[i];
			if(nodes[node2].component == -1)
			{
				nodes[node2].component = label;
				nodes[node2].depth = nodes[s].depth+1;
				DFS(numNodes, nodes, node2,label);
			}
		}
		
	}//end_DFS()
	//----------------------------------------------------------------------------------------------------------------------
	public static Node[] loadCompTesting(int numNodes2, int tempA[], Node nodes[])
	{
		// relevant nodes to be in integer array tempA[].  Use nodes[] to determine what is connected to what.  
		// Populate nodes2[] accordingly and return 
		
		Node nodes2[] = new Node[maxNodes];
		for(int i=0;i<numNodes2;++i)
		{	nodes2[i] = new Node();		
			int node = tempA[i];
			nodes2[i].weight = nodes[node].weight;
			int next = 0;
			for(int j=0;j< nodes[node].degAll;++j)
			{
				int node2 = nodes[node].nBors[j];
				for(int k=0;k<numNodes2;++k)
					if(tempA[k]==node2)
					{
						nodes2[i].nBors[next] = k;
						++next;
						nodes2[i].degAll = next;
						break;
					}
				
			}
		}
		
		return(nodes2);
		
	}/*end_loadCompTesting()*/
	//------------------------------------------------------------------------------------------------
	public static Node[] loadCompTesting2(int numNodes2, int tempA[], Node nodes[], Arc arcs[])
	{
		// relevant nodes to be in integer array tempA[].  Use nodes[] to determine what is connected to what.  
		// Populate nodes2[] accordingly and return 
		
		Node nodes2[] = new Node[maxNodes];
		for(int i=0;i<numNodes2;++i)
		{	nodes2[i] = new Node();		
			int node = tempA[i];
			nodes2[i].weight = nodes[node].weight;  // note: node2[] organized not in terms of original node numbers - those are lost in node2[]
			int next = 0;
			for(int j=0;j< nodes[node].degAll;++j)
			{
				int node2 = nodes[node].nBors[j];
				int ilink = nodes[node].incidentArcs[j];
				if(arcs[ilink].type == 1)  // inter-network link
				{
					for(int k=0;k<numNodes2;++k)
						if(tempA[k]==node2)
						{
							nodes2[i].nBors[next] = k;
							++next;
							nodes2[i].degAll = next;
							break;
						}
				}
				
			}
		}
		
		return(nodes2);
		
	}/*end_loadCompTesting()*/
	//------------------------------------------------------------------------------------------------

}
