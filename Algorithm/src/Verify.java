// verify solution 
// 1) remaining network a connected component in each case
// 2) end-state condition satisfied
public class Verify {
	
	static public int maxNodes = 500;
	static public double EPS = 0.001;
	
	//---------------------------------------------------------------------------------------
	public static int verifyConnectivity(double solnX[], int nsize1, int nsize2, int nsize3, Node nodes[])
	{
		int[] tempA = new int[maxNodes];
		Node[] nodes2 = new Node[maxNodes];
		int ncomp;
		int next;
		int flag = 1;
		
		next = 0;
		for(int i=0;i<nsize1;++i)
			if(Math.abs(solnX[i]-0.0) < EPS)
			  {
				tempA[next] = i;
			    ++next;
			  }
		nodes2 = Component.loadCompTesting(next,tempA,nodes);
		if((ncomp = Component.find_components(next, nodes2)) != 1)
		{
		  System.out.println(" remanant of network1 not connected; no. of components =  "+ncomp);
		  flag = 0;
		}
			
		next = 0;
		for(int i=nsize1;i<nsize1+nsize2;++i)
			if(Math.abs(solnX[i]-0.0) < EPS)
			  {
				tempA[next] = i;
			    ++next;
			  }
		nodes2 = Component.loadCompTesting(next,tempA,nodes);
		if((ncomp = Component.find_components(next, nodes2)) != 1)
		{
		  System.out.println(" remenant of network2 not connected; no. of components =  "+ncomp);
		  flag = 0;
		}
		
		//System.out.println(" network3 ");
		next = 0;
		for(int i=nsize1+nsize2;i<nsize1+nsize2+nsize3;++i)
			if(Math.abs(solnX[i]-0.0) < EPS)
			  {
				tempA[next] = i;
				//System.out.println(" i "+next+" node "+i);
			    ++next;
			  }
		nodes2 = Component.loadCompTesting(next,tempA,nodes);
		if((ncomp = Component.find_components(next, nodes2)) != 1)
		{
		  System.out.println(" network3 not connected; no. of components =  "+ncomp);
		  flag = 0;
		}
			
		
		return(flag);
		
	}// end_find_components()
	//-------------------------------------------------------------------------------------------------------
	public static int verifyEndState(int numNodes, int numArcs, double solnX[], Node nodes[], Arc arcs[])
	{
		int flag = 1;
		
		for(int i=0;i<numNodes;++i)
			if(Math.abs(solnX[i]-1.0) < EPS)  // disabled node
			{
				for(int j=0;j<nodes[i].degAll;++j)  // go thu neighbors on inter-network links
				{
					int ilink = nodes[i].incidentArcs[j];
					if(arcs[ilink].type == 1)  // inter-network link
					{
						int node2 = nodes[i].nBors[j];
						if(Math.abs(solnX[node2]-0.0) < EPS)
						{
							// violates end-state condition
							flag = 0;
							break;
						}
					}
				}
				if(flag == 0)
					break;
			}
		
		return(flag);	

	}//end_verifyEndState()
	//------------------------------------------------------------------------------------------------------

}//end-class
