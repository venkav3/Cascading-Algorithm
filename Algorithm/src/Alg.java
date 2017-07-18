// Solves the interdependent networks CNA Problem directly as a minimization problem (with quadratic constraints)
//import ilog.concert.IloException;
// full optimization
import ilog.concert.IloException;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.*;




public class Alg {

	static final double EPS1 = 0.001;
	
	static double lambda1, lambda2, lambda3;
	static public int maxArcs = 1000;
	static public int maxNodes = 500;
	static final int maxDegree = 20;
	static public int MAX = 60000;
	
	static public int nsize2;
	static public int nsize3;
	static public int nsize1;
	static public int numArcs;
	static public int numNodes;
	
	static Arc arcs[] = new Arc[maxArcs];
	static Node nodes[] = new Node[maxNodes];
	
	
		
	/*------------------------------------------------------------------------------*/
	public static void main(String[] args) throws IloException, IOException  {
		
	    		
         
			
			Alg driver = new Alg();

			long startTimeMs = System.currentTimeMillis( );
			
			lambda1 = Double.parseDouble(args[0]);
			lambda2 = Double.parseDouble(args[1]);
			lambda3 = Double.parseDouble(args[2]);
			
			System.out.println("\n Demand attenuation factor, lambda1, lambda2, lambda3: "
					+ args[0]+", "+args[1]+", "+args[2]);
			//driver.read_input(args[3],args[4],nsize1,nsize2,nsize3,numArcs,nodes,arcs);
			try {
				
				driver.read_input(args[3],args[4]);
				System.out.println("nsize1 "+nsize1+" nsize2 "+nsize2+" nsize3 "+nsize3);
				//Diagnostic.dumpNodes(nsize1,nsize2,nsize3,numArcs,nodes,arcs);
				} catch (NumberFormatException | IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				}
			double[] x = new double[MAX];
			x = Optimize.optimize(lambda1, lambda2, lambda3, nsize1, nsize2, nsize3, numNodes, numArcs, nodes, arcs);
			// count # of intra-network links
			int nIntra = 0;
			for(int i=0;i<numArcs;++i)
				if(arcs[i].type == 0)
					++nIntra;
			
			// count # of inter-network links
			int nInter = 0;
			for(int i=0;i<numArcs;++i)
				if(arcs[i].type == 1)
					++nInter;
			
			int nvars = (int) x[0];
			int ncons = (int) x[1];
			
			
			double[] soln = new double[MAX];
			for(int i = 0;i<2*numNodes;++i)  // pass only x-values and attacked nodes info; not other ancillary variables
				soln[i] = x[2+i];
			
			double[] slacks = new double[MAX];
			
			/******************************  
			for(int i=0;i<ncons;++i)
				slacks[i] = x[2+nvars+i];
				***********************************/
			Diagnostic.output(nvars,ncons,soln,slacks,nsize1,nsize2,nsize3,numNodes,numArcs,nodes,arcs, startTimeMs);
		
	}/* end_main() */

	/*----------------------------------------------------------------------------------------*/
public void read_input(String path1, String path2) throws NumberFormatException,
			IOException {
		DataInputStream istream = new DataInputStream(new FileInputStream(path1));
		BufferedReader buffer = new BufferedReader(new InputStreamReader(
				istream));
		String str;

		// String field[];
		str = buffer.readLine(); /* first line in Bienstock data */
		String field[] = str.split("\t");
		nsize1 = Integer.parseInt(field[0]);
		nsize2 = Integer.parseInt(field[1]);
		nsize3 = Integer.parseInt(field[2]);
		numArcs = Integer.parseInt(field[3]);
		numNodes = nsize1+nsize2+nsize3;
		// System.out.println("nodes "+numNodes+" arcs "+numArcs+" generators "+numGenerators+" demands "+numDemands);
		if (numNodes > maxNodes)
			Diagnostic.verror(1, numNodes, maxNodes);
		if (numArcs > maxArcs)
			Diagnostic.verror(2, numArcs, maxArcs);
		

		for (int i = 0; i < numNodes; ++i)
		{
			nodes[i] = new Node();
			if(i < nsize1)
				nodes[i].network = 0;
			else if (i < (nsize1+nsize2))
				nodes[i].network = 1;
			else
				nodes[i].network = 2;
		}
		for (int i = 0; i < numArcs; ++i)
			arcs[i] = new Arc();
		
		int seqIntra = 0;
		int seqInter = 0;
		for (int i = 0; i < numArcs; ++i) /* read arc info */
		{
			str = buffer.readLine();

			int node1, node2;
			field = str.split("\t");
			arcs[i].from = node1 = Integer.parseInt(field[0]);
			arcs[i].to = node2 = Integer.parseInt(field[1]);
			if(nodes[node1].network == nodes[node2].network)
			{
				arcs[i].type = 0;
				arcs[i].seqIntra = seqIntra;
				++seqIntra;
			}
			else
			{
				arcs[i].type = 1;
				arcs[i].seqInter = seqInter;
				++seqInter;
			}
			load(node1,i,nodes,arcs);
			load(node2,i,nodes,arcs);
		}
		System.out.println("numNodes "+numNodes+" numArcs "+numArcs);
		//Diagnostic.dumpNodes( nsize1,  nsize2,  nsize3,  numArcs,  nodes,  arcs);
		//Diagnostic.dumpArcs( numArcs,  arcs);
		
		// verify if each input network is connected
		
		int[] tempA = new int[maxNodes];
		Node[] nodes2 = new Node[maxNodes];
		int ncomp;
		
		for(int i=0;i<nsize1;++i)
			tempA[i] = i;
		nodes2 = Component.loadCompTesting(nsize1,tempA,nodes);
		if((ncomp = Component.find_components(nsize1, nodes2)) != 1)
		{
			System.out.println(" network1 not connected; no. of components =  "+ncomp);
			System.exit(0);
		}
		for(int i=nsize1;i<nsize1+nsize2;++i)
			tempA[i] = i;
		nodes2 = Component.loadCompTesting(nsize2,tempA,nodes);
		if((ncomp = Component.find_components(nsize2, nodes2)) != 1)
		{
			System.out.println(" network2 not connected; no. of components =  "+ncomp);
			System.exit(0);
		}
		for(int i=nsize1+nsize2;i<numNodes;++i)
			tempA[i] = i;
		nodes2 = Component.loadCompTesting(nsize3,tempA,nodes);
		if((ncomp = Component.find_components(nsize3, nodes2)) != 1)
		{
			System.out.println(" network3 not connected; no. of components =  "+ncomp);
			System.exit(0);
		}
		System.out.println("input networks tested and found to be connected.");
		// read node weights
		
		istream = new DataInputStream(new FileInputStream(path2));
	    buffer = new BufferedReader(new InputStreamReader(istream));
		
	    Random r = new Random();
	    
        for(int i=0;i<numNodes;++i)
        {
		  str = buffer.readLine(); 
		  field = str.split("\t");
		  nodes[i].weight = Double.parseDouble(field[0]);
		  // perturb to make weights unique
		  /*************************
		  double x;
			  do{
				  x = r.nextGaussian()*EPS1;  
			  }while(Math.abs(x) > 3* EPS1);
		  nodes[i].weight += x;
		  *****************************/
		}
		
	} /* end_read_input() */
	/*----------------------------------------------------------------------------------------------*/
	public void load(int inode, int ilink, Node nodes[], Arc arcs[])
	{
		/* load ilink into set of incident links on inode
		 * load 
		 */
		int ind = nodes[inode].degAll;
		if(ind >= maxDegree)
			Diagnostic.verror(3, inode, maxDegree);
		// determine neighbor from ilink appropriately
		int nbor = arcs[ilink].from;
		if(nbor == inode)
			nbor = arcs[ilink].to;
		// check if this is a repeat link
		for(int i = 0; i<ind;++i)
			if(nodes[inode].nBors[i]==nbor)
				Diagnostic.verror(4, inode, nbor);
		nodes[inode].nBors[ind] = nbor;
		nodes[inode].incidentArcs[ind] = ilink;
		++nodes[inode].degAll;
		if(arcs[ilink].type == 0)
			++nodes[inode].degIntra;
				
	}/*end_load()*/
	/*----------------------------------------------------------------------------------------------*/
	
	
}/* end_class */

