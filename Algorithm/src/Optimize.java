import java.io.IOException;
import java.text.DecimalFormat;

import ilog.concert.*;
import ilog.cplex.IloCplex;
public class Optimize {	
	/* using sparse matrix techniques of filling constraint matrix */
  static DecimalFormat decfor = new DecimalFormat("0.0000");
  static DecimalFormat decfor2 = new DecimalFormat("0.00");
  static int numcols;
  static final double EPS3 = 0.1;
  static final double MARGIN = 0.3;  // to enforce degradation is withing stipulated levels
  static public int MAX = 60000;
  static Arc mstArcs[] = new Arc[1000];
			
	public static double[] optimize(double lambda1, double lambda2, double lambda3, int nsize1, int nsize2, int nsize3, int numNodes, int numArcs, Node nodes[], Arc arcs[])throws IloException, IOException
	{
		IloCplex cplex = new IloCplex();
		cplex.setOut(null);

		IloNumVar[][] var = new IloNumVar[1][];
		IloRange[][] rng = new IloRange[1][];
			
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
		System.out.println(" nInter "+nInter+" nIntra "+nIntra);
		// compute # of variables including those required for Martins' MST formulation
		int numNodes2 = numNodes+1;
		int numArcs2 = nInter+numNodes;
		numcols = 6*numNodes+numArcs2+(2*numArcs2*numNodes2)+2*nIntra; 
		System.out.println(" numcol (# of vars. in the problem) "+numcols);
		
		double[] x = new double[numcols];
		double[] y = new double[MAX];
		
		Optimize.ByRow2(numcols, cplex, var, rng, lambda1, lambda2, lambda3, nsize1, nsize2, nsize3, numNodes, numArcs, nodes, arcs);
		
		if (cplex.solve()) {
			System.out.print("cplex.solve completed \n");
			System.out.println(" Solution status = " + cplex.getStatus());
			
			double objVal = cplex.getObjValue();
			
			System.out.println("objval achieved  = " + decfor.format(objVal));
			
			/*********
			System.out.println("no. of binary vars "+cplex.getNbinVars());
			System.out.println("no. of integer vars "+cplex.getNintVars());
			System.out.println("no. of semi continuous vars "+cplex.getNsemiContVars());
			System.out.println("no. of vars "+cplex.getNcols());
			*********/
			
			x = cplex.getValues(var[0]);
			
			
			
			// verify if remaining networks are connected
			if(Verify.verifyConnectivity(x, nsize1, nsize2, nsize3, nodes) == 1)
				System.out.println("remaining networks verified for connectivity");
			
			//verify end-state condition met
			if(Verify.verifyEndState(numNodes, numArcs, x, nodes, arcs) == 1)
				System.out.println("end-state condition satisfied");
			
			int offset = 6*numNodes+nInter;
			
			
			/*************
			for(int i = 0;i<numNodes;++i)
			{if(Math.abs(x[offset+i]-1.0) < 0.001)
				System.out.println("z_uv for "+i+" = "+x[offset+i]+" x = "+x[i]+" alpha "+x[3*numNodes+i]+" weight "+nodes[i].weight);
			}
			*****************/
			
			//double[] slack = cplex.getSlacks(rng[0]);
			int nvars = cplex.getNcols();
			int nlinear = cplex.getNrows();
			int nqcs = cplex.getNQCs();
			int ncons = nlinear+nqcs;
			//System.out.println("nvar "+nvars+" nlinear "+nlinear+" nqcs "+nqcs+" ncons "+ncons);
			
			y[0] = nvars;
			y[1] = ncons;
			for(int i = 0;i<numNodes;++i)  // RETRIEVE ONLY x-values AND attacked node info; not other large number of ancillary variables
				y[2+i] = x[i];
			for(int i = 0;i<numNodes;++i)
			{
				if(Math.abs(x[3*numNodes+i]) < 0.001)
				 y[2+numNodes+i] = x[offset+i];
				else
					y[2+numNodes+i] = 0.0;
			}
			
			/****************
			for(int i = 0;i<ncons;++i)
				y[2+nvars+i]= slack[i];
			*********************************/
		}
		else
			{System.out.println("cplex didn't solve ");System.out.println("Solution status = " + cplex.getStatus());
			System.out.println("Attack on any node will break-up all 3 networks through cascading. ");
			System.exit(0);}
		cplex.end();
		return(y);
	}/*end_solveWA()*/
	/*---------------------------------------------------------------------------------------------------------------------------*/
	static void ByRow2(int numvar, IloMPModeler model, IloNumVar[][] var, IloRange[][] rng, double lambda1, double lambda2, double lambda3, int nsize1, int nsize2, int nsize3, int numNodes, int numArcs,  Node nodes[], Arc arcs[]) throws IloException {
		
		IloLPMatrix matrix = model.LPMatrix();
		
		//double[] dummy = new double[numvar];
		double[] zlb = new double[numvar];
		double[] zub = new double[numvar];
		IloNumVarType[] zt = new IloNumVarType[numvar];
		double[] objvals = new double[numvar];
		//System.out.println("numvar in ByRow2 "+numvar);
		 
		 for (int i = 0; i < numArcs; ++i)
				mstArcs[i] = new Arc();
		 
		// determine suitable value for BigM
		int BigM = 1;
		for(int i=0;i<numNodes;++i)
			if(nodes[i].degIntra > BigM) BigM = nodes[i].degIntra;
		
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
		
		// populate mstArcs for MST construction.  These inter-network arcs and arcs from the super source to each of the nodes.  
		// The number of such arcs must be numNodes2 = (n Inter+numNodes)
		int j1 = 0;
		for(int i=0;i<numArcs;++i)
			if(arcs[i].type == 1)
			{
				mstArcs[j1].from = arcs[i].from;
				mstArcs[j1].to = arcs[i].to;
				mstArcs[j1].type = 1;
				++j1;
			}
		for(int i=0;i<numNodes;++i) // make fictitious arcs from super source to each of the nodes
		{
			mstArcs[j1].from = numNodes; // the nod number of the fictitious super source
			mstArcs[j1].to = i;
			mstArcs[j1].type = 2; // new type introduced to mark fictitious arcs
			++j1;
		}
		
		int numArcs2 = numNodes+nInter;  // this captures the total # of inter-network links, including the virtual links from the super source
		int numNodes2 = numNodes+1;     /* includes fictitious super source */
		
		int xstart = 0;
		int gstart = xstart+numNodes;
		int tstart = gstart+numNodes;
		int alphstart = tstart+numNodes;
		int cstart = alphstart+numNodes;
		int sstart = cstart+numNodes;
		int zstart = sstart+numNodes;
		int ystart = zstart+numArcs2;
		int hstart = ystart+2*numArcs2*numNodes2;
		int fstart = hstart+nIntra;
			
		int next = 0;
		
		for(int i=0;i<numNodes;++i)  /* x_u */
		{
		 zlb[i] = 0;
		 zub[i] = 1;
		 zt[i] = IloNumVarType.Bool;
		 objvals[i] = 0.0;
		 ++next;
		}
		
		for(int i=0;i<numNodes;++i)  /* g_u; a complement of x_u */
		{
		 zlb[next] = 0;
		 zub[next] = 1;
		 zt[next] = IloNumVarType.Bool;
		 objvals[next] = 0.00;
		 ++next;
		}
		
		for(int i=0;i<numNodes;++i)  /* t_u */
		{
		 zlb[next] = 0;
		 zub[next] = 1;
		 zt[next] = IloNumVarType.Bool;
		 objvals[next] = 0.00;
		 ++next;
		}
		
		for(int i=0;i<numNodes;++i)  /* alpha_u */
		{
			 zlb[next] = 0;
			 zub[next] = 1;
			 zt[next] = IloNumVarType.Bool;
			 objvals[next] = -nodes[i].weight;
			 ++next;
		} 	
		
		for(int i=0;i<numNodes;++i)  /* c_u */
		{
			 zlb[next] = 0;
			 zub[next] = Integer.MAX_VALUE;;
			 zt[next] = IloNumVarType.Int;
			 objvals[next] = 0.0;
			 ++next;
		} 	
		
		for(int i=0;i<numNodes;++i)  /* s_u */
		{
			 zlb[next] = 0.0;
			 zub[next] = Double.MAX_VALUE;
			 zt[next] = IloNumVarType.Float;
			 objvals[next] = 0.00;
			 ++next;
		} 
		//System.out.println("check zstart "+zstart+" next "+next);
		for(int i=0;i<numArcs2;++i)  /* z_uv; NOTE the last numNodes of this set denotes the (v*,u) links, v* = super source */
		{
			 zlb[next] = 0;
			 zub[next] = 1;
			 zt[next] = IloNumVarType.Bool;
			 objvals[next] = 0.00;
			 if(i >= nInter)
				 objvals[next] = nodes[i-nInter].weight;
			 ++next;
		} 
		
		int lim = 2*numArcs2*numNodes2;
		for(int i=0;i<lim;++i)  /* y_{uv}^{k}  AND y_{vu}^{k}  */
		{
			 zlb[next] = 0;
			 zub[next] = 1;
			 zt[next] = IloNumVarType.Bool;
			 objvals[next] = 0.00;
			 ++next;
		} 
		
		for(int i=0;i<nIntra;++i)  /* h_uv */
		{
			 zlb[next] = -numNodes;
			 zub[next] = numNodes;
			 zt[next] = IloNumVarType.Float;
			 objvals[next] = 0.00;
			 ++next;
		} 
		
		for(int i=0;i<nIntra;++i)  /* f_uv */
		{
			 zlb[next] = 0;
			 zub[next] = 1;
			 zt[next] = IloNumVarType.Bool;
			 objvals[next] = 0.00;
			 ++next;
		} 	
		
	   	if(numvar != next)
			Diagnostic.verror(5,next,numvar);
		
		IloNumVar[] z = model.numVarArray(numvar, zlb, zub, zt);
		
		
		for(int i=0;i<numvar;++i)
		       matrix.addCols( new IloNumVar[] { z[i] }); 
		var[0] = z;
		model.addMinimize(model.scalProd(z, objvals));
		
		/* add constraints */
		
		IloLinearNumExpr lhs; 
		
		int ncons = 0;
		lhs = model.linearNumExpr();  // constraint (1) 
		for(int u = 0;u < numNodes;++u)
		  lhs.addTerm(z[xstart+u], -1.0);	
		for(int i = 0;i<numArcs2;++i)
			  lhs.addTerm(z[zstart+i], 1.0);	
		matrix.addRow(model.eq(lhs, 0.0));
		++ncons;
		
		 // constraint (2)
	    for(int i = 0;i<numArcs2;++i)  
		  {
		    int u = mstArcs[i].from;
		    int v = mstArcs[i].to;
		    for(int k = 0;k<numNodes2;++k) 
		      if((k != u) && (k != v))
		      {		
		    	lhs = model.linearNumExpr();
				lhs.addTerm(z[zstart+i], -1.0);
				int ind = encode(numNodes,numArcs,nInter,arcs,ystart,u,v,k);
				lhs.addTerm(z[ind], 1.0);
				ind = encode(numNodes,numArcs,nInter,arcs,ystart,v,u,k);
				lhs.addTerm(z[ind], 1.0);
				matrix.addRow(model.eq(lhs, 0.0));
				++ncons;
		      }
		  }
		
		// constraint (3a) next
		for(int i = 0;i<nInter;++i)  // inter-network 
		  {
		    int u = mstArcs[i].from;
		    int v = mstArcs[i].to;
		    lhs = model.linearNumExpr();
		    lhs.addTerm(z[xstart+u], -0.5);
			lhs.addTerm(z[xstart+v], -0.5);
			lhs.addTerm(z[zstart+i], 1.0);
		    for(int j = 0;j<nodes[u].degAll;++j) 
		    {
		      int k = nodes[u].nBors[j];
		      if(nodes[u].network != nodes[k].network) // so an inter-network link, a valid link in the MST problem
		    	  if(k != v)  // don't pick the same inter-arc 
		          {
	     	         int yukv = encode(numNodes,numArcs,nInter,arcs,ystart,u,k,v); // i.e., y_{uk}^v
			         lhs.addTerm(z[yukv], 1.0);
		          }
		    }
		    // next do y_{uv*}^v
		    int yukv = encode(numNodes,numArcs,nInter,arcs,ystart,u,numNodes,v);
		    lhs.addTerm(z[yukv], 1.0);
			matrix.addRow(model.eq(lhs, 0.0));
			++ncons;
		  }
				
		// constraint (3b) next
		for(int u = 0;u<numNodes;++u)  // fictitious arc from super source v* 
		  {
			lhs = model.linearNumExpr();
			lhs.addTerm(z[xstart+u], -1.0);
			lhs.addTerm(z[zstart+nInter+u], 1.0);
			for(int k=0;k<numNodes;++k)
				if(k != u)
				{
				  int yvstarku = encode(numNodes,numArcs,nInter,arcs,ystart,numNodes,k,u); // i.e., y_{uk}^v
		          lhs.addTerm(z[yvstarku], 1.0);
				}
			matrix.addRow(model.eq(lhs, 0.0));
			++ncons;
	      }	
			
			
	// constraint (4)
	for(int u = 0;u<numNodes;++u)  
	{
		lhs = model.linearNumExpr();
	    lhs.addTerm(z[sstart+u], 1.0);
		for(int j=0;j<nodes[u].degAll;++j)
		{
			int node = nodes[u].nBors[j];
			if(nodes[node].network == nodes[u].network)  //intra-link
			  lhs.addTerm(z[xstart+node], 1.0);
		}
		matrix.addRow(model.eq(lhs, nodes[u].degIntra));
		++ncons;
	}
				
	// constraint (5)
	for(int u = 0;u<numNodes;++u)  
	{
		lhs = model.linearNumExpr();
		lhs.addTerm(z[sstart+u], 1.0/BigM);
		lhs.addTerm(z[alphstart+u], 1.0);
		matrix.addRow(model.le(lhs, 1.0));
		++ncons;
	}
			
	// constraint (6)
	for(int u = 0;u<numNodes;++u)  
	{
		lhs = model.linearNumExpr();
		lhs.addTerm(z[zstart+nInter+u], -1.0); 
		lhs.addTerm(z[alphstart+u], 1.0);
		matrix.addRow(model.le(lhs, 0.0));
		++ncons;
	}
			
	// new constraint (6b); to ensure attacked nodes are nodes with x_{u} = 1
		for(int u = 0;u<numNodes;++u)  
		{
			lhs = model.linearNumExpr();
			lhs.addTerm(z[zstart+nInter+u], 1.0); 
			lhs.addTerm(z[u], -1.0);
			matrix.addRow(model.le(lhs, 0.0));
			++ncons;
		}		
			
	// constraint (7); add 3 constraints enforcing degradation objectives, (1) in report
	lhs = model.linearNumExpr();
	for(int u=0;u<nsize1;++u)
		lhs.addTerm(z[xstart+u], 1.0);
	matrix.addRow(model.ge(lhs, lambda1*nsize1));
	++ncons;
	
	lhs = model.linearNumExpr();
	for(int u=nsize1;u<nsize1+nsize2;++u)
		lhs.addTerm(z[xstart+u], 1.0);
	matrix.addRow(model.ge(lhs, lambda2*nsize2));
	++ncons;
	
	lhs = model.linearNumExpr();
	for(int u=nsize1+nsize2;u<nsize1+nsize2+nsize3;++u)
		lhs.addTerm(z[xstart+u], 1.0);
	matrix.addRow(model.ge(lhs, lambda3*nsize3));
	++ncons;
	
	// constraint (7b); we are adding these to ensure degradation not excessive, for then the cost can be driven down to zero (see report)
		lhs = model.linearNumExpr();
		for(int u=0;u<nsize1;++u)
			lhs.addTerm(z[xstart+u], 1.0);
		matrix.addRow(model.le(lhs, lambda1*nsize1*(1+MARGIN)));
		++ncons;
		
		lhs = model.linearNumExpr();
		for(int u=nsize1;u<nsize1+nsize2;++u)
			lhs.addTerm(z[xstart+u], 1.0);
		matrix.addRow(model.le(lhs, lambda2*nsize2*(1+MARGIN)));
		++ncons;
		
		lhs = model.linearNumExpr();
		for(int u=nsize1+nsize2;u<nsize1+nsize2+nsize3;++u)
			lhs.addTerm(z[xstart+u], 1.0);
		matrix.addRow(model.le(lhs, lambda3*nsize3*(1+MARGIN)));
		++ncons;
	// artificial constraints only to enable cplex to process quadratic constraints; USE g_u INSTEAD !!!! So this is now constraint (9)
	for(int i = 0;i < numNodes;++i)
	{
		lhs = model.linearNumExpr();
		lhs.addTerm(z[xstart+i], 1.0);
		lhs.addTerm(z[gstart+i], 1.0);
		matrix.addRow(model.eq(lhs, 1.0));
		++ncons;
	}
	
	// add terminating configuration constraints (8) and (8b) 
	for(int i = 0;i<numArcs;++i)	
	  if(arcs[i].type == 1)  // inter-network link
		{
			int u = arcs[i].from;
			int v = arcs[i].to;
			matrix.addRow(model.addLe(model.prod(z[xstart+u], z[gstart+v]),EPS3));
			++ncons;
			matrix.addRow(model.addLe(model.prod(z[xstart+v], z[gstart+u]),EPS3));
			++ncons;
		}	
	
	
				
	// constraint (10) 
	for(int i = 0;i < numArcs;++i)
		if(arcs[i].type == 0)  // intra-link; between nodes in same network
	    {
		  lhs = model.linearNumExpr();
		  int u = arcs[i].from;
		  int v = arcs[i].to;
		  int seq = arcs[i].seqIntra;
		  lhs.addTerm(z[fstart+seq], 1.0);
		  lhs.addTerm(z[gstart+u], -0.5);
		  lhs.addTerm(z[gstart+v], -0.5);
		  matrix.addRow(model.le(lhs, 0.0));
		  ++ncons;
	  }
	
	// constraint (11)
	
	for(int i = 0;i < numArcs;++i)
		if(arcs[i].type == 0)  // intra-link; between nodes in same network
	    {
		  lhs = model.linearNumExpr();
		  int u = arcs[i].from;
		  int v = arcs[i].to;
		  int networkSize= 0; // initialized with dummy value
		  switch(nodes[u].network)
		  {
		    case 0:
		    	networkSize = nsize1;
		    	break;
		    case 1:
		    	networkSize = nsize2;
		    	break;
		    case 2:
		    	networkSize = nsize3;
		    	break;
		    default:
		    	Diagnostic.verror(7,u,nodes[u].network);
		   }
		  int seq = arcs[i].seqIntra;
		  lhs.addTerm(z[fstart+seq], -networkSize);
		  lhs.addTerm(z[hstart+seq], 1.0);
		  matrix.addRow(model.le(lhs, 0.0));
		  ++ncons;
	  }
	
	// constraint (12)
	
	for(int i = 0;i < numArcs;++i)
		if(arcs[i].type == 0)  // intra-link; between nodes in same network
	    {
		  lhs = model.linearNumExpr();
		  int u = arcs[i].from;
		  int v = arcs[i].to;
		  int networkSize = 0; // initialize with dummy value
		  switch(nodes[u].network)
		  {
		    case 0:
		    	networkSize = nsize1;
		    	break;
		    case 1:
		    	networkSize = nsize2;
		    	break;
		    case 2:
		    	networkSize = nsize3;
		    	break;
		    default:
		    	Diagnostic.verror(7,u,nodes[u].network);
		   }
		  int seq = arcs[i].seqIntra;
		  lhs.addTerm(z[fstart+seq], networkSize);
		  lhs.addTerm(z[hstart+seq], 1.0);
		  matrix.addRow(model.ge(lhs, 0.0));
		  ++ncons;
	  }
	 	
	// constraint (13)
	for(int i = 0;i < numNodes;++i)
	{
		lhs = model.linearNumExpr();
		lhs.addTerm(z[gstart+i], -1.0);
		lhs.addTerm(z[tstart+i], 1.0);
		matrix.addRow(model.le(lhs, 0.0));
		++ncons;
	}
	
	// constraint (14)
	for(int j=0;j<3;++j)
	{
		lhs = model.linearNumExpr();
		for(int u = 0; u<numNodes;++u)
			if(nodes[u].network == j)
			  lhs.addTerm(z[tstart+u], 1.0);
		matrix.addRow(model.eq(lhs, 1.0));
		++ncons;
	}
	
	// constraint (15)
		for(int u = 0;u < numNodes;++u)
		{
			int networkSize = 0; // initialize w/ dummy value
			switch(nodes[u].network)
			  {
			    case 0:
			    	networkSize = nsize1;
			    	break;
			    case 1:
			    	networkSize = nsize2;
			    	break;
			    case 2:
			    	networkSize = nsize3;
			    	break;
			    default:
			    	Diagnostic.verror(7,u,nodes[u].network);
			   }
			lhs = model.linearNumExpr();
			lhs.addTerm(z[cstart+u], 1.0);
			lhs.addTerm(z[tstart+u], -networkSize);
			matrix.addRow(model.le(lhs, 0.0));
			++ncons;
		}
	
	// constraint (16)
	for(int u = 0;u < numNodes;++u)
	{
		lhs = model.linearNumExpr();
		lhs.addTerm(z[cstart+u], 1.0);
		lhs.addTerm(z[gstart+u], -1.0);
		for(int i = 0; i < nodes[u].degAll;++i)
		{
			int v = nodes[u].nBors[i];
			if(nodes[v].network == nodes[u].network) // intra-link then
			{
				int ilink = nodes[u].incidentArcs[i];
				int seq = arcs[ilink].seqIntra;
				if(arcs[ilink].from == v)
					lhs.addTerm(z[hstart+seq], 1.0);
				else
					lhs.addTerm(z[hstart+seq], -1.0);
			}
		}
		matrix.addRow(model.eq(lhs, 0.0));
		++ncons;
	}
			
	//System.out.println("(2) ends at "+(ncons-1));	
	
	model.add(matrix);	
	
	IloRange rng2[] =  matrix.getRanges();
	rng[0] = rng2;
	
	System.out.println("model set up for cplex");
	}/*end_ByRows2()*/
	/*-------------------------------------------------------------------------------------------------------*/
	public static int encode(int numNodes, int numArcs, int nInter, Arc arcs[], int ystart, int u, int v, int k)
	{
		// finds the position of y_{uv}^{k}
		// k must be different from u, v
		// (u,v) must be a valid inter-network link or a link from the super source (s,u)
		
		int val = -1;
		
		int numNodes2 = numNodes+1;
		int numArcs2 = numNodes+nInter;
		if(u == v)
		{
			System.out.println("**** in encode() u "+u+" equals v "+v);
			System.exit(0);
		}
		if(k == u)
		{
			System.out.println("**** in encode() k "+k+" equals u "+u);
			System.exit(0);
		}
		if(k == v)
		{
			System.out.println("**** in encode() k "+k+" equals v "+v);
			System.exit(0);
		}
		if((u >= numNodes2) || (u < 0))
		{
			System.out.println("**** in encode() u "+u+" out of bounds");
			System.exit(0);
		}
		if((v >= numNodes2) || (v < 0))
		{
			System.out.println("**** in encode() v "+v+" out of bounds");
			System.exit(0);
		}
		if((k >= numNodes2) || (k < 0))
		{
			System.out.println("**** in encode() k "+k+" out of bounds");
			System.exit(0);
		}
		
		
		int found = 0;
		int debugVar = -1;
		for(int i=0;i<numArcs;++i)  // what about fictitious arcs?
			if(arcs[i].type == 1)
			{
				int node1 = arcs[i].from;
				int node2 = arcs[i].to;
				if((node1 == u) && (node2 == v))
				{
					val = ystart + 2*numArcs2*k + arcs[i].seqInter;
					//System.out.println("val "+val);
					debugVar = arcs[i].seqInter;
					//System.out.println(" seq no "+arcs[i].seqInter);
					found = 1;
					break;
				}
				if((node1 == v) && (node2 == u))
				{
					val = ystart + 2*numArcs2*k + numArcs2 + arcs[i].seqInter;
					//System.out.println("val "+val);
					debugVar = arcs[i].seqInter;
					found = 1;
					break;
				}
			}
		if(found==0)
		{ // handling fictitious arcs here
		  if((u != numNodes) && (v != numNodes))
					{
				      //System.out.println("**** in encode() u v k "+u+" "+v+" "+k+" do not represent a valid inter-network link or virtual link from super source");
				      return(-1);
					}
		
			if(u == numNodes)  // so y_{v*v}^k
				val = ystart+2*numArcs2*k+nInter+v;
			else
				val = ystart+2*numArcs2*k+numArcs2+nInter+u;
		}
		
		/****************
		if(val > 11140)
		{
			System.out.println("val "+val+" u "+u+" v "+v+" k "+k+" ystart "+ystart+" numArcs2 "+numArcs2);
			System.out.println("debugVar "+debugVar);
		}
		*****************/
		
		
		return(val);
		
	}// end_encode()
	//-------------------------------------------------------------------------------------------------------------------------
	public static int[] decode(int numNodes, int numArcs, int nInter, Arc arcs[], int ystart, int val)
	{
		// given position in y_{uv}^{k} array find u,v,k
		
		int[] uvk = new int[4]; // contains solution: values of "found",u,v,k
		
		int numNodes2 = numNodes+1;
		int numArcs2 = nInter+numNodes;
		int lim1 = ystart;
		int lim2 = ystart+2*numArcs2*numNodes2;
		
		if((val < lim1) || (val >= lim2))
		{
			System.out.println("**** in decode() val out of bounds ["+lim1+", "+lim2+" )");
			System.exit(0);
		}
		
		int uv = (val-ystart) % (2*numArcs2);
		int k = (val-ystart)/(2*numArcs2);
		int u=-1;
		int v=-1;
		int found = 0;
		if(uv < numArcs2)
		{
			if(uv < nInter)
			{
				for(int i=0;i<numArcs;++i)
					if(arcs[i].type == 1)
					{
						if(arcs[i].seqInter == uv)
						{
							u = arcs[i].from;
							v = arcs[i].to;
							found = 1;
							break;
						}
					}
			}
			else
			{
				u = (numNodes+1);
				v = uv-nInter;
				found = 1;
			}
		}
		if(uv >= nInter)
		{
			uv -= nInter;
			if(uv < nInter)
			{
				for(int i=0;i<numArcs;++i)
					if(arcs[i].type == 1)
					{
						if(arcs[i].seqInter == uv)
						{
							v = arcs[i].from;
							u = arcs[i].to;
							found = 1;
							break;
						}
					}
			}
			else
			{
				v = (numNodes+1);
				u = uv-nInter;
				found = 1;
			}
			
		}
		if(found==0)
		  {
		      System.out.println("**** in decode() unable to find u,v,k ");
		      System.exit(0);;
		  }
		
		uvk[0] = found;
		uvk[1] = u;
		uvk[2] = v;
		uvk[3] = k;
		return(uvk);
		
	}// end_decode()
	//-------------------------------------------------------------------

}
