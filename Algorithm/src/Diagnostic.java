import java.io.BufferedReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Date;



public class Diagnostic {
	static DecimalFormat decfor = new DecimalFormat("0.00");
	static DecimalFormat decfor2 = new DecimalFormat("0.0000");
	static final double EPS2 = 0.01;
	static public int maxNodes = 500;
	static String towrite;
   	static String outfile = "//cob-file.scheller.gatech.edu/profiles/vvenkate3/Documents/Desktop/Computational_Tests/output.txt";
	/*------------------------------------------------------------------------------------*/
    
   public static void verror(int icase, int val,int max)
    {
    	switch(icase){
    	case 1:
    		System.out.println("numNodes "+val+" > max "+max+" aborting\n");
    		break;
    	case 2:
    		System.out.println("numArcs "+val+" > max "+max+" aborting\n");
    		break;
    	case 3:
    		System.out.println("degree of node "+val+" > max "+max+" aborting\n");
    		break;
    	case 4:
    		System.out.println("node  "+val+" has neighbor "+max+" repeated in multiple links\n");
    		break;
    	case 5:
    		System.out.println("variables in constraint matrix "+val+" not equal to predicted, numvar "+max+"; aborting\n");
    		break;
    	case 6:
    		System.out.println("for node "+val+" discrepency between no.intra network neighbors measured and predicted; delta "+max+"; aborting\n");
    		break;
    	case 7:
    		System.out.println("for node "+val+" invalid network number "+max+"; aborting\n");
    		break;
    	default:
    	}
    	System.exit(0);
    }/*end_verror()*/
 /*-------------------------------------------------------------------*/
   public static void dumpNodes(int nsize1, int nsize2, int nsize3, int numArcs, Node nodes[], Arc arcs[]) 
   {
	 System.out.println("enter dumpNodes()");
	 int numNodes = nsize1+nsize2+nsize3;
	 for(int i=0;i<numNodes;++i)
	 {
		 System.out.println("\n node "+i+" weight "+decfor.format(nodes[i].weight)+" degAll "+nodes[i].degAll+" degree intra "+nodes[i].degIntra+" neighbors:");
		 for(int j=0;j<nodes[i].degAll;++j)
			 System.out.println("\n node "+nodes[i].nBors[j]+" type "+arcs[nodes[i].incidentArcs[j]].type);
	 }
   }/*end_dumpNodes()*/
   /*--------------------------------------------------------------------------------------*/
   public static void dumpArcs( int numArcs, Arc arcs[]) 
   {
	 System.out.println("enter dumpArcs()");
	 
	 for(int i=0;i<numArcs;++i)
	 {
		 
			 System.out.print("\n link "+i+" ("+arcs[i].from+" , "+arcs[i].to+" )"+" type "+arcs[i].type);
			 if(arcs[i].type == 0)
				 System.out.print("Intraseq "+arcs[i].seqIntra);
			 else
				 System.out.print("Intraseq "+arcs[i].seqInter);
			 System.out.print("\n");
				 
	 }
   }/*end_dumpNodes()*/
   /*--------------------------------------------------------------------------------------*/
   public static void output(int nvars, int ncons, double x[], double slacks[],int nsize1, int nsize2, int nsize3, int numNodes, int numArcs, Node nodes[], Arc arcs[],  long startTimeMs) throws IOException
   {
	   
	   OutPut.initFile(outfile);
	// count # of inter-network links
		
	   
	   int nattack = 0;
	   System.out.println("Nodes attacked:");
	   
	   for(int i = 0;i < numNodes ;++i)
	   {
		   
		   if(Math.abs(x[numNodes+i]-1.0)< EPS2)
		   {
			   System.out.print("node "+i+" weight "+nodes[i].weight+"\n");
		       ++nattack;
		   }
		   
	   }
	   if(nattack == 0)
		   System.out.println(" None \n");
	   
	   
	   
	   
	   
	   int[] tempA = new int[maxNodes];
	   System.out.println("Nodes disabled:");
	   int count = 0;
	   
	   for(int i = 0;i < numNodes ;++i)
		   if(Math.abs(x[i]-1.0)< EPS2)
		   {
			   System.out.print(i+" ");
			   tempA[count] = i;  // tempA[] contains list of nodes disabled
			   ++count;
			   if((count%15) == 0)
				   System.out.print("\n");
		   }
	   System.out.println(" count "+count);
	   Node[] nodes2 = new Node[maxNodes];
	   nodes2 = Component.loadCompTesting2(count,tempA,nodes,arcs);
	   int ncomp = Component.find_components(count, nodes2);
	   //System.out.println("\nncomp "+ncomp);
	   double xobj = 0;
	   int xnode = -1;
	   for(int i=0;i<ncomp;++i)
	   {
		   double xmin = 999;
		   int num = 0;
		   for(int j=0;j<count;++j)
		   {
			   if(nodes2[j].component == i)	
			   {
				   ++num;
				   
			     if(nodes2[j].weight < xmin)
			     {
			    	 xmin = nodes2[j].weight;
			    	 xnode = tempA[j];
			     }
			   }
		   }
		   
		   if(num == 1)
		   {
			   // Singleton! Add its weight to obj function only if it is directly connected to a x= 0 node
			   int flag = 0;
			   for(int k=0;k<nodes[xnode].degAll;++k)
			   {
				  int ilink = nodes[xnode].incidentArcs[k];
				  if(arcs[ilink].type == 0)
				  {
					  int node2 = nodes[xnode].nBors[k];
					  if(Math.abs(x[node2]-0.0)< EPS2)
						  flag = 1;
				  }
			   }
			   if(flag == 0)
				 {
					 // no, this node will disable when attacked nodes fail, set xmin to 0 - no contribution to obj. value
					 xmin = 0;
				 }
			   
		    }
			xobj += xmin;    	 
	   }
	   //System.out.println("\nObj value computed by second method (approx. method) "+xobj);
	   
	   System.out.println("\nDegradation effected:");
	   double frac1 = 0;
	   double frac2 = 0;
	   double frac3 = 0;
	   for(int i = 0;i < numNodes ;++i)
	   {
		   if(i < nsize1)
		   {
		    if(Math.abs(x[i]-1.0)< EPS2)
		    	frac1 += 1.0;
		   }
		    else if(i < nsize1+nsize2)
		    {
			    if(Math.abs(x[i]-1.0)< EPS2)
			    	frac2 += 1.0;
			}
		    else
		    {
			    if(Math.abs(x[i]-1.0)< EPS2)
			    	frac3 += 1.0;
			 }
	   }
	   System.out.println(decfor.format(frac1/nsize1)+", "+decfor.format(frac2/nsize2)+", "+decfor.format(frac3/nsize3));
	   long taskTimeMs  = System.currentTimeMillis( ) - startTimeMs;
	   System.out.println("elapsed time = "+taskTimeMs+" (ms)");
   }/*end_output()*/
   /*-----------------------------------------------------------------------------------------*/
   
   public static void dump_x(double x[],int numNodes,int node)
   {
   //System.out.println(" Dump of x[] \n");
   //for(int i = 0;i < numNodes ;++i)
	   //System.out.println(" node "+i+ " x "+x[i]+"\n");
   System.out.println(" node "+node+ " x "+x[node]+"\n");
   //System.out.println("-----------------------\n");
   }/*end_dump_y()*/
   /*-----------------------------------------------------------------------------------------*/
   public static void dump_z(double x[],int numNodes, int node)
   {
   System.out.println(" Dump of z[]; numnodes = "+numNodes+ "\n");
   for(int i = 0;i < numNodes ;++i)
	   System.out.println(" node "+i+ " z "+decfor2.format(x[numNodes+i])+"\n");
	   System.out.println(" node "+node+ " z "+decfor2.format(x[numNodes+node])+"\n");
   System.out.println("-----------------------\n");
   }/*end_dump_y()*/
   /*-----------------------------------------------------------------------------------------*/
   public static void dump_s(double x[],int numNodes, int node)
   {
   //System.out.println(" Dump of s[] \n");
   //for(int i = 0;i < numNodes ;++i)
	   //System.out.println(" node "+i+ " s "+x[2*numNodes+i]+"\n");
	   System.out.println(" node "+node+ " s "+x[2*numNodes+node]+"\n");
   //System.out.println("-----------------------\n");
   }/*end_dump_y()*/
   /*-----------------------------------------------------------------------------------------*/
   public static void dump_y(double x[],int numNodes, int node)
   {
   //System.out.println(" Dump of y[] \n");
   //for(int i = 0;i < numNodes ;++i)
	   //System.out.println(" node "+i+ " y "+x[3*numNodes+i]+"\n");
	   System.out.println(" node "+node+ " y "+x[3*numNodes+node]+"\n");
   //System.out.println("-----------------------\n");
   }/*end_dump_y()*/
   /*-----------------------------------------------------------------------------------------*/
   public static void dump_all(double x[],int numNodes) throws IOException
   {
	   
   	
       //System.out.println(" ------------------------------ soln -----------------------------");
       
   System.out.println(" Dump of z[] \n");
   for(int i = 0;i < 5*numNodes ;++i)
   {
	   //System.out.println(" x["+i+ "] = "+decfor2.format(x[i])+"\n");
	   towrite = " x["+i+ "] = "+decfor2.format(x[i])+"\r\n";
       OutPut.writeOut(outfile, towrite);
   }
	   //System.out.println(" node "+node+ " z "+decfor2.format(x[numNodes+node])+"\n");
   System.out.println("-----------------------\n");
   }/*end_dump_y()*/
   /*-----------------------------------------------------------------------------------------*/
   public static void dump_slacks(double slacks[],int ncons) throws IOException
   {
   //System.out.println(" ----------------------- \n Dump of slacks \n");
   for(int i = 0;i < ncons ;++i)
   {
	   towrite = " slack["+i+ "] = "+decfor2.format(slacks[i])+"\r\n";
       OutPut.writeOut(outfile, towrite);
	}
	   //System.out.println(" node "+node+ " z "+decfor2.format(x[numNodes+node])+"\n");
   //System.out.println("-----------------------\n");
   }/*end_dump_y()*/
   /*-----------------------------------------------------------------------------------------*/
}/*end-class*/
