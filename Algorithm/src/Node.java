
public class Node {
	static final int maxDegree = 20;
	int nnum; /* network number 0,1,2 */
	int disFlag; /* =1 if disabled, 0 otherwise */
	int attackFlag; /* =1 if disabled, 0 otherwise */
	double wt; /* node weight 1-5, but algorithm requires weights to be distinct, so perturb data accordingly */
	int degAll; /* degree counting both inter- and intra- links */
	int degIntra; /* degree based only on intra-network links */
	int nBors[] = new int[maxDegree]; // neighbor list for finding connected components
	int incidentArcs[] = new int[maxDegree]; 
	int xval; /* = 1 if disabled */
	int yval; /* =1 if attacked */
	int network; /* 0,1, or 2 */
	double weight;
	int component;
	int depth;
	Node(){}
}
