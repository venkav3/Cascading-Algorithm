
public class Arc {
int from, to;
int type; // type = 0 for intra, type = 1 for inter-network link, type = 2 for fictious arcs from super source to each of the nodes 
int seqIntra; // sequence number in set of Intra links if appropriate
int seqInter;
Arc(){};
}
