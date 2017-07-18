The algorithm will find the minimum weight set of nodes to target in a system of interconnected networks so that a desired level of system degradation is achieved allowing for damage propagation through cascading.

1. The algorithm requires ILOG CPLEX to solve the resulting MIPs.
2. The folder Data Files contains the problems used for testing the algorithm. 
3. The data files contain problems each of which is an interconnected system of 3 networks.  As an example, prob1.txt in the folder 40nodes folder contains the topology for a system of 3 interconnected networks, each with 40 nodes.
	The weights.txt file contains the random weights for the 40 nodes; the same random weights were used for all problems in the folder. 
4. lambda_1, lambda_2, lambda_3 are the degradations sought on the three networks.  All problems were tested with lambda_i = 0.15.
5. Th ecommand line arguments are:
     lambda_1 lambda_2 lambda_3  prob1.txt weights.txt 
6. output is printed to the scree - sample included.
