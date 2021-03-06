# CommDet

This project contains C++ implementation of Louvain algorithm for fast community detection.
This algorithm is described in following paper:
Fast unfolding of communities in large networks 
by Vincent D. Blondel, Jean-Loup Guillaume, Renaud Lambiotte and Etienne Lefebvre.

This is a heuristic method that is based on modularity optimization. 
It is shown to outperform all other known community detection method in terms of computation time. 
Moreover, the quality of the communities detected is very good, as measured by the modularity.

I have implemented this and tested its performance on different sized datasets.
Specifically, I have used dblp and youtube datasets publicly available from SNAP:
https://snap.stanford.edu/data/com-DBLP.html
https://snap.stanford.edu/data/com-Youtube.html
The resultant communities for both datasets are presented in corresponding output files. 
A detailed timing report during runs on a machine with 3GB of memory running on Intel i3 CPU (32-bit OS) is included.