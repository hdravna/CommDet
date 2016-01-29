/*
//Author: Anvardh Nanduri
//This is a C++ implementation of Louvain algorithm.
//Reference: Fast unfolding of communities in large networks Blondel et.al
//Department of Mathematical Engineering, Universit´e catholique de Louvain.
//This code has been developed as part of Social Network Analysis course work
//at George Mason University in Mar 2015 solely for academic purpose.
//No part of this code can be used without author's prior permission
//email: anvardh@gmail.com
*/

#ifndef _NETWORK_H_
#define _NETWORK_H_

#include <math.h>
#include <list>
#include <fstream>
#include <vector>
#include <sstream>
#include <random>
#include <hash_map>
#include <unordered_map>
#include <map>
#include <iostream>

using namespace std;
class Network
{
private:
   int m_num_vertices;
   int m_init_num_vertices;
   int m_num_clusters;
   //# of vertices in the graph. Either hardcoded or known from # of lines in input file
   int m_num_edges;
   string m_input_file_path; //path of input file which may contain adjacency matrix

   double m_modularity;
   //multimap can have multiple values with same key. They will be stored in one bucket.
   //We can have fromNode as the key with multiple toNodes. 
   typedef std::unordered_multimap<int, int> EdgeMap;
   EdgeMap m_edge_map;
   hash_map<int, int> m_node_hash;
   int maxNodeId = 0;
   int addNode(int);
   int findIfNodePresent(int);
   vector<bool> m_is_vertex_present;
   bool unit_weights_b;

   int m_latest_cluster_count;
   double m_latest_mod;
   vector<int> nodeId;
   vector<int> nodeIndex;
   vector<int> oldNewClusterMappings;
   vector<int> newOldClusterMappings; //give newclusterId as input and get old clusterId
   vector<int> currentClusterId;
   vector<int> globalCurrentClusterId;
   vector<int> selfLoops;

   typedef unordered_multimap<int, int> adjListType;
   typedef unordered_multimap<int, pair<int,int>> clusterAdjListType; //extended version of clusterAssocs
   adjListType adjacencyList;
   
   typedef map<pair<int, int>, int> assocsType;
   assocsType edgeAssocs;
   map<pair<int, int>, int> clusterAssocs;
   vector <list <int> > componentNodes;
   vector <list <int> > globalComponentNodes;

   vector<int> totalNodes;
   vector<int> nodeDegree_Ki;
   vector<int> SumWeightsInsideC; //Sigma(in)
   vector<int> totalLinksIncidentSum; //Sigma(tot)
   vector<int> m_nodes_array;
   vector<int> m_degrees_array;
   vector<int> m_node_ptr_array;
   vector<int> neighboringClusters;
   vector<int> nodeToNeighborsWeights;
   vector<int> edgeWeights;
   int ki_in(int nodeId, int clusterId);
   Network* new_obj;
   void linearizeNodes(int, int);

public:
   Network(string);
   Network();
   Network(const Network& rhs);
   Network& operator=(const Network& rhs);
   ~Network();

   void readGraph();
   void printClusters();
   bool getEdge(int, int, EdgeMap::iterator&); //returns the Edge of the input fromNode and toNode
   void resetAllDataStructures();
   void computeModularityGains(); //phase I
   void createNewNetwork(int, Network*);      //phase II
   double calculateGain(int, int);
   double computeNetworkModularity();
   void regroupNodes(int);
   void identifyNewNeighbors(Network* nw, Network* parent);
   int findNumClusters();
   void outputToFile(string);
   unordered_multimap<int, int>* getAdjList();
   Network* getNewObj();
   void setNewObj();
   void setLatestMod(double);
   int getNumvertices();
   double getModularity();
};


#endif
