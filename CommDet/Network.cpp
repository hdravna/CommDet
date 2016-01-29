#include <chrono>
#include "Network.h"

using namespace std;
using namespace std::chrono;


Network* Network::getNewObj()
{
   return new_obj;
}

double Network::getModularity()
{
   return m_modularity;
}

Network::~Network()
{
   m_edge_map.clear();
}

int Network::addNode(int NodeId)
{
   int newId;
   typedef pair<int, int> intPair;
   pair<hash_map<int, int>::iterator, bool> returnPair;
   hash_map<int, int>::iterator iter;
   iter = m_node_hash.find(NodeId);
   if (iter == m_node_hash.end()) //nodeId not present- so insert
   {
      returnPair = m_node_hash.insert(intPair(NodeId, maxNodeId));
      newId = maxNodeId;
      maxNodeId++;
   }
   else //already present
   {
      newId = iter->second;
   }
   return newId;
}

Network& Network::operator=(const Network& rhs)
{
   this->m_num_vertices = rhs.m_num_vertices;
   this->componentNodes = rhs.componentNodes;
   this->m_edge_map = rhs.m_edge_map;
   this->totalLinksIncidentSum = rhs.totalLinksIncidentSum;
   this->currentClusterId = rhs.currentClusterId;
   this->SumWeightsInsideC = rhs.SumWeightsInsideC;

   return *this;
}

Network::Network(const Network& rhs)
{
   this->m_num_vertices = rhs.m_num_vertices;
   this->componentNodes = rhs.componentNodes;
   this->m_edge_map = rhs.m_edge_map;
   this->totalLinksIncidentSum = rhs.totalLinksIncidentSum;
   this->currentClusterId = rhs.currentClusterId;
   this->SumWeightsInsideC = rhs.SumWeightsInsideC;
}


Network::Network(string file_path)
{
   m_input_file_path = file_path;
   fstream input_file_stream;
   string single_line = "";
   input_file_stream.open(m_input_file_path, std::fstream::in);
   if (input_file_stream.is_open())
   {
      getline(input_file_stream, single_line);
      getline(input_file_stream, single_line);
      //This contains # of nodes
      getline(input_file_stream, single_line);
      std::stringstream mystream(single_line);
      string num_nodes;
      mystream >> num_nodes;
      mystream >> num_nodes;
      //This gives # nodes
      mystream >> num_nodes;
      m_num_vertices = atoi(num_nodes.c_str());
      m_init_num_vertices = m_num_vertices;
      string num_edges;
      mystream >> num_edges;
      mystream >> num_edges; //this has num of edges
      m_num_edges = atoi(num_edges.c_str());
      m_num_edges = 2 * m_num_edges; //since undirected graph we have a corresponding reverse edge (not given in file)
      m_is_vertex_present.resize(m_num_vertices);
      selfLoops.resize(m_num_vertices);
      for (int i = 0; i < m_num_vertices; i++)
      {
         m_is_vertex_present[i] = false;
         selfLoops[i] = 0;
      }

      componentNodes.resize(m_num_vertices);
      globalComponentNodes.resize(m_num_vertices);
      nodeId.resize(m_num_vertices);
      nodeIndex.resize(m_num_vertices);
      currentClusterId.resize(m_num_vertices);
      globalCurrentClusterId.resize(m_num_vertices);
      totalNodes.resize(m_num_vertices);
      nodeDegree_Ki.resize(m_num_vertices); //degree of node
      SumWeightsInsideC.resize(m_num_vertices); //Sigma(in) sum of weights inside a cluster
      totalLinksIncidentSum.resize(m_num_vertices); //Sigma(tot) degree of cluster (sum of degree of all component nodes)
      oldNewClusterMappings.resize(m_num_vertices);
      newOldClusterMappings.resize(m_num_vertices);
      getline(input_file_stream, single_line);

      int fromNodeId;
      int toNodeId;
      int fromNodeIndex = 0;
      int toNodeIndex = 0;
      int i = 1;
      int adj_len = 0;
      int adj_index = 0;
      while (getline(input_file_stream, single_line)) /*!input_file_stream.eof()*/
      {
         stringstream nodeValuesStream(single_line);
         nodeValuesStream >> fromNodeId >> toNodeId;

         fromNodeIndex = addNode(fromNodeId); //nodeIndex returned by hash
         if (!m_is_vertex_present[fromNodeIndex])
         {
            nodeId[fromNodeIndex] = fromNodeId;
            nodeIndex[fromNodeIndex] = fromNodeIndex;
            currentClusterId[fromNodeIndex] = fromNodeIndex;
            newOldClusterMappings[fromNodeIndex] = fromNodeIndex;
            componentNodes[fromNodeIndex].push_back(fromNodeIndex);
            globalComponentNodes[fromNodeIndex].push_back(fromNodeIndex);
            SumWeightsInsideC[fromNodeIndex] = 0;
            totalNodes[fromNodeIndex] = 1;
         }

         toNodeIndex = addNode(toNodeId);
         if (!m_is_vertex_present[toNodeIndex])
         {
            nodeIndex[toNodeIndex] = toNodeIndex;
            nodeId[toNodeIndex] = toNodeId;
            currentClusterId[toNodeIndex] = toNodeIndex;
            newOldClusterMappings[toNodeIndex] = toNodeIndex;
            componentNodes[toNodeIndex].push_back(toNodeIndex);
            globalComponentNodes[toNodeIndex].push_back(toNodeIndex);
            SumWeightsInsideC[toNodeIndex] = 0;
            totalNodes[toNodeIndex] = 1;
         }

         bool fromNode_modified_b = false;
         bool toNode_modified_b = false;

         if (m_is_vertex_present[fromNodeIndex] && !m_is_vertex_present[toNodeIndex])
         {
            adjacencyList.insert({ fromNodeIndex, toNodeIndex });
            nodeDegree_Ki[fromNodeIndex] = adjacencyList.count(fromNodeIndex);
            totalLinksIncidentSum[fromNodeIndex] = nodeDegree_Ki[fromNodeIndex];
            adjacencyList.insert({ toNodeIndex, fromNodeIndex });

            nodeDegree_Ki[toNodeIndex] = adjacencyList.count(toNodeIndex);
            totalLinksIncidentSum[toNodeIndex] = nodeDegree_Ki[toNodeIndex];
            fromNode_modified_b = true;
         }

         else if (m_is_vertex_present[fromNodeIndex] && m_is_vertex_present[toNodeIndex])
         {
            adjacencyList.insert({ fromNodeIndex, toNodeIndex });
            nodeDegree_Ki[fromNodeIndex] = adjacencyList.count(fromNodeIndex);
            totalLinksIncidentSum[fromNodeIndex] = nodeDegree_Ki[fromNodeIndex];
            fromNode_modified_b = true;

            adjacencyList.insert({ toNodeIndex, fromNodeIndex });
            nodeDegree_Ki[toNodeIndex] = adjacencyList.count(toNodeIndex);
            totalLinksIncidentSum[toNodeIndex] = nodeDegree_Ki[toNodeIndex];
            toNode_modified_b = true;
         }

         else if (!m_is_vertex_present[fromNodeIndex] && m_is_vertex_present[toNodeIndex])
         {
            adjacencyList.insert({ toNodeIndex, fromNodeIndex });
            nodeDegree_Ki[toNodeIndex] = adjacencyList.count(toNodeIndex);
            totalLinksIncidentSum[toNodeIndex] = nodeDegree_Ki[toNodeIndex];
            toNode_modified_b = true;

            adjacencyList.insert({ fromNodeIndex, toNodeIndex });
            nodeDegree_Ki[fromNodeIndex] = adjacencyList.count(fromNodeIndex);
            totalLinksIncidentSum[fromNodeIndex] = nodeDegree_Ki[fromNodeIndex];
         }
         else if (!m_is_vertex_present[fromNodeIndex] && !m_is_vertex_present[toNodeIndex])
         {
            adjacencyList.insert({ fromNodeIndex, toNodeIndex });
            nodeDegree_Ki[fromNodeIndex] = adjacencyList.count(fromNodeIndex);
            totalLinksIncidentSum[fromNodeIndex] = nodeDegree_Ki[fromNodeIndex];

            adjacencyList.insert({ toNodeIndex, fromNodeIndex });
            nodeDegree_Ki[toNodeIndex] = adjacencyList.count(toNodeIndex);
            totalLinksIncidentSum[toNodeIndex] = nodeDegree_Ki[toNodeIndex];
         }

         if (!m_is_vertex_present[fromNodeIndex] || fromNode_modified_b)
         {
            m_is_vertex_present[fromNodeIndex] = true;
         }

         if (!m_is_vertex_present[toNodeIndex] || toNode_modified_b)
         {
            m_is_vertex_present[toNodeIndex] = true;
         }

         m_edge_map.insert({ fromNodeIndex, toNodeIndex });
         edgeAssocs.insert({ make_pair(fromNodeIndex, toNodeIndex), 1 });
         //edgeWeights.push_back(1);

         m_edge_map.insert({ toNodeIndex, fromNodeIndex });
         edgeAssocs.insert({ make_pair(toNodeIndex, fromNodeIndex), 1 });
         //edgeWeights.push_back(1);
         nodeValuesStream.clear();
      }
      input_file_stream.close();
      linearizeNodes(m_num_edges,m_num_vertices);
      unit_weights_b = true;
      mystream.clear();
   }
   else
   {
      cerr << "\nError! Input file Cannot be opened/Not found!\n";
   }
}

void Network::linearizeNodes(int edges, int vertices)
{
   m_nodes_array.resize(edges);
   m_degrees_array.resize(vertices);
   m_node_ptr_array.resize(vertices);
   EdgeMap::iterator begin;
   pair<EdgeMap::iterator, EdgeMap::iterator> bounds;
   int i = 0;
   int j = 0;
   int start_j = 0;
   int degree = 0;
   for (int i = 0; i < vertices; i++)
   {
      bounds = m_edge_map.equal_range(i);
      begin = bounds.first;
      degree = 0;
      start_j = j;
      while (begin != bounds.second)
      {
         m_nodes_array[j] = begin->second;
         ++degree;
         ++j;
         ++begin;
      }
      m_degrees_array[i] = degree;
      m_node_ptr_array[i] = start_j;
   }
}

int Network::getNumvertices()
{
   return m_num_vertices;
}

void Network::setLatestMod(double mod)
{
   m_latest_mod = mod;
}
double Network::calculateGain(int current_node, int neighbor_cluster_id)
{
   double gain = 0.0;
   int sigma_in = 0;
   int sigma_total = 0;
   int k_i = 0; //degree of node
   int k_i_in = nodeToNeighborsWeights[neighbor_cluster_id];
   //int k_i_in = ki_in(current_node, neighbor_cluster_id);
   sigma_in = SumWeightsInsideC[neighbor_cluster_id];
   sigma_total = totalLinksIncidentSum[neighbor_cluster_id];
   bool is_edge_present_b = false;
   k_i = nodeDegree_Ki[current_node];  //sum of weights of links incident to current node
   gain = (k_i_in - (double(sigma_total*k_i) / m_num_edges));
   return gain;
}


int Network::findNumClusters()
{
   int num_clusters = 0; //# of clusters
   for (int i = 0; i < m_num_vertices; i++)
   {
      if (totalNodes[i] > 0)
      {
         num_clusters++;
      }
   }
   m_num_clusters = num_clusters;
   return m_num_clusters;
}

void Network::computeModularityGains() //Phase I
{
   bool improvement_possible_b = false;
   double old_gain = 0.0;
   double new_gain = 0.0;
   double local_max = 0.0;
   int j = 0;
   int currentNode;
   int current_cluster_id;
   int target_cluster_id;
   pair<adjListType::iterator, adjListType::iterator> bounds;
   adjListType::iterator local_iter;
   pair<adjListType::iterator, adjListType::iterator> adjListBounds;
   neighboringClusters.resize(m_num_vertices);
   int lastestDegree=0;
   nodeToNeighborsWeights.resize(m_num_vertices);
   do
   {
      improvement_possible_b = false;
      for (int i = 0; i < m_num_vertices; i++)
      {
         currentNode = i; //currentNode is int starts from 0 to m_num_vertices
         current_cluster_id = currentClusterId[currentNode];
         target_cluster_id = currentClusterId[currentNode];
         //bounds = adjacencyList.equal_range(currentNode);
         //local_iter = bounds.first;
         int start_pos = m_node_ptr_array[i];
         int neighbor = 0;

         for (int k = 0; k < lastestDegree; k++)
            nodeToNeighborsWeights[neighboringClusters[k]] = 0;

         lastestDegree = m_degrees_array[i];
         for (int d = 0; d < lastestDegree; d++)
         {
            neighbor = m_nodes_array[start_pos + d];
            neighboringClusters[d] = currentClusterId[neighbor];
            if (unit_weights_b)
               nodeToNeighborsWeights[currentClusterId[neighbor]]++;
            else
               nodeToNeighborsWeights[currentClusterId[neighbor]] = nodeToNeighborsWeights[currentClusterId[neighbor]] + edgeWeights[start_pos + d];
         }


         //remove from current cluster
         totalNodes[current_cluster_id] = totalNodes[current_cluster_id] - 1;
         int edgeWeight_contribution = 2 * nodeToNeighborsWeights[current_cluster_id] + selfLoops[currentNode];
         //int edgeWeight_contribution = 2 * ki_in(currentNode,current_cluster_id) + selfLoops[currentNode];
         totalLinksIncidentSum[current_cluster_id] = totalLinksIncidentSum[current_cluster_id] - nodeDegree_Ki[currentNode];
         SumWeightsInsideC[current_cluster_id] = SumWeightsInsideC[current_cluster_id] - edgeWeight_contribution;
         componentNodes[current_cluster_id].remove(currentNode);

         int neighClusterId = 0;
         old_gain = 0.0;
         //while (local_iter != bounds.second)
         for (int n = 0; n < m_degrees_array[i]; n++)
         {
            neighClusterId = neighboringClusters[n];
            //neighClusterId = currentClusterId[local_iter->second];
            new_gain = calculateGain(currentNode, neighClusterId);
            if (new_gain > 0 && (new_gain > old_gain))
            {
               //target_cluster_id = currentClusterId[neighClusterId];
               target_cluster_id = neighClusterId;
               old_gain = new_gain;
               if (new_gain > local_max)
               {
                  local_max = new_gain;
                  improvement_possible_b = true;
               }
            }
            //++local_iter;
         }

      //move to better cluster
      if (target_cluster_id != current_cluster_id)
      {
         int edgeWeight_to_tgt = 2 * nodeToNeighborsWeights[target_cluster_id] + selfLoops[currentNode]; //since weight of every link inside C is 2 
         //int edgeWeight_to_tgt = 2 * ki_in(currentNode,target_cluster_id) + selfLoops[currentNode]; //since weight of every link inside C is 2 
         SumWeightsInsideC[target_cluster_id] = SumWeightsInsideC[target_cluster_id] + edgeWeight_to_tgt;
         totalLinksIncidentSum[target_cluster_id] = totalLinksIncidentSum[target_cluster_id] + nodeDegree_Ki[currentNode];
         componentNodes[target_cluster_id].push_back(currentNode);
         totalNodes[target_cluster_id] = totalNodes[target_cluster_id] + 1;
         currentClusterId[currentNode] = target_cluster_id;
      }
      else
      {
         //move back to same old cluster
         currentClusterId[currentNode] = current_cluster_id;
         componentNodes[current_cluster_id].push_back(currentNode);
         totalNodes[current_cluster_id] = totalNodes[current_cluster_id] + 1;
         totalLinksIncidentSum[current_cluster_id] = totalLinksIncidentSum[current_cluster_id] + nodeDegree_Ki[currentNode];
         SumWeightsInsideC[current_cluster_id] = SumWeightsInsideC[current_cluster_id] + edgeWeight_contribution;
      }
   }
   j++;
} 
//while (improvement_possible_b);
while (j<3);
std::cout << "\nStill improvement possible? " << improvement_possible_b << "\n";
}

//finds number of links present from currentNode to given Cluster
//For this we need to check how many of the current nodes neighbors belong to given cluster
int Network::ki_in(int currentNodeId, int clusterId)
{
   pair<EdgeMap::iterator, EdgeMap::iterator> bounds = m_edge_map.equal_range(currentNodeId);
   EdgeMap::local_iterator begin_iter = bounds.first;
   EdgeMap::local_iterator end_iter = bounds.second;
   assocsType::iterator edgeAssocIterValue;
   int sum = 0;
   int weight = 0;
   while (begin_iter != end_iter)
   {
      weight = 0;
      if (currentClusterId[begin_iter->second] == clusterId)
      {
         edgeAssocIterValue = edgeAssocs.find(make_pair(currentNodeId, begin_iter->second));
         weight = edgeAssocIterValue->second;
         sum = sum + weight;
      }
      ++begin_iter;
   }
   return sum;
}

Network::Network()
{
}

void Network::identifyNewNeighbors(Network* nw, Network* parent)
{
   EdgeMap::iterator edgeIterBegin = parent->m_edge_map.begin();
   EdgeMap::iterator edgeIterEnd = parent->m_edge_map.end();
   map<pair<int, int>, int>::iterator assocsIterEnd = nw->clusterAssocs.end();
   map<pair<int, int>, int>::iterator assocsValueIter;
   while (edgeIterBegin != edgeIterEnd)
   {
      //both ends of the edge belong to same cluster- increment number of self loops for clusterId
      //Otherwise increment the edge weight in clusterAssocs between two clusters
      int fromNode = edgeIterBegin->first;
      int toNode = edgeIterBegin->second;
      int fromNodeCluster = parent->currentClusterId[fromNode];
      int toNodeCluster = parent->currentClusterId[toNode];
      int newFromNode = parent->oldNewClusterMappings[fromNodeCluster];
      int newToNode = parent->oldNewClusterMappings[toNodeCluster];
      int new_val = 0;
      pair<int, int> p;
      if (newFromNode == newToNode)
      {
         nw->selfLoops[newFromNode]++;
      }
      else
      {
         assocsValueIter = nw->clusterAssocs.find(make_pair(newFromNode, newToNode));
         if (assocsValueIter == assocsIterEnd) //fresh link between clusters
         {
            nw->clusterAssocs.insert({ make_pair(newFromNode, newToNode), 1 });
            nw->edgeAssocs.insert({ make_pair(newFromNode, newToNode), 1 });
            nw->m_edge_map.insert({ newFromNode, newToNode});
            nw->adjacencyList.insert(make_pair(newFromNode, newToNode));
         }
         else //link between clusters already exists- increment the weight
         {
            p = make_pair(newFromNode, newToNode);
            new_val = assocsValueIter->second;
            new_val++;
            nw->clusterAssocs.erase(p);
            nw->edgeAssocs.erase(p);
            nw->clusterAssocs.insert({ p, new_val });
            nw->edgeAssocs.insert({ p, new_val });
         }
      }
      ++edgeIterBegin;
   }
}

void Network::createNewNetwork(int clusters, Network* parent) //Phase II
{ //in this Index means Index starting from 0 to clusters-1
   new_obj = new Network();
   int newNodeId;
   int newNodeIndex;
   maxNodeId = 0;
   adjacencyList.clear();
   oldNewClusterMappings.clear();
   oldNewClusterMappings.resize(m_num_vertices);
   newOldClusterMappings.clear();
   newOldClusterMappings.resize(m_num_vertices);
   new_obj->nodeIndex.resize(clusters);
   new_obj->nodeId.resize(clusters);
   new_obj->totalNodes.resize(clusters);
   new_obj->currentClusterId.resize(clusters);
   new_obj->nodeDegree_Ki.resize(clusters);
   new_obj->selfLoops.resize(clusters);
   new_obj->totalLinksIncidentSum.resize(clusters);
   new_obj->m_is_vertex_present.resize(clusters);
   new_obj->componentNodes.resize(clusters); //IMPORTANT **now this holds the previous iterations summary**
   new_obj->SumWeightsInsideC.resize(clusters);
   for (int i = 0; i < m_num_vertices; i++)
   {
      if (totalNodes[i] > 0)
      {
         newNodeId = nodeId[i];
         newNodeIndex = new_obj->addNode(newNodeId);
         new_obj->nodeId[newNodeIndex] = newNodeId;
         new_obj->totalNodes[newNodeIndex] = totalNodes[i];
         //copy all the components from parent node (IMPORTANT **now this holds the previous iterations summary**)
         std::copy(componentNodes[i].begin(), componentNodes[i].end(),
            std::back_insert_iterator<std::list<int> >(new_obj->componentNodes[newNodeIndex]));
         new_obj->SumWeightsInsideC[newNodeIndex] = SumWeightsInsideC[i];
         new_obj->nodeDegree_Ki[newNodeIndex] = totalLinksIncidentSum[i];
         new_obj->currentClusterId[newNodeIndex] = newNodeIndex; //these by themselves are clusters initially
         new_obj->totalLinksIncidentSum[newNodeIndex] = new_obj->nodeDegree_Ki[newNodeIndex];
         oldNewClusterMappings[i] = newNodeIndex;
         componentNodes[i].clear(); //free this for next iteration
      }
   }

   clusterAssocs.clear();
   new_obj->identifyNewNeighbors(new_obj, parent);
   swap(edgeAssocs, new_obj->edgeAssocs);
   swap(m_edge_map, new_obj->m_edge_map);
   assocsType::iterator begin_iter = edgeAssocs.begin();
   assocsType::iterator end_iter = edgeAssocs.end();
   edgeWeights.clear();
   edgeWeights.resize(m_edge_map.size());
   int i = 0;
   while (begin_iter != end_iter)
   {
      edgeWeights[i] = begin_iter->second;
      i++;
      ++begin_iter;
   }
   unit_weights_b = false;
   linearizeNodes(m_edge_map.size(), clusters);

   //By the time links amongst clusters are identified- we have 
   // 1. Degrees of each node (cluster) nodeDegree_Ki[] = sum of all the counts for this node
   // 2. For clusters to be formed from these we shall have totalLinksIncidentSum[]- initially degree of cluster = nodeDegree_Ki
   // 3. SumWeightsInsideC already known. Equals SumWeightsInsideC from last iteration
   // 4. Ki_in have to be computed on the fly using the function

   int total = 0;
   nodeDegree_Ki.clear();
   nodeDegree_Ki.resize(clusters);
   totalLinksIncidentSum.clear();
   totalLinksIncidentSum.resize(clusters);
   totalNodes.clear();
   currentClusterId.clear();
   totalNodes.resize(clusters);
   componentNodes.resize(clusters);
   currentClusterId.resize(clusters);

   for (int i = 0; i < clusters; i++)
   {
      totalNodes[i] = 1;
      currentClusterId[i] = i;
      componentNodes[i].clear();
      componentNodes[i].push_back(i);
   }

   swap(SumWeightsInsideC, new_obj->SumWeightsInsideC);
   swap(selfLoops, new_obj->selfLoops);
   swap(nodeDegree_Ki, new_obj->nodeDegree_Ki);
   swap(adjacencyList, new_obj->adjacencyList);
   swap(totalLinksIncidentSum, new_obj->totalLinksIncidentSum);
   swap(clusterAssocs, new_obj->clusterAssocs);
   this->m_num_vertices = clusters;
}

void Network::regroupNodes(int clusters)
{
   for (int i = 0; i < clusters; i++)
   {
      list<int>::iterator c_iter_begin = new_obj->componentNodes[i].begin();
      list<int>::iterator c_iter_end = new_obj->componentNodes[i].end();
      list<int>::iterator c_this_list_begin;
      list<int>::iterator c_this_list_end;
      list<int>::iterator c_new_list_begin;
      list<int>::iterator c_new_list_end;
      while (c_iter_begin != c_iter_end)
      {
         size_t this_size = globalComponentNodes[i].size();
         size_t new_size = globalComponentNodes[*c_iter_begin].size();
         c_this_list_begin = globalComponentNodes[i].begin();
         c_this_list_end = globalComponentNodes[i].end();
         globalComponentNodes[i].splice(c_this_list_end, globalComponentNodes[*c_iter_begin]); //constant time :)
         ++c_iter_begin;
      }
   }
   cerr << "\n# clusters after regrouping = " << clusters << "\n";
}

double Network::computeNetworkModularity()
{
   double ls = 0.0;
   double ds = 0.0;
   double q = 0.0;
   for (int i = 0; i < m_num_vertices; i++)
   {
      if (totalNodes[i]>0)
      {
         ls = (double(SumWeightsInsideC[i]) / m_num_edges); //m_num_edges here also counts reverse edges so its = 2*|m|
         ds = (double(totalLinksIncidentSum[i]) / m_num_edges);
         q += (ls - (ds*ds));
      }
   }
   m_modularity = q;
   std::cout << "\nmodularity of N/W = " << m_modularity << "\n";
   std::cout << "\nnumber of clusters formed = " << m_num_clusters;
   return q;
}
void Network::outputToFile(string opf)
{
   ofstream opfile(opf);
   if (opfile.is_open())
   {
      opfile << "\n--Louvian's Community Detection Algorithm --\n";
      opfile << "\n--# of Communities detected = "<< m_latest_cluster_count <<"\n";
      opfile << "\n--Modularity of network = " << m_latest_mod << "\n";
      opfile << "NodeId \t ClusterId\n";
      for (int i = 0; i < m_num_vertices; i++)
      {
         list<int>::iterator c_iter_begin = globalComponentNodes[i].begin();
         list<int>::iterator c_iter_end = globalComponentNodes[i].end();
         while (c_iter_begin != c_iter_end)
         {
            opfile << nodeId[*c_iter_begin] << "\t" << i << "\n";
            ++c_iter_begin;
         }
      }
      opfile.close();
   }
   else
   {
      cout << "Unable to open output file!\n";
   }
}
int main(int argc, char** argv)
{
   string input_file_name = "";
   string output_file_name = "";
   string num_clusters = "";
   if (argc != 4)
   {
      cerr << "Please provide Filename ip-file #-clusters op-file\n";
   }
   for (int i = 0; i < argc; i++)
   {
      if (i == 1)
         input_file_name = argv[i];
      else if (i == 2)
         num_clusters = argv[i];
      else if (i == 3)
         output_file_name = argv[i];
   }

   /*input_file_name = "D:/com-dblp.ungraph.txt";*/
   //input_file_name = "D:/com-youtube.ungraph.txt";
   //output_file_name = "D:/youtube-op3.txt";
   //input_file_name = "D:/Sample2.txt";
   /*output_file_name = "D:/sample-op.txt";*/
   //output_file_name = "D:/dblp-op.txt";
   cerr << "\nStarting file read... \n";
   high_resolution_clock::time_point t1 = high_resolution_clock::now();
   Network* network_obj = new Network(input_file_name);
   high_resolution_clock::time_point t2 = high_resolution_clock::now();
   auto duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
   cerr << "\nReading done! Time taken: " << duration << " seconds\n";
   int clusters = network_obj->getNumvertices();
   int pass = 1;
   double oldModularity = 0.0;
   double newModularity = 0.0;
   bool continue_b = false;
   do
   {
      cerr << "\n---------- Starting pass # " << pass <<"--------------\n";
      cerr << "\nStarted Phase#1... \n";
      t1 = high_resolution_clock::now();
      network_obj->computeModularityGains(); //Phase I
      int clusters = network_obj->findNumClusters();
      t2 = high_resolution_clock::now();

      duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
      cerr << "\nPhase I done! Time taken: " << duration << " seconds\n";
      oldModularity = newModularity;
      newModularity = network_obj->computeNetworkModularity();
      if (newModularity < oldModularity)
      {
         network_obj->setLatestMod(oldModularity);
         cerr << "\nNo improvement in modularity. Ignoring this pass! " << "\n";
         cerr << "\nLatest best modularity: " << oldModularity << "\n";
      }
      continue_b = newModularity > oldModularity;
      if (continue_b)
      {
         cerr << "\nNumber of clusters: " << clusters << "\n";
         cerr << "\nStarted Phase#2... \n";
         t1 = high_resolution_clock::now();
         network_obj->createNewNetwork(clusters, network_obj); //Phase II
         network_obj->regroupNodes(clusters);
         t2 = high_resolution_clock::now();
         duration = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1).count();
         cerr << "\nPhase II success! Time taken: " << duration << " seconds\n\n";
         pass++;
      }
   } while (continue_b);

   cerr << "\n\n-------Writing to Output File--------\n";
   network_obj->outputToFile(output_file_name);
   delete network_obj->getNewObj();
   delete network_obj;
   int i = 0;
   cout << "\n\nPress any key to exit...\n";
   cin >> i;
}
