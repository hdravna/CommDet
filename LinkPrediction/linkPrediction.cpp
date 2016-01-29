//Author: Anvardh, Nanduri
//(June 2015, GMU)
//This file contains the code (all interfaces and definitions) for reading the input network, 
//analyzing it and generating various static and temporal features for network analysis 
//over various evolving snapshots of the given network.
//Makes use of C++ SNAP library (Stanford University) for analysis
//Also makes use of Microsoft Office Objects for redirecting output to excel sheets.
//This code cannot be reused in any form without prior permission from author.
//email: anvardh@gmail.com

//microsoft office objects
#import \
   "C:\Program Files\Common Files\microsoft shared\OFFICE12\mso.dll" \
   rename("DocumentProperties", "DocumentPropertiesXL")\
   rename("RGB", "RGBXL")

//MS VBA objects
#import \
   "C:\Program Files\Common Files\microsoft shared\VBA\VBA6\VBE6EXT.OLB"

//Excel application objects
#import \
   "C:\Program Files\Microsoft Office\Office12\EXCEL.EXE"\
   rename("DialogBox", "DialogBoxXL")\
   rename("RGB", "RGBXL")\
   rename("DocumentProperties", "DocumentPropertiesXL")\
   rename("ReplaceText", "ReplaceTextXL")\
   rename("CopyFile", "CopyFileXL")\
   exclude("IFont", "IPicture") no_dual_interfaces
/*************************************************************/

#include<iostream>
#include <Snap.h>
//#include <graph.h>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <ctime>
#include <map>
#include <random>
#include <hash_map>
#include <unordered_map>

using namespace std;
using namespace TSnap;

class GraphU
{
public:
   GraphU(int, char*, char*);
   ~GraphU();
   GraphU(const GraphU&);
   const GraphU& operator=(const GraphU&);
   void extractFeatureSet();
   void calculateMovingAverages(int, int, int, bool);
   void setGraph(PTimeNet);
   PTimeNet getGraph(void);
   vector<vector<double>>* getNonTemporalFeatureSet();
   vector<vector<double>>* getTemporalFeatureSet();

   void setNonTemporalFeatureSet(vector<vector<double>>*);
   void setTemporalFeatureSet(vector<vector<double>>*);

   const vector<bool>* getTargetConceptVec();
   void setTargetConceptVec(const vector<bool>*);
   void analyzeTemporalCharacteristics();
   void generateTemoralGraphModels();

   const TVec<PNGraph>* getSnapShotVec();

   const vector<int>* getNodesInTimeSlot();
   const vector<int>* getEdgesInTimeSlot();
   const vector<int>* getLinkTimeBucketKeyVec();
   const vector<double>* getEffDiamInTimeSlot();
   const vector<double>* getFFMEffecDiamV1();
   const vector<double>* getFFMEffecDiamV2();
   const vector<double>* getdegreeInTimeSlot();
   const vector<double>* getMaxSccNodesInTimeSlot();
   const vector<pair<int, int>>* getFFMDensificationV1();
   const vector<pair<int, int>>* getFFMDensificationV2();

private:
   int maxNodeId;
   int m_num_edges;
   PTimeNet m_graph;
   PNGraph m_msg_graph;
   PNGraph m_ungraph;
   int m_num_vertices;
   int clustering_coeff;
   PNGraph forestFireGraph1;
   PNGraph forestFireGraph2;
   bool m_is_directed_b = false;

   vector<int> nodesInTimeSlot;
   vector<int> edgesInTimeSlot;
   vector<pair<int,int>> edgeTimeBucketVec;
   vector<double> FFMEffecDiamV1;
   vector<double> FFMEffecDiamV2;
   vector<double> degreeInTimeSlot;
   vector<double> effecDiamInTimeSlot;
   vector<double> maxSccNodesFraction;
   vector<bool> m_target_concept_vec;
   vector<PNGraph> m_modelGraphsVector;
   vector<vector<double>> m_temporal_fv;
   vector<vector<double>> m_placeholder_fv;
   vector<vector<double>> m_non_temporal_fv;


   vector<pair<int, int>> FFMDensificationV1;
   vector<pair<int, int>> FFMDensificationV2;

   vector<vector<pair<int, int>>> testSamplesVectorsVec;
   vector<vector<pair<int, int>>> positiveSamplesVectorsVec;
   vector<vector<pair<int, int>>> negativeSamplesVectorsVec;
   vector<vector<pair<int, int>>> hardNegativeSamplesVectorsVec;


   PNGraph graphSnapShotPt;
   TNGraph graphSnapShot;
   TVec<PNGraph> snapShotVec;

   int addNode(int);
   int getNodeIndexGivenNodeId(int);
   TVec<TIntV> messageTimeBucketNodeVec;
   TVec<TIntPrV> messageTimeBucketEdgeVec;
   vector<int> messageTimeBucketKeyVec;
   map<int, int> keyMessageTimeMap;

   TVec<TIntV> linkTimeBucketNodeVec;
   TVec<TIntPrV> linkTimeBucketEdgeVec;
   vector<int> linkTimeBucketKeyVec;

   vector<int> nodeId;
   vector<int> nodeIndex;
   vector<int> nodeDegree_Ki;
   vector<bool> m_is_vertex_present;
   int selfies;

   typedef std::unordered_multimap<int, int> adjListType;
   adjListType adjacencyList;

   typedef map<pair<int, int>, int> assocsType;
   assocsType edgeAssocs;
   hash_map<int, int> m_node_hash;
   THash <TInt, TSecTm> m_birthTime_hash;
   multimap<int, int> m_recent_activity_hash;
   hash_map<int, int> m_node_month_hash;
   unordered_multimap<int, pair<int, int>> messagesEdgeTimeMap;
   unordered_multimap<int, pair<int, int>> linksEdgeTimeMap;

   typedef THash<TPair <TInt, TInt>, TSecTm>  timeMap;
   timeMap messagesFullTimeMap;
   timeMap linksFullTimeMap;

   //multimap can have multiple values with same key. They will be stored in one bucket.
   //We can have fromNode as the key with multiple toNodes. 
   typedef std::unordered_multimap<int, int> EdgeMap;
   EdgeMap m_edge_map;

   TTimeNet::TTmBucketV TmBucketV;
   int calculateRecency(int, int);
   int numEdgesAmongNeighbors(PNGraph, int, int);
   int numEdgesAmongActiveNeighbors(PNGraph,int, int, int);
   bool calculateFeatures(int, int, int, bool, bool);
   void calculateTemporalFeatures(int, int, int, bool);
   void getInnerGraphEdges(int, int, TIntPrV&);
   void calculateMessageTimeBucketVector();
   void calculateLinksTimeBucketVector();
   int calculateActiveDegree(int, TIntV&);
   bool isNodeActive(int, int);
   int calculateActiveCN(int, TIntV&);
   
public:
   void readGraph(string, string);
   int activityThreshold;
};

const vector<int>* GraphU::getNodesInTimeSlot()
{
   return &nodesInTimeSlot;
}

const vector<int>* GraphU::getEdgesInTimeSlot()
{
   return &edgesInTimeSlot;
}

const vector<int>* GraphU::getLinkTimeBucketKeyVec()
{
   return &linkTimeBucketKeyVec;
}

const vector<double>* GraphU::getFFMEffecDiamV1()
{
   return &FFMEffecDiamV1;
}

const TVec<PNGraph>* GraphU::getSnapShotVec()
{
   return &snapShotVec;
}

const vector<double>* GraphU::getFFMEffecDiamV2()
{
   return &FFMEffecDiamV2;
}

const vector<pair<int, int>>* GraphU::getFFMDensificationV1()
{
   return &FFMDensificationV1;
}

const vector<pair<int, int>>* GraphU::getFFMDensificationV2()
{
   return &FFMDensificationV2;
}

GraphU::GraphU(int actThreshold, char* inputPath1, char* inputPath2)
{
   //const TStr tFileInput = inputPath1;
   //m_graph = TSnap::LoadEdgeList<PUNGraph>(inputPath, 0, 1);
   //m_graph = TTimeNet::LoadArxiv(inputPath1, inputPath2);
   //TGnuPlot::GnuPlotPath = "C:/gnuplot/bin";
   //m_graph->PlotCCfOverTm("linkPredOP1", "", TTmUnit::tmuMonth, 20);
   //readGraph(inputPath1);
   activityThreshold = actThreshold;
}

void GraphU::calculateMessageTimeBucketVector()
{
   TVec<TIntPr> timeBucket;
   std::sort(messageTimeBucketKeyVec.begin(), messageTimeBucketKeyVec.end());
   TIntV nodeTimeBucket;
   for (int i = 0; i < messageTimeBucketKeyVec.size(); i++)
   {
      int key = messageTimeBucketKeyVec[i];
      keyMessageTimeMap.insert({ key, i });
      typedef unordered_multimap<int, pair<int, int>>::iterator Iter;
      pair<Iter, Iter> bounds;
      bounds = messagesEdgeTimeMap.equal_range(key);
      Iter begin = bounds.first;
      nodeTimeBucket.Clr();
      while (begin != bounds.second)
      {
         nodeTimeBucket.AddUnique(begin->second.first);
         nodeTimeBucket.AddUnique(begin->second.second);
         ++begin;
      }
      messageTimeBucketEdgeVec.AddUnique(timeBucket);
      messageTimeBucketNodeVec.AddUnique(nodeTimeBucket);
   }
}

void GraphU::calculateLinksTimeBucketVector()
{
   TVec<TIntPr> edgeTimeBucket;
   std::sort(linkTimeBucketKeyVec.begin(), linkTimeBucketKeyVec.end());
   TIntV nodeTimeBucket;
   for (int i = 1; i < linkTimeBucketKeyVec.size(); i++) //i = 0 will have months with -1 (link formation date unknown)
   {
      int key = linkTimeBucketKeyVec[i];
      typedef unordered_multimap<int, pair<int, int>>::iterator Iter;
      pair<Iter, Iter> bounds;
      bounds = linksEdgeTimeMap.equal_range(key);
      Iter begin = bounds.first;
      nodeTimeBucket.Clr();
      while (begin != bounds.second)
      {
         int first = begin->second.first;
         int second = begin->second.second;
         const TInt v1 = begin->second.first;
         const TInt v2 = begin->second.second;
         TPair<TInt, TInt> pr(v1, v2);
         edgeTimeBucket.Add(pr);
         nodeTimeBucket.AddUnique(begin->second.first);
         nodeTimeBucket.AddUnique(begin->second.second);
         ++begin;
      }
      linkTimeBucketEdgeVec.AddUnique(edgeTimeBucket);
      linkTimeBucketNodeVec.AddUnique(nodeTimeBucket);
   }
}
void GraphU::analyzeTemporalCharacteristics()
{
   int NodesBucket = 20;
   TIntV NodeIdV;
   TIntV degSeqV;
   TIntV degDistributionV;
   TIntPrV degV;
   //TTmUnit TmUnit = TTmUnit::tmuNodes;
   //TTmUnit TmUnit = TTmUnit::tmuMonth;
   TTmUnit TmUnit = TTmUnit::tmuMonth;
   TFltPrV DegToCCfV, CcfV, OpClV, OpV, DegToCntV;
   int XVal = 0;
   TTuple<TFlt, 4> Tuple;
   TVec<TTuple<TFlt, 4> > OpenClsV;
   TExeTm ExeTm;
   PNGraph NGraph = TSnap::ConvertGraph<PNGraph>(m_graph);
   m_graph->GetTmBuckets(TmUnit, TmBucketV);
   PGStatVec growthStatsVector = m_graph->TimeGrowth(tmuMonth, TGStat::AllStat(), TmBucketV.BegI()->BegTm);
   growthStatsVector->SaveTxt("IPGraphGrowthProp", "Temporal Graph Growth Properties");
   TIntFltKdV DistNbrsV;
   m_graph->PlotEffDiam("EffecDiam", "Effective Diameter for temporal graph", tmuMonth, TmBucketV[0].BegTm, 5, false, false);
   for (int t = 0; t < TmBucketV.Len(); t++)
   {
      std::printf("\r  %d/%d: ", t + 1, TmBucketV.Len());
      NodeIdV.AddV(TmBucketV[t].NIdV); // edges up to time T
      int numNodesInTimeslot = NodeIdV.Len();
      nodesInTimeSlot.push_back(numNodesInTimeslot);

      int64 Open = 0, Close = 0;
      const PNGraph Graph = TSnap::GetSubGraph(NGraph, NodeIdV);
      edgesInTimeSlot.push_back(Graph->GetEdges());
      //int effectiveDiameter = TSnap::GetBfsEffDiam(Graph, Graph->GetNodes(), false);
      TSnap::GetAnf(Graph, DistNbrsV, -1, false, 32);
      double effectiveDiamD = TSnap::TSnapDetail::CalcEffDiam(DistNbrsV, 0.9);
      effecDiamInTimeSlot.push_back(effectiveDiamD);
      double maxSccSize = TSnap::GetMxSccSz(Graph);
      double maxSccFrac = double(maxSccSize / numNodesInTimeslot);
      maxSccNodesFraction.push_back(maxSccFrac);
      TSnap::GetDegSeqV(Graph, degSeqV); //degSeqV contains degrees of all nodes in the graph till now
      double averageDegree = 0.0;
      for (int i = 0; i < degSeqV.Len(); i++)
      {
         averageDegree += degSeqV[i]();
      }
      averageDegree /= degSeqV.Len();
      degreeInTimeSlot.push_back(averageDegree);

      TSnap::GetDegCnt(Graph, DegToCntV);

      const double CCf = TSnap::GetClustCf(Graph, DegToCCfV, Open, Close);
      if (TmUnit == tmuNodes)
      {
         XVal = Graph->GetNodes();
      }
      else
      {
         XVal = TmBucketV[t].BegTm.GetInUnits(TmUnit);
      }
      CcfV.Add(TFltPr(XVal, CCf));
      OpClV.Add(TFltPr(XVal, Open + Close == 0 ? 0 : Close / (Open + Close)));
      OpV.Add(TFltPr(XVal, Open == 0 ? 0 : Close / Open));
      Tuple[0] = Graph->GetNodes();
      Tuple[1] = Graph->GetEdges();
      Tuple[2] = Close;  Tuple[3] = Open;
      OpenClsV.Add(Tuple);
   }
}

void GraphU::generateTemoralGraphModels()
{
   //generates a forest fire model with forward probability = 0.37
   //and backward probablity = 0.32 and 0.33

   int numNodes;
   int numEdges;
   double effDiam;
   vector<PNGraph> forestFireGraphV;

   //********************************************
   TFfGGen::GenFFGraphs(0.38, 0.35, "ForestFire.plt");
   //********************************************

   //for (int i = 0; i < 9; i++)
   //{
   //   PNGraph modelGraphSnapShot = TSnap::GenForestFire(10 + i * 10, 0.37, 0.33);
   //   //forestFireGraphV.push_back(modelGraphSnapShot);
   //   numNodes = modelGraphSnapShot->GetNodes();
   //   numEdges = modelGraphSnapShot->GetEdges();
   //   effDiam = TSnap::GetBfsEffDiam(modelGraphSnapShot, numNodes, false);
   //   FFMEffecDiamV1.push_back(effDiam);
   //   FFMDensificationV1.push_back(make_pair(numNodes, numEdges));
   //}

   //for (int i = 1; i < 5; i++)
   //{
   //   PNGraph modelGraphSnapShot = TSnap::GenForestFire(100 + i * 200, 0.37, 0.33);
   //   //forestFireGraphV.push_back(modelGraphSnapShot);
   //   numNodes = modelGraphSnapShot->GetNodes();
   //   numEdges = modelGraphSnapShot->GetEdges();
   //   effDiam = TSnap::GetBfsEffDiam(modelGraphSnapShot, numNodes, false);
   //   FFMEffecDiamV1.push_back(effDiam);
   //   FFMDensificationV1.push_back(make_pair(numNodes, numEdges));
   //}

   //for (int i = 1; i < 10; i++)
   //{
   //   PNGraph modelGraphSnapShot = TSnap::GenForestFire(1000 + i*1000, 0.37, 0.33);
   //   //forestFireGraphV.push_back(modelGraphSnapShot);
   //   numNodes = modelGraphSnapShot->GetNodes();
   //   numEdges = modelGraphSnapShot->GetEdges();
   //   effDiam = TSnap::GetBfsEffDiam(modelGraphSnapShot, numNodes, false);
   //   FFMEffecDiamV1.push_back(effDiam);
   //   FFMDensificationV1.push_back(make_pair(numNodes, numEdges));
   //}

   ////PNGraph forestFireGraph2 = TSnap::GenForestFire(10000, 0.38, 0.35);
   //for (int i = 0; i < 9; i++)
   //{
   //   PNGraph modelGraphSnapShot = TSnap::GenForestFire(10 + i * 10, 0.38, 0.35);
   //   //forestFireGraphV.push_back(modelGraphSnapShot);
   //   numNodes = modelGraphSnapShot->GetNodes();
   //   numEdges = modelGraphSnapShot->GetEdges();
   //   effDiam = TSnap::GetBfsEffDiam(modelGraphSnapShot, numNodes, false);
   //   FFMEffecDiamV2.push_back(effDiam);
   //   FFMDensificationV2.push_back(make_pair(numNodes, numEdges));
   //}

   //for (int i = 1; i < 5; i++)
   //{
   //   PNGraph modelGraphSnapShot = TSnap::GenForestFire(100 + i * 200, 0.38, 0.35);
   //   //forestFireGraphV.push_back(modelGraphSnapShot);
   //   numNodes = modelGraphSnapShot->GetNodes();
   //   numEdges = modelGraphSnapShot->GetEdges();
   //   effDiam = TSnap::GetBfsEffDiam(modelGraphSnapShot, numNodes, false);
   //   FFMEffecDiamV2.push_back(effDiam);
   //   FFMDensificationV2.push_back(make_pair(numNodes, numEdges));
   //}

   //for (int i = 1; i < 10; i++)
   //{
   //   PNGraph modelGraphSnapShot = TSnap::GenForestFire(1000 + i * 1000, 0.38, 0.35);
   //   //forestFireGraphV.push_back(modelGraphSnapShot);
   //   numNodes = modelGraphSnapShot->GetNodes();
   //   numEdges = modelGraphSnapShot->GetEdges();
   //   effDiam = TSnap::GetBfsEffDiam(modelGraphSnapShot, numNodes, false);
   //   FFMEffecDiamV2.push_back(effDiam);
   //   FFMDensificationV2.push_back(make_pair(numNodes, numEdges));
   //}
   //m_modelGraphsVector.push_back(forestFireGraph1);
   //m_modelGraphsVector.push_back(forestFireGraph2);
}

void GraphU::setGraph(PTimeNet incomingG)
{
   m_graph = incomingG;
}

PTimeNet GraphU::getGraph()
{
   return m_graph;
}

const vector<bool>* GraphU::getTargetConceptVec()
{
   return &m_target_concept_vec;
}

const vector<double>* GraphU::getEffDiamInTimeSlot()
{
   return &effecDiamInTimeSlot;
}

const vector<double>* GraphU::getdegreeInTimeSlot()
{
   return &degreeInTimeSlot;
}


const vector<double>* GraphU::getMaxSccNodesInTimeSlot()
{
   return &maxSccNodesFraction;
}

void GraphU::setTargetConceptVec(const vector<bool>* tcv)
{
   m_target_concept_vec = *tcv;
}

vector<vector<double>>* GraphU::getNonTemporalFeatureSet()
{
   return &m_non_temporal_fv;
}

vector<vector<double>>* GraphU::getTemporalFeatureSet()
{
   return &m_temporal_fv;
}

void GraphU::setNonTemporalFeatureSet(vector<vector<double>>* fs)
{
   m_non_temporal_fv = *fs;
}

int GraphU::addNode(int NodeId)
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

//int GraphU::getNodeIndexGivenNodeId(int id)
//{
//   int index = 0;
//   hash_map<int, int>::iterator iter;
//   iter = m_node_hash.find(id);
//   if (iter == m_node_hash.end()) //nodeId not present return -1
//   {
//      index = -1;
//   }
//   else
//   {
//      iter = m_node_hash.find(id);
//      index = iter->second;
//   }
//   return index;
//}

//implemented for more control and flexibility 
void GraphU::readGraph(string linksFile, string messagesFile)
{
   fstream input_file_stream;
   string single_line = "";
   input_file_stream.open(messagesFile, std::fstream::in);
   m_ungraph = PNGraph::New();
   m_graph = PTimeNet::New();
   m_msg_graph = PNGraph::New();
   bool isFromPresent = false;
   bool isToPresent = false;
   int Day = 0;
   int Hour = 0;
   int Min = 0;
   int Sec = 0;

   if (input_file_stream.is_open())
   {
      int fromNodeId;
      int toNodeId;
      int weight;
      time_t ts;
      string timeStamp;
      string ts_s;
      int fromNodeIndex = 0;
      int toNodeIndex = 0;
      int i = 1;
      int adj_len = 0;
      int adj_index = 0;
      string ts_mmmyyyy;
      string month;
      int Year = 0;
      int Month = 0;
      int YYYYMM = 0;
      char MonthStr[4];
      while (getline(input_file_stream, single_line)) /*!input_file_stream.eof()*/
      {
         stringstream nodeValuesStream(single_line);
         if (single_line.back() != 'N')
         {
            nodeValuesStream >> fromNodeId >> toNodeId >> ts;
            timeStamp = std::ctime(&ts);
            Year = atoi(timeStamp.substr(20, 4).c_str());
            month = timeStamp.substr(4, 3); //gives month mmm
            MonthStr[0] = toupper(month.at(0));  MonthStr[1] = tolower(month.at(1));
            MonthStr[2] = tolower(month.at(2));  MonthStr[3] = 0;
            Month = TTmInfo::GetMonthN(MonthStr, lUs);
            YYYYMM = Year * 100 + Month;
            Day = atoi(timeStamp.substr(8, 2).c_str());
            Hour = atoi(timeStamp.substr(11, 2).c_str());
            Min = atoi(timeStamp.substr(14, 2).c_str());
            Sec = atoi(timeStamp.substr(17, 2).c_str());
            isFromPresent = m_graph->IsNode(fromNodeId);
            //if (fromNodeId != toNodeId)
            {
               if (!isFromPresent)
               {
                  m_msg_graph->AddNode(fromNodeId);
                  m_graph->AddNode(fromNodeId, TSecTm(Year, Month, Day, Hour, Min, Sec));
                  //log the initial activity- this is also the approx. time of joining the network
                  m_birthTime_hash.AddDat(fromNodeId, TSecTm(Year, Month, Day, Hour, Min, Sec));
               }
               //log the nodes latest activity every time
               m_recent_activity_hash.insert({ fromNodeId, YYYYMM });
               isToPresent = m_graph->IsNode(toNodeId);
               if (!isToPresent)
               {
                  m_msg_graph->AddNode(toNodeId);
                  m_graph->AddNode(toNodeId, TSecTm(Year, Month, Day, Hour, Min, Sec));
                  m_birthTime_hash.AddDat(toNodeId, TSecTm(Year, Month, Day, Hour, Min, Sec));
               }
               m_recent_activity_hash.insert({ toNodeId, YYYYMM });
               if (!m_graph->IsEdge(fromNodeId, toNodeId))
               {
                  m_msg_graph->AddEdge(fromNodeId, toNodeId);
                  m_graph->AddEdge(fromNodeId, toNodeId);
               }
               //ts_mmmyyyy.clear();
               //ts_mmmyyyy.append(timeStamp.substr(4, 3));
               //ts_mmmyyyy.append(timeStamp.substr(20, 4));
               //TStr mmmyyyy = TStr(ts_mmmyyyy.c_str());
               if (std::find(messageTimeBucketKeyVec.begin(), messageTimeBucketKeyVec.end(), YYYYMM) == messageTimeBucketKeyVec.end())
               {
                  messageTimeBucketKeyVec.push_back(YYYYMM);
               }
               messagesEdgeTimeMap.insert({ YYYYMM, make_pair(fromNodeId, toNodeId) });
            }
         }
         nodeValuesStream.clear();
      }
      input_file_stream.close();
      int numEdges = m_graph->GetEdges();
      int numNodes = m_graph->GetNodes();
      std::printf("\r  %d # messages in edge time map \n", messagesEdgeTimeMap.size());
      std::printf("\r  %d unique nodes in messages\n", numNodes);
      std::printf("\r  %d messages read\n", numEdges);
   }
   input_file_stream.open(linksFile, std::fstream::in);
   if (input_file_stream.is_open())
   {
      int fromNodeId;
      int toNodeId;
      time_t ts;
      string timeStamp;
      string ts_s;
      single_line.clear();
      int Year;
      string month;
      int Month = 0;
      int YYYYMM = 0;
      char MonthStr[4];

      while (getline(input_file_stream, single_line)) /*!input_file_stream.eof()*/
      {
         stringstream nodeValuesStream(single_line);
         if (single_line.back() != 'N')
         {
            nodeValuesStream >> fromNodeId >> toNodeId >> ts;
            timeStamp = std::ctime(&ts);
            Year = atoi(timeStamp.substr(20, 4).c_str());
            month = timeStamp.substr(4, 3); //gives month mmm
            MonthStr[0] = toupper(month.at(0));  MonthStr[1] = tolower(month.at(1));
            MonthStr[2] = tolower(month.at(2));  MonthStr[3] = 0;
            Month = TTmInfo::GetMonthN(MonthStr, lUs);
            YYYYMM = Year * 100 + Month;
            Day = atoi(timeStamp.substr(8, 2).c_str());
            Hour = atoi(timeStamp.substr(11, 2).c_str());
            Min = atoi(timeStamp.substr(14, 2).c_str());
            Sec = atoi(timeStamp.substr(17, 2).c_str());
         }
         else
         {
            nodeValuesStream >> fromNodeId >> toNodeId;
            YYYYMM = -1;
         }
         isFromPresent = m_ungraph->IsNode(fromNodeId);
         isToPresent = m_ungraph->IsNode(toNodeId);
         if (!isFromPresent)
         {
            m_ungraph->AddNode(fromNodeId);
         }
         if (!isToPresent)
         {
            m_ungraph->AddNode(toNodeId);
         }
         if (!m_ungraph->IsEdge(fromNodeId, toNodeId))
         {
            m_ungraph->AddEdge(fromNodeId, toNodeId);
            TPair<TInt, TInt> pr;
            pr.Val1 = fromNodeId;
            pr.Val2 = toNodeId;
            linksFullTimeMap.AddDat(pr, TSecTm(Year, Month, Day, Hour, Min, Sec));
         }

         if (find(linkTimeBucketKeyVec.begin(), linkTimeBucketKeyVec.end(), YYYYMM) == linkTimeBucketKeyVec.end())
         {
            linkTimeBucketKeyVec.push_back(YYYYMM);
         }
         linksEdgeTimeMap.insert({ YYYYMM, make_pair(fromNodeId, toNodeId) });
         nodeValuesStream.clear();
      }
      input_file_stream.close();
      int numEdges = m_ungraph->GetEdges();
      int numNodes = m_ungraph->GetNodes();
      std::printf("\n=================================\n");
      std::printf("\r  %d unique nodes in Network..\n", numNodes);
      std::printf("\r  %d friendships recorded\n", numEdges);
   }
   else
   {
      cout << "\nError! Input file(s) cannot be opened or it's not found!\n";
   }
}

int GraphU::calculateActiveDegree(int currentSnapshot, TIntV& neighborhoodVec)
{
   int neighborhoddSize = neighborhoodVec.Len();
   int activeDeg = 0;
   for (int i = 0; i < neighborhoddSize; i++) //search for all the neighbors' participation in wall posts
   {
      for (int n = 0; n < activityThreshold; n++) //search in previous n snapshots
      {
         if (currentSnapshot >= n)
         {
            if (messageTimeBucketNodeVec[currentSnapshot - n].Count(neighborhoodVec[i]) > 0)
            {
               activeDeg++;
               break; //Node participation confirmed in one of the past n snapshots, proceed to next neighbor
            }
         }
      }
   }
   return activeDeg;
}

int GraphU::calculateActiveCN(int i, TIntV& commonNeighborhoodVec) //returns number of active common neighbors
{
   int activeCommonNeigh = calculateActiveDegree(i, commonNeighborhoodVec);
   return activeCommonNeigh;
}

bool GraphU::isNodeActive(int currentSnapshot, int nodeId)
{
   bool isActive = false;
   for (int n = 0; n < activityThreshold; n++) //search in previous n snapshots
   {
      if (currentSnapshot >= n)
      {
         if (messageTimeBucketNodeVec[currentSnapshot - n].Count(nodeId) > 0)
         {
            isActive = true;
            return isActive; //Node participation confirmed in one of the past n snapshots, proceed to next neighbor
         }
      }
   }
   return isActive;
}

void GraphU::extractFeatureSet()
{
   //m_graph->GetTmBuckets(tmuMonth, TmBucketV);
   //PNGraph NGraph = TSnap::ConvertGraph<PTimeNet>(m_graph);

   TVec<TInt> NodeIdV;
   TVec<TIntPr> EdgeIdV;
   calculateMessageTimeBucketVector();
   calculateLinksTimeBucketVector();
   vector<double> featureVec;
   vector<pair<int, int>> posSamplesDyadsVec;
   vector<pair<int, int>> negSamplesDyadsVec;
   vector<pair<int, int>> hardNegSamplesDyadsVec;
   vector<pair<int, int>> testSamplesDyadsVec;
   for (int t = 0; t < linkTimeBucketNodeVec.Len(); t++)
   {
      //printf("\r %d/%d: ", t + 1, timeBucketNodeVec.Len());
      //NodeIdV.AddV(TmBucketV[t].NIdV); // edges up to time T
      //NodeIdV.AddVMerged(linkTimeBucketNodeVec[t]); // edges up to time T
      //NodeIdV.AddVMerged(messageTimeBucketNodeVec[t]); // edges up to time T
      //NodeIdV.Clr();
      NodeIdV.AddVMerged(linkTimeBucketNodeVec[t]); // edges up to time T
      EdgeIdV.AddVMerged(linkTimeBucketEdgeVec[t]);
      //NodeIdV.AddVMerged(messageTimeBucketNodeVec[t]); // edges up to time T

      int num_nodes_added_now = linkTimeBucketNodeVec[t].Len();
      if (t > 0)
      {
         num_nodes_added_now = linkTimeBucketNodeVec[t].Len() - linkTimeBucketNodeVec[t - 1].Len();
      }
      std::printf("\n curr bucket.. %d, nodes added now.. %d, all nodes( %d )", t, num_nodes_added_now, NodeIdV.Len());

      //The graph Snapshot at this point in time
      //We start from initial snapshot t=1 and move till last snapshot t=n
      //We progressively analyze the graph structure as it evolves in time
      //Idea is to take equal number of positive and negative examples from each snapshot
      //Forming set- If a pair of nodes has no link in current snapshot but has link in next snapshot then its +ve example
      //Unconn set- If a pair of nodes has no link in this snapshot and also no link in next snapshot then its -ve example
      //For each pair of nodes, calculate the temporal and non-temporal features:
      //graphSnapShot = TSnap::GetSubGraph(m_ungraph, NodeIdV);
      graphSnapShotPt = PNGraph::New();
      TNGraph& graphSnapShot = *graphSnapShotPt;
      graphSnapShot.Reserve(NodeIdV.Len(), EdgeIdV.Len());
      for (int nodes = 0; nodes < NodeIdV.Len(); nodes++)
      {
         if (!graphSnapShot.IsNode(NodeIdV[nodes]))
         {
            graphSnapShot.AddNode(NodeIdV[nodes]);
         }
      }
      for (int edges = 0; edges < EdgeIdV.Len(); edges++)
      {
         TIntPr pr = EdgeIdV[edges];
         if (!graphSnapShot.IsEdge(pr.Val1, pr.Val2))
         {
            graphSnapShot.AddEdge(pr.Val1,pr.Val2);
         }
      }
      //graphSnapShot = TSnap::GetSubGraph(m_msg_graph, NodeIdV);
      snapShotVec.Add(graphSnapShotPt);
      std::printf("\n # nodes in this snapshot= %d", graphSnapShotPt->GetNodes());
      std::printf("\n # edges in this snapshot= %d", graphSnapShotPt->GetEdges());
   }
   for (int i = 2; i < snapShotVec.Len() - 1; i++)
   {
      PNGraph currentSnapShot = snapShotVec[i];
      PNGraph previousSnapShot = snapShotVec[i - 1];
      PNGraph SecondPreviousSnapShot = snapShotVec[i - 2];
      int numPosExamples = 0;
      int numNegExamples = 0;
      PNGraph nextSnapShot = snapShotVec[i + 1];

      TNGraph::TEdgeI nxt_beginEdgeI = nextSnapShot->BegEI();
      TNGraph::TEdgeI nxt_endEdgeI = nextSnapShot->EndEI();

      //while (numPosExamples < 1000)
      //{
      //   int fromNode = currentSnapShot->GetRndNId();
      //   int toNode = currentSnapShot->GetRndNId();
      //   if (currentSnapShot->IsEdge(fromNode, toNode))
      //   {
      //      posSamplesDyadsVec.push_back(make_pair(fromNode, toNode));
      //      numPosExamples++;
      //   }
      //}

      TNGraph::TEdgeI edgeIter = previousSnapShot->BegEI();
      TNGraph::TEdgeI endEdgeI = previousSnapShot->EndEI();

      while (edgeIter != endEdgeI)
      {
         int fromNode = edgeIter.GetSrcNId();
         int toNode = edgeIter.GetDstNId();

         if (!SecondPreviousSnapShot->IsEdge(fromNode, toNode) &&
            previousSnapShot->IsEdge(fromNode, toNode))
         {
            posSamplesDyadsVec.push_back(make_pair(fromNode, toNode));
            numPosExamples++;
         }
         edgeIter++;
      }

      edgeIter = currentSnapShot->BegEI();
      endEdgeI = currentSnapShot->EndEI();

      while (edgeIter != endEdgeI)
      {
         int fromNode = edgeIter.GetSrcNId();
         int toNode = edgeIter.GetDstNId();

         if (!previousSnapShot->IsEdge(fromNode, toNode) && 
             currentSnapShot->IsEdge(fromNode, toNode))
         {
            posSamplesDyadsVec.push_back(make_pair(fromNode, toNode));
            numPosExamples++;
         }
         edgeIter++;
      }

      while (numNegExamples < numPosExamples)
      {
         const int& fromNode = currentSnapShot->GetRndNId();
         const int& toNode = currentSnapShot->GetRndNId();
         if (!currentSnapShot->IsEdge(fromNode, toNode)) //there is no edge in current snapshot
         {
            negSamplesDyadsVec.push_back(make_pair(fromNode, toNode));
            numNegExamples++;
         }
      }
      numNegExamples = 0;
      while (numNegExamples < numPosExamples)
      {
         TIntV neighborhoodVec;
         const int& fromNode = currentSnapShot->GetRndNId();
         int fromNode2HopNeighbors = TSnap::GetNodesAtHop(currentSnapShot, fromNode, 2, neighborhoodVec, true);
         if (neighborhoodVec.Len() > 0 && !currentSnapShot->IsEdge(fromNode, neighborhoodVec[0])) //there is no edge in current snapshot
         {
            hardNegSamplesDyadsVec.push_back(make_pair(fromNode, neighborhoodVec[0]));
            numNegExamples++;
         }
      }
      /*while (nxt_beginEdgeI != nxt_endEdgeI)
      {
         int fromNode = beginEdgeI.GetSrcNId();
         int toNode = beginEdgeI.GetDstNId();
         testSamplesDyadsVec.push_back(make_pair(fromNode, toNode));
         nxt_beginEdgeI++;
      }*/
      positiveSamplesVectorsVec.push_back(posSamplesDyadsVec);
      negativeSamplesVectorsVec.push_back(negSamplesDyadsVec);
      //testSamplesVectorsVec.push_back(testSamplesDyadsVec);
      hardNegativeSamplesVectorsVec.push_back(hardNegSamplesDyadsVec);
      posSamplesDyadsVec.clear();
      negSamplesDyadsVec.clear();
      testSamplesDyadsVec.clear();
      hardNegSamplesDyadsVec.clear();
   }


   //==================CALCULATE TEMPORAL + NONTEMPORAL FEATURES==================//

   //for (int i = 11; i < positiveSamplesVectorsVec.size(); i++) //i denotes the time
   for (int i = 22; i < 23; i++) //i denotes the time
   {
      for (int j = 0; j < positiveSamplesVectorsVec[i].size(); j++) //j denotes the jth sample dyad (pair) at ith snap shot
      {
         int fromNode = positiveSamplesVectorsVec[i][j].first;
         int toNode = positiveSamplesVectorsVec[i][j].second;
         calculateFeatures(i, fromNode, toNode, true, false);
      }
      //for (int j = 0; j < negativeSamplesVectorsVec[i].size(); j++) //j denotes the jth sample dyad (pair) at ith snap shot
      //{
      //   int fromNode = negativeSamplesVectorsVec[i][j].first;
      //   int toNode = negativeSamplesVectorsVec[i][j].second;
      //   calculateFeatures(i, fromNode, toNode, false, false);
      //}
         for (int j = 0; j < hardNegativeSamplesVectorsVec[i].size(); j++) //j denotes the jth sample dyad (pair) at ith snap shot
         {
            int fromNode = hardNegativeSamplesVectorsVec[i][j].first;
            int toNode = hardNegativeSamplesVectorsVec[i][j].second;
            calculateFeatures(i, fromNode, toNode, false, false);
         }
   }

   //temporal metrics
   //for (int i = 11; i < positiveSamplesVectorsVec.size(); i++) //i denotes the time
   for (int i = 22; i < 23; i++) //i denotes the time
   {
      for (int j = 0; j < positiveSamplesVectorsVec[i].size(); j++) //j denotes the jth sample dyad (pair) at ith snap shot
      {
         int fromNode = positiveSamplesVectorsVec[i][j].first;
         int toNode = positiveSamplesVectorsVec[i][j].second;
         calculateMovingAverages(i, fromNode, toNode, true);
      }
      //for (int j = 0; j < negativeSamplesVectorsVec[i].size(); j++) //j denotes the jth sample dyad (pair) at ith snap shot
      //{
      //   int fromNode = negativeSamplesVectorsVec[i][j].first;
      //   int toNode = negativeSamplesVectorsVec[i][j].second;
      //   calculateMovingAverages(i, fromNode, toNode, false);
      //}
         for (int j = 0; j < hardNegativeSamplesVectorsVec[i].size(); j++) //j denotes the jth sample dyad (pair) at ith snap shot
         {
            int fromNode = hardNegativeSamplesVectorsVec[i][j].first;
            int toNode = hardNegativeSamplesVectorsVec[i][j].second;
            calculateMovingAverages(i, fromNode, toNode, false);
         }
   }
}


void GraphU::calculateMovingAverages(int i, int fromNode, int toNode, bool isPositiveExample)
{
   vector<double> movingAvg2V;
   vector<double> movingAvg5V;
   vector<double> movingAvg10V;
   vector<double> temporalFeatureVector;
   //moving average of window size 10, 5, 2
   bool is_calc_b = false;
   for (int k = 1; k <= 2; k++)
   {
      is_calc_b = calculateFeatures(i - k, fromNode, toNode, isPositiveExample, true);
   }

   if (is_calc_b)
   {
      int size = m_placeholder_fv[0].size();
      movingAvg2V.resize(size);
      movingAvg5V.resize(size);
      movingAvg10V.resize(size);
      // m_placeholder_fv has metrics for t = i, i - 1 and i - 2 snapshots
      for (int l = 0; l < m_placeholder_fv.size(); l++)
      {
         vector<double>* featureVec = &m_placeholder_fv[l];
         for (int m = 0; m < featureVec->size(); m++)
         {
            movingAvg2V[m] = movingAvg2V[m] + featureVec->at(m);
         }
      }
      for (int m = 0; m < movingAvg2V.size(); m++)
      {
         if (movingAvg2V[m] > 0)
         {
            movingAvg2V[m] = double(movingAvg2V[m] / 2.0);
         }
         else
         {
            movingAvg2V[m] = 0.0;
         }
         temporalFeatureVector.push_back(movingAvg2V[m]);
      }
   }
   else //if features were not calculated then insert sequence of 0.0s
   {
      for (int m = 0; m < movingAvg2V.size(); m++)
      {
         temporalFeatureVector.push_back(0.0);
      }
   }

   is_calc_b = false;
   for (int k = 3; k <= 5; k++)
   {
      is_calc_b = calculateFeatures(i - k, fromNode, toNode, isPositiveExample, true);
   }

   if (is_calc_b)
   {
      // m_placeholder_fv has metrics for t = i, (i - 1).... (i - 5) snapshots
      for (int l = 0; l < m_placeholder_fv.size(); l++)
      {
         vector<double>* featureVec = &m_placeholder_fv[l];
         for (int m = 0; m < featureVec->size(); m++)
         {
            movingAvg5V[m] = movingAvg5V[m] + featureVec->at(m);
         }
      }
      for (int m = 0; m < movingAvg5V.size(); m++)
      {
         if (movingAvg5V[m] > 0)
         {
            movingAvg5V[m] = double(movingAvg5V[m] / 5.0);
         }
         else
         {
            movingAvg5V[m] = 0.0;
         }
         temporalFeatureVector.push_back(movingAvg5V[m]);
      }
   }
   else
   {
      for (int m = 0; m < movingAvg5V.size(); m++)
      {
         temporalFeatureVector.push_back(0.0);
      }
   }

   is_calc_b = false;
   for (int k = 6; k <= 10; k++)
   {
      is_calc_b = calculateFeatures(i - k, fromNode, toNode, isPositiveExample, true);
   }

   if (is_calc_b)
   {
      // m_placeholder_fv has metrics for t = i, (i - 1).... (i - 10) snapshots
      for (int l = 0; l < m_placeholder_fv.size(); l++)
      {
         vector<double>* featureVec = &m_placeholder_fv[l];
         for (int m = 0; m < featureVec->size(); m++)
         {
            movingAvg10V[m] = movingAvg10V[m] + featureVec->at(m);
         }
      }
      for (int m = 0; m < movingAvg10V.size(); m++)
      {
         if (movingAvg10V[m] > 0)
         {
            movingAvg10V[m] = double(movingAvg10V[m] / m_placeholder_fv.size());
         }
         else
         {
            movingAvg10V[m] = 0.0;
         }
         temporalFeatureVector.push_back(movingAvg10V[m]);
      }
   }
   else
   {
      for (int m = 0; m < movingAvg10V.size(); m++)
      {
         temporalFeatureVector.push_back(0.0);
      }
   }

   //19. Link delay measure d(u,v) = T(u,v) - max(b(u), b(v))
   TPair<TInt, TInt> pr;
   pr.Val1 = fromNode;
   pr.Val2 = toNode;
   TSecTm time = linksFullTimeMap.GetDat(pr);
   int linkCreationTime = time.GetYearN() + time.GetMonthN() + time.GetDayN();

   time = m_birthTime_hash.GetDat(fromNode);
   int fromNodeBirth = time.GetYearN() + time.GetMonthN() + time.GetDayN();

   time = m_birthTime_hash.GetDat(toNode);
   int toNodeBirth = time.GetYearN() + time.GetMonthN() + time.GetDayN();

   int linkDelay = linkCreationTime - max(fromNodeBirth, toNodeBirth);
   if (linkDelay < 0)
   {
      temporalFeatureVector.push_back(-1);
   }
   else
   {
      temporalFeatureVector.push_back(linkDelay);
   }

   m_temporal_fv.push_back(temporalFeatureVector);
   temporalFeatureVector.clear();
   m_placeholder_fv.clear(); //clear the placeholder vector and make it available for new sample pair
}

bool GraphU::calculateFeatures(int i, int fromNode, int toNode, bool isPositiveInstance, bool isTemporal)
{
   TVec<TInt> fromNodeNeighVec;
   TVec<TInt> toNodeNeighVec;
   TVec<TInt> neighborsIntrsVec;
   vector<double> featureVector;
   bool has_calculated_b = false;

   //============== CALCULATE THE MONADIC FEATURES ============//
   // 1. Degree of nodes
   // 2. Recency of nodes

   if (snapShotVec[i]->IsNode(fromNode) && snapShotVec[i]->IsNode(toNode))
   {
      int inDegreeFromNode = snapShotVec[i]->GetNI(fromNode).GetInDeg();
      int outDegreeFromNode = snapShotVec[i]->GetNI(fromNode).GetOutDeg();
      int inDegreeToNode = snapShotVec[i]->GetNI(toNode).GetInDeg();
      int outDegreeToNode = snapShotVec[i]->GetNI(toNode).GetOutDeg();
      int recencyFromNode = calculateRecency(i, fromNode);
      int recencyToNode = calculateRecency(i, toNode);

      if (!isTemporal)
      {
         featureVector.push_back(i);
         featureVector.push_back(fromNode);
         featureVector.push_back(toNode);
      }

      featureVector.push_back(inDegreeFromNode);
      featureVector.push_back(outDegreeFromNode);
      featureVector.push_back(inDegreeToNode);
      featureVector.push_back(outDegreeToNode);
      featureVector.push_back(recencyFromNode);
      featureVector.push_back(recencyToNode);

      //============== CALCULATE THE DYADIC FEATURES ============//
      TIntV neighborhoodVec;
      TIntV fromNodeNeighborVec;
      TIntV toNodeNeighborVec;
      // 10. Common Neighbors
      int numCommonNeigh = TSnap::GetCmnNbrs(m_graph, fromNode, toNode, neighborhoodVec);
      featureVector.push_back(numCommonNeigh);
      int activeCommonNeigh = calculateActiveCN(i, neighborhoodVec);

      // 11. Size2F
      int fromNodeSize2F = TSnap::GetNodesAtHop(snapShotVec[i], fromNode, 2, neighborhoodVec, true);
      int toNodeSize2F = TSnap::GetNodesAtHop(snapShotVec[i], toNode, 2, neighborhoodVec, true);
      featureVector.push_back(fromNodeSize2F);
      featureVector.push_back(toNodeSize2F);

      // 12. Approximate Katz coeff
      int approxKatzMeasure = numEdgesAmongNeighbors(snapShotVec[i], fromNode, toNode);
      featureVector.push_back(approxKatzMeasure);

      //13. Active Approximate Katz Coeff (Active Friend's measure)
      int activeApproxKatz = numEdgesAmongActiveNeighbors(snapShotVec[i], i, fromNode, toNode);

      //14. preferential attachment score-
      //estimates how "rich" the two nodes are by multiplying # neighbors (total degrees) of each node
      int prefAttachScore = (inDegreeFromNode + outDegreeFromNode)*(inDegreeToNode + outDegreeToNode);
      featureVector.push_back(prefAttachScore);
      
      TSnap::GetNodesAtHop(snapShotVec[i], fromNode, 1, fromNodeNeighborVec);
      TSnap::GetNodesAtHop(snapShotVec[i], toNode, 1, toNodeNeighborVec);
      int activeFromDegree = calculateActiveDegree(i, fromNodeNeighborVec);
      int activeToDegree = calculateActiveDegree(i, toNodeNeighborVec);
      //15. Active PA Score
      int activePAScore = activeFromDegree * activeToDegree;

      //populate fromNodeNeighVec
      fromNodeNeighVec.Clr();
      const TNGraph::TNodeI fromNodeIter = snapShotVec[i]->GetNI(fromNode);
      for (int edge = 0; edge < fromNodeIter.GetDeg(); edge++)
      {
         const int nodeId = fromNodeIter.GetNbrNId(edge);
         fromNodeNeighVec.Add(nodeId);
      }

      //populate toNodeNeighVec
      toNodeNeighVec.Clr();
      const TNGraph::TNodeI toNodeIter = snapShotVec[i]->GetNI(toNode);
      for (int edge = 0; edge < toNodeIter.GetDeg(); edge++)
      {
         const int nodeId = toNodeIter.GetNbrNId(edge);
         toNodeNeighVec.Add(nodeId);
      }

      //16. Adamic/Adar measure
      fromNodeNeighVec.Intrs(toNodeNeighVec, neighborsIntrsVec);
      double aaMeasure = 0.0;
      int deg = 0;
      double temp = 0.0;
      for (int k = 0; k < neighborsIntrsVec.Len(); k++)
      {
         TNGraph::TNodeI nodeIter = snapShotVec[i]->GetNI(neighborsIntrsVec[k]);
         deg = nodeIter.GetDeg();
         temp = 1 / double(log(deg));
         aaMeasure += temp;
      }
      featureVector.push_back(aaMeasure);
      //17. Active AA measure
      temp = 0.0;
      int activeDeg = 0;
      double activeAAMeasure = 0.0;
      for (int k = 0; k < neighborsIntrsVec.Len(); k++)
      {
         neighborhoodVec.Clr();
         int node = neighborsIntrsVec[k];
         TSnap::GetNodesAtHop(snapShotVec[i], node, 1, neighborhoodVec);
         activeDeg = calculateActiveDegree(i, neighborhoodVec);
         temp = 1 / double(log(activeDeg));
         activeAAMeasure += temp;
      }

      // 18. Jackard's coeff
      //Higher value of JC denotes stronger tie between two nodes
      //Ratio between # of common friends to # of total friends
      double sumFromTodegree = inDegreeFromNode + outDegreeFromNode + inDegreeToNode + outDegreeToNode;
      double jackardsCoeff = 0.0;
      if (sumFromTodegree > 0)
      {
         jackardsCoeff = numCommonNeigh / sumFromTodegree;
      }
      featureVector.push_back(jackardsCoeff);

      //19. Active Jackard's coeff
      double sumActiveDegree = activeFromDegree + activeToDegree;
      double activeJackardsCoeff = 0.0;
      if (sumActiveDegree > 0)
      {
         activeJackardsCoeff = activeCommonNeigh / sumActiveDegree;
      }

      //20. Num of triads From Node participates in
      int numTriadsFromN = TSnap::GetNodeTriads(m_ungraph, fromNode);
      featureVector.push_back(numTriadsFromN);

      //21. Num of triads To Node participates in
      int numTriadsToN = TSnap::GetNodeTriads(m_ungraph, toNode);
      featureVector.push_back(numTriadsToN);

      featureVector.push_back(activeFromDegree);
      featureVector.push_back(activeToDegree);
      featureVector.push_back(activeAAMeasure);
      featureVector.push_back(activeApproxKatz);
      featureVector.push_back(activeCommonNeigh);
      featureVector.push_back(activeJackardsCoeff);

      if (!isTemporal)
      {
         m_target_concept_vec.push_back(isPositiveInstance);
         m_non_temporal_fv.push_back(featureVector);
         featureVector.clear();
      }
      else //if the metrics are calculated in temporal fashion store them in placeholder fv
         //the caller will process the contents after certain rounds 
         //and populate the temporal fv with resultant metrics and clear the placeholder fv
      {
         m_placeholder_fv.push_back(featureVector);
      }
      has_calculated_b = true;
   }
   return has_calculated_b;
}

int GraphU::calculateRecency(int currentSnapShot, int nodeId)
{
   int rec = 0;
   pair<multimap<int, int>::iterator, multimap<int, int>::iterator> bounds = m_recent_activity_hash.equal_range(nodeId);
   //200810, 200811,...all the snapshots in which the node was active
   multimap<int, int>::iterator begin = bounds.first;
   vector<int> activityTimesTillNow;
   while (begin != bounds.second)
   {
      if (keyMessageTimeMap[begin->second] <= currentSnapShot)
      {
         activityTimesTillNow.push_back(begin->second);
      }
      ++begin;
   }

   if (activityTimesTillNow.size() > 0)
   {
      rec = activityTimesTillNow[activityTimesTillNow.size() - 1];
   }
   else
   {
      return 0; //no activity till now
   }
   int recentActivitySnapShot = keyMessageTimeMap[rec]; //gives i
   rec = (currentSnapShot - recentActivitySnapShot) + 1;
   return rec;
   //0 = recentActivitySnapShot is in next time step
   // -ve means recentActivitySnapShot is in future time step or |rec| steps ahead
}


//void GraphU::extractNonTemporalFeatureSet(PTimeNet Graph) 
//{
//   TGStat* statistics = new TGStat(m_graph, (TSecTm)0);
//   int num_nodes = statistics->GetNodes();
//   PNGraph NGraph = TSnap::ConvertGraph<PNGraph>(Graph, true);
//   double avg_cc = TSnap::GetClustCf(m_graph,num_nodes);
//   TVec<TInt> degSequenceVec;
//   TVec <TPair<TInt, TInt>> degDistVec;
//   TVec <TPair<TInt, TInt>> nodeInDegreeVec;
//   TVec<TInt> commonNeighVec;
//   TSnap::GetDegCnt(m_graph, degDistVec);
//
//   //Computes the in-degree for every node in Graph. 
//   //The result is stored in NodeIdInDegreeVec, 
//   //a vector of pairs (node id, node in degree).
//   TSnap::GetNodeInDegV(m_graph, nodeInDegreeVec);
//
//   //Computes the degree sequence vector for nodes in Graph. 
//   //The degree sequence vector is stored in DegV.
//   //A vector containing the degree for each node in the graph.
//   TSnap::GetDegSeqV(m_graph, degSequenceVec);
//
//   int numCommonNeigh = commonNeighVec.Len();
//   cout << "\nstat op: num nodes= " << num_nodes;
//   cout << "\nstat op: avg cc= " << avg_cc;
//   int maxDegreeNodeId = TSnap::GetMxInDegNId(m_graph);
//   cout << "\Max degree= " << maxDegreeNodeId;
//   int numEdges = Graph->GetEdges();
//   int numNodes = Graph->GetNodes();
//   TVec<TInt> neighVec;
//   TVec<TInt> nodeIDVec;
//   TVec<TInt> fromNodeNeighVec;
//   TVec<TInt> toNodeNeighVec;
//   TVec<TInt> NeighborIntersectionVec;
//   PTimeNet fromNodeInducedSubGraph;
//   PTimeNet fromNodeInducedSubGraph_withSelf;
//   PTimeNet toNodeInducedSubGraph;
//   PTimeNet toNodeInducedSubGraph_withSelf;
//   PTimeNet linkSubGraph;
//   PTimeNet linkSubGraph_withSelf;
//   PTimeNet innerSubGraph;
//   TCnComV CommunityVec;
//   TCnComV CNMCommunityVec;
//   vector<int> featureVector;
//
//   int totalExamples = 0;
//   int numNegEx = 0;
//   bool enoughNegEx = false;
//   while (totalExamples < 6000)
//   {
//      if (numNegEx >= 3000)
//      {
//         enoughNegEx = true;
//      }
//      const int& fromNodeId = Graph->GetRndNId();
//      const int& toNodeId = Graph->GetRndNId();
//      //int fromNodeId = 2;
//      //int toNodeId = 7;
//      bool isEdgePresent = Graph->IsEdge(fromNodeId, toNodeId);
//      if ((enoughNegEx && isEdgePresent) ||
//         !enoughNegEx)
//      {
//         featureVector.clear();
//         featureVector.push_back(fromNodeId);
//         featureVector.push_back(toNodeId);
//
//         //Feature 1- number of common neighbors
//         int numCommonNeigh = TSnap::GetCmnNbrs(m_graph, fromNodeId, toNodeId, neighVec);
//         featureVector.push_back(numCommonNeigh);
//
//         //Feature 2- degrees of both nodes
//         int fromNodeDeg = Graph->GetNI(fromNodeId).GetDeg();
//         int toNodeDeg = Graph->GetNI(toNodeId).GetDeg();
//         featureVector.push_back(fromNodeDeg);
//         featureVector.push_back(toNodeDeg);
//
//         //Feature 3- Jackard's Coefficient*100
//         //Higher value of JC denotes stronger tie between two nodes
//         //Ratio between # of common friends to # of total friends
//         if ((fromNodeDeg + toNodeDeg) > 0)
//         {
//            double jackardsCoeff = double(numCommonNeigh) / (fromNodeDeg + toNodeDeg);
//            jackardsCoeff *= 100;
//            featureVector.push_back(jackardsCoeff);
//         }
//
//         //Feature 4- preferential attachment score-
//         //estimates how "rich" the two nodes are by multiplying # neighbors of each node
//         int prefAttachScore = fromNodeDeg*toNodeDeg;
//         featureVector.push_back(prefAttachScore);
//
//         //Feature 5- Approximated Katz Measure
//         // Original Katz measure is cubic complexity
//         //Number of edges among neighbors of u and v of an edge (u,v)
//         int approxKatzMeasure = numEdgesAmongNeighbors(TSnap::ConvertGraph<PNGraph, PTimeNet>(m_graph), fromNodeId, toNodeId);
//         featureVector.push_back(approxKatzMeasure);
//
//         //Feature 6- Same Community TODO
//
//         //Feature 7- Subgraphs: used by WJ Cukierski et al.
//         fromNodeNeighVec.Clr();
//         const TTimeNet::TNodeI fromNodeIter = Graph->GetNI(fromNodeId);
//         for (int edge = 0; edge < fromNodeIter.GetOutDeg(); edge++)
//         {
//            const int OutNId = fromNodeIter.GetOutNId(edge);
//            fromNodeNeighVec.Add(OutNId);
//         }
//         fromNodeInducedSubGraph = TSnap::GetSubGraph(m_graph, fromNodeNeighVec);
//         fromNodeNeighVec.AddUnique(fromNodeId);
//         fromNodeInducedSubGraph_withSelf = TSnap::GetSubGraph(m_graph, fromNodeNeighVec);
//         int fromNodeInducedSubgraphLinkNumber = fromNodeInducedSubGraph->GetEdges();
//         int fromNodeInducedSubgraphLinkNumber_withSelf = fromNodeInducedSubGraph_withSelf->GetEdges();
//
//         featureVector.push_back(fromNodeInducedSubgraphLinkNumber);
//         featureVector.push_back(fromNodeInducedSubgraphLinkNumber_withSelf);
//
//         toNodeNeighVec.Clr();
//         const TTimeNet::TNodeI toNodeIter = Graph->GetNI(toNodeId);
//         for (int edge = 0; edge < toNodeIter.GetOutDeg(); edge++)
//         {
//            const int OutNId = toNodeIter.GetOutNId(edge);
//            toNodeNeighVec.Add(OutNId);
//         }
//         toNodeInducedSubGraph = TSnap::GetSubGraph(m_graph, toNodeNeighVec);
//         toNodeNeighVec.AddUnique(toNodeId);
//         toNodeInducedSubGraph_withSelf = TSnap::GetSubGraph(m_graph, toNodeNeighVec);
//         int toNodeInducedSubGraphLinkNumber = toNodeInducedSubGraph->GetEdges();
//         int toNodeInducedSubGraphLinkNumber_withSelf = toNodeInducedSubGraph_withSelf->GetEdges();
//
//         featureVector.push_back(toNodeInducedSubGraphLinkNumber);
//         featureVector.push_back(toNodeInducedSubGraphLinkNumber_withSelf);
//
//         //Link Subgraph Features
//         TVec<TInt> allNodes;
//         //neighVec.AddV(fromNodeNeighVec);
//         TVec<TInt> uniqueNeighborsVec;
//         TVec<TInt> neighborsIntrsVec;
//         TVec<TInt> unionNeighborsVec;
//         toNodeNeighVec.Clr();
//         fromNodeNeighVec.Clr();
//         for (int edge = 0; edge < toNodeIter.GetOutDeg(); edge++)
//         {
//            const int OutNId = toNodeIter.GetOutNId(edge);
//            toNodeNeighVec.AddUnique(OutNId);
//         }
//         for (int edge = 0; edge < fromNodeIter.GetOutDeg(); edge++)
//         {
//            const int OutNId = fromNodeIter.GetOutNId(edge);
//            fromNodeNeighVec.AddUnique(OutNId);
//         }
//
//         for (int n = 0; n < toNodeIter.GetOutDeg(); n++)
//         {
//            uniqueNeighborsVec.AddUnique(toNodeNeighVec[n]);
//         }
//         for (int n = 0; n < fromNodeIter.GetOutDeg(); n++)
//         {
//            uniqueNeighborsVec.AddUnique(fromNodeNeighVec[n]);
//         }
//         vector<int> nv;
//         for (int z = 0; z < uniqueNeighborsVec.Len(); z++)
//         {
//            nv.push_back(uniqueNeighborsVec[z]);
//         }
//
//         linkSubGraph = TSnap::GetSubGraph(m_graph, uniqueNeighborsVec);
//         int linkSubGraphLinkNumber = linkSubGraph->GetEdges();
//
//         featureVector.push_back(linkSubGraphLinkNumber);
//
//         TSnap::GetSccs(linkSubGraph, CommunityVec);
//         int SCCsInLinkSubGraph = CommunityVec.Len();
//         TSnap::GetWccs(linkSubGraph, CommunityVec);
//         int WCCsInLinkSubGraph = CommunityVec.Len();
//
//         featureVector.push_back(SCCsInLinkSubGraph);
//         featureVector.push_back(WCCsInLinkSubGraph);
//
//         uniqueNeighborsVec.AddUnique(fromNodeId);
//         uniqueNeighborsVec.AddUnique(toNodeId);
//         vector<int> v;
//         for (int z = 0; z < uniqueNeighborsVec.Len(); z++)
//         {
//            v.push_back(uniqueNeighborsVec[z]);
//         }
//         linkSubGraph_withSelf = TSnap::GetSubGraph(m_graph, uniqueNeighborsVec);
//         int linkSubGraphLinkNumber_withSelf = linkSubGraph_withSelf->GetEdges();
//         featureVector.push_back(linkSubGraphLinkNumber_withSelf);
//
//         TSnap::GetSccs(linkSubGraph_withSelf, CommunityVec);
//         int SCCsInLinkSubGraph_withSelf = CommunityVec.Len();
//         featureVector.push_back(SCCsInLinkSubGraph_withSelf);
//
//         TIntPrV* innerGraphEdges = new TIntPrV();
//         getInnerGraphEdges(fromNodeId, toNodeId, *innerGraphEdges);
//         innerSubGraph = TSnap::GetESubGraph(m_graph, *innerGraphEdges);
//         int innerSubGraphLinkNumber = innerSubGraph->GetEdges();
//         featureVector.push_back(innerSubGraphLinkNumber);
//
//         TSnap::GetSccs(innerSubGraph, CommunityVec);
//         int SCCsInInnerSubGraph = CommunityVec.Len();
//         featureVector.push_back(SCCsInInnerSubGraph);
//
//         TSnap::GetWccs(innerSubGraph, CommunityVec);
//         int WCCsInInnerSubGraph = CommunityVec.Len();
//         featureVector.push_back(WCCsInInnerSubGraph);
//
//         //shortest path Length
//         int shortestPathLen = TSnap::GetShortPath(m_graph, fromNodeId, toNodeId, m_is_directed_b);
//         featureVector.push_back(shortestPathLen);
//
//         //Adamic/Adar measure
//         fromNodeNeighVec.Intrs(toNodeNeighVec, neighborsIntrsVec);
//         double aaMeasure = 0.0;
//         int deg = 0;
//         double temp = 0.0;
//         for (int i = 0; i < neighborsIntrsVec.Len(); i++)
//         {
//            TTimeNet::TNodeI nodeIter = Graph->GetNI(neighborsIntrsVec[i]);
//            deg = nodeIter.GetDeg();
//            temp = 1/double(log(deg));
//            aaMeasure += temp;
//         }
//         featureVector.push_back(aaMeasure*100);
//
//         //TODO: Have betweeness measure as a feature
//
//         bool is_pos_example = Graph->IsEdge(fromNodeId, toNodeId);
//         featureVector.push_back(is_pos_example);
//         if (!is_pos_example)
//         {
//            numNegEx++;
//         }
//         totalExamples++;
//         m_target_concept_vec.push_back(is_pos_example);
//         m_non_temporal_fv.push_back(featureVector);
//      }
//   }
//}

void GraphU::getInnerGraphEdges(int fromNodeId, int toNodeId, TIntPrV& edgeIdVec)
{
   vector<int> fromNodeNeigh;
   vector<int> toNodeNeigh;

   TTimeNet::TNodeI fromNodeIter = m_graph->GetNI(fromNodeId);
   TTimeNet::TNodeI toNodeIter = m_graph->GetNI(toNodeId);
   for (int edge = 0; edge < fromNodeIter.GetOutDeg(); edge++)
   {
      const int OutNId = fromNodeIter.GetOutNId(edge);
      fromNodeNeigh.push_back(OutNId);
   }

   for (int edge = 0; edge < toNodeIter.GetOutDeg(); edge++)
   {
      const int OutNId = toNodeIter.GetOutNId(edge);
      toNodeNeigh.push_back(OutNId);
   }
   int size = 0;
   edgeIdVec.Gen(0);
   for (int n1 = 0; n1 < fromNodeNeigh.size(); n1++)
   {
      for (int n2 = 0; n2 < toNodeNeigh.size(); n2++)
      {
         if (m_graph->IsEdge(fromNodeNeigh[n1], toNodeNeigh[n2]))
         {
            const TIntPr edge = TIntPr(fromNodeNeigh[n1], toNodeNeigh[n2]);
            edgeIdVec.Add(edge);
         }
      }
   }
}

int GraphU::numEdgesAmongNeighbors(PNGraph graph, int fromNodeId, int toNodeId)
{
   TVec<TInt> fromNodeNeigh;
   TVec<TInt> toNodeNeigh;

   TNGraph::TNodeI fromNodeIter = graph->GetNI(fromNodeId);
   TNGraph::TNodeI toNodeIter = graph->GetNI(toNodeId);
   int e = 0;
   for (int edge = 0; edge < fromNodeIter.GetOutDeg(); edge++)
   {
      const int OutNId = fromNodeIter.GetOutNId(edge);
      fromNodeNeigh.Add(OutNId);
   }

   for (int edge = 0; edge < toNodeIter.GetOutDeg(); edge++)
   {
      const int OutNId = toNodeIter.GetOutNId(edge);
      toNodeNeigh.Add(OutNId);
   }

   //todo- see how to make use of node iters
   for (int i = 0; i < fromNodeNeigh.Len(); i++)
   {
      for (int j = 0; j < toNodeNeigh.Len(); j++)
      {
         bool isEdge = graph->IsEdge(fromNodeNeigh[i], toNodeNeigh[j]);
         if (isEdge)
         {
            e++;
         }
      }
   }
   return e / 2; //TODO- check if e/2 or e?
}

int GraphU::numEdgesAmongActiveNeighbors(PNGraph graph, int currentSnapshot, int fromNodeId, int toNodeId)
{
   TVec<TInt> activeFromNodeNeigh;
   TVec<TInt> activeToNodeNeigh;

   TNGraph::TNodeI fromNodeIter = graph->GetNI(fromNodeId);
   TNGraph::TNodeI toNodeIter = graph->GetNI(toNodeId);
   int e = 0;
   for (int edge = 0; edge < fromNodeIter.GetOutDeg(); edge++)
   {
      const int OutNId = fromNodeIter.GetOutNId(edge);
      if (isNodeActive(currentSnapshot, OutNId))
      {
         activeFromNodeNeigh.Add(OutNId);
      }
   }

   for (int edge = 0; edge < toNodeIter.GetOutDeg(); edge++)
   {
      const int OutNId = toNodeIter.GetOutNId(edge);
      if (isNodeActive(currentSnapshot, OutNId))
      {
         activeToNodeNeigh.Add(OutNId);
      }
   }

   //todo- see how to make use of node iters
   for (int i = 0; i < activeFromNodeNeigh.Len(); i++)
   {
      for (int j = 0; j < activeToNodeNeigh.Len(); j++)
      {
         bool isEdge = graph->IsEdge(activeFromNodeNeigh[i], activeToNodeNeigh[j]);
         if (isEdge)
         {
            e++;
         }
      }
   }
   return e / 2; //TODO- check if e/2 or e?
}

//int GraphU::numEdgesAmongNeighbors(int fromNodeId, int toNodeId)
//{
//int fromNode = getNodeIndexGivenNodeId(fromNodeId);
//int toNode = getNodeIndexGivenNodeId(toNodeId);

//pair<EdgeMap::iterator, EdgeMap::iterator> fromBounds = adjacencyList.equal_range(fromNode);
//EdgeMap::local_iterator from_begin_iter = fromBounds.first;
//EdgeMap::local_iterator from_end_iter = fromBounds.second;

//pair<EdgeMap::iterator, EdgeMap::iterator> toBounds = adjacencyList.equal_range(fromNode);
//EdgeMap::local_iterator to_begin_iter = toBounds.first;
//EdgeMap::local_iterator to_end_iter = toBounds.second;

//vector<int> toNodeNeighbors;
//vector<int> fromNodeNeighbors;
//assocsType::iterator edgeAssocIterValue;

//int e = 0; //number of edges among neighbors of u and v
//while (from_begin_iter != from_end_iter)
//{
//   int n = from_begin_iter->second;
//   fromNodeNeighbors.push_back(n);
//   ++from_begin_iter;
//}
//while (to_begin_iter != to_end_iter)
//{
//   int n = to_begin_iter->second;
//   toNodeNeighbors.push_back(n);
//   ++to_begin_iter;
//}

//for (unsigned int i = 0; i < fromNodeNeighbors.size(); i++)
//{
//   int n1 = fromNodeNeighbors[i];
//   for (unsigned int j = 0; j < toNodeNeighbors.size(); j++)
//   {
//      int n2 = toNodeNeighbors[j];
//      if (edgeAssocs.count(make_pair(n1, n2)))
//      {
//         e++;
//      }
//   }
//}
//return e / 2; //return # of links
//}

int main(int argc, char** argv)
{
   int dummy;
   char* fileIp1 = "D:/ACADS/CS695_SNA/Project/Datasets/facebook-links.txt";
   char* fileIp2 = "D:/ACADS/CS695_SNA/Project/Datasets/facebook-wall.txt";
   //char* fileIp = "D:/ACADS/CS695_SNA/Project/sample.txt";
   //char* fileIp1 = "D:/ACADS/CS695_SNA/Project/Datasets/Cit-HepTh-dates.txt";
   //char* fileIp2 = "D:/ACADS/CS695_SNA/Project/Datasets/Cit-HepTh.txt";
   GraphU* g = new GraphU(5, "", "");
   g->readGraph(fileIp1, fileIp2);
   g->extractFeatureSet();
   //g->extractTemporalFeatureSet();
   //g->analyzeTemporalCharacteristics();
   //g->generateTemoralGraphModels();

   const vector<int>* nodesTimeSlot = g->getNodesInTimeSlot();
   const vector<int>* edgesTimeSlot = g->getEdgesInTimeSlot();
   const vector<int>* linkTimeBucketKeyV = g->getLinkTimeBucketKeyVec();
   const vector<double>* effDiamTimeSlot = g->getEffDiamInTimeSlot();
   const vector<double>* avgDegTimeSlot = g->getdegreeInTimeSlot();
   const vector<double>* maxSccNodesTimeSlot = g->getMaxSccNodesInTimeSlot();
   const TVec<PNGraph>* snapShotVec = g->getSnapShotVec();

   const vector<double>* FFMEffDiam1 = g->getFFMEffecDiamV1();
   const vector<double>* FFMEffDiam2 = g->getFFMEffecDiamV2();

   const vector<pair<int, int>>* FFMDensV1 = g->getFFMDensificationV1();
   const vector<pair<int, int>>* FFMDensV2 = g->getFFMDensificationV2();

   Excel::_ApplicationPtr XL_ptr;
   CoInitialize(NULL);
   XL_ptr.CreateInstance(L"Excel.Application");
   XL_ptr->Visible = true;

   //std::wstring title1 = L"From Node Id";
   //BSTR bs_title1 = SysAllocStringLen(title1.data(), (UINT)title1.size());
   //range_of_cells_ptr->Item[1][1] = bs_title1;

   //std::wstring title2 = L"To Node Id";
   //bs_title1 = SysAllocStringLen(title1.data(), (UINT)title1.size());
   //range_of_cells_ptr->Item[1][2] = bs_title1;

   //std::wstring title0 = L"Num common Neighbors";
   //BSTR bs_title0 = SysAllocStringLen(title0.data(), (UINT)title0.size());
   //range_of_cells_ptr->Item[1][3] = bs_title0;

   //title1 = L"FromNode Degree";
   //bs_title1 = SysAllocStringLen(title1.data(), (UINT)title1.size());
   //range_of_cells_ptr->Item[1][4] = bs_title1;

   //title2 = L"toNode Degree";
   //BSTR bs_title2 = SysAllocStringLen(title2.data(), (UINT)title2.size());
   //range_of_cells_ptr->Item[1][5] = bs_title2;

   //std::wstring title3 = L"Jackard*100";
   //BSTR bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
   //range_of_cells_ptr->Item[1][6] = bs_title3;

   //std::wstring title4 = L"Pref attachment score";
   //BSTR bs_title4 = SysAllocStringLen(title4.data(), (UINT)title4.size());
   //range_of_cells_ptr->Item[1][7] = bs_title4;

   //std::wstring title5 = L"Katz";
   //BSTR bs_title5 = SysAllocStringLen(title5.data(), (UINT)title5.size());
   //range_of_cells_ptr->Item[1][8] = bs_title5;

   //std::wstring title6 = L"fromNodeInducedSubgraphLinkNumber";
   //BSTR bs_title6 = SysAllocStringLen(title6.data(), (UINT)title6.size());
   //range_of_cells_ptr->Item[1][9] = bs_title6;

   //std::wstring title7 = L"fromNodeInducedSubgraphLinkNumber_withSelf";
   //BSTR bs_title7 = SysAllocStringLen(title7.data(), (UINT)title7.size());
   //range_of_cells_ptr->Item[1][10] = bs_title7;

   //std::wstring title8 = L"toNodeInducedSubGraphLinkNumber";
   //BSTR bs_title8 = SysAllocStringLen(title8.data(), (UINT)title8.size());
   //range_of_cells_ptr->Item[1][11] = bs_title8;

   //std::wstring title9 = L"toNodeInducedSubgraphLinkNumber_withSelf";
   //BSTR bs_title9 = SysAllocStringLen(title9.data(), (UINT)title9.size());
   //range_of_cells_ptr->Item[1][12] = bs_title9;

   //std::wstring title10 = L"linkSubGraphLinkNumber";
   //BSTR bs_title10 = SysAllocStringLen(title10.data(), (UINT)title10.size());
   //range_of_cells_ptr->Item[1][13] = bs_title10;

   //std::wstring title11 = L"SCCsInLinkSubGraph";
   //BSTR bs_title11 = SysAllocStringLen(title11.data(), (UINT)title11.size());
   //range_of_cells_ptr->Item[1][14] = bs_title11;

   //std::wstring title12 = L"WCCsInLinkSubGraph";
   //BSTR bs_title12 = SysAllocStringLen(title12.data(), (UINT)title12.size());
   //range_of_cells_ptr->Item[1][15] = bs_title12;

   //std::wstring title13 = L"linkSubGraphLinkNumber_withSelf";
   //BSTR bs_title13 = SysAllocStringLen(title13.data(), (UINT)title13.size());
   //range_of_cells_ptr->Item[1][16] = bs_title13;

   //std::wstring title14 = L"SCCsInLinkSubGraph_withSelf";
   //BSTR bs_title14 = SysAllocStringLen(title14.data(), (UINT)title14.size());
   //range_of_cells_ptr->Item[1][17] = bs_title14;

   //std::wstring title15 = L"innerSubGraphLinkNumber";
   //BSTR bs_title15 = SysAllocStringLen(title15.data(), (UINT)title15.size());
   //range_of_cells_ptr->Item[1][18] = bs_title15;

   //std::wstring title16 = L"SCCsInInnerSubGraph";
   //BSTR bs_title16 = SysAllocStringLen(title16.data(), (UINT)title16.size());
   //range_of_cells_ptr->Item[1][19] = bs_title16;

   //std::wstring title17 = L"WCCsInInnerSubGraph";
   //BSTR bs_title17 = SysAllocStringLen(title17.data(), (UINT)title17.size());
   //range_of_cells_ptr->Item[1][20] = bs_title17;

   //std::wstring title18 = L"shortest Path Len";
   //BSTR bs_title18 = SysAllocStringLen(title18.data(), (UINT)title18.size());
   //range_of_cells_ptr->Item[1][21] = bs_title18;

   //std::wstring title19 = L"Adamic-Adar*100";
   //BSTR bs_title19 = SysAllocStringLen(title19.data(), (UINT)title19.size());
   //range_of_cells_ptr->Item[1][22] = bs_title19;

   //std::wstring title45 = L"Link?";
   //BSTR bs_title45 = SysAllocStringLen(title45.data(), (UINT)title45.size());
   //range_of_cells_ptr->Item[1][23] = bs_title45;

   vector<vector<double>>* ntfv = g->getNonTemporalFeatureSet();
   vector<vector<double>>* tfv = g->getTemporalFeatureSet();
   vector<double> NTfeaturesV;
   vector<double> TfeaturesV;
   const vector<bool>* targetConceptV = g->getTargetConceptVec();
   int j = 0;
   int k = 0;
   std::wstring title0;
   BSTR bs_title0;

   BSTR bs_title1;
   std::wstring title1;
   BSTR bs_title2;
   std::wstring title2;
   BSTR bs_title3;
   std::wstring title3;
   BSTR bs_title4;
   std::wstring title4;
   BSTR bs_title5;
   std::wstring title5;
   Excel::_WorkbookPtr workbook_ptr1 = XL_ptr->Workbooks->Add(Excel::xlWorksheet);
   Excel::_WorksheetPtr worksheet_ptr1 = XL_ptr->ActiveWorkbook->Worksheets->Add();
   Excel::RangePtr range_of_cells_ptr = worksheet_ptr1->Cells;
   int previousSnapShot = 0;
   int z = 2;
   for (int i = 0; i < ntfv->size(); i++)
   {
      NTfeaturesV = ntfv->at(i);
      TfeaturesV = tfv->at(i);
      int currentSnapShot = NTfeaturesV.at(0);
      if (i > 0)
      {
         previousSnapShot = ntfv->at(i - 1).at(0);
      }
      if (previousSnapShot != currentSnapShot)
      {
         z = 2;
         worksheet_ptr1 = XL_ptr->ActiveWorkbook->Worksheets->Add();
         range_of_cells_ptr = worksheet_ptr1->Cells;
         title0 = L"SnapShot";
         bs_title0 = SysAllocStringLen(title0.data(), (UINT)title0.size());
         range_of_cells_ptr->Item[1][1] = bs_title0;

         title1 = L"From Node Id";
         bs_title1 = SysAllocStringLen(title1.data(), (UINT)title1.size());
         range_of_cells_ptr->Item[1][2] = bs_title1;

         title2 = L"To Node Id";
         bs_title2 = SysAllocStringLen(title2.data(), (UINT)title2.size());
         range_of_cells_ptr->Item[1][3] = bs_title2;

         title1 = L"InDegreeF";
         bs_title1 = SysAllocStringLen(title1.data(), (UINT)title1.size());
         range_of_cells_ptr->Item[1][4] = bs_title1;

         title2 = L"OutDegreeF";
         bs_title2 = SysAllocStringLen(title2.data(), (UINT)title2.size());
         range_of_cells_ptr->Item[1][5] = bs_title2;

         title1 = L"InDegreeTo";
         bs_title1 = SysAllocStringLen(title1.data(), (UINT)title1.size());
         range_of_cells_ptr->Item[1][6] = bs_title1;

         title2 = L"OutDegreeTo";
         bs_title2 = SysAllocStringLen(title2.data(), (UINT)title2.size());
         range_of_cells_ptr->Item[1][7] = bs_title2;

         title2 = L"RecencyFrom";
         bs_title2 = SysAllocStringLen(title2.data(), (UINT)title2.size());
         range_of_cells_ptr->Item[1][8] = bs_title2;

         title2 = L"RecencyTo";
         bs_title2 = SysAllocStringLen(title2.data(), (UINT)title2.size());
         range_of_cells_ptr->Item[1][9] = bs_title2;

         title0 = L"CN";
         bs_title0 = SysAllocStringLen(title0.data(), (UINT)title0.size());
         range_of_cells_ptr->Item[1][10] = bs_title0;

         title0 = L"Size2From";
         bs_title0 = SysAllocStringLen(title0.data(), (UINT)title0.size());
         range_of_cells_ptr->Item[1][11] = bs_title0;

         title0 = L"Size2To";
         bs_title0 = SysAllocStringLen(title0.data(), (UINT)title0.size());
         range_of_cells_ptr->Item[1][12] = bs_title0;

         title5 = L"Katz";
         bs_title5 = SysAllocStringLen(title5.data(), (UINT)title5.size());
         range_of_cells_ptr->Item[1][13] = bs_title5;

         title4 = L"PA";
         bs_title4 = SysAllocStringLen(title4.data(), (UINT)title4.size());
         range_of_cells_ptr->Item[1][14] = bs_title4;

         title4 = L"AA";
         bs_title4 = SysAllocStringLen(title4.data(), (UINT)title4.size());
         range_of_cells_ptr->Item[1][15] = bs_title4;

         //title4 = L"Distance";
         //bs_title4 = SysAllocStringLen(title4.data(), (UINT)title4.size());
         //range_of_cells_ptr->Item[1][16] = bs_title4;

         title3 = L"JC";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][16] = bs_title3;

         title3 = L"FNTriads";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][17] = bs_title3;

         title3 = L"TNTriads";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][18] = bs_title3;

         title3 = L"actFromDeg";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][64] = bs_title3;

         title3 = L"actToDeg";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][65] = bs_title3;

         title3 = L"actAA";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][66] = bs_title3;

         title3 = L"actApproxKatz";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][67] = bs_title3;

         title3 = L"actCN";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][68] = bs_title3;

         title3 = L"actJC";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][69] = bs_title3;

         /*=========================================================================*/

         title1 = L"InDegreeFA2";
         bs_title1 = SysAllocStringLen(title1.data(), (UINT)title1.size());
         range_of_cells_ptr->Item[1][19] = bs_title1;

         title2 = L"OutDegreeFA2";
         bs_title2 = SysAllocStringLen(title2.data(), (UINT)title2.size());
         range_of_cells_ptr->Item[1][20] = bs_title2;

         title1 = L"InDegreeToA2";
         bs_title1 = SysAllocStringLen(title1.data(), (UINT)title1.size());
         range_of_cells_ptr->Item[1][21] = bs_title1;

         title2 = L"OutDegreeToA2";
         bs_title2 = SysAllocStringLen(title2.data(), (UINT)title2.size());
         range_of_cells_ptr->Item[1][22] = bs_title2;

         title2 = L"RecencyFromA2";
         bs_title2 = SysAllocStringLen(title2.data(), (UINT)title2.size());
         range_of_cells_ptr->Item[1][23] = bs_title2;

         title2 = L"RecencyToA2";
         bs_title2 = SysAllocStringLen(title2.data(), (UINT)title2.size());
         range_of_cells_ptr->Item[1][24] = bs_title2;

         title0 = L"CNA2";
         bs_title0 = SysAllocStringLen(title0.data(), (UINT)title0.size());
         range_of_cells_ptr->Item[1][25] = bs_title0;

         title0 = L"Size2FromA2";
         bs_title0 = SysAllocStringLen(title0.data(), (UINT)title0.size());
         range_of_cells_ptr->Item[1][26] = bs_title0;

         title0 = L"Size2ToA2";
         bs_title0 = SysAllocStringLen(title0.data(), (UINT)title0.size());
         range_of_cells_ptr->Item[1][27] = bs_title0;

         title5 = L"KatzA2";
         bs_title5 = SysAllocStringLen(title5.data(), (UINT)title5.size());
         range_of_cells_ptr->Item[1][28] = bs_title5;

         title4 = L"PAA2";
         bs_title4 = SysAllocStringLen(title4.data(), (UINT)title4.size());
         range_of_cells_ptr->Item[1][29] = bs_title4;

         title4 = L"AAA2";
         bs_title4 = SysAllocStringLen(title4.data(), (UINT)title4.size());
         range_of_cells_ptr->Item[1][30] = bs_title4;

         //title4 = L"DistanceA2";
         //bs_title4 = SysAllocStringLen(title4.data(), (UINT)title4.size());
         //range_of_cells_ptr->Item[1][32] = bs_title4;

         title3 = L"JCA2";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][31] = bs_title3;

         title3 = L"FNTriadsA2";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][32] = bs_title3;

         title3 = L"TNTriadsA2";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][33] = bs_title3;

         title3 = L"actFromDeg";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][34] = bs_title3;

         title3 = L"actToDeg";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][35] = bs_title3;

         title3 = L"actAA";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][36] = bs_title3;

         title3 = L"actApproxKatz";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][37] = bs_title3;

         title3 = L"actCN";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][38] = bs_title3;

         title3 = L"actJC";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][39] = bs_title3;

         /*===========================================*/
         title1 = L"InDegreeFA5";
         bs_title1 = SysAllocStringLen(title1.data(), (UINT)title1.size());
         range_of_cells_ptr->Item[1][40] = bs_title1;

         title2 = L"OutDegreeFA5";
         bs_title2 = SysAllocStringLen(title2.data(), (UINT)title2.size());
         range_of_cells_ptr->Item[1][41] = bs_title2;

         title1 = L"InDegreeToA5";
         bs_title1 = SysAllocStringLen(title1.data(), (UINT)title1.size());
         range_of_cells_ptr->Item[1][42] = bs_title1;

         title2 = L"OutDegreeToA5";
         bs_title2 = SysAllocStringLen(title2.data(), (UINT)title2.size());
         range_of_cells_ptr->Item[1][43] = bs_title2;

         title2 = L"RecencyFromA5";
         bs_title2 = SysAllocStringLen(title2.data(), (UINT)title2.size());
         range_of_cells_ptr->Item[1][44] = bs_title2;

         title2 = L"RecencyToA5";
         bs_title2 = SysAllocStringLen(title2.data(), (UINT)title2.size());
         range_of_cells_ptr->Item[1][45] = bs_title2;

         title0 = L"CNA5";
         bs_title0 = SysAllocStringLen(title0.data(), (UINT)title0.size());
         range_of_cells_ptr->Item[1][46] = bs_title0;

         title0 = L"Size2FromA5";
         bs_title0 = SysAllocStringLen(title0.data(), (UINT)title0.size());
         range_of_cells_ptr->Item[1][47] = bs_title0;

         title0 = L"Size2ToA5";
         bs_title0 = SysAllocStringLen(title0.data(), (UINT)title0.size());
         range_of_cells_ptr->Item[1][48] = bs_title0;

         title5 = L"KatzA5";
         bs_title5 = SysAllocStringLen(title5.data(), (UINT)title5.size());
         range_of_cells_ptr->Item[1][49] = bs_title5;

         title4 = L"PAA5";
         bs_title4 = SysAllocStringLen(title4.data(), (UINT)title4.size());
         range_of_cells_ptr->Item[1][50] = bs_title4;

         title4 = L"AAA5";
         bs_title4 = SysAllocStringLen(title4.data(), (UINT)title4.size());
         range_of_cells_ptr->Item[1][51] = bs_title4;

         //title4 = L"DistanceA5";
         //bs_title4 = SysAllocStringLen(title4.data(), (UINT)title4.size());
         //range_of_cells_ptr->Item[1][48] = bs_title4;

         title3 = L"JCA5";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][52] = bs_title3;

         title3 = L"FNTriadsA5";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][53] = bs_title3;

         title3 = L"TNTriadsA5";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][54] = bs_title3;

         title3 = L"actFromDeg5";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][55] = bs_title3;

         title3 = L"actToDeg5";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][56] = bs_title3;

         title3 = L"actAA5";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][57] = bs_title3;

         title3 = L"actApproxKatz5";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][58] = bs_title3;

         title3 = L"actCN5";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][59] = bs_title3;

         title3 = L"actJC5";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][60] = bs_title3;

         /*=====================================================*/
         title1 = L"InDegreeFA10";
         bs_title1 = SysAllocStringLen(title1.data(), (UINT)title1.size());
         range_of_cells_ptr->Item[1][61] = bs_title1;

         title2 = L"OutDegreeFA10";
         bs_title2 = SysAllocStringLen(title2.data(), (UINT)title2.size());
         range_of_cells_ptr->Item[1][62] = bs_title2;

         title1 = L"InDegreeToA10";
         bs_title1 = SysAllocStringLen(title1.data(), (UINT)title1.size());
         range_of_cells_ptr->Item[1][63] = bs_title1;

         title2 = L"OutDegreeToA10";
         bs_title2 = SysAllocStringLen(title2.data(), (UINT)title2.size());
         range_of_cells_ptr->Item[1][64] = bs_title2;

         title2 = L"RecencyFromA10";
         bs_title2 = SysAllocStringLen(title2.data(), (UINT)title2.size());
         range_of_cells_ptr->Item[1][65] = bs_title2;

         title2 = L"RecencyToA10";
         bs_title2 = SysAllocStringLen(title2.data(), (UINT)title2.size());
         range_of_cells_ptr->Item[1][66] = bs_title2;

         title0 = L"CNA10";
         bs_title0 = SysAllocStringLen(title0.data(), (UINT)title0.size());
         range_of_cells_ptr->Item[1][67] = bs_title0;

         title0 = L"Size2FromA10";
         bs_title0 = SysAllocStringLen(title0.data(), (UINT)title0.size());
         range_of_cells_ptr->Item[1][68] = bs_title0;

         title0 = L"Size2ToA10";
         bs_title0 = SysAllocStringLen(title0.data(), (UINT)title0.size());
         range_of_cells_ptr->Item[1][69] = bs_title0;

         title5 = L"KatzA10";
         bs_title5 = SysAllocStringLen(title5.data(), (UINT)title5.size());
         range_of_cells_ptr->Item[1][70] = bs_title5;

         title4 = L"PAA10";
         bs_title4 = SysAllocStringLen(title4.data(), (UINT)title4.size());
         range_of_cells_ptr->Item[1][71] = bs_title4;

         title4 = L"AAA10";
         bs_title4 = SysAllocStringLen(title4.data(), (UINT)title4.size());
         range_of_cells_ptr->Item[1][72] = bs_title4;

         //title4 = L"DistanceA10";
         //bs_title4 = SysAllocStringLen(title4.data(), (UINT)title4.size());
         //range_of_cells_ptr->Item[1][64] = bs_title4;

         title3 = L"JCA10";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][73] = bs_title3;

         title3 = L"FNTriadsA10";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][74] = bs_title3;

         title3 = L"TNTriadsA10";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][75] = bs_title3;

         title3 = L"actFromDeg10";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][76] = bs_title3;

         title3 = L"actToDeg10";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][77] = bs_title3;

         title3 = L"actAA10";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][78] = bs_title3;

         title3 = L"actApproxKatz10";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][79] = bs_title3;

         title3 = L"actCN10";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][80] = bs_title3;

         title3 = L"actJC10";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][81] = bs_title3;

         title3 = L"LD";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][82] = bs_title3;

         /*=====================================================*/
         title3 = L"IsPos?";
         bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
         range_of_cells_ptr->Item[1][83] = bs_title3;
      }
      else
      {
         z++;
      }

      for (j = 0; j < NTfeaturesV.size(); j++)
      {
         range_of_cells_ptr->Item[z][j + 1] = NTfeaturesV.at(j);
      }
      for (k = 0; k < TfeaturesV.size(); k++)
      {
         range_of_cells_ptr->Item[z][k + j + 1] = TfeaturesV.at(k);
      }
      //range_of_cells_ptr->Item[i + 2][k + j+ 1] = targetConceptV->at(i);
      range_of_cells_ptr->Item[z][83] = targetConceptV->at(i); //hardcoded as of now. TODO- Remove this hack!
   }
   //***********************************************
   Excel::_WorksheetPtr worksheet_ptr2 = XL_ptr->ActiveWorkbook->Worksheets->Add();
   Excel::RangePtr range_of_cells_ptr2 = worksheet_ptr2->Cells;

   title1 = L"Nodes";
   bs_title1 = SysAllocStringLen(title1.data(), (UINT)title1.size());
   range_of_cells_ptr2->Item[1][1] = bs_title1;

   title2 = L"Edges";
   bs_title2 = SysAllocStringLen(title2.data(), (UINT)title2.size());
   range_of_cells_ptr2->Item[1][2] = bs_title2;

   title3 = L"Eff Diam";
   bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
   range_of_cells_ptr2->Item[1][3] = bs_title3;

   title4 = L"Avg degree";
   bs_title4 = SysAllocStringLen(title4.data(), (UINT)title4.size());
   range_of_cells_ptr2->Item[1][4] = bs_title4;

   title5 = L"% Nodes in MaxSCC";
   bs_title5 = SysAllocStringLen(title5.data(), (UINT)title5.size());
   range_of_cells_ptr2->Item[1][5] = bs_title5;

   for (int i = 0; i < nodesTimeSlot->size(); i++)
   {
      range_of_cells_ptr2->Item[i + 2][1] = nodesTimeSlot->at(i);
      range_of_cells_ptr2->Item[i + 2][2] = edgesTimeSlot->at(i);
      range_of_cells_ptr2->Item[i + 2][3] = effDiamTimeSlot->at(i);
      range_of_cells_ptr2->Item[i + 2][4] = avgDegTimeSlot->at(i);
      range_of_cells_ptr2->Item[i + 2][5] = maxSccNodesTimeSlot->at(i);
   }

   //*****************Forest Fire Model Graphs****************************
   Excel::_WorksheetPtr worksheet_ptr3 = XL_ptr->ActiveWorkbook->Worksheets->Add();
   Excel::RangePtr range_of_cells_ptr3 = worksheet_ptr3->Cells;

   title1 = L"EffDiameter1";
   bs_title1 = SysAllocStringLen(title1.data(), (UINT)title1.size());
   range_of_cells_ptr3->Item[1][1] = bs_title1;

   title2 = L"Nodes1";
   bs_title2 = SysAllocStringLen(title2.data(), (UINT)title2.size());
   range_of_cells_ptr3->Item[1][2] = bs_title2;

   title3 = L"EffDiameter2";
   bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
   range_of_cells_ptr3->Item[1][4] = bs_title3;

   title4 = L"Nodes2";
   bs_title4 = SysAllocStringLen(title4.data(), (UINT)title4.size());
   range_of_cells_ptr3->Item[1][5] = bs_title4;

   for (int i = 0; i < FFMEffDiam1->size(); i++)
   {
      range_of_cells_ptr3->Item[i + 2][1] = FFMEffDiam1->at(i);
      range_of_cells_ptr3->Item[i + 2][2] = FFMDensV1->at(i).first;
      range_of_cells_ptr3->Item[i + 2][3] = FFMDensV1->at(i).second;

      range_of_cells_ptr3->Item[i + 2][4] = FFMEffDiam2->at(i);
      range_of_cells_ptr3->Item[i + 2][5] = FFMDensV2->at(i).first;
      range_of_cells_ptr3->Item[i + 2][6] = FFMDensV2->at(i).second;
   }

   //basic graph stats
   Excel::_WorksheetPtr worksheet_ptr4 = XL_ptr->ActiveWorkbook->Worksheets->Add();
   Excel::RangePtr range_of_cells_ptr4 = worksheet_ptr4->Cells;

   title1 = L"TimePeriod";
   bs_title1 = SysAllocStringLen(title1.data(), (UINT)title1.size());
   range_of_cells_ptr4->Item[1][1] = bs_title1;

   title2 = L"NumNodes";
   bs_title2 = SysAllocStringLen(title2.data(), (UINT)title2.size());
   range_of_cells_ptr4->Item[1][2] = bs_title2;

   title3 = L"NumEdges";
   bs_title3 = SysAllocStringLen(title3.data(), (UINT)title3.size());
   range_of_cells_ptr4->Item[1][3] = bs_title3;

   for (int i = 0; i < snapShotVec->Len(); i++)
   {
      PNGraph currentSnapShot = snapShotVec->GetVal(i);
      range_of_cells_ptr4->Item[i + 2][1] = linkTimeBucketKeyV->at(i);
      range_of_cells_ptr4->Item[i + 2][2] = currentSnapShot->GetNodes();
      range_of_cells_ptr4->Item[i + 2][3] = currentSnapShot->GetEdges();
   }
   // XL_ptr->Quit();
   // CoUninitialize();
   cout << "\n Press any key to exit!!\n";
   cin >> dummy;
   return 0;
}



