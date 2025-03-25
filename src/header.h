#ifndef _HEADER_H_
#define _HEADER_H_
#include <bits/stdc++.h>
using namespace std;

#define BITSET_WIDTH 4096
//For random generation ---------------------------------------------
extern mt19937_64 rand_generator;
inline uint64_t rand_ull(uint64_t l,uint64_t r){
    std::uniform_int_distribution<uint64_t> distribution{l,r};
    return distribution(rand_generator);
}
inline void shuffle_vector(vector<int> &vec, int from=0){
    shuffle(vec.begin()+from,vec.end(),rand_generator);
}
inline bitset<BITSET_WIDTH> rand_bits() {
    std::uniform_int_distribution<uint64_t> distribution{0,std::numeric_limits<uint64_t>::max()};
    bitset<BITSET_WIDTH> res;
    for(size_t i=0;i<BITSET_WIDTH;i+=64){
        res<<=64;
        res|=distribution(rand_generator);
    }
    return res;
}
//-----------------------------------------------------------------

//Global
class Global{
public:
    string output_filename;
    vector<vector<string>> PO_vectors;
    vector<double> PO_weights;
    string ed_mode="";
    chrono::time_point<chrono::system_clock> start_time;
    int time_limit;
    int runtime;
    void setup_PO_vec_and_PO_weights(string ed_mode);
    bool check_if_time_out();
    void show_time();
    void setTimeLimit(int);
    void setOutputFilename(string filename);
    Global();
}; 



// Circuit Parser
enum NodeType{Dot, DotXI, DotYI, DotZI, DotXYI, DotYZI, DotXZI, DotXYZI
            , DotOI, DotXOI, DotYOI, DotZOI, DotXYOI, DotYZOI, DotXZOI, DotXYZOI
            , PI, PO_BUF, PO_NOT, CONST1, CONST0};

inline unordered_map <NodeType, string> typeName = {
    {Dot,"Dot"},{DotXI,"DotXI"},{DotYI,"DotYI"},{DotZI,"DotZI"},
    {DotXYI,"DotXYI"},{DotYZI,"DotYZI"},{DotXZI,"DotXZI"},{DotXYZI,"DotXYZI"},
    {DotOI,"DotOI"},{DotXOI,"DotXOI"},{DotYOI,"DotYOI"},{DotZOI,"DotZOI"},
    {DotXYOI,"DotXYOI"},{DotYZOI,"DotYZOI"},{DotXZOI,"DotXZOI"},{DotXYZOI,"DotXYZOI"},
};

inline bitset<4> getInvMapByNodeType(NodeType ntype){
    //bit0: x, bit1: y, bit2: z, bit3: o
    //0: no inversion, 1: inversion
    switch(ntype){
        case Dot: return bitset<4>(0);
        case DotXI: return bitset<4>(1);
        case DotYI: return bitset<4>(2);
        case DotZI: return bitset<4>(4);    
        case DotXYI: return bitset<4>(3);
        case DotYZI: return bitset<4>(6);
        case DotXZI: return bitset<4>(5);
        case DotXYZI: return bitset<4>(7);
        case DotOI: return bitset<4>(8);
        case DotXOI: return bitset<4>(9);
        case DotYOI: return bitset<4>(10);
        case DotZOI: return bitset<4>(12);
        case DotXYOI: return bitset<4>(11);
        case DotYZOI: return bitset<4>(14);
        case DotXZOI: return bitset<4>(13);
        case DotXYZOI: return bitset<4>(15);
        default:
            cerr << "error in get_inv_mapping" << endl;
            exit(0);
    }
}
inline NodeType getNodeTypeByInvMap(bitset<4> invMap){
    switch(invMap.to_ulong()){
        case 0: return Dot;
        case 1: return DotXI;
        case 2: return DotYI;
        case 4: return DotZI;
        case 3: return DotXYI;
        case 6: return DotYZI;
        case 5: return DotXZI;
        case 7: return DotXYZI;
        case 8: return DotOI;
        case 9: return DotXOI;
        case 10: return DotYOI;
        case 12: return DotZOI;
        case 11: return DotXYOI;
        case 14: return DotYZOI;
        case 13: return DotXZOI;
        case 15: return DotXYZOI;
        default:
            cerr << "error in getNodeType" << endl;
            exit(0);
    }
}
class Node{
public:
    string name;
    NodeType type;
    vector<Node*> fi;
    vector<Node*> fo;
    Node(string _name, NodeType _type):name(_name),type(_type){};
};
class Circuit{
public:
    string modelName;
    vector<Node*> allNodes;
    vector<Node*> PIs;
    vector<Node*> POs;    
    Node* constant0;
    Node* constant1;
    unordered_map<string, Node*> nodeName_2_node;
    void read_dig(string filename);
    void write_dig(string filename);
    vector<Node*> getNodesInTopoOrder();    //[constant0, constant1, PI1, PI2, ...]
    vector<Node*> getNodesInTopoOrder_POsAtEnd();
    void replace_nt_with_ns(Node* nt, Node* ns);
    void dfsFromPOstoCleanNetwork(); //All PIs will be kept.
    int getSize();
    int getLevel();
};

// GA engine
class Chromosome{
public:
   vector<vector<int>> genes; //{type, id, x, y, z}
   vector<int> activeGenes; //index of active genes
   double error_rate;
   double error_rate2;
   int size=INT_MAX;
   int depth;
   double area_rate;
   double diff_rate;
   double depth_rate;
   double fitness_val;
   double fitness_val2;
   bool data_avail=false;
   void print_info();
   Chromosome gen_new(){
        Chromosome newCand;
        newCand=*this;
        newCand.data_avail=false;
        return newCand;
   }
};
class GA_Engine{
    Circuit* golden_circuit;    //for checking error rate
    vector<Node*> golden_topo;  //for checking error rate
    int golden_size;
    int orig_size;
    int orig_depth;
    string output_filename;

    Circuit* orig_circuit;
    vector<Node*> orig_topo;
    unordered_map<Node*, int> getID;        //Node* -> id
    unordered_map<int, Node*> getNode;      //id -> Node*
    unordered_map<Node*,int> getGeneIndex;  //Node* -> gene index (position in chromo.genes)

    double ERROR_RATE_THRESHOLD;
    Chromosome orig;
    Chromosome bestChromo;
    bool terminate_signal=false;
    int terminate_counter=0;
    Global* global;

public:
    GA_Engine(Global& global,Circuit* circuit, Circuit* golden_circuit, string output_blif, double error_rate_threshold);
    double get_difference_rate(Chromosome& chromo1, Chromosome& chromo2);
    void simulate_chromosome_by_bitset(Chromosome& chromo, unordered_map<Node*,bitset<BITSET_WIDTH>> &vals);
    double check_error_rate_of(Chromosome& chromo, int batch_count=30);
    double check_error_dist_of(Chromosome& chromo, int batch_count=30);
    set<Node*> get_active_Nodes(Chromosome& chromo);
    int check_size_of(Chromosome& chromo);
    int check_depth_of(const Chromosome& chromo);
    void gen_blif(Chromosome& chromo, string filename);
    Chromosome encode_orig_circuit();
    vector<Chromosome> gen_initial_population(Chromosome& orig, int size);
    void mutate(vector<Chromosome>& cands,int ptr1, int ptr2, int count_for_each, int step_size);
    void crossover_onepoint(vector<Chromosome>& cands, int bound,int count_for_each);
    void crossover_universal(vector<Chromosome>& cands, int bound,int count_for_each);
    int learn_from_mutation_result(vector<Chromosome>& cands, int mutate_ptr1, int mutate_ptr2, int mutate_step_size);
    void learn_from_survivor(vector<Chromosome>& cands);
    void eval(vector<Chromosome> & cands);
    void select(vector<Chromosome>& cands, int trim_size1, int trim_size2);
    int start(int total_generation, int round_indicator);
    int get_merged(Chromosome& chromo, int nt);
    unordered_map<int,int> getMA(bool,Chromosome&,int,int,unordered_map<int,vector<int>>&, unordered_map<int, bitset<65536>>&,unordered_map<int, bitset<262144>>&);
    int check_subs(int, Chromosome& chromo, unordered_map<int,int>& ma0, unordered_map<int,int>& ma1, unordered_map<int,bitset<65536>>&);
    void trim_more(Chromosome &chromo);
};

// Simulator - Circuit
const unordered_map<Node*,int>& constructMFFC(Circuit* orig);
void simulate_by_bitset(vector<Node*> &nodesInTopoOrder, unordered_map<Node*,bitset<BITSET_WIDTH>> &vals);
double check_error_rate(Circuit* orig, Circuit* modified, int batch_count=30);
double check_error_distance(Circuit* orig, Circuit* modified, Global& global, int batch_count=30);
void nodes_to_constant(Circuit* orig, double alpha, int batch_count=30);
int nodes_to_constant_ver2(Circuit* orig, double alpha, int batch_count=30, int numberOfNodesToCheck=1);
extern vector<unordered_map<Node*, bitset<BITSET_WIDTH>>> goldens;
void build_golden_output_vals_for_error_eval(Circuit* orig, int batch_count=30);

//pre
void preprocess_and_postprocess(Global& global);

#endif