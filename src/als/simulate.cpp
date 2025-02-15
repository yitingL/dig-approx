#include <bits/stdc++.h>
#include "../header.h"
using namespace std;

vector<unordered_map<Node*, bitset<BITSET_WIDTH>>> goldens;
const unordered_map<Node*,int>& constructMFFC(Circuit* orig){  //maximum transitive fanout-free cone
    unordered_map<Node*,int> mffc_size;
    for(Node* nd:orig->allNodes){
        unordered_set<Node*> mffc;
        queue<Node*> q;         //BFS
        q.push(nd);
        mffc.insert(nd);
        while(!q.empty()){
            Node* cur=q.front();
            q.pop();
            if(!mffc.count(cur)){
                bool all_fanout_in_mffc=true;
                for(Node* fanout:cur->fo){
                    if(!mffc.count(fanout)){
                        all_fanout_in_mffc=false;
                        break;
                    }
                }
                if(all_fanout_in_mffc){
                    mffc.insert(cur);
                    for(Node* fanin:cur->fi) q.push(fanin);
                }
            }
        }
        mffc_size[nd]=mffc.size();
    }
    return mffc_size;
}
void simulate_by_bitset(vector<Node*> &nodesInTopoOrder, unordered_map<Node*,bitset<BITSET_WIDTH>> &vals){
    for(Node* nd : nodesInTopoOrder){
        switch(nd->type){
            case PI:{
                continue;
            }
            case CONST0:{
                vals[nd].reset();
                break;
            }
            case CONST1:{
                vals[nd].set();
                break;
            }
            case PO_BUF:{
                vals[nd]=vals[nd->fi[0]];
                break;
            }
            case PO_NOT:{
                vals[nd]=~vals[nd->fi[0]];
                break;
            }
            case Dot:{
                auto x = vals[nd->fi[0]];
                auto y = vals[nd->fi[1]];
                auto z = vals[nd->fi[2]];
                vals[nd] = x^(z|(x&y));
                break;
            }
            case DotXI:{
                auto x = vals[nd->fi[0]];
                auto y = vals[nd->fi[1]];
                auto z = vals[nd->fi[2]];
                x=~x;
                vals[nd] = x^(z|(x&y));
                break;
            }
           case DotYI:{
                auto x = vals[nd->fi[0]];
                auto y = vals[nd->fi[1]];
                auto z = vals[nd->fi[2]];
                y = ~y;
                vals[nd] = x^(z|(x&y));
                break;
            }
            case DotZI:{
                auto x = vals[nd->fi[0]];
                auto y = vals[nd->fi[1]];
                auto z = vals[nd->fi[2]];
                z=~z;
                vals[nd] = x^(z|(x&y));
                break;
            }
            case DotXYI:{
                auto x = vals[nd->fi[0]];
                auto y = vals[nd->fi[1]];
                auto z = vals[nd->fi[2]];
                x=~x; y=~y; 
                vals[nd] = x^(z|(x&y));
                break;
            }
            case DotYZI:{
                auto x = vals[nd->fi[0]];
                auto y = vals[nd->fi[1]];
                auto z = vals[nd->fi[2]];
                y=~y; z=~z;
                vals[nd] = x^(z|(x&y));
                break;
            }
            case DotXZI:{
                auto x = vals[nd->fi[0]];
                auto y = vals[nd->fi[1]];
                auto z = vals[nd->fi[2]];
                x=~x; z=~z;
                vals[nd] = x^(z|(x&y));
                break;
            }
            case DotXYZI:{
                auto x = vals[nd->fi[0]];
                auto y = vals[nd->fi[1]];
                auto z = vals[nd->fi[2]];
                x=~x; y=~y; z=~z;
                vals[nd] = x^(z|(x&y));
                break;
            }
            case DotOI:{
                auto x = vals[nd->fi[0]];
                auto y = vals[nd->fi[1]];
                auto z = vals[nd->fi[2]];
                vals[nd] = ~( x^(z|(x&y)));
                break;
            }
            case DotXOI:{
                auto x = vals[nd->fi[0]];
                auto y = vals[nd->fi[1]];
                auto z = vals[nd->fi[2]];
                x=~x;
                vals[nd] = ~( x^(z|(x&y)));
                break;
            }
            case DotYOI:{
                auto x = vals[nd->fi[0]];
                auto y = vals[nd->fi[1]];
                auto z = vals[nd->fi[2]];
                y = ~y;
                vals[nd] = ~( x^(z|(x&y)));
                break;
            }
            case DotZOI:{
                auto x = vals[nd->fi[0]];
                auto y = vals[nd->fi[1]];
                auto z = vals[nd->fi[2]];
                z=~z;
                vals[nd] = ~( x^(z|(x&y)));
                break;
            }
            case DotXYOI:{
                auto x = vals[nd->fi[0]];
                auto y = vals[nd->fi[1]];
                auto z = vals[nd->fi[2]];
                x=~x; y=~y; 
                vals[nd] = ~( x^(z|(x&y)));
                break;
            }
            case DotYZOI:{
                auto x = vals[nd->fi[0]];
                auto y = vals[nd->fi[1]];
                auto z = vals[nd->fi[2]];
                y=~y; z=~z;
                vals[nd] = ~( x^(z|(x&y)));
                break;
            }
            case DotXZOI:{
                auto x = vals[nd->fi[0]];
                auto y = vals[nd->fi[1]];
                auto z = vals[nd->fi[2]];
                x=~x; z=~z;
                vals[nd] = ~( x^(z|(x&y)));
                break;
            }
            case DotXYZOI:{
                auto x = vals[nd->fi[0]];
                auto y = vals[nd->fi[1]];
                auto z = vals[nd->fi[2]];
                x=~x; y=~y; z=~z;
                vals[nd] = ~( x^(z|(x&y)));
                break;
            }
            default:{
                cerr << "error in simulation: " << nd->name << ", " << nd->type << endl;
                exit(0);
                break;
            }
        }
    }
}
void build_golden_output_vals_for_error_eval(Circuit* orig, int batch_count){
    goldens.clear();
    for(int i=0; i<batch_count; i++){
        unordered_map<Node*,bitset<BITSET_WIDTH>> vals;
        for(Node* nd:orig->PIs){
            vals[nd]=rand_bits();
        }
        vector<Node*> nodesInTopoOrder=orig->getNodesInTopoOrder();
        simulate_by_bitset(nodesInTopoOrder,vals);
        goldens.push_back(vals);
    }
}
double check_error_rate(Circuit* gol, Circuit* modified, int batch_count){
    assert(gol->PIs.size()==modified->PIs.size());
    assert(gol->POs.size()==modified->POs.size());
    assert(goldens.size()>=batch_count);
    int error_count=0;
    //vector<Node*> orig_nodes_topo=orig->getNodesInTopoOrder();
    vector<Node*> modified_nodes_topo=modified->getNodesInTopoOrder();
    for(int i=0; i<batch_count; i++){
        //unordered_map<Node*,bitset<BITSET_WIDTH>> goldens;
        unordered_map<Node*,bitset<BITSET_WIDTH>> vals;
        for(Node* nd1: gol->PIs){
            Node* nd2 = modified->nodeName_2_node[nd1->name];
            //goldens[nd1]=rand_bits();
            vals[nd2]=goldens[i][nd1];
        }
        //simulate_by_bitset(orig_nodes_topo, goldens);
        simulate_by_bitset(modified_nodes_topo, vals);
        bitset<BITSET_WIDTH> temp1;
        temp1.reset();
        for(Node* nd1: gol->POs){
            Node* nd2 = modified->nodeName_2_node[nd1->name];
            temp1 |= goldens[i][nd1] ^ vals[nd2];
        }
        error_count += temp1.count();
    }
    return (double) error_count / (BITSET_WIDTH*batch_count);
}
double check_error_distance(Circuit* gol, Circuit* modified, Global& global,int batch_count){
    assert(gol->PIs.size()==modified->PIs.size());
    assert(gol->POs.size()==modified->POs.size());
    assert(global.ed_mode!="");
    assert(goldens.size()>=batch_count);
    assert(batch_count>0);
    double error_dist=0;
    //vector<Node*> orig_nodes_topo=orig->getNodesInTopoOrder();
    vector<Node*> modified_nodes_topo=modified->getNodesInTopoOrder();
    for(int i=0; i<batch_count; i++){
        //unordered_map<Node*,bitset<BITSET_WIDTH>> goldens;
        unordered_map<Node*,bitset<BITSET_WIDTH>> vals;
        for(Node* nd1: gol->PIs){
            Node* nd2 = modified->nodeName_2_node[nd1->name];
            //goldens[nd1]=rand_bits();
            vals[nd2]=goldens[i][nd1];
        }
        //simulate_by_bitset(orig_nodes_topo, goldens);
        simulate_by_bitset(modified_nodes_topo, vals);
        for(auto& vec : global.PO_vectors){
            for(int k=0; k<vec.size(); k++){
                Node* nd1 = gol->nodeName_2_node[vec[k]];
                Node* nd2 = modified->nodeName_2_node[vec[k]];
                int diff = (goldens[i][nd1] ^ vals[nd2]).count();
                error_dist += ((double)diff * global.PO_weights[k]);
            }
        }
    }
    return error_dist / (BITSET_WIDTH*batch_count);
}
void nodes_to_constant(Circuit* orig, double alpha, int batch_count){
    unordered_map<Node*,int> times_of_1;
    for(Node* nd : orig->allNodes) times_of_1[nd]=0;
    for(int i=0; i<batch_count; i++){
        unordered_map<Node*,bitset<BITSET_WIDTH>> vals;
        for(Node* nd:orig->PIs){
            vals[nd]=rand_bits();
        }
        vector<Node*> nodesInTopoOrder=orig->getNodesInTopoOrder();
        simulate_by_bitset(nodesInTopoOrder,vals);
        for(Node* nd:orig->allNodes){
            if(nd->type==PI||nd->type==PO_BUF||nd->type==PO_NOT||nd->type==CONST0||nd->type==CONST1) continue;
            else times_of_1[nd]+=vals[nd].count();
        }
    }
    for(Node* nd:orig->allNodes){
        if(nd->type==PI||nd->type==PO_BUF||nd->type==PO_NOT||nd->type==CONST0||nd->type==CONST1) continue;
        else{
            Node* ns=NULL;
            assert(times_of_1[nd]<=batch_count*BITSET_WIDTH);
            if(times_of_1[nd]>=batch_count*BITSET_WIDTH*alpha){
                ns=orig->constant1;
            }else if(times_of_1[nd]<=batch_count*BITSET_WIDTH*(1-alpha)){
                ns=orig->constant0;
            }
            if(ns){
                orig->replace_nt_with_ns(nd,ns);
            }
        }
    }
    orig->dfsFromPOstoCleanNetwork();
}
int nodes_to_constant_ver2(Circuit* orig, double alpha, int batch_count, int numberOfNodesToCheck){
    //this version removes one node, instead of several nodes, and then re-simulate in each iteration.
    unordered_map<Node*,int> mffc_size=constructMFFC(orig);
    vector<Node*> allNodesMFFCsorted=orig->allNodes;
    sort(allNodesMFFCsorted.begin(),allNodesMFFCsorted.end(),[&](Node* a, Node* b){
        return mffc_size[a]>mffc_size[b];
    });
    unordered_map<Node*,int> times_of_1;
    for(Node* nd : orig->allNodes) times_of_1[nd]=0;
    for(int i=0; i<batch_count; i++){
        unordered_map<Node*,bitset<BITSET_WIDTH>> vals;
        for(Node* nd:orig->PIs){
            vals[nd]=rand_bits();
        }
        vector<Node*> nodesInTopoOrder=orig->getNodesInTopoOrder();
        simulate_by_bitset(nodesInTopoOrder,vals);
        for(Node* nd:orig->allNodes){
            if(nd->type==PI||nd->type==PO_BUF||nd->type==PO_NOT||nd->type==CONST0||nd->type==CONST1) continue;
            else times_of_1[nd]+=vals[nd].count();
        }
    }
    int counter=0;
    for(Node* nd:allNodesMFFCsorted){
        if(nd->type==PI||nd->type==PO_BUF||nd->type==PO_NOT||nd->type==CONST0||nd->type==CONST1) continue;
        else{
            Node* ns=NULL;
            assert(times_of_1[nd]<=batch_count*BITSET_WIDTH);
            if(times_of_1[nd]>=batch_count*BITSET_WIDTH*alpha){
                ns=orig->constant1;
            }else if(times_of_1[nd]<=batch_count*BITSET_WIDTH*(1-alpha)){
                ns=orig->constant0;
            }
            if(ns){
                orig->replace_nt_with_ns(nd,ns);
                counter++;
                if(counter>=numberOfNodesToCheck) break;
            }
        }
    }
    orig->dfsFromPOstoCleanNetwork();
    return counter;
}
