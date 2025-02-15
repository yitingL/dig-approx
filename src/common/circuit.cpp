#include <bits/stdc++.h>
#include "../header.h"

using namespace std;
int Circuit::getSize(){
    int size=0;
    for(Node* n:allNodes){
        if(n->type==PI||n->type==PO_BUF||n->type==PO_NOT||n->type==CONST0||n->type==CONST1) continue;
        else{
            size++;
        }
    }
    return size;
}
int Circuit::getLevel(){
    vector<Node*> nodesInTopoOrder=getNodesInTopoOrder();
    unordered_map<Node*,int> level;
    for(Node* n:nodesInTopoOrder) level[n]=0;
    for(Node* n:nodesInTopoOrder){
        if(n->type==PI||n->type==CONST0||n->type==CONST1||n->type==PO_BUF||n->type==PO_NOT){
            for(Node* nfo:n->fo){
                level[nfo]=level[n];
            }
        }else{
            for(Node* nfo:n->fo){
                level[nfo]=max(level[nfo],level[n]+1);
            }
        }
    }
    int max_level=0;
    for(Node* n:nodesInTopoOrder) max_level=max(max_level,level[n]);
    return max_level;
}
void Circuit::replace_nt_with_ns(Node* nt, Node* ns){
    if(nt==ns) return;
    //move all nt's fanout to ns.
    for(int i=0; i<nt->fo.size(); i++){
        ns->fo.push_back(nt->fo[i]);
        //change fanin variable of nt's fanout.
        for(int j=0; j<nt->fo[i]->fi.size(); j++){
            if(nt->fo[i]->fi[j] == nt){
                nt->fo[i]->fi[j] = ns;
                break;
            }
        }
    }
    nt->fo.clear();
}
void Circuit::dfsFromPOstoCleanNetwork(){   //All PIs will be kept.
    set<Node*> visited;
    stack<Node*> stk;
    for(Node* po:POs) stk.push(po);
    for(Node* pi:PIs) visited.insert(pi);
    visited.insert(constant0);
    visited.insert(constant1);
    while(!stk.empty()){
        Node* n=stk.top();
        stk.pop();
        if(!visited.count(n)){
            visited.insert(n);
            for(Node* nfi:n->fi) stk.push(nfi);
        }
    }
    //remove all nodes not visited, using remove_if only needs one scan of the vector.
    auto new_end = std::remove_if(allNodes.begin(), allNodes.end(),
        [&visited](Node* n) {return !visited.count(n);});
    allNodes.erase(new_end, allNodes.end());
}
vector<Node*> Circuit::getNodesInTopoOrder_POsAtEnd(){
    vector<Node*> ans;
    queue<Node*> Q;
    unordered_map<Node*, int> in_degree;
    for(Node* n: allNodes){
        in_degree[n]=n->fi.size();
    }
    if(constant0) Q.push(constant0);
    if(constant1) Q.push(constant1);
    for(Node* pi:PIs) Q.push(pi);
    while(!Q.empty()){
        Node* n=Q.front();
        Q.pop();
        ans.push_back(n);
        for(Node* nfo:n->fo){
            --in_degree[nfo];
            if(in_degree[nfo]==0 && nfo->type!=PO_BUF && nfo->type!=PO_NOT) Q.push(nfo);
        }
    }
    for(Node* n:allNodes)
        if(n->type==PO_BUF || n->type==PO_NOT) ans.push_back(n);
    assert(ans.size()==allNodes.size());
    return ans;
}
vector<Node*> Circuit::getNodesInTopoOrder(){
    vector<Node*> ans;
    queue<Node*> Q;
    unordered_map<Node*, int> in_degree;
    for(Node* n: allNodes){
        in_degree[n]=n->fi.size();
    }
    if(constant0) Q.push(constant0);
    if(constant1) Q.push(constant1);
    for(Node* pi:PIs) Q.push(pi);
    while(!Q.empty()){
        Node* n=Q.front();
        Q.pop();
        ans.push_back(n);
        for(Node* nfo:n->fo){
            --in_degree[nfo];
            if(in_degree[nfo]==0) Q.push(nfo);
        }
    }
    assert(ans.size()==allNodes.size());
    return ans;
}
void Circuit::read_dig(string filename){
    ifstream file;
    file.open(filename);

    Node* current_node;
    vector<string> planes;
    string line_str;
    while(getline(file, line_str)){
        stringstream line_ss(line_str);
        string token;
        line_ss >> token;
        if(token == ".model"){
            line_ss >> token;
            modelName = token;
        }else if(token == ".inputs"){
            while(line_ss >> token){
                Node* n = new Node(token,PI);
                nodeName_2_node[token]=n;
                allNodes.push_back(n);
                PIs.push_back(n);
            }
        }else if(token == ".outputs"){
            while(line_ss >> token){
                Node* n = new Node(token,PO_BUF);
                nodeName_2_node[token]=n;
                allNodes.push_back(n);
                POs.push_back(n);
            }
        }else if(token == ".names" || token == ".end"){
            if(!planes.empty()){
                if(!current_node){
                    cerr << "error: line of unknown syntax: " << line_str << endl;
                    exit(0);
                }
                if(current_node->fi.size()==0 && planes.size()==1){
                    if(planes[0]=="0"){
                        current_node->type=CONST0;
                        constant0=current_node;
                    }else if(planes[0]=="1"){
                        current_node->type=CONST1;
                        constant1=current_node;
                    }else{
                        cerr << "error in current_node with fi.size = 0" << endl;
                        exit(0);
                    }
                }else if(current_node->type==PO_BUF && current_node->fi.size()==1 && planes.size()==1){
                    if(planes[0]=="0 1"){
                        current_node->type=PO_NOT;
                    }else if(planes[0]=="1 1"){
                        current_node->type=PO_BUF;
                    }else{
                        cerr << "error in reading a PO" << endl;
                        exit(0);
                    }
                }else if(current_node->fi.size()==3 && planes.size()==2){
                    // A Dot gate without INV should be:
                    //      xyz   out
                    //      100   1
                    //      0-1   1
                    bool x_flip;
                    if(planes[0][0]=='1' && planes[1][0]=='0'){
                        x_flip=false;
                    }else if(planes[0][0]=='0' && planes[1][0]=='1'){
                        x_flip=true;
                    }else{
                        cerr << "the node's supposed not a DOT gate? x port errors: " << current_node->name << endl;
                        exit(0);
                    }
                    bool y_flip;
                    if(planes[0][1]=='0' && planes[1][1]=='-'){
                        y_flip=false;
                    }else if(planes[0][1]=='1' && planes[1][1]=='-'){
                        y_flip=true;
                    }else{
                        cerr << "the node's supposed not a DOT gate? y port errors: " << current_node->name << endl;
                        exit(0);
                    }
                    bool z_flip;
                    if(planes[0][2]=='0' && planes[1][2]=='1'){
                        z_flip=false;
                    }else if(planes[0][2]=='1' && planes[1][2]=='0'){
                        z_flip=true;
                    }else{
                        cerr << "the node's supposed not a DOT gate? z port errors: " << current_node->name << endl;
                        exit(0);
                    }
                    bool o_flip;
                    if(planes[0][4]=='1' && planes[1][4]=='1'){
                        o_flip=false;
                    }else if(planes[0][4]=='0' && planes[1][4]=='0'){
                        o_flip=true;
                    }
                    if(!x_flip && !y_flip && !z_flip && !o_flip){
                        current_node->type=Dot;
                    }else if(x_flip && !y_flip && !z_flip && !o_flip){
                        current_node->type=DotXI;
                    }else if(!x_flip && y_flip && !z_flip && !o_flip){
                        current_node->type=DotYI;
                    }else if(!x_flip && !y_flip && z_flip && !o_flip){
                        current_node->type=DotZI;
                    }else if(x_flip && y_flip && !z_flip && !o_flip){
                        current_node->type=DotXYI;
                    }else if(x_flip && !y_flip && z_flip && !o_flip){
                        current_node->type=DotXZI;
                    }else if(!x_flip && y_flip && z_flip && !o_flip){
                        current_node->type=DotYZI;
                    }else if(x_flip && y_flip && z_flip && !o_flip){
                        current_node->type=DotXYZI;
                    }else if(!x_flip && !y_flip && !z_flip && o_flip){
                        current_node->type=DotOI;
                    }else if(x_flip && !y_flip && !z_flip && o_flip){
                        current_node->type=DotXOI;
                    }else if(!x_flip && y_flip && !z_flip && o_flip){
                        current_node->type=DotYOI;
                    }else if(!x_flip && !y_flip && z_flip && o_flip){
                        current_node->type=DotZOI;
                    }else if(x_flip && y_flip && !z_flip && o_flip){
                        current_node->type=DotXYOI;
                    }else if(x_flip && !y_flip && z_flip && o_flip){
                        current_node->type=DotXZOI;
                    }else if(!x_flip && y_flip && z_flip && o_flip){
                        current_node->type=DotYZOI;
                    }else if(x_flip && y_flip && z_flip && o_flip){
                        current_node->type=DotXYZOI;
                    }
                }else{
                    cerr << "error in identifying a node: " << current_node->name << endl;
                    exit(0);
                }
                planes.clear();
            }
            if(token==".names"){
                vector<Node*> signal_list;
                while(line_ss >> token){
                    Node* n;
                    if(!nodeName_2_node.count(token)){
                        n = new Node(token, Dot);
                        nodeName_2_node[token]=n;
                        allNodes.push_back(n);
                    }else{
                        n = nodeName_2_node[token];
                    }
                    signal_list.push_back(n);
                }
                current_node=signal_list[signal_list.size()-1];
                for(Node* n:signal_list){
                    if(n!=current_node){
                        current_node->fi.push_back(n);
                        n->fo.push_back(current_node);
                    }
                }
            }
        }else{
            planes.push_back(line_str);
        }
    }
    if(!constant0 || !constant1){
        cerr<<"the circuit has no constant nodes." << endl;
        exit(0);
    }
}
void Circuit::write_dig(string filename){
    ofstream file;
    file.open(filename);

    file << ".model" << " " << modelName << endl;
    file << ".inputs";
    for(Node* n:PIs){
        file << " " << n->name;
    }
    file << endl;
    file << ".outputs";
    for(Node* n:POs){
        file << " " << n->name;
    }
    file << endl;
    file << ".names" << " " << constant0->name << endl;
    file << "0" << endl;
    file << ".names" << " " << constant1->name << endl;
    file << "1" << endl;
    for(Node* n:allNodes){
        if(n->type==PI || n->type==PO_BUF || n->type==PO_NOT || n->type==CONST0 || n->type==CONST1) continue;
        file << ".names";
        for(Node* fi:n->fi){
            file << " " << fi->name;
        }
        file << " " << n->name << endl;
        // A Dot gate without INV should be:
        //      xyz   out
        //      100   1
        //      0-1   1
        string out1="100 1";
        string out2="0-1 1";
        bool x_flip=(n->type==DotXI||n->type==DotXYI||n->type==DotXZI||n->type==DotXYZI||n->type==DotXOI||n->type==DotXYOI||n->type==DotXZOI||n->type==DotXYZOI);
        bool y_flip=(n->type==DotYI||n->type==DotXYI||n->type==DotYZI||n->type==DotXYZI||n->type==DotYOI||n->type==DotXYOI||n->type==DotYZOI||n->type==DotXYZOI);
        bool z_flip=(n->type==DotZI||n->type==DotYZI||n->type==DotXZI||n->type==DotXYZI||n->type==DotZOI||n->type==DotYZOI||n->type==DotXZOI||n->type==DotXYZOI);
        bool o_flip=(n->type==DotOI||n->type==DotXOI||n->type==DotYOI||n->type==DotZOI||n->type==DotXYOI||n->type==DotYZOI||n->type==DotXZOI||n->type==DotXYZOI);
        if(x_flip){ out1[0]='0'; out2[0]='1'; }
        if(y_flip){ out1[1]='1'; out2[1]='-'; }
        if(z_flip){ out1[2]='1'; out2[2]='0'; }
        if(o_flip){ out1[4]='0'; out2[4]='0'; }
        file << out1 << endl;
        file << out2 << endl;               
    }
    for(Node* n:POs){
        file << ".names" << " " << n->fi[0]->name << " " << n->name << endl;
        if(n->type==PO_BUF){
            file << "1 1" << endl;
        }else if(n->type==PO_NOT){
            file << "0 1" << endl;
        }else{
            cerr << "error in writing a PO" << endl;
            exit(0);
        }
    }
    file << ".end" << endl;
    file.close();
}