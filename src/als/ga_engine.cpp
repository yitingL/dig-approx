#include <bits/stdc++.h>
#include "../header.h"
using namespace std;

//*/
#define debug(x) do{cout<<#x<<": "<<x<<endl;}while(0)
/*/
#define debug(x)
//*/ 

inline vector<int> shuffle_fanins(const vector<int>& fis){
    vector<int> ans=fis;
    assert(fis.size()==3);
    random_shuffle(ans.begin(),ans.end());
    return ans;
}

inline void print_genes(const vector<vector<int>>& v){
    for(auto i:v){
        cout<<"(";
        for(auto j:i){
            cout<< j << ",";
        }
        cout <<")";
    }
    cout<<endl;
}

double GA_Engine::get_difference_rate(Chromosome& chromo1, Chromosome& chromo2){
    int diff=0;
    for(int i=0; i<chromo1.genes.size(); i++){
        for(int j=0; j<chromo1.genes[i].size(); j++){
            if(chromo1.genes[i][j]!=chromo2.genes[i][j]){
                diff++;
            }
        }
    }
    return (double) diff / (double) (chromo1.genes.size()*chromo1.genes[0].size());
}   
void GA_Engine::simulate_chromosome_by_bitset(Chromosome& chromo, unordered_map<Node*,bitset<BITSET_WIDTH>> &vals){
    vals[orig_circuit->constant0].reset();
    vals[orig_circuit->constant1].set();
    for(vector<int> gene : chromo.genes){
        NodeType type= (NodeType) gene[0];
        if(type==PI||type==CONST0||type==CONST1||type==PO_BUF||type==PO_NOT){
            cerr << "error in simulate chromosome by bitset" << endl;
            exit(0);
        }
        Node* nd=getNode[gene[1]];
        auto x=vals[getNode[gene[2]]];
        auto y=vals[getNode[gene[3]]];
        auto z=vals[getNode[gene[4]]];
        switch(type){
            case (int)Dot:{
                vals[nd] = x^(z|(x&y));
                break;
            }
            case (int)DotXI:{
                x=~x;
                vals[nd] = x^(z|(x&y));
                break;
            }
            case (int)DotYI:{
                y = ~y;
                vals[nd] = x^(z|(x&y));
                break;
            }
            case (int)DotZI:{
                z=~z;
                vals[nd] = x^(z|(x&y));
                break;
            }
            case (int)DotXYI:{
                x=~x; y=~y; 
                vals[nd] = x^(z|(x&y));
                break;
            }
            case (int)DotYZI:{
                y=~y; z=~z;
                vals[nd] = x^(z|(x&y));
                break;
            }
            case (int)DotXZI:{
                x=~x; z=~z;
                vals[nd] = x^(z|(x&y));
                break;
            }
            case (int)DotXYZI:{
                x=~x; y=~y; z=~z;
                vals[nd] = x^(z|(x&y));
                break;
            }
            case (int)DotOI:{
                vals[nd] = ~( x^(z|(x&y)));
                break;
            }
            case (int)DotXOI:{
                x=~x;
                vals[nd] = ~( x^(z|(x&y)));
                break;
            }
            case (int)DotYOI:{
                y = ~y;
                vals[nd] = ~( x^(z|(x&y)));
                break;
            }
            case (int)DotZOI:{
                z=~z;
                vals[nd] = ~( x^(z|(x&y)));
                break;
            }
            case (int)DotXYOI:{
                x=~x; y=~y; 
                vals[nd] = ~( x^(z|(x&y)));
                break;
            }
            case (int)DotYZOI:{
                y=~y; z=~z;
                vals[nd] = ~( x^(z|(x&y)));
                break;
            }
            case (int)DotXZOI:{
                x=~x; z=~z;
                vals[nd] = ~( x^(z|(x&y)));
                break;
            }
            case (int)DotXYZOI:{
                x=~x; y=~y; z=~z;
                vals[nd] = ~( x^(z|(x&y)));
                break;                    
            }
            default:{
                cerr << "error in simulation (chromosome): " << nd->name << ", " << nd->type << endl;
                exit(0);
                break;                    
            }
        }
    }
    //get PO's value
    for(Node* po:orig_circuit->POs){
        switch(po->type){
            case PO_BUF:{
                vals[po]=vals[po->fi[0]];
                break;
            }
            case PO_NOT:{
                vals[po]=~vals[po->fi[0]];
                break;
            }
            default:{
                cerr << "error in check_error_rate_of()" << endl;
                exit(0);
                break;
            }
        }
    }
}
double GA_Engine::check_error_rate_of(Chromosome& chromo, int batch_count){
    assert(batch_count <= goldens.size());
    int error_count=0;
    int error_count2=0;
    for(int i=0; i<batch_count; i++){
        //unordered_map<Node*,bitset<BITSET_WIDTH>> goldens;
        unordered_map<Node*,bitset<BITSET_WIDTH>> vals;
        for(Node* nd: orig_circuit->PIs){
            Node* nd2=golden_circuit->nodeName_2_node[nd->name];
            //goldens[nd2]=rand_bits();
            vals[nd]=goldens[i][nd2];
        }
        //simulate_by_bitset(golden_topo, goldens);
        simulate_chromosome_by_bitset(chromo, vals);
        bitset<BITSET_WIDTH> temp1;
        temp1.reset();
        for(Node* nd: orig_circuit->POs){
            Node* nd2=golden_circuit->nodeName_2_node[nd->name];
            temp1 |= goldens[i][nd2] ^ vals[nd];
            error_count2 += (goldens[i][nd2] ^ vals[nd]).count();
        }
        error_count += temp1.count();
    }
    chromo.error_rate2=(double) error_count2 / (double)(BITSET_WIDTH*batch_count*orig_circuit->POs.size());
    return (double) error_count / (double)(BITSET_WIDTH*batch_count);
}

double GA_Engine::check_error_dist_of(Chromosome& chromo, int batch_count){
    assert(batch_count <= goldens.size());
    double error_dist=0;
    for(int i=0; i<batch_count; i++){
        //unordered_map<Node*,bitset<BITSET_WIDTH>> goldens;
        unordered_map<Node*,bitset<BITSET_WIDTH>> vals;
        for(Node* nd1: orig_circuit->PIs){
            Node* nd2=golden_circuit->nodeName_2_node[nd1->name];
            //goldens[nd2]=rand_bits();
            vals[nd1]=goldens[i][nd2];
        }
        //simulate_by_bitset(golden_topo, goldens);
        simulate_chromosome_by_bitset(chromo, vals);
        for(auto& vec : global->PO_vectors){
            for(int k=0; k<vec.size(); k++){
                Node* nd1 = orig_circuit->nodeName_2_node[vec[k]];      //vals is on this circuit
                Node* nd2 = golden_circuit->nodeName_2_node[vec[k]];
                int diff = (goldens[i][nd2] ^ vals[nd1]).count();
                error_dist += ((double) diff * global->PO_weights[k]);
            }
        }
    }
    chromo.error_rate2=0.0;
    return error_dist / (BITSET_WIDTH*batch_count);
}

set<Node*> GA_Engine::get_active_Nodes(Chromosome& chromo){
    set<Node*> visited;
    set<Node*> activeNodes;
    stack<Node*> stk;
    for(Node* po:orig_circuit->POs) stk.push(po);
    while(!stk.empty()){
        Node* n=stk.top();
        stk.pop();
        if(!visited.count(n)){
            visited.insert(n);
            if(n->type==PO_BUF||n->type==PO_NOT){
                for(Node* fi:n->fi){
                    stk.push(fi);
                }
            }else if(n->type==PI||n->type==CONST0||n->type==CONST1){
                continue;
            }else{
                activeNodes.insert(n);
                int geneIndex=getGeneIndex[n];
                for(int i=2;i<=4;i++){
                    Node* fi=getNode[chromo.genes[geneIndex][i]];
                    stk.push(fi);
                }
            }
        }
    }
    return activeNodes;      
}
int GA_Engine::check_depth_of(const Chromosome& chromo){
    unordered_map<Node*,int> depth;
    for(Node* nd:orig_circuit->PIs) depth[nd]=0;
    depth[orig_circuit->constant0]=0;
    depth[orig_circuit->constant1]=0;
    for(auto& gene : chromo.genes){
        NodeType type=(NodeType) gene[0];
        if(type==PI||type==CONST0||type==CONST1||type==PO_BUF||type==PO_NOT){
            cerr << "error in check depth of" << endl;
            exit(0);
        }
        Node* nd=getNode[gene[1]];
        auto x=depth[getNode[gene[2]]];
        auto y=depth[getNode[gene[3]]];
        auto z=depth[getNode[gene[4]]];
        depth[nd]=max({x,y,z})+1;
    }
    int ans=0;
    for(Node* po:orig_circuit->POs){
        NodeType type=po->type;
        if(type==PO_BUF||type==PO_NOT){
            ans=max(ans,depth[po->fi[0]]);
        }else{
            cerr << "error in check depth of" << endl;
            exit(0);
        }
    }
    return ans;
}
int GA_Engine::check_size_of(Chromosome& chromo){
    chromo.activeGenes.clear();
    set<Node*> visited;
    set<Node*> activeNodes;
    stack<Node*> stk;
    for(Node* po:orig_circuit->POs) stk.push(po);
    while(!stk.empty()){
        Node* n=stk.top();
        stk.pop();
        if(!visited.count(n)){
            visited.insert(n);
            if(n->type==PO_BUF||n->type==PO_NOT){
                for(Node* fi:n->fi){
                    stk.push(fi);
                }
            }else if(n->type==PI||n->type==CONST0||n->type==CONST1){
                continue;
            }else{
                activeNodes.insert(n);
                chromo.activeGenes.push_back(getGeneIndex[n]);
                int geneIndex=getGeneIndex[n];
                for(int i=2;i<=4;i++){
                    Node* fi=getNode[chromo.genes[geneIndex][i]];
                    stk.push(fi);
                }
            }
        }
    }
    return activeNodes.size();
}
void GA_Engine::gen_blif(Chromosome& chromo, string filename){
    set<Node*> activeNodes=get_active_Nodes(chromo);
    ofstream file;
    file.open(filename);
    file << ".model " << orig_circuit->modelName << endl;
    file << ".inputs";
    for(Node* pi:orig_circuit->PIs){
        file << " " << pi->name;
    }
    file << endl;
    file << ".outputs";
    for(Node* po:orig_circuit->POs){
        file << " " << po->name;
    }
    file << endl;
    file << ".names " << orig_circuit->constant0->name << endl;
    file << "0" << endl;
    file << ".names " << orig_circuit->constant1->name << endl;
    file << "1" << endl;
    for(auto gene:chromo.genes){
        Node* n=getNode[gene[1]];
        if(activeNodes.count(n)){
            file << ".names";
            for(int i=2; i<=4; i++){
                file << " " << getNode[gene[i]]->name;
            }
            file << " " << n->name << endl;
            NodeType type= (NodeType) gene[0];
            // A Dot gate without INV should be:
            //      xyz   out
            //      100   1
            //      0-1   1
            string out1="100 1";
            string out2="0-1 1";
            bool x_flip=(type==DotXI||type==DotXYI||type==DotXZI||type==DotXYZI||type==DotXOI||type==DotXYOI||type==DotXZOI||type==DotXYZOI);
            bool y_flip=(type==DotYI||type==DotXYI||type==DotYZI||type==DotXYZI||type==DotYOI||type==DotXYOI||type==DotYZOI||type==DotXYZOI);
            bool z_flip=(type==DotZI||type==DotYZI||type==DotXZI||type==DotXYZI||type==DotZOI||type==DotYZOI||type==DotXZOI||type==DotXYZOI);
            bool o_flip=(type==DotOI||type==DotXOI||type==DotYOI||type==DotZOI||type==DotXYOI||type==DotYZOI||type==DotXZOI||type==DotXYZOI);
            if(x_flip){ out1[0]='0'; out2[0]='1'; }
            if(y_flip){ out1[1]='1'; out2[1]='-'; }
            if(z_flip){ out1[2]='1'; out2[2]='0'; }
            if(o_flip){ out1[4]='0'; out2[4]='0'; }
            file << out1 << endl;
            file << out2 << endl;
        }
    }
    for(Node* n:orig_circuit->POs){
        file << ".names " << n->fi[0]->name << " " << n->name << endl;
        if(n->type==PO_BUF){
            file << "1 1" << endl;
        }else{
            file << "0 1" << endl;
        }
    }
    file << ".end" << endl;
    file.close();
}
Chromosome GA_Engine::encode_orig_circuit(){
    Chromosome chromo;
    int id=0;
    for(Node* nd : orig_topo){
        getID[nd]=id;
        getNode[id]=nd;
        id++;
        if(nd->type==PI||nd->type==PO_BUF||nd->type==PO_NOT||nd->type==CONST0||nd->type==CONST1){
            continue;
        }else{
            vector<int> gene;
            gene.push_back(nd->type);
            gene.push_back(getID[nd]);
            for(Node* nfi:nd->fi){
                gene.push_back(getID[nfi]);
            }
            getGeneIndex[nd]=chromo.genes.size();
            chromo.genes.push_back(gene);
        }
    }
    return chromo;
}
vector<Chromosome> GA_Engine::gen_initial_population(Chromosome& orig, int size){
    vector<Chromosome> cands;
    cands.push_back(orig);
    for(int i=1; i<size; i++){
        Chromosome cand=orig.gen_new();
        int locus=rand_ull(0, cand.genes.size()-1);
        bool mode=rand_ull(0,1);
        if(mode){
            for(int g=locus; g>=0; g--){
                vector<int>& gene=cand.genes[g];
                NodeType type_new=(NodeType) rand_ull(Dot, DotXYZOI);
                gene[0]=type_new;
                //shuffle_vector(gene,2);
                int id=gene[1];
                for(int j=2; j<gene.size(); j++){
                    gene[j]=rand_ull(0, id-1);
                }
            }
        }else{
            for(int g=locus; g<cand.genes.size(); g++){
                vector<int>& gene=cand.genes[g];
                NodeType type_new=(NodeType) rand_ull(Dot, DotXYZOI);
                gene[0]=type_new;
                //shuffle_vector(gene,2);
                int id=gene[1];
                for(int j=2; j<gene.size(); j++){
                    gene[j]=rand_ull(0, id-1);
                }
            }
        }
        cands.push_back(cand);
    }
    return cands;
}
void GA_Engine::mutate(vector<Chromosome>& cands,int ptr1, int ptr2, int count_for_each, int step_size){
    for(int i=ptr1; i<=ptr2; i++){
        for(int j=0; j<count_for_each; j++){
            Chromosome new_cand=cands[i].gen_new();
            for(int k=0; k<step_size; k++){
                int mut_type=rand_ull(0, 1);
                if(mut_type==0){    //REMOVE
                    int target_gene=new_cand.activeGenes[rand_ull(0,new_cand.activeGenes.size()-1)];
                    //int target_gene=rand_ull(0,new_cand.genes.size()-1);
                    //int select_fi=rand_ull(0,2);
                    for(int g=target_gene+1;g<new_cand.genes.size()-1;g++){
                        int id=new_cand.genes[g][1];
                        for(int gi=2; gi<=4; gi++){
                            if(new_cand.genes[g][gi]==new_cand.genes[target_gene][1]){
                                //new_cand.genes[g][gi]=new_cand.genes[target_gene][select_fi+2];
                                new_cand.genes[g][gi]=rand_ull(0, id-1);
                            }
                        }
                    }
                }else{              //CHANGE
                    int mut=rand_ull(0, (new_cand.genes.size()-1)*4);
                    if(mut%4==0){   //change gate type
                        new_cand.genes[mut/4][0]=(NodeType) rand_ull(Dot, DotXYZOI);
                    }else{          //change fi
                        int id=new_cand.genes[mut/4][1];
                        new_cand.genes[mut/4][mut%4+1]=rand_ull(0, id-1);
                    }
                }
            }
            cands.push_back(new_cand);
        }
    }
    // for(int i=0; i<iteration; i++){
    //     int mut_type=rand_ull(0, 1);
    //     Chromosome &targetChr=cands[rand_ull(0, bound-1)];
    //     if(mut_type==0){    //REMOVE
    //         int idx=rand_ull(0,targetChr.genes.size()-1);
    //         auto &targetGene=targetChr.genes[idx];
    //         int select_fi=rand_ull(0,2);
    //         for(int i=idx+1; i<targetChr.genes.size()-1; i++){
    //             for(int j=2; j<=4; j++){
    //                 if(targetChr.genes[i][j]==targetGene[1]){
    //                     targetChr.genes[i][j]=targetGene[select_fi+2];
    //                 }
    //             }
    //         }
    //     }else if(mut_type==1){  //CHANGE
    //         int mut=rand_ull(0, (targetChr.genes.size()-1)*4);
    //         if(mut%4==0){   //change gate type
    //             targetChr.genes[mut/4][0]=(NodeType) rand_ull(Dot, DotXYZOI);
    //         }else{          //change fi
    //             int id=targetChr.genes[mut/4][1];
    //             targetChr.genes[mut/4][mut%4+1]=rand_ull(0, id-1);
    //         }
    //     }
    // }
}
void GA_Engine::crossover_onepoint(vector<Chromosome>& cands, int bound,int count_for_each){
    for(int i=1;i<=bound; i++){
        for(int time=0; time<count_for_each; time++){
            Chromosome new_cand;
            int j=rand_ull(0, i-1);
            int cut=rand_ull(0, cands[i].genes.size()-1);
            for(int g=0; g<cut; g++){
                new_cand.genes.push_back(cands[i].genes[g]);
            }
            for(int g=cut; g<cands[j].genes.size(); g++){
                new_cand.genes.push_back(cands[j].genes[g]);
            }
            cands.push_back(new_cand);
        }
    }
}
void GA_Engine::crossover_universal(vector<Chromosome>& cands, int bound,int count_for_each){
    for(int i=1; i<=bound; i++){
        for(int time=0; time<count_for_each; time++){
            Chromosome new_cand;
            int j=rand_ull(0, i-1);
            for(int g=0; g<cands[i].genes.size(); g++){
                if(rand_ull(0,1)){
                    new_cand.genes.push_back(cands[i].genes[g]);
                }else{
                    new_cand.genes.push_back(cands[j].genes[g]);
                }
            }
            cands.push_back(new_cand);
        }
    }
    // for(int i=bound-1; i>0; i--){
    //     for(int time=0; time<count_for_each; time++){
    //         Chromosome cand1=cands[i].gen_new();
    //         //Chromosome cand2=cands[rand_ull(i+1, bound-1)].gen_new();
    //         Chromosome cand2=cands[rand_ull(0, i-1)].gen_new();
    //         for(int g=0; g<cand1.genes.size(); g++){
    //             if(rand_ull(0,1)){
    //                 swap(cand1.genes[g], cand2.genes[g]);
    //             }
    //         }
    //         cands.push_back(cand1);
    //         cands.push_back(cand2);
    //     }
    // }
}
int GA_Engine::learn_from_mutation_result(vector<Chromosome>& cands, int mutate_ptr1, int mutate_ptr2, int mutate_step_size){
    double acc=0;
    double avg_mutate_error_rate;
    for(int i=mutate_ptr1; i<=mutate_ptr2; i++) acc+=cands[i].error_rate;
    avg_mutate_error_rate=acc/(mutate_ptr2-mutate_ptr1+1);
    //for(int i=0; i<cands.size(); i++) acc+=cands[i].error_rate;
    //avg_mutate_error_rate=acc/cands.size();
    //cout << "<stats> avg_mutate_error_rate: " << avg_mutate_error_rate << endl;
    int ans;
    if(avg_mutate_error_rate < ERROR_RATE_THRESHOLD){
        ans=(mutate_step_size+1)*1.3;
    }else if(avg_mutate_error_rate > ERROR_RATE_THRESHOLD*2){
        ans=(mutate_step_size-1)*0.7;
    }else{
        ans=mutate_step_size;
    }
    if(ans<=0) return 1;
    else return min(ans, orig_size/2);
}
void GA_Engine::learn_from_survivor(vector<Chromosome>& cands){
    // bool all_illegal=true;
    // for(Chromosome& cand:cands){
    //     if(cand.error_rate<=ERROR_RATE_THRESHOLD){
    //         all_illegal=false;
    //         break;
    //     }
    // }
    // if(all_illegal){
    //     cands[cands.size()-1]=bestChromo;
    // }
    bool best_not_inside=true;
    for(Chromosome& cand:cands){
        if(cand.error_rate==bestChromo.error_rate && cand.size==bestChromo.size){
            best_not_inside=false;
            break;
        }
    }
    if(best_not_inside){
        cands[cands.size()-1]=bestChromo;
    }
}
void GA_Engine::eval(vector<Chromosome> & cands){
    if(global->ed_mode==""){
        for(Chromosome& cand:cands){
            if(global->check_if_time_out()) return;
            if(cand.data_avail){
                //cand.fitness_val*=0.95;   //aging
            }else{
                cand.error_rate=check_error_rate_of(cand,3);
                cand.size=check_size_of(cand);
                cand.diff_rate=get_difference_rate(cand,orig);
                cand.depth=check_depth_of(cand);
                cand.area_rate=(double) cand.size / (double) orig_size;
                cand.depth_rate=(double) cand.depth / (double) orig_depth;
                double penalty;
                if(cand.error_rate<=ERROR_RATE_THRESHOLD){
                    penalty=1;
                }else{
                    penalty=exp(1000*(cand.error_rate/ERROR_RATE_THRESHOLD));
                }
                if(cand.size < bestChromo.size && cand.error_rate <= ERROR_RATE_THRESHOLD){
                    cand.error_rate=check_error_rate_of(cand,30);
                    if(cand.error_rate <= ERROR_RATE_THRESHOLD){
                        bestChromo=cand;
                        terminate_counter=0;
                        //cout << "# output size: " << bestChromo.size << endl;
                        //cout << "# error rate: " << bestChromo.error_rate << endl;
                        gen_blif(bestChromo, output_filename);
                    }
                }else if(cand.size==bestChromo.size && cand.error_rate < bestChromo.error_rate && cand.error_rate <= ERROR_RATE_THRESHOLD){
                    cand.error_rate=check_error_rate_of(cand,30);
                    if(cand.error_rate < bestChromo.error_rate && cand.error_rate < ERROR_RATE_THRESHOLD){
                        bestChromo=cand;
                        //cout << "# output size: " << bestChromo.size << endl;
                        //cout << "# error rate: " << bestChromo.error_rate << endl;
                        gen_blif(bestChromo, output_filename);
                    }
                }
                cand.fitness_val=1/(1*penalty*(cand.error_rate + cand.error_rate2) + 100*cand.area_rate + 0*cand.depth_rate);
                cand.fitness_val2=1/(1*penalty*(cand.error_rate + cand.error_rate2) + 1*cand.area_rate + 0*cand.depth_rate);
                // cand.fitness_val=1/(1*penalty*(cand.error_rate + cand.error_rate2) + 1*cand.area_rate + 100*cand.depth_rate);
                // cand.fitness_val2=1/(1*penalty*(cand.error_rate + cand.error_rate2) + 1*cand.area_rate + 1*cand.depth_rate);
                // if(cand.fitness_val > bestChromo.fitness_val && cand.error_rate <= ERROR_RATE_THRESHOLD){
                //     cand.error_rate=check_error_rate_of(cand,30);
                //     if(cand.error_rate <= ERROR_RATE_THRESHOLD){
                //         bestChromo=cand;
                //         terminate_counter=0;
                //         cout << "# output size: " << bestChromo.size << endl;
                //         cout << "# error rate: " << bestChromo.error_rate << endl;
                //         gen_blif(bestChromo, output_filename);
                //     }
                // }
                // cand.fitness_val=1/(1*penalty*(cand.error_rate + cand.error_rate2) + 99*cand.area_rate + 1*cand.depth_rate);
                // cand.fitness_val2=1/(1*penalty*(cand.error_rate + cand.error_rate2) + 1*cand.area_rate + 1*cand.depth_rate);
                
                cand.data_avail=true;
            }
        }
    }else{
        for(Chromosome& cand:cands){
            if(global->check_if_time_out()) return;
            cand.error_rate=check_error_dist_of(cand,3);
            cand.size=check_size_of(cand);
            cand.diff_rate=get_difference_rate(cand,orig);
            cand.depth=check_depth_of(cand);
            cand.area_rate=(double) cand.size / (double) orig_size;
            cand.depth_rate=(double) cand.depth / (double) orig_depth;
            double penalty;
            if(cand.error_rate<=ERROR_RATE_THRESHOLD){
                penalty=1;
            }else{
                penalty=exp(100*(cand.error_rate/ERROR_RATE_THRESHOLD));
                //penalty=1;
            }
            if(cand.size < bestChromo.size && cand.error_rate <= ERROR_RATE_THRESHOLD){
                cand.error_rate=check_error_dist_of(cand,30);
                if(cand.error_rate <= ERROR_RATE_THRESHOLD){
                    bestChromo=cand;
                    terminate_counter=0;
                    cout << "# output size: " << bestChromo.size << endl;
                    cout << "# error rate: " << bestChromo.error_rate << endl;
                    gen_blif(bestChromo, output_filename);
                }
            }else if(cand.size==bestChromo.size && cand.error_rate < bestChromo.error_rate && cand.error_rate <= ERROR_RATE_THRESHOLD){
                cand.error_rate=check_error_dist_of(cand,30);
                if(cand.error_rate < bestChromo.error_rate && cand.error_rate < ERROR_RATE_THRESHOLD){
                    bestChromo=cand;
                    cout << "# output size: " << bestChromo.size << endl;
                    cout << "# error rate: " << bestChromo.error_rate << endl;
                    gen_blif(bestChromo, output_filename);
                }
            }
            cand.fitness_val=1/(1*penalty*(cand.error_rate) + 100*cand.area_rate);
            cand.fitness_val2=1/(1*penalty*(cand.error_rate) + 1*cand.area_rate);
            cand.data_avail=true;
        }
    }
}
void GA_Engine::select(vector<Chromosome>& cands, int trim_size1, int trim_size2){
    //cout << "selecting... from " << cands.size() << endl;
    // sort(cands.begin(), cands.end(), [](Chromosome& a, Chromosome& b){
    //     return a.fitness_val > b.fitness_val;
    // });
    // if(cands.size()>size){
    //     cands.resize(size);
    // }
    vector<Chromosome> cands2=cands;
    sort(cands.begin(), cands.end(), [](const Chromosome& a, const Chromosome& b){
        return a.fitness_val > b.fitness_val;
    });
    if(cands.size()>trim_size1){
        cands.resize(trim_size1);
    }
    sort(cands2.begin(), cands2.end(), [](const Chromosome& a, const Chromosome& b){
        return a.fitness_val2 > b.fitness_val2;
    });
    cands.insert(cands.end(), cands2.begin(), cands2.begin()+trim_size2);
    assert(cands.size()==trim_size1+trim_size2);
}
GA_Engine::GA_Engine(Global& global, Circuit* circuit, Circuit* golden_circuit, string output_blif, double error_rate_threshold){
    this->orig_circuit = circuit;
    orig_topo=orig_circuit->getNodesInTopoOrder_POsAtEnd();
    orig_size=orig_circuit->getSize();
    orig_depth=orig_circuit->getLevel();
    
    this->golden_circuit = golden_circuit;
    golden_topo=golden_circuit->getNodesInTopoOrder_POsAtEnd();
    golden_size=golden_circuit->getSize();
    
    bestChromo.size=INT_MAX;
    bestChromo.fitness_val=0;
    terminate_signal=false;

    output_filename=output_blif;
    this->ERROR_RATE_THRESHOLD=error_rate_threshold; //for safety
    this->global=&global;
}
void Chromosome::print_info(){
    //cout << fixed << setprecision(5);
    cout<<"error_rate: "<<error_rate
    //<<"\terror_rate2: "<<error_rate2
    <<"\tdepth: " <<depth
    <<"\tsize: "<<size
    //<<"\tarea_rate: " <<area_rate
    //<<"\tfitness_val: "<<fitness_val
    //<<"\tdiff_rate: "<<diff_rate
    <<endl;
}
int GA_Engine::start(int total_generation, int round_indicator){
    orig=encode_orig_circuit();
    vector<Chromosome> cands=gen_initial_population(orig,500);
    eval(cands);
    if(global->check_if_time_out()) return 1;
    int setSize1=15, setSize2=15;
    select(cands, setSize1,setSize2);
    //for(auto cand:cands) cand.print_info();
    int generation_count=0;
    int mutate_step_size1=1;
    int mutate_step_size2=1;
    while(terminate_counter<total_generation){
        //cout << "-----------------" << endl;
        //cout << "GA Generation: " << round_indicator << "." << ++generation_count << endl;
        //cout << "orig_size: " << orig_size << endl;
        mutate(cands,0,setSize1-1,2,mutate_step_size1);
        mutate(cands,setSize1,setSize1+setSize2-1,2,mutate_step_size2);
        crossover_universal(cands,setSize1+setSize2-1,1);
        //crossover_onepoint(cands,24,1);
        eval(cands);
        if(global->check_if_time_out()) return 1;
        mutate_step_size1=learn_from_mutation_result(cands, setSize1+setSize2, 3*setSize1+setSize2-1, mutate_step_size1);
        mutate_step_size2=learn_from_mutation_result(cands, 3*setSize1+setSize2, 3*setSize1+3*setSize2-1, mutate_step_size2);
        //cout << "<stats>... mutation_step_size1: " << mutate_step_size1 << endl;
        //cout << "<stats>... mutation_step_size2: " << mutate_step_size2 << endl;
        select(cands,setSize1,setSize2);
        
        //for(auto cand:cands) cand.print_info();
        //if((++generation_count % 10)==0) cout << ">>> " << endl;
        cout << "[" << round_indicator << ", " << generation_count << ", " << terminate_counter << "] Completed." << endl;
        //bestChromo.print_info();
        //cout << "terminate_counter: " << terminate_counter << endl;
        terminate_counter++;
        if(bestChromo.size<=orig_size*0.5){
            return 0;
        }
    }
    if(bestChromo.size>orig_size){   //GA took no effect
        return 2;
    }
    return 1;
}