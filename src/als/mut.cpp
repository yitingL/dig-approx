#include "../header.h"
#include <bits/stdc++.h>
using namespace std;

unordered_map<int,int> GA_Engine::getMA(bool sa, Chromosome& chromo, int ntID, int domID, unordered_map<int, vector<int>>& fos_of, unordered_map<int, bitset<65536>>& dom, unordered_map<int, bitset<262144>>& domSI){
    bool a=domSI[domID][3*domID];
    bool b=domSI[domID][3*domID+1]; 
    bool c=domSI[domID][3*domID+2];
    unordered_map<int,int> ans;
    unordered_map<int, bitset<4>> assign_bit;
    unordered_map<int, bitset<4>> assign_known;
    bitset<4> assign_valid;     //to record if the assignment is valid: 0 means conflict
    for(auto g : chromo.genes){
        assign_known[g[1]].reset();
    }
    auto g = chromo.genes[getGeneIndex[getNode[domID]]];
    if(a && !b && !c){
        assign_bit[g[3]].reset();
        assign_bit[g[4]].reset();
        assign_bit[g[3]][0]=0; assign_bit[g[4]][0]=0;
        assign_bit[g[3]][1]=0; assign_bit[g[4]][1]=1;
        assign_bit[g[3]][2]=1; assign_bit[g[4]][2]=1;
        assign_known[g[3]].reset();
        assign_known[g[4]].reset();
        assign_known[g[3]][0]=1; assign_known[g[4]][0]=1;
        assign_known[g[3]][1]=1; assign_known[g[4]][1]=1;
        assign_known[g[3]][2]=1; assign_known[g[4]][2]=1;
        assign_valid.reset();
        assign_valid[0]=1;
        assign_valid[1]=1;
        assign_valid[2]=1;
    }else if(!a && b && !c){
        assign_bit[g[2]].reset();
        assign_bit[g[4]].reset();
        assign_bit[g[2]][0]=1; assign_bit[g[4]][0]=0;
        assign_known[g[2]].reset();
        assign_known[g[4]].reset();
        assign_known[g[2]][0]=1; assign_known[g[4]][0]=1;
        assign_valid.reset();
        assign_valid[0]=1;
    }else if(!a && !b && c){
        assign_bit[g[3]].reset();
        assign_bit[g[4]].reset();
        assign_bit[g[3]][0]=0; assign_bit[g[4]][0]=0;
        assign_bit[g[3]][1]=0; assign_bit[g[4]][1]=1;
        assign_bit[g[3]][2]=1; assign_bit[g[4]][2]=1;
        assign_known[g[3]].reset();
        assign_known[g[4]].reset();
        assign_known[g[3]][0]=1; assign_known[g[4]][0]=1;
        assign_known[g[3]][1]=1; assign_known[g[4]][1]=1;
        assign_known[g[3]][2]=1; assign_known[g[4]][2]=1;
        assign_valid.reset();
        assign_valid[0]=1;
        assign_valid[1]=1;
        assign_valid[2]=1;
    }else if(a && b && !c){
        assign_bit[g[2]].reset();
        assign_bit[g[3]].reset();
        assign_bit[g[4]].reset();
        assign_bit[g[2]][0]=(sa?1:0); assign_bit[g[3]][0]=(sa?1:0); assign_bit[g[4]][0]=1;
        assign_bit[g[2]][1]=(sa?0:1); assign_bit[g[3]][1]=(sa?0:1); assign_bit[g[4]][1]=1;
        assign_bit[g[2]][2]=(sa?1:0); assign_bit[g[3]][2]=(sa?0:1);
        assign_bit[g[2]][3]=(sa?0:1); assign_bit[g[3]][3]=(sa?1:0);
        assign_known[g[2]].reset();
        assign_known[g[3]].reset();
        assign_known[g[4]].reset();
        assign_known[g[2]][0]=1; assign_known[g[3]][0]=1; assign_known[g[4]][0]=1;
        assign_known[g[2]][1]=1; assign_known[g[3]][1]=1; assign_known[g[4]][1]=1;
        assign_known[g[2]][2]=1; assign_known[g[3]][2]=1;
        assign_known[g[2]][3]=1; assign_known[g[3]][3]=1;
        assign_valid.reset();
        assign_valid[0]=1;
        assign_valid[1]=1;
        assign_valid[2]=1;
        assign_valid[3]=1;
    }else if(!a && b && c){
        assign_bit[g[2]].reset();
        assign_bit[g[3]].reset();
        assign_bit[g[4]].reset();
        assign_bit[g[2]][0]=0; assign_bit[g[3]][0]=(sa?1:0); assign_bit[g[4]][0]=(sa?0:1);
        assign_bit[g[2]][1]=0; assign_bit[g[3]][1]=(sa?0:1); assign_bit[g[4]][1]=(sa?1:0);
        assign_bit[g[3]][2]=(sa?1:0); assign_bit[g[4]][2]=(sa?1:0);
        assign_bit[g[3]][3]=(sa?0:1); assign_bit[g[4]][3]=(sa?0:1);
        assign_known[g[2]].reset();
        assign_known[g[3]].reset();
        assign_known[g[4]].reset();
        assign_known[g[2]][0]=1; assign_known[g[3]][0]=1; assign_known[g[4]][0]=1;
        assign_known[g[2]][1]=1; assign_known[g[3]][1]=1; assign_known[g[4]][1]=1;
        assign_known[g[3]][2]=1; assign_known[g[4]][2]=1;
        assign_known[g[3]][3]=1; assign_known[g[4]][3]=1;
        assign_valid.reset();
        assign_valid[0]=1;
        assign_valid[1]=1;
        assign_valid[2]=1;
        assign_valid[3]=1;
    }else if(a && !b && c){
        assign_bit[g[2]].reset();
        assign_bit[g[3]].reset();
        assign_bit[g[4]].reset();
        assign_bit[g[2]][0]=(sa?1:0); assign_bit[g[3]][0]=1; assign_bit[g[4]][0]=(sa?0:1);
        assign_bit[g[2]][1]=(sa?0:1); assign_bit[g[3]][1]=0; assign_bit[g[4]][1]=(sa?1:0);
        assign_known[g[2]].reset();
        assign_known[g[3]].reset();
        assign_known[g[4]].reset();
        assign_known[g[2]][0]=1; assign_known[g[3]][0]=1; assign_known[g[4]][0]=1;
        assign_known[g[2]][1]=1; assign_known[g[3]][1]=1; assign_known[g[4]][1]=1;
        assign_valid.reset();
        assign_valid[0]=1;
        assign_valid[1]=1;
    }else if(a && b && c){
        assign_bit[g[2]].reset();
        assign_bit[g[3]].reset();
        assign_bit[g[4]].reset();
        assign_bit[g[2]][0]=(sa?1:0); assign_bit[g[3]][0]=(sa?1:0); assign_bit[g[4]][0]=(sa?0:1);
        assign_bit[g[2]][0]=(sa?0:1); assign_bit[g[3]][0]=(sa?0:1); assign_bit[g[4]][0]=(sa?1:0);
        assign_known[g[2]].reset();
        assign_known[g[3]].reset();
        assign_known[g[4]].reset();
        assign_known[g[2]][0]=1; assign_known[g[3]][0]=1; assign_known[g[4]][0]=1;
        assign_known[g[2]][1]=1; assign_known[g[3]][1]=1; assign_known[g[4]][1]=1;
        assign_valid.reset();
        assign_valid[0]=1;
        assign_valid[1]=1;
    }else{
        ans[-1]=0;
        return ans;
    }
    if(sa){
        assign_bit[ntID].reset();
        assign_known[ntID].set();    
    }else{
        assign_bit[ntID].set();
        assign_known[ntID].set();            
    }
    assign_bit[getID[orig_circuit->constant1]].set();
    assign_known[getID[orig_circuit->constant1]].set();
    assign_bit[getID[orig_circuit->constant0]].reset();  
    assign_known[getID[orig_circuit->constant0]].set();    
    /*logic implication by bitset
    The forward implications are:
        x   y   z  | Dot
        ---------------- 
        --  11  10 | 10
        11  11  -- | 10
        10  --  10 | 10             
        10  --  11 | 11  
        11  --  11 | 10
        11  10  10 | 11
        else       | 0-

    The backward implications are:
        Dot y   z  | x          Dot x   z  | y      Dot x   y  | z
        ---------------         ---------------     ---------------       
        11  11  -- | 10         11  11  -- | 10     11  11  -- | 10
        11  --  10 | 11         11  --  10 | 10     11  --  11 | 11         
        11  --  11 | 10         10  11  10 | 11     10  10  -- | 10
        10  --  11 | 11         else       | 0-     11  10  -- | 11  
        10  10  10 | 10                             10  11  10 | 11
        else       | 0-                             else       | 0-

    Check-Then-Assign: For any node, we have to check implicated value with the original value before
    assigning the value to the node. If the values conflicts, we assert the valid bit of the assginments to 0.

        orig imp | new  valid
        ---------------------
        0-   0-  | 0-   1
        0-   10  | 10   1
        0-   11  | 11   1
        10   0-  | 10   1
        11   0-  | 11   1
        10   10  | 10   1
        11   11  | 11   1
        10   11  | --   0
        11   10  | --   0
    */
    bitset<65536> implied_nodes;
    implied_nodes.reset(); 
    list<int> Q;
    for(int k : fos_of[getID[orig_circuit->constant1]]){
        Q.push_back(k);
    }
    for(int k : fos_of[getID[orig_circuit->constant0]]){
        Q.push_back(k);
    }
    for(int k : fos_of[ntID]){
        Q.push_back(k);
    }
    Q.push_back(domID);
    Q.push_back(g[2]);
    Q.push_back(g[3]);
    Q.push_back(g[4]);
    Q.push_back(ntID);
    while(implied_nodes.count()<300 && !Q.empty()){
        int nd=Q.front();
        Q.pop_front();
        implied_nodes.set(nd);
        auto g = chromo.genes[getGeneIndex[getNode[nd]]];
        auto ob = assign_bit[g[1]];
        auto ok = assign_known[g[1]];
        auto xb = assign_bit[g[2]];
        auto xk = assign_known[g[2]];
        auto yb = assign_bit[g[3]];
        auto yk = assign_known[g[3]];
        auto zb = assign_bit[g[4]];
        auto zk = assign_known[g[4]];
        switch((NodeType) g[0]){
            case Dot:{
                break;
            }
            case DotXI:{
                xb.flip();
                break;
            }
            case DotYI:{
                yb.flip();
                break;
            }
            case DotZI:{
                zb.flip();
                break;
            }
            case DotXYI:{
                xb.flip(); yb.flip();
                break;
            }
            case DotYZI:{
                yb.flip(); zb.flip();
                break;
            }
            case DotXZI:{
                xb.flip(); zb.flip();
                break;
            }
            case DotXYZI:{
                xb.flip(); yb.flip(); zb.flip();
                break;
            }
            case DotOI:{
                ob.flip();
                break;
            }
            case DotXOI:{
                ob.flip(); xb.flip();
                break;
            }
            case DotYOI:{
                ob.flip(); yb.flip();
                break; 
            }
            case DotZOI:{
                ob.flip(); zb.flip();
                break;
            }
            case DotXYOI:{
                ob.flip(); xb.flip(); yb.flip();
                break;
            }
            case DotYZOI:{
                ob.flip(); yb.flip(); zb.flip();
                break;
            }
            case DotXZOI:{
                ob.flip(); xb.flip(); zb.flip();
                break;
            }
            case DotXYZOI:{
                ob.flip(); xb.flip(); yb.flip(); zb.flip();
                break;                    
            }
            default:{
                cerr << "error in getMA (0): " << g[0] << endl;
                exit(0);
                break;
            }
        }
        {   //forward
            auto origk = ok;
            auto origb = ob;
            auto impk = yk&yb&zk&~zb|xk&xb&yk&yb|xk&~xb&zk|xk&xb&zk&zb|xk&xb&yk&~yb&zk&~zb;
            auto impb = xk&~xb&zk&zb |xk&xb&yk&~yb&zk&~zb;
            ok = origk | impk;
            ob = origk&origb | impk&impb;
            assign_valid = assign_valid & ~((origk & origb & impk & ~impb) | (origk & ~impb & impk & impb));
            // we AND the assign_valid because we want to prevent the case that the assignment is already invalid.
            auto changed = assign_valid & ~origk & impk;
            // if valid==1 && origk == 0 && new1 == 1, we have to push the fanout to the list to do implication.
            if(changed.any()){
                for(int k : fos_of[g[1]]) Q.push_back(k);
            }
        }
        {   //backward - x
            auto origk = xk;
            auto origb = xb;
            auto impk = ok&ob&yk&yb | ok&ob&zk&~zb | ok&zk&zb | ok&~ob&yk&~yb&zk&~zb;
            auto impb = ~ob&zk&zb | ob&zk&~zb;
            xk = origk | impk;
            xb = origk&origb | impk&impb;
            assign_valid = assign_valid & ~((origk & origb & impk & ~impb) | (origk & ~impb & impk & impb));
            auto changed = assign_valid & ~origk & impk;
            if(changed.any()){
                Q.push_back(g[2]);
                for(int k : fos_of[g[2]])
                    if(k!=g[1]) Q.push_back(k);
            }
        }
        {   //backward - y
            auto origk = yk;
            auto origb = yb;
            auto impk = ok&ob&xk&xb | ok&ob&zk&~zb | ok&~ob&xk&xb&zk&~zb;
            auto impb = ~ob&xk&xb&zk&~zb;
            yk = origk | impk;
            yb = origk&origb | impk&impb;
            assign_valid = assign_valid & ~((origk & origb & impk & ~impb) | (origk & ~impb & impk & impb));  
            auto changed = assign_valid & ~origk & impk; 
            if(changed.any()){
                Q.push_back(g[3]);
                for(int k : fos_of[g[3]])
                    if(k!=g[1]) Q.push_back(k);
            }
        }
        {   //backward - z
            auto origk = zk;
            auto origb = zb;
            auto impk = ok&ob&xk&xb | ok&ob&yk&yb | ok&xk&~xb | ok&~ob&xk&xb&yk&~yb;
            auto impb = ob&yk&yb    | ob&xk&~xb |~ob&xk&xb&yk&~yb;
            zk = origk | impk;
            zb = origk&origb | impk&impb;
            assign_valid = assign_valid & ~((origk & origb & impk & ~impb) | (origk & ~impb & impk & impb));
            auto changed = assign_valid & ~origk & impk; 
            if(changed.any()){
                Q.push_back(g[4]);
                for(int k : fos_of[g[4]])
                    if(k!=g[1]) Q.push_back(k);
            }
        }
        switch((NodeType) g[0]){
            case Dot:{
                break;
            }
            case DotXI:{
                xb.flip();
                break;
            }
            case DotYI:{
                yb.flip();
                break;
            }
            case DotZI:{
                zb.flip();
                break;
            }
            case DotXYI:{
                xb.flip(); yb.flip();
                break;
            }
            case DotYZI:{
                yb.flip(); zb.flip();
                break;
            }
            case DotXZI:{
                xb.flip(); zb.flip();
                break;
            }
            case DotXYZI:{
                xb.flip(); yb.flip(); zb.flip();
                break;
            }
            case DotOI:{
                ob.flip();
                break;
            }
            case DotXOI:{
                ob.flip(); xb.flip();
                break;
            }
            case DotYOI:{
                ob.flip(); yb.flip();
                break; 
            }
            case DotZOI:{
                ob.flip(); zb.flip();
                break;
            }
            case DotXYOI:{
                ob.flip(); xb.flip(); yb.flip();
                break;
            }
            case DotYZOI:{
                ob.flip(); yb.flip(); zb.flip();
                break;
            }
            case DotXZOI:{
                ob.flip(); xb.flip(); zb.flip();
                break;
            }
            case DotXYZOI:{
                ob.flip(); xb.flip(); yb.flip(); zb.flip();
                break;                    
            }
            default:{
                cerr << "error in getMA (1): " << g[0] << endl;
                exit(0);
                break;
            }
        }
        assign_bit[g[1]] = ob;
        assign_known[g[1]] = ok;
        assign_bit[g[2]] = xb;
        assign_known[g[2]] = xk;
        assign_bit[g[3]] = yb;
        assign_known[g[3]] = yk;
        assign_bit[g[4]] = zb;
        assign_known[g[4]] = zk;
    }
    if(assign_valid.none()){
        ans[-1]=1;
        return ans;
    }
    /*
    To find the intersection of bit1 and bit0, we can use the following method:
        Define: <1> valid_bits,   <2> assign_known,   <3> assign_bit
        Condition [1]: <1> & <2> & <3> all zero --> there is no assignment-(1,1,1) in these three bitsets.
        Condition [2]: <1> & <2> & ~<3> all zero --> there is no assignment-(1,1,0) in these three bitsets.
        Condition [3]: <1> & ~<2> all zero --> there is no assignment-(1,0,-) in these three bitsets.
    */
    for(auto g : chromo.genes){
        auto cond1 = assign_valid & assign_known[g[1]] & assign_bit[g[1]];
        auto cond2 = assign_valid & assign_known[g[1]] & ~assign_bit[g[1]];
        auto cond3 = assign_valid & ~assign_known[g[1]];
        if(!cond1.none() && cond2.none() && cond3.none()){
            ans[g[1]] = 1;
        }else if(cond1.none() && !cond2.none() && cond3.none()){
            ans[g[1]] = 0;
        }
    }
    return ans;
}
int GA_Engine::check_subs(int nt, Chromosome& chromo, unordered_map<int,int>& ma0, unordered_map<int,int>& ma1, unordered_map<int, bitset<65536>>& tfo){
    for(auto& g : chromo.genes){
        if(g[1]==nt) continue;
        if(tfo[nt].test(g[1])) continue;
        if(ma0[g[1]]==1 && ma1[g[1]]==0){
            return g[1];
        }
        if(ma0[g[1]]==0 && ma1[g[1]]==1){
            switch((NodeType) g[0]){
                case Dot:{
                    g[0]=DotOI;
                    break;
                }
                case DotXI:{
                    g[0]=DotXOI;
                    break;
                }
                case DotYI:{
                    g[0]=DotYOI;
                    break;
                }
                case DotZI:{
                    g[0]=DotZOI;
                    break;
                }
                case DotXYI:{
                    g[0]=DotXYOI;
                    break;
                }
                case DotYZI:{
                    g[0]=DotYZOI;
                    break;
                }
                case DotXZI:{
                    g[0]=DotXZOI;
                    break;
                }
                case DotXYZI:{
                    g[0]=DotXYZOI;
                    break;
                }
                case DotOI:{
                    g[0]=Dot;
                    break;
                }
                case DotXOI:{
                    g[0]=DotXI;
                    break;
                }
                case DotYOI:{
                    g[0]=DotYI;
                    break;
                }
                case DotZOI:{
                    g[0]=DotZI;
                    break;
                }
                case DotXYOI:{
                    g[0]=DotXYI;
                    break;
                }
                case DotYZOI:{
                    g[0]=DotYZI;
                    break;
                }
                case DotXZOI:{
                    g[0]=DotXZI;
                    break;
                }
                case DotXYZOI:{
                    g[0]=DotXYZI;
                    break;
                }
                default:{
                    cerr << "error in check_subs" << endl;
                    exit(0);
                    break;
                }
            }
            return g[1];
        }
    }
    return -1;
}
int GA_Engine::get_merged(Chromosome& chromo, int nt){ 
    /*
    domSI:  001  fault from input x
            010  fault from input y
            100  fault from input z
    dom[nd][id] : 1 if id is in the dominators set of nd
    */
    assert(chromo.genes.size() < 65536);
    unordered_map<int, int> level;
    for(auto nd:orig_circuit->PIs) level[getID[nd]]=0;
    level[getID[orig_circuit->constant0]]=0;
    level[getID[orig_circuit->constant1]]=0;
    for(auto g : chromo.genes){
        NodeType type=(NodeType) g[0];
        if(type==PI||type==CONST0||type==CONST1||type==PO_BUF||type==PO_NOT){
            cerr << "error in merge_node_with_a_merger" << endl;
            exit(0);
        }
        auto x=level[g[2]];
        auto y=level[g[3]];
        auto z=level[g[4]];
        level[g[1]]=max({x,y,z})+1; 
    }
    unordered_map<int, vector<int>> fos_of;
    for(auto g : chromo.genes){
        fos_of[g[2]].push_back(g[1]);
        fos_of[g[3]].push_back(g[1]);
        fos_of[g[4]].push_back(g[1]);
    }
    unordered_map<int, bitset<65536>> dom; 
    unordered_map<int, bitset<262144>> domSI;
    unordered_map<int, bitset<65536>> tfo; 
    for(int i=chromo.genes.size()-1; i>=0; i--){
        auto g=chromo.genes[i];
        NodeType type=(NodeType) g[0];
        if(type==PI||type==CONST0||type==CONST1||type==PO_BUF||type==PO_NOT){
            cerr << "error in merge_node_with_a_merger" << endl;
            exit(0);
        }
        if(fos_of[g[1]].size()==0){
            dom[g[1]].reset();
            domSI[g[1]].reset();
            dom[g[1]].set(g[1]);
        }else{
            dom[g[1]]=dom[fos_of[g[1]][0]];
            domSI[g[1]]=domSI[fos_of[g[1]][0]];
            for(int fo : fos_of[g[1]]){
                dom[g[1]] &= dom[fo];
                domSI[g[1]] |= domSI[fo];
                auto foGene = chromo.genes[getGeneIndex[getNode[fo]]];
                if(foGene[2]==g[1]){
                    domSI[fo].set(3*fo);
                }else if(foGene[3]==g[1]){
                    domSI[fo].set(3*fo+1);
                }else if(foGene[4]==g[1]){
                    domSI[fo].set(3*fo+2);
                }
            }
            dom[g[1]].set(g[1]);
        }
        if(g[1]==nt) break;
    }
    for(int i=chromo.genes.size()-1; i>=0; i--){
        auto g=chromo.genes[i];
        NodeType type=(NodeType) g[0];
        if(type==PI||type==CONST0||type==CONST1||type==PO_BUF||type==PO_NOT){
            cerr << "error in merge_node_with_a_merger" << endl;
            exit(0);
        }
        if(fos_of[g[1]].size()==0){
            tfo[g[1]].reset();
            tfo[g[1]].set(g[1]);
        }else{
            tfo[g[1]]=tfo[fos_of[g[1]][0]];
            for(int fo : fos_of[g[1]]){
                tfo[g[1]] |= tfo[fo];
            }
            tfo[g[1]].set(g[1]);
        }
        if(g[1]==nt) break;
    }
    int cnt=0;
    for(auto g : chromo.genes){
        if(dom[nt].test(g[1])){
            if((g[1] != nt) && level[g[1]] <= level[nt]+5){
                cnt++;
                unordered_map<int,int> ma0 = getMA(0,chromo,nt,g[1],fos_of,dom,domSI);
                if(ma0[-1]==1){
                    return getID[orig_circuit->constant0];
                }else{
                    unordered_map<int,int> ma1 = getMA(1,chromo,nt,g[1],fos_of,dom,domSI);
                    if(ma1[-1]==1){
                        return getID[orig_circuit->constant1];
                    }else if(ma1[-1]!=0 && ma0[-1]!=0){
                        return check_subs(nt,chromo,ma0,ma1,tfo);
                    }
                }
            }
        }
    }
    assert(cnt <= 5);
    return -1;
}