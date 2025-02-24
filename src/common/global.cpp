#include <bits/stdc++.h>
#include "../header.h"
using namespace std;

void Global::setup_PO_vec_and_PO_weights(string ed_mode){
    /* 0: none
       1: adder
            { f[0:127], cOut } <--> { po0-po127, po128 }
       2: div
            quotient[0:63] <--> po0-po63
            remainder[0:63] <--> po64-po127
       3: log2
            result[0:31] <--> po0-po31
       4: multiplier
            f[0:127] <--> po0-po127
       5: sqrt
            asqrt[0:63] <--> po0-po63
       6: square 
            asquared[0:127] <--> po0-po127
       7: bar
            result[0:127] <--> po0-po127
       8: sin

       9: max
            result[0:127] <--> po0-po127
            address[0:1]   <--> po128-po129

    */
    this->ed_mode=ed_mode;
    if(ed_mode=="adder"){ //adder
        double bitWeight = 1.0 / (pow(2.0, 129) - 1.0);
        int po_counter=0;
        PO_vectors.resize(1);
        for(int i=0; i<129; i++){
            PO_vectors[0].push_back("po"+to_string(po_counter++));
            PO_weights.push_back(bitWeight);
            bitWeight*=2;
        }
    }else if(ed_mode=="div"){   //div
        double bitWeight = 1.0 / (pow(2.0, 64) - 1.0);
        int po_counter=0;
        PO_vectors.resize(2);
        for(int i=0; i<64; i++){
            PO_vectors[0].push_back("po"+to_string(po_counter++));
            PO_weights.push_back(bitWeight);
            bitWeight*=2;
        }
        for(int i=0; i<64; i++){
            PO_vectors[1].push_back("po"+to_string(po_counter++));
        }
    }else if(ed_mode=="log2"){   //log2
        double bitWeight = 1.0 / (pow(2.0, 64) - 1.0);
        int po_counter=0;
        PO_vectors.resize(1);
        for(int i=0; i<32; i++){
            PO_vectors[0].push_back("po"+to_string(po_counter++));
            PO_weights.push_back(bitWeight);
            bitWeight*=2;
        }
    }else if(ed_mode=="multiplier"){   //multiplier
        double bitWeight = 1.0 / (pow(2.0, 128) - 1.0);
        int po_counter=0;
        PO_vectors.resize(1);
        for(int i=0; i<128; i++){
            PO_vectors[0].push_back("po"+to_string(po_counter++));
            PO_weights.push_back(bitWeight);
            bitWeight*=2;
        }
    }else if(ed_mode=="sqrt"){   //sqrt
        double bitWeight = 1.0 / (pow(2.0, 64) - 1.0);
        int po_counter=0;
        PO_vectors.resize(1);
        for(int i=0; i<64; i++){
            PO_vectors[0].push_back("po"+to_string(po_counter++));
            PO_weights.push_back(bitWeight);
            bitWeight*=2;
        }
    }else if(ed_mode=="square"){   //square
        double bitWeight = 1.0 / (pow(2.0, 128) - 1.0);
        int po_counter=0;
        PO_vectors.resize(1);
        for(int i=0; i<128; i++){
            PO_vectors[0].push_back("po"+to_string(po_counter++));
            PO_weights.push_back(bitWeight);
            bitWeight*=2;
        }
    }else if(ed_mode=="bar"){
        double bitWeight = 1.0 / (pow(2.0, 128) - 1.0);
        int po_counter=0;
        PO_vectors.resize(1);
        for(int i=0; i<128; i++){
            PO_vectors[0].push_back("po"+to_string(po_counter++));
            PO_weights.push_back(bitWeight);
            bitWeight*=2;
        }
    }else if(ed_mode=="max"){
        double bitWeight = 1.0 / (pow(2.0, 128) - 1.0);
        int po_counter=0;
        PO_vectors.resize(2);
        for(int i=0; i<128; i++){
            PO_vectors[0].push_back("po"+to_string(po_counter++));
            PO_weights.push_back(bitWeight);
            bitWeight*=2;
        }
        PO_vectors[1].push_back("po"+to_string(po_counter++));
        PO_vectors[1].push_back("po"+to_string(po_counter++));
    }else{
        this->ed_mode="";
    }
}
void Global::show_time(){
    auto current = chrono::system_clock::now();
    chrono::seconds elapsed = chrono::duration_cast<chrono::seconds>(current-start_time);
    cout << "\033[101m Current Runtime (s): " << elapsed.count() << "\033[0m" << endl;
}
bool Global::check_if_time_out(){
    auto current = chrono::system_clock::now();
    chrono::seconds elapsed = chrono::duration_cast<chrono::seconds>(current-start_time);
    if(elapsed.count() % 60==0 && elapsed.count()!=runtime){
        show_time();
        runtime=elapsed.count();
    }
    if(elapsed.count()>=time_limit){
        return true;
    }else{
        return false;
    }
}
void Global::setTimeLimit(int limit){
    time_limit=limit;
}
void Global::setOutputFilename(string filename){
    output_filename=filename;
}
Global::Global(){
    start_time = std::chrono::system_clock::now();
    //cout << "Time: 0s" << endl;
    runtime=0;
}