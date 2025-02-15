#include <bits/stdc++.h>
#include "../header.h"
using namespace std;

void preprocess_and_postprocess(Global& global){
    //./ext/abc -c "source ext/abc.rc; read_blif $1; strash; resyn2; write_verilog temp/aig.v; quit;"
    string filename = global.output_filename;
    string command;
    command="./ext/abc -c \"source ./ext/abc.rc; read_blif "+filename+"; strash; resyn2; write_verilog temp/aig.v; quit;\" > /dev/null 2>&1";
    if(system(command.c_str())){
        cerr << "Error in preprocess_and_postprocess (1)" << endl;
        exit(0);
    }
    //./ext/tig dig temp/aig.v temp/test1.blif
    command="./ext/tig dig temp/aig.v temp/test1.blif > /dev/null 2>&1";
    if(system(command.c_str())){
        cerr << "Error in preprocess_and_postprocess (2)" << endl;
        exit(0);
    }
    Circuit c1, c2;
    c1.read_dig(filename);
    c2.read_dig("temp/test1.blif");
    if(c1.getSize() <= c2.getSize()){
        //c1.write_dig(blif2);
    }else{
        c2.write_dig(filename);
    }
}