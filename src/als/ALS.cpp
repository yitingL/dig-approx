#include <bits/stdc++.h>
#include "../header.h"
using namespace std;

mt19937_64 rand_generator;

int main(int argc, char* argv[]){
    string orig_blif="";       //circuit before optimization
    string output_blif="";     //circuit after optimization
    string golden_blif="";     //circuit with error rate=0
    string ed_mode="";
    double ERROR_RATE_THRES=0.0;
    int timelimit=3600;
    int opt;
    while((opt=getopt(argc,argv,"t:g:c:o:e:d:"))!=-1){
        switch(opt){
            case 't':
                timelimit=stoi(optarg);
                break;
            case 'g':
                golden_blif=optarg;
                break;
            case 'c':
                orig_blif=optarg;
                break;
            case 'o':
                output_blif=optarg;
                break;
            case 'e':
                ERROR_RATE_THRES=stod(optarg)/100.0;
                break;
            case 'd':
                ed_mode=optarg;
                break;
            default:
                break;
        }
    }
    if(golden_blif==""|| output_blif=="" || ERROR_RATE_THRES==0.0 ){
        cerr << "Usage: ./main ..." << endl;
        cerr << "\t-g <golden_circuit>" << endl;
        cerr << "\t-o <output_circuit>" << endl;
        cerr << "\t-e <error_rate in %>" << endl;        
        cerr << "\t[-t <timeout> = 3600]" << endl;
        cerr << "\t[-c <original_circuit>]" << endl;
        cerr << "\t[-d <ed_mode>]" << endl;
        exit(0);
    }
    if(orig_blif=="") orig_blif=golden_blif;
    cout << "golden: " << golden_blif << endl;
    cout << "output: " << output_blif << endl;
    cout << "error rate: " << ERROR_RATE_THRES << endl;
    cout << "timeout: " << timelimit << endl;
    cout << "ed mode: " << ed_mode << endl;
    // if(argc==5){
    //     golden_blif=argv[1];
    //     orig_blif=golden_blif;
    //     output_blif=argv[2];
    //     ERROR_RATE_THRES=stod(argv[3])/100.0;
    // }else if (argc==6){
    //     golden_blif=argv[1];
    //     orig_blif=argv[2];
    //     output_blif=argv[3];
    //     ERROR_RATE_THRES=stod(argv[4])/100.0;
    // }else if(argc==7){
    //     golden_blif=argv[1];
    //     orig_blif=golden_blif;
    //     output_blif=argv[2];
    //     ERROR_RATE_THRES=stod(argv[3])/100.0;
    //     assert(string(argv[4])=="ED");
    //     ed_mode=argv[5];
    // }else{
    //     cout << "Usage: ./ga <timeout> <golden_circuit> <output_circuit> <error_rate>" << endl;
    //     cout << "Usage: ./ga <timeout> <golden_circuit> <original_circuit> <output_circuit> <error_rate>" << endl;
    //     cout << "Usage: ./ga <timeout> <golden_circuit> <output_circuit> <error_rate> <ed_mode>" << endl;
    //     return 0;
    // }
    
    // Random seed generation
    auto seed=time(0);
    rand_generator.seed(seed);
    cout << "seed: " << seed << endl;
    ofstream fout("../temp/seed.txt");
    fout << seed << endl;

    //-----------------
    // read golden circuit (dont-touch)
    //-----------------
    Circuit golden_crt;
    golden_crt.read_dig(golden_blif);
    build_golden_output_vals_for_error_eval(&golden_crt);
    Global global;
    global.setup_PO_vec_and_PO_weights(ed_mode);
    global.setOutputFilename(output_blif);
    global.setTimeLimit(timelimit);
    cout << "----- golden circuit -----" << endl;
    cout << "#PIs: " << golden_crt.PIs.size() << endl;
    cout << "#POs: " << golden_crt.POs.size() << endl;
    cout << "#Nodes: " << golden_crt.getSize() << endl;
    //cout << scientific << setprecision(5);

    Circuit* orig_crt=new Circuit();
    orig_crt->read_dig(orig_blif);
    orig_crt->write_dig(output_blif);
    preprocess_and_postprocess(global);

    double target_error_rate=ERROR_RATE_THRES;
    int consecutive_gen=50;
    int round=0;
    bool goGA=true;
    bool goNTC=true;
    while(true){
        //-------------------
        // node to constant
        //-------------------
        Circuit* crt;
        if(goNTC){
            double node_to_constant_error_rate;
            double alpha=1.0-target_error_rate*0.1;
            int counter=1;
            double node_to_constant_error_rate_thres;
            node_to_constant_error_rate_thres=target_error_rate;
            do{
                crt=new Circuit();
                crt->read_dig(output_blif);
                nodes_to_constant(crt,alpha);
                if(ed_mode==""){
                    node_to_constant_error_rate=check_error_rate(&golden_crt, crt);
                }else{
                    node_to_constant_error_rate=check_error_distance(&golden_crt, crt, global);
                }
                //cout << "Node To Constant with Error Rate Constraint: " << node_to_constant_error_rate_thres << endl;
                //cout << "after node to constant with alpha = " << alpha << endl;
                //cout << "golden circuit's size: " << golden_crt.getSize() << endl;
                //cout << "modified circuit's size: " << crt->getSize() << endl;
                //cout << "error rate: " << node_to_constant_error_rate*100 << "%"<< endl;
                // if(counter==2){
                //     //alpha=1.0;
                //     alpha=alpha*0.1+0.9;
                // }else{
                //     alpha=alpha*0.1+0.9;
                // }
            }while(--counter && node_to_constant_error_rate>node_to_constant_error_rate_thres);
            
            if(node_to_constant_error_rate<=node_to_constant_error_rate_thres){
                crt->write_dig(output_blif);
                //cout << "# output size: " << crt->getSize() << endl;
                //cout << "# error rate: " << check_error_rate(&golden_crt, crt) << endl;
            }
        }

        //-----------------
        // genetic algorithm
        //-----------------
        if(goGA){
            crt=new Circuit();
            crt->read_dig(output_blif);
            GA_Engine ga(global, crt, &golden_crt, output_blif, target_error_rate);
            int ga_return_signal=ga.start(consecutive_gen,++round);
            orig_blif=output_blif;
            if(ga_return_signal==0){            // give up some nodes, and repeat the process again.
                goGA=true;
                goNTC=false;
            }else if(ga_return_signal==1){      // go to the find round of NTC
                goGA=false;
                goNTC=true;
            }else if(ga_return_signal==2){      // terminate the process
                break;
            }
        }else{
            break;
        }
    }
    preprocess_and_postprocess(global);
    Circuit* crt=new Circuit();
    crt->read_dig(output_blif);
    cout << "----- final circuit -----" << endl;
    cout << "### Area: " << crt->getSize() << endl;
    cout << "### Depth: " << crt->getLevel() << endl;
    if(ed_mode=="")
        cout << "### Error Rate: " << check_error_rate(&golden_crt, crt) << endl;
    else
        cout << "### Error Distance: " << check_error_distance(&golden_crt, crt, global) << endl;

    return 0;
}