#include "aqueous_sweep.h"
#include <exception>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
int main(int argc,char** argv) {
    if(argc<7 || argc>8) {
        std::cerr<<"Usage:\n"
                 <<"  aqueous_sweep MODEL temperature START_K STOP_K STEPS OUTPUT.csv\n"
                 <<"  aqueous_sweep MODEL total COMPONENT START STOP STEPS OUTPUT.csv\n";
        return 2;
    }
    try {
        const std::string parameter=argv[2];
        if((parameter=="temperature" && argc!=7) ||
           (parameter=="total" && argc!=8))
            throw std::invalid_argument("Wrong number of parameters for sweep mode");
        const std::string component=parameter=="total"?argv[3]:"";
        const int offset=parameter=="total"?1:0;
        std::size_t a=0,b=0,n=0;
        const std::string startStr=argv[3+offset],stopStr=argv[4+offset],stepsStr=argv[5+offset];
        const double start=std::stod(startStr,&a),stop=std::stod(stopStr,&b);
        const int steps=std::stoi(stepsStr,&n);
        if(a!=startStr.size()||b!=stopStr.size()||n!=stepsStr.size())
            throw std::invalid_argument("Invalid numerical sweep parameter");
        const auto model=aqueous::readModel(argv[1]);
        const auto series=aqueous::runSweep(model,parameter,component,start,stop,steps);
        const std::string csv=aqueous::sweepCsv(series,parameter);
        std::ofstream out(argv[6+offset]);
        if(!out) throw std::runtime_error("Cannot write sweep output");
        out<<csv;
        if(!out) throw std::runtime_error("Sweep output failed to write");
        std::cout<<"Wrote "<<series.size()<<" equilibrium sweep points\n";
        return 0;
    }catch(const std::exception& err) {
        std::cerr<<"Sweep failed: "<<err.what()<<"\n";
        return 1;
    }
}
