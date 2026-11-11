#include "aqueous_engine.h"
#include <exception>
#include <fstream>
#include <iostream>
#include <string>
int main(int argc,char** argv) {
    if(argc<2 || argc>4) {
        std::cerr<<"Usage: aqueous_cli MODEL.txt [RESULT.json] [SPECIES.csv]\n";
        return 2;
    }
    try {
        const auto model=aqueous::readModel(argv[1]);
        const auto result=aqueous::solve(model);
        const auto json=aqueous::toJson(result);
        if(argc>=3) {
            std::ofstream file(argv[2]);
            if(!file) throw std::runtime_error("Cannot write JSON output");
            file<<json;
        } else std::cout<<json;
        if(argc>=4) aqueous::writeCsv(result,argv[3]);
        return 0;
    } catch(const std::exception& error) {
        std::cerr<<"Equilibrium calculation failed: "<<error.what()<<"\n";
        return 1;
    }
}
