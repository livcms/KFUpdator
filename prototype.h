#ifndef PROTOTYPE
#define PROTOTYPE

#include <vector> 



struct output{ 
  std::vector<double> finalStateVector = std::vector<double>(5);
  std::vector<double> finalStateError = std::vector<double>(15);
};


output testfunction(std::vector<double> xin, std::vector<double> Cin, std::vector<double> rin, std::vector<double> Vin, std::vector<double> VMeasin );
#endif 
