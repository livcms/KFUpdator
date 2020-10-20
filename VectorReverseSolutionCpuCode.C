/**
 * Document: MaxCompiler Tutorial (maxcompiler-tutorial.pdf)
 * Chapter: 6      Exercise Solution: 1      Name: Vector reverse solution
 * MaxFile name: VectorReverseSolution
 * Summary:
 * 	   Reverse a vector.
 */
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <MaxSLiCInterface.h>
#include "VectorReverseSolution.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>


void printvec(std::vector<double> invec, std::string name){
	for( auto i = invec.begin(); i != invec.end(); ++i){
		std::cout<<name<<" at  is "<<*i<<std::endl;
	}
}


//Old, change
int check(int size, double *dataOut, double *expected)
{
	int status = 0;
	for (int i = 0; i < size; i++) {
		if (dataOut[i] != expected[i]) {
			fprintf(
				stderr, "Output data @ %f = %f (expected %f)\n",
				i, dataOut[i], expected[i]);static int vectorSize3 = 3;
			status = 1;
		}
	}
	return status;
}

//old, change to a proper CPU solution
void VectorReverseSolutionCPU(int streamSize, double *dataIn, double *dataOut)
{
	const int vectorSize = VectorReverseSolution_vectorSize;
	for (int i = 0; i < streamSize; i++) {
		for (int j = 0; j < vectorSize; j++) {
			dataOut[i * vectorSize + j] = dataIn[i * vectorSize + (vectorSize - 1 - j)];
		}
	}
}

template <class T>
char* as_bytes(T& x) {
    return &reinterpret_cast<char&>(x);
    // or:
    // return reinterpret_cast<char*>(std::addressof(x));
}


int main()
{
	//Reading in elements into vectors. There's probably a much better way of doing this...
	std::vector<double>x;
	std::vector<double>C;
	std::vector<double>r;
	std::vector<double>rMeas;
	std::vector<double>V;
	std::vector<double>VMeas;
	std::ifstream ifs{"/home/maxeler/workspace/KFPrototype-CPU/src/testprint.dat", std::ios_base::binary};
	if (!ifs) std::cout<<"can't open file"<<std::endl;
	int count = 0;
	int xcount = 0;
	int overallcount = 0;
	for(double in; ifs.read(as_bytes(in), sizeof(double));){
		overallcount++;
		if(count <5){
			x.push_back(in);
			count++;
			xcount++;
			//std::cout<<" I put "<<in<<" in x "<<std::endl;
			} else if(count>4 && count <20){
				C.push_back(in);
				count++;
			} else if(count>19 && count<22){
				r.push_back(in);
				count++;
			} else if(count>21 && count<24){
				rMeas.push_back(in);
				count++;
			} else if(count>23  && count<27){
				V.push_back(in);
				count++;
			} else if(count>26 && count <30){
				VMeas.push_back(in);
				count++;
			// there are 30 elements in each vector, so after 30 it repeats
			} else if(count==30){
				count = 0;
			}
	}



	std::vector<double>testIn;

	//This needs changing - do more performance modelling
	const int streamSize = 16;
	const int vectorSize = VectorReverseSolution_vectorSize;
	size_t sizeBytes = streamSize * vectorSize * sizeof(double);
	double *dataIn = (double *)malloc(sizeBytes);
	double *expected = (double *)malloc(sizeBytes);
	for (int i = 0; i < (streamSize * vectorSize); i++) {
		dataIn[i] = i*1.0;
		testIn.push_back(i*1.0);
	}

	double *dataOut = (double *)malloc(sizeBytes/2*25);
	for ( int i =0; i < streamSize * 25; i++) {
		dataOut[i] = 0*1.0;
	}

	VectorReverseSolutionCPU(streamSize, dataIn, expected);
	//VectorReverseSolutionCPU(streamSize, r.data(), expected);

	printf("Running DFE.\n");
	//The order of inputs has been jumbled, use advanced interface for slac
	VectorReverseSolution(streamSize, C.data(), V.data(), VMeas.data(), r.data(), x.data(),  dataOut);

	for( int i = 0; i<streamSize*25/2; i++){
		std::cout<<"data out at "<<i<<" is "<<dataOut[i]<<std::endl;
	}

	//printvec(VMeas, "VMeas");
	//int status = check(streamSize * vectorSize, dataOut, expected);
	//if (status)
	//	printf("Test failed.\n");
	//else
	//	printf("Test passed OK!\n");

	return 0;
}

