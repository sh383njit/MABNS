
//fist do: $export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/salman/extern/bin/glnxa64/:matlabroot/sys/os/glnxa64
//then compile: $g++ -o main -std=c++11 -I /home/salman/extern/include/ -L /home/salman/extern/bin/glnxa64/ -pthread matlabtest.cpp -lMatlabDataArray -lMatlabEngine
//the matlabroot is /home/salman
//more info: https://www.mathworks.com/help/matlab/matlab_external/call-matlab-functions-from-c-1.html, https://www.mathworks.com/help/matlab/matlab_external/build-c-engine-programs.html

#include "MatlabEngine.hpp"
#include "MatlabDataArray.hpp"
#include <iostream>

using namespace std;

void callFevalgcd() {

	float x,y;

    // Pass vector containing MATLAB data array scalar
    using namespace matlab::engine;
    // Start MATLAB engine synchronously
    std::unique_ptr<MATLABEngine> matlabPtr = startMATLAB();
    // Create MATLAB data array factory
    matlab::data::ArrayFactory factory;

    // Pass vector containing 2 scalar args in vector    
   // std::vector<matlab::data::Array> args({factory.createScalar<int16_t>(30),factory.createScalar<int16_t>(56) });

	y=0.25;

	std::vector<matlab::data::Array> args({factory.createScalar<float>(1),factory.createScalar<float>(y)});

    // Call MATLAB function and return result
	
	for(int i=0;i<100;i++) {
    	matlab::data::TypedArray<float> result = matlabPtr->feval(u"ncx2rnd", args);
    	std::cout << "Result: " << result[0] << std::endl;
	}
}


int main() {

	callFevalgcd();

}
