#include "matrix.h"
#include "iomatrix.h"


int main() {
  
    float dataA[2*3] = {1,2,3,
                        3,2,1,};

    std::vector<float> dataB
    {4,5,
     6,7,
     8,9};

    Matrix<float> A(2, 3, dataA);
    Matrix<float> B(3, 2, dataB);
    Matrix<float> C{A * B};

    std::cout << C;
    
    std::vector<float> dataD
    { 1,2,2};

    std::vector<float> dataE
    { 1,2,
      1,3,
      1,4};



    Matrix<float> D(1, 3, dataD);
    Matrix<float> E(3, 2, dataE);
    Matrix<float> F = D*E;


  

    return 0;
}

