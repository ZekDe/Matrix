#include "matrix.h"
#include "iomatrix.h"




using namespace std;


int main() {
  
    float dataA[3*2] = {1,2,3,
                        3,2,1};

    std::vector<float> dataB = 
    {4,5,
     6,7,
     8,9};

    auto A = Matrix<float>::makeMatrix(2, 3, dataA);
    auto B = Matrix<float>::makeMatrix(3, 2, dataB);
    auto C = A * B;
    std::cout << C;
    

   
  

    return 0;
}

