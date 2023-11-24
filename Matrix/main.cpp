#include "matrix.h"
#include "iomatrix.h"




using namespace std;


int main() {
  
    int dataA[3*2] = {1,2,3,
                      3,2,1};

    std::vector<float> dataB = 
    {4,5,
     6,7,
     8,9};

    auto A = Matrix<int>::makeMatrix(2, 3, dataA);
    auto B = Matrix<float>::makeMatrix(3, 2, dataB);
    auto C = Matrix<float>::makeMatrix(2, 2);
   
    C =  2.2 * A * B;
    C *= 2;
    std::cout << C;

    auto D = Matrix<float>::makeLinSpace(1.0, 10.0, 15);

    std::cout << D;
   


    return 0;
}

