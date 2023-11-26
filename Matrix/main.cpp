#include "matrix.h"
#include "iomatrix.h"
#include <sstream>


using namespace std;


int main() {
  
    std::vector<int> dataA = 
    {1,2,3,
     3,2,1};

    std::vector<float> dataB = 
    {8.1,5.2,
     6.3,7.44,
     8.5,9.6};

    std::vector<double> dataD =
    {4.100001,5.2,
     6.3,7.44,
     8.5,9.6 };

    std::vector<int> dataE =
    {1, 2,
     2, 3,
     3, 4};

    auto A = Matrix<int>::makeMatrix(2, 3, dataA);
    auto B = Matrix<float>::makeMatrix(3, 2, dataB);
    auto C = Matrix<double>::makeMatrix(2, 2);
    auto D = Matrix<double>::makeMatrix(3, 2, dataD);
    auto E = Matrix<int>::makeMatrix(3, 2, dataE);
    auto F = Matrix<double>::makeMatrix(1, 3, 1.0);


    auto I = Matrix<double>::makeMatrix(1, 4, std::vector<double>{1,2,3,4});
    auto J = Matrix<double>::makeMatrix(4, 1, std::vector<double>{
                                                                  1, 
                                                                  2, 
                                                                  3,
                                                                  4});
    //auto H = B * D;

    cout << (B / A);
    
    return 0;
}

