#include "matrix.h"
#include "iomatrix.h"
#include <sstream>


using namespace std;


int main() {
  
    std::vector<int> dataA = 
    {1,2,3,
     3,2,1};

    std::vector<double> dataB = 
    {4.1,5.2,
     6.3,7.44,
     8.5,9.6};

    auto A = Matrix<int>::makeMatrix(2, 3, dataA);
    auto B = Matrix<double>::makeMatrix(3, 2, dataB);
    auto C = Matrix<double>::makeMatrix(2, 2);
    auto iB = Matrix<int>(B);
    C =  A * B * 2;
    C *= 2;

    cout << B << C;

    /*auto D = Matrix<double>::makeMatrix(2,2);
    cin >> D;*/

    //cout << D << C;

   
        
    
   


    return 0;
}

