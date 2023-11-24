#include "matrix.h"
#include "iomatrix.h"




using namespace std;


int main() {
  
    std::vector<int> dataA = 
    {1,2,3,
     3,2,1};

    std::vector<float> dataB = 
    {4.1,5.2,
     6.3,7.4,
     8.5,9.6};

    auto A = Matrix<int>::makeMatrix(2, 3, dataA);
    auto B = Matrix<float>::makeMatrix(3, 2, dataB);
    auto C = Matrix<float>::makeMatrix(2, 2);
    auto iB = Matrix<int>(B);
    C =  A * B * 2.2F;
    C *= 2;




    for (int i = 0; i < A.size().first; ++i)
    {
        for (int j = 0; j < A.size().second; ++j)
            std::cout << A(i, j) << " ";
        std::cout << "\n";
    }
        
    
    


    auto D = Matrix<float>::makeLinSpace(1.0, 10.0, 20);

    //std::cout << D;
   


    return 0;
}

