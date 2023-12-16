#include "matrix.h"
#include "matrixoperations.h"
#include "iomatrix.h"
#include <sstream>
#include <chrono>
#include <typeinfo>


using namespace std;
using namespace MathLab;

template<typename T>
void printMatrix(const Matrix<T> &A)
{
    auto [rows, cols] = A.size();
    for (size_t i{}; i < rows; ++i){
        for (size_t j{}; j < cols; ++j)
            cout << A(i, j) << " ";
        cout << "\n";
    }
}; 



#define PI 3.1415926535

int main() 
{
    
    auto B = Matrix<int>::makeMatrix(4, 4, { 1,2,4,-80,
                                               -5,2,1,-10,
                                                1,10,3,20,
                                                1,2, 2,3 });

    auto C = Matrix<float>::makeMatrix(3, 3, std::vector<float>
        {1,2,3,4,5,6,7,8,9});

    auto D = Matrix<>::makeMatrix(3, 3, {1,2,3, 3,4,5, 6,7,8});
    auto J = Matrix<>::makeRandomMatrix(4, 3,0,20);

    auto I = fMatrix::makeEyeMatrix(4);
    cout << I;

    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    //cout << duration.count();

  



    



    return 0;
}

