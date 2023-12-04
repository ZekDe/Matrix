#include "matrix.h"
#include "matrixoperations.h"
#include "iomatrix.h"
#include <sstream>
#include <chrono>



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
    float arr[] = { 1, 2, 4, -80,
                    -5, 2, 1, -10,
                     1, 10, 3, 20,
                     1, 2, 2, 3 };

    auto A = Matrix<float>::makeMatrix(4, 4, arr);
  

    auto B = Matrix<int>::makeMatrix(4, 4, { 1,2,4,-80,
                                               -5,2,1,-10,
                                                1,10,3,20,
                                                1,2, 2,3 });

    auto C = Matrix<float>::makeMatrix(3, 2, std::vector<float>
        {8.1,5.2,
         6.3,7.444,
         8.5,9.6 });


    auto D = Matrix<float>::makeRandomMatrix(1, 3);
    auto E = Matrix<float>::makeRandomMatrix(3, 1);
  
    auto[row, col] = A(10);
    cout << "A(" << row << "," << col << ")";

    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    //cout << duration.count();
    


    return 0;
}

