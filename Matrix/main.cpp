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
  

    auto A = Matrix<int>::makeMatrix(4, 4, { 1,2,4,-80,
                                               -5,2,0,-10,
                                                1,10,3,20,
                                                1,2, 2,3 });

    auto B = Matrix<float>::makeMatrix(3, 2, std::vector<float>
        {8.1,5.2,
         6.3,7.444,
         8.5,9.6 });

    auto C = Matrix<double>::makeMatrix(3, 2, std::vector<double>
        {4.1, 5.2,
         9.6, 7.44,
         4.1, 9.6 });

    auto D = Matrix<int>::makeMatrix(3, 2, std::vector<int>
        {1, 2,
         2, 3,
         3, 4});


    cout << det(A);

    auto start = std::chrono::high_resolution_clock::now();
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    //cout << duration.count();
    


    return 0;
}

