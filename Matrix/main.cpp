#include "matrix.h"
#include "iomatrix.h"
#include <sstream>




using namespace std;

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
int main() {
  
    std::vector<int> dataA = 
    {1,2,3,4,5,6,7,8,9,
     9,7,8,6,5,4,3,2,1};

    std::vector<float> dataB = 
    {8.1,5.2,
     6.3,7.444,
     8.5,9.6};

    std::vector<double> dataD =
    {4.1,5.2,
     9.6,7.44,
     4.1,9.6 };

    std::vector<int> dataE =
    {1, 2,
     2, 3,
     3, 4};

    int a[] = { 1,2,3,4,5,6,7,8,9,
               9,8,7,6,5,4,3,2,1 };

    auto A = Matrix<int>::makeMatrix(2, 9, a);

    auto B = Matrix<float>::makeMatrix(3, 2, move(dataB));
    auto C = Matrix<double>::makeMatrix(2, 2);
    auto D = Matrix<double>::makeMatrix(3, 2, move(dataD));
    auto E = Matrix<int>::makeMatrix(3, 2, move(dataE));
    auto F = Matrix<double>::makeMatrix(3, 1, 1.0);
    auto G = Matrix<float>::makeLinSpace(1, 5.0, 5);
    auto H = Matrix<double>::makeLinInc(-1, -0.1, -2);
    auto I = Matrix<double>::makeMatrix(2, 2, { 1,2,3,4 });

    cout << I;


   
    return 0;
}

