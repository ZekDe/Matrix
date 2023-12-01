#include "matrix.h"
#include "iomatrix.h"
#include <sstream>




using namespace std;

template<typename T>
void printMatrix(const Matrix<T> &A)
{
    auto [rows, cols] = A.size();
    for (int i{}; i < rows; ++i){
        for (int j{}; j < cols; ++j)
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

    auto A = Matrix<int>::makeMatrix(2, 9, move(dataA));

    auto B = Matrix<float>::makeMatrix(3, 2, move(dataB));
    auto C = Matrix<double>::makeMatrix(2, 2);
    auto D = Matrix<double>::makeMatrix(3, 2, move(dataD));
    auto E = Matrix<int>::makeMatrix(3, 2, move(dataE));
    auto F = Matrix<double>::makeMatrix(3, 1, 1.0);

    auto K = Matrix<float>::makeLinSpace(1, PI, 1);
    auto L = Matrix<double>::makeLinSpace(1, PI, 20);
    auto M = Matrix<int>::makeLinSpace(1, PI, 20);
    auto N = Matrix<int>::makeMatrix(1, 3, 4);
    auto O = Matrix<int>::makeMatrix(1, 3, 5);
  


    return 0;
}

