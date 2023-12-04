#pragma once
#include "matrix.h"

namespace MathLab
{
    template<typename T>
    void swapRows(Matrix<T>& A, int row1, int row2)
    {
        auto [_, cols] = A.size();
        for (int j = 0; j < cols; ++j)
            std::swap(A(row1, j), A(row2, j));
    }
  
    template<typename T>
    inline double det(const Matrix<T>& A)
    {
        // LUP decomposition and pivot
        auto C{Matrix<double>(A)};
        auto [rows, _] = A.size();
        std::vector<int> P(A.size().first);

        
        for (int i = 0; i < rows; ++i) {
            P[i] = i;

            // Pivot
            int pivot_row = i;

            for (int j = i + 1; j < rows; ++j) 
                if (std::abs(C(j, i)) > std::abs(C(pivot_row, i)))
                    pivot_row = j;

 
            if (pivot_row != i) {
                swapRows(C, i, pivot_row);
                std::swap(P[i], P[pivot_row]);
            }

            // decompose
            for (int j = i + 1; j < rows; ++j) {
                C(j, i) /= C(i, i);
                for (int k = i + 1; k < rows; ++k) 
                    C(j, k) -= C(j, i) * C(i, k);
            }
        }

        // calculate determinant
        double det = 1.0;

        for (int i = 0; i < rows; i++)
            det *= C(i, i);

        int sign = -1;
        for (int i = 0; i < rows; i++) {
            if (P[i] != i)
                sign = -sign; 
        }

        return det * sign;
    }



}