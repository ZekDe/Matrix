#pragma once
#include "matrix.h"

namespace MathLab
{
    constexpr double tol = 1E-7;
    constexpr size_t max_iter = 100;

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
        auto C{ Matrix<double>(A) };
        auto [rows, _] = A.size();
        std::vector<int> P(A.size().first);


        for (int i = 0; i < rows; ++i) {
            P[i] = i;

            // Pivot
            int pivot_row = i;

            double tmp = 0.0;
            for (int j = i + 1; j < rows; ++j)
                if (double x = std::abs(C(j, i)) > tmp) {
                    tmp = x;
                    pivot_row = j;
                }


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

    template<typename T>
    inline Matrix<double> inv(const Matrix<T>& A)
    {
        auto [rows, cols] = A.size();

        if (rows != cols)
            return Matrix<double>::makeMatrix(0, 0); //todo: exception runtime_error("Square Matrix needed.");

        auto I{ Matrix<double>::makeEyeMatrix(rows, cols) };
        auto C{ Matrix<double>(A) };

        // Gauss-Jordan 
        for (int i{}; i < rows; ++i) {

            if (C(i, i) == 0.0)
                return Matrix<double>::makeMatrix(0, 0); // todo:exception throw runtime_error("Zero pivot.");

            // make diagonal element to 1
            double pivot = C(i, i);
            for (int j = 0; j < cols; ++j) {
                C(i, j) /= pivot;
                I(i, j) /= pivot;
            }

            // make 0 the elements in other rows
            for (int k{}; k < rows; ++k) {
                if (k != i) {
                    double factor = C(k, i);
                    for (int j = 0; j < cols; ++j) {
                        C(k, j) -= factor * C(i, j);
                        I(k, j) -= factor * I(i, j);
                    }
                }
            }
        }

        return I;

    }

    template<typename T>
    inline Matrix<T> transpose(const Matrix<T>& A)
    {
        auto [rows, cols] = A.size();
        auto C{ A };

        for (size_t i{}; i < rows; ++i)
            for (size_t j{}; j < cols; ++j)
                C(j, i) = A(i, j);

        return C;
    }


}