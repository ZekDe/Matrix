#pragma once
#include "matrix.h"

//namespace MathLab
//{
//    template <typename T>
//    class MatrixOperations 
//    {
//    public:
//        static Matrix<T> det(Matrix<T> A)
//        {
//           auto C =  Matrix<T>::makeMatrix(2, 2, 2);
//           C.m_data[0] = 1;
//           return C;
//        }
//
//    };
//    
//
//}

namespace MathLab
{
    template <typename T>
    inline Matrix<T> det(Matrix<T> A)
    {
        auto C = Matrix<T>::makeMatrix(2, 2, 2);
        C(0, 0) = 1;
        //C.m_data[0] = 1;
        return C;
    }
}