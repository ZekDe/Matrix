#pragma once
#include <iostream>
#include "matrix.h"

template <typename U>
std::ostream& operator<<(std::ostream& os, const Matrix<U>& A) // inserter 
{
    os << A.m_rows << "x" << A.m_cols << "\n";
    for (size_t i{}; i < A.m_rows; i++)
    {
        for (size_t j{}; j < A.m_cols; j++)
            os << A.m_data[i + A.m_rows * j] << " ";

        os << "\n";
    }

    return os;
}


