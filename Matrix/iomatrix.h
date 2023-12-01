#pragma once
#include <iostream>
#include "matrix.h"
#include <string>
#include <sstream>

template <typename U>
std::ostream& operator<<(std::ostream& os, const Matrix<U>& A) // inserter 
{
    auto [rows, cols] = A.size();

    if (!rows || !cols) {
        os << rows << "x" << cols << "\n" ;
        return os;
    }
        

    int line_length_max = 0;
    for (int i{}; i < rows; ++i){

        int line_length = 0;

        for (size_t j{}; j < cols; ++j) {
            std::stringstream ss;
            ss << A(i, j);
            line_length += ss.str().length();
        }
        line_length_max = std::max(line_length_max, line_length);
        line_length = 0;
    }
    // get max length of line
    line_length_max += (cols - 1);

    std::stringstream ss;
    ss << rows << cols;
    int rows_x_cols_length = ss.str().length() + 1;
    rows_x_cols_length = std::abs(line_length_max - rows_x_cols_length);
    
    os << std::string(line_length_max, '-') << '\n';

    os  << std::string(rows_x_cols_length / 2, ' ')
        << rows << "x" << cols 
        << std::string(rows_x_cols_length / 2, ' ')
        << '\n';

    os << std::string(line_length_max, '-') << '\n';
   

    for (size_t i{}; i < rows; i++){
        for (size_t j{}; j < cols; j++)
            os << A(i,j) << " ";

        os << "\n";
    }

    os << std::string(line_length_max, '-') << '\n';

    return os;
}

template <typename U>
std::istream& operator>>(std::istream& is, Matrix<U>& A) // extractor
{
    auto [rows, cols] = A.size();

    for (int i{}; i < rows; ++i)
        for (int j{}; j < cols; ++j)
            is >> A(i, j);

    return is;
}


