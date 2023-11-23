#pragma once
#include <memory>
#include <vector>


template <typename T>
class Matrix 
{
public:
    
    static Matrix makeMatrix(size_t rows, size_t cols, std::vector<T> tvec)
    {
        return Matrix(rows, cols, tvec);
    }
    
    static Matrix makeMatrix(size_t rows, size_t cols, const T* data)
    {
        return Matrix(rows, cols, data);
    }

    static Matrix makeMatrix(size_t rows, size_t cols, int val = 0)
    {
        std::vector<T> data(rows * cols, val);
        return Matrix(rows, cols, data);
    }
 
    
    std::pair<size_t, size_t> size() const;

    template <typename U>
    friend Matrix<U> operator*(const Matrix<U>& A, const Matrix<U>& B);
  
    template <typename U>
    friend std::ostream& operator<<(std::ostream& os, const Matrix<U>& A);

private:
    Matrix(size_t rows, size_t cols, std::vector<T> tvec);
    /* This ctor must be used carefully, in case of mismatched rows, cols
    * and data is undefined-behaviour  */
    Matrix(size_t rows, size_t cols, const T* data);

    static void align(Matrix &);


    size_t m_rows, m_cols;
    std::vector<T> m_data;
};




template <typename T>
Matrix<T>::Matrix(size_t rows, size_t cols, std::vector<T> tvec) :
    m_rows(rows), m_cols(cols), m_data(tvec)
{
    if (tvec.empty() || rows * cols != tvec.size()) {
        m_rows = m_cols = 0;
        return;
    }
    
    align(*this);
}


template <typename T>
Matrix<T>::Matrix(size_t rows, size_t cols, const T *data) : 
Matrix(rows, cols, std::vector<T>(data, data + rows * cols))
{} 


template <typename T>
std::pair<size_t, size_t> Matrix<T>::size() const {
    return { m_rows, m_cols };
}

template <typename U>
Matrix<U> operator*(const Matrix<U>& A, const Matrix<U>& B)
{
    size_t m, i, inner, n, j, coffset, boffset, b_i, k, aoffset;
    U temp;

    if (A.m_cols != B.m_rows)
    {
        return Matrix<U>::makeMatrix(0, 0);
    }

    auto C = Matrix<U>::makeMatrix(A.m_rows, B.m_cols);

    m = A.m_rows;

    if (A.m_cols == 1 || B.m_rows == 1)
    {
        for (i = 0; i < m; i++)
        {
            inner = B.m_cols;
            for (n = 0; n < inner; n++)
            {
                C.m_data[i + C.m_rows * n] = 0.0F;
                j = A.m_cols;
                for (coffset = 0; coffset < j; coffset++)
                {
                    C.m_data[i + C.m_rows * n] += A.m_data[i + A.m_rows * coffset] *
                        B.m_data[coffset + B.m_rows * n];
                }
            }
        }
    }
    else {

        inner = A.m_cols;
        n = B.m_cols;

        for (j = 0; j < n; j++)
        {
            coffset = j * m;
            boffset = j * inner;
            for (b_i = 0; b_i < m; b_i++)
                C.m_data[coffset + b_i] = 0.0F;

            for (k = 0; k < inner; k++)
            {
                aoffset = k * m;
                temp = B.m_data[boffset + k];
                for (b_i = 0; b_i < m; b_i++)
                {
                    i = coffset + b_i;
                    C.m_data[i] += temp * A.m_data[aoffset + b_i];
                }
            }
        }
    }
    return C;
}

//template <typename T>
//Matrix<T> Matrix<T>::operator*(const Matrix<T>& B)
//{
//
//    size_t m,i,inner,n,j,coffset,boffset,b_i,k,aoffset;
//    T temp;
//
//    if (m_cols != B.m_rows) 
//    {
//        return Matrix<T>(0, 0, std::vector<T>());
//    }
//
//    Matrix<T> C(m_rows, B.m_cols);
//
//    m = m_rows;
//
//    if (m_cols == 1 || B.m_rows == 1)
//    {
//        for (i = 0; i < m; i++)
//        {
//            inner = B.m_cols;
//            for (n = 0; n < inner; n++)
//            {
//                C.m_data[i + C.m_rows * n] = 0.0F;
//                j = m_cols;
//                for (coffset = 0; coffset < j; coffset++)
//                {
//                    C.m_data[i + C.m_rows * n] += m_data[i + m_rows * coffset] *
//                        B.m_data[coffset + B.m_rows * n];
//                }
//            }
//        }
//    }
//    else {
//
//        inner = m_cols;
//        n = B.m_cols;
//
//        for (j = 0; j < n; j++)
//        {
//            coffset = j * m;
//            boffset = j * inner;
//            for (b_i = 0; b_i < m; b_i++)
//                C.m_data[coffset + b_i] = 0.0F;
//
//            for (k = 0; k < inner; k++)
//            {
//                aoffset = k * m;
//                temp = B.m_data[boffset + k];
//                for (b_i = 0; b_i < m; b_i++)
//                {
//                    i = coffset + b_i;
//                    C.m_data[i] += temp * m_data[aoffset + b_i];
//                }
//            }
//        }
//    }
//    return C;
//}


template <typename T>
void Matrix<T>::align(Matrix &A)
{
    size_t size = A.m_rows * A.m_cols;
    auto tvec{std::vector<T>(size)};
    
    tvec = A.m_data;

    for (size_t i{}, k{}; i < A.m_rows; ++i)
        for (size_t j{}; j < A.m_cols; ++j)
            A.m_data[i + A.m_rows * j] = tvec[k++];
            
}







