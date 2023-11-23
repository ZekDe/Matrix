#pragma once
#include <memory>
#include <vector>


template <typename T>
class Matrix 
{

public:
    Matrix(size_t rows, size_t cols, const T *data);
    Matrix(size_t rows, size_t cols, std::vector<T> tvec);
    Matrix(size_t rows, size_t cols);
    
    std::pair<size_t, size_t> size() const;

    //inserter
    template <typename U>
    friend std::ostream& operator<<(std::ostream& os, const Matrix<U>& A);

    Matrix<T> operator*(const Matrix<T>& B);

private:
    void align();


    const size_t m_rows, m_cols;
    std::vector<T> m_data;
};




template <typename T>
Matrix<T>::Matrix(size_t rows, size_t cols, const T *data) : m_rows(rows), m_cols(cols),
m_data(rows * cols)
{
    if (!data)
        return;
 
    std::copy(data, data + rows * cols, m_data.begin());

    align();
} 

template <typename T>
Matrix<T>::Matrix(size_t rows, size_t cols, std::vector<T> tvec) :
    Matrix(rows, cols, tvec.data())
{}

template <typename T>
Matrix<T>::Matrix(size_t rows, size_t cols) :
    Matrix(rows, cols, nullptr)
{}

template <typename T>
std::pair<size_t, size_t> Matrix<T>::size() const {
    return { m_rows, m_cols };
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& B)
{

    size_t m,i,inner,n,j,coffset,boffset,b_i,k,aoffset;
    T temp;

    if (m_cols != B.m_rows) 
    {
        return Matrix<T>(0, 0);
    }

    Matrix<T> C(m_rows, B.m_cols);

    m = m_rows;

    if (m_cols == 1 || B.m_rows == 1)
    {
        for (i = 0; i < m; i++)
        {
            inner = B.m_cols;
            for (n = 0; n < inner; n++)
            {
                C.m_data[i + C.m_rows * n] = 0.0F;
                j = m_cols;
                for (coffset = 0; coffset < j; coffset++)
                {
                    C.m_data[i + C.m_rows * n] += m_data[i + m_rows * coffset] *
                        B.m_data[coffset + B.m_rows * n];
                }
            }
        }
    }
    else {

        inner = m_cols;
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
                    C.m_data[i] += temp * m_data[aoffset + b_i];
                }
            }
        }
    }
    return C;
}


template <typename T>
void Matrix<T>::align()
{
    size_t size = m_rows * m_cols;
    auto tvec{std::vector<T>(size)};
    
    tvec = m_data;

    for (size_t i{}, k{}; i < m_rows; ++i)
        for (size_t j{}; j < m_cols; ++j)
            m_data[i + m_rows * j] = tvec[k++];
            
}








