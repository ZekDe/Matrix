#pragma once
#include <memory>
#include <vector>
#include <cmath>




template <typename T>
class Matrix 
{
public:

    template <typename U, typename T>
    using ut = Matrix<std::common_type_t<U, T>>;

    static Matrix<T> makeMatrix(size_t rows, size_t cols, std::vector<T> tvec)
    {
        if (tvec.empty() || rows * cols != tvec.size())
            return Matrix(0,0);

        align(rows, cols, tvec);
        return Matrix(rows, cols, tvec);
    }
    
    /* This must be used carefully, in case of mismatched rows, cols
    * and data is undefined-behaviour  */
    static Matrix makeMatrix(size_t rows, size_t cols, const T* data)
    {
        return makeMatrix(rows, cols, std::vector<T>(data, data + rows * cols));
    }

    static Matrix makeMatrix(size_t rows, size_t cols, T val = 0)
    {
        std::vector<T> data(rows * cols, val);
        return Matrix(rows, cols, data);
    }

    static Matrix makeLinSpace(T begin, T end, size_t n);


    std::pair<size_t, size_t> size() const;

    template<typename U>
    explicit operator Matrix<U>() const; //conversion

    const T& operator()(size_t, size_t) const; // get index
    T& operator()(size_t, size_t); // set index

    template <typename U>
    Matrix& operator*=(U);

    template<typename U>
    Matrix& operator/=(U); 


    template<typename U>
    Matrix& operator+=(const Matrix<U>&); 

    template<typename U>
    Matrix& operator+=(U);

    template<typename U>
    Matrix& operator-=(const Matrix<U>&);

    template<typename U>
    Matrix& operator-=(U);

    Matrix operator-() const;
    Matrix operator+() const;


    // friends
    friend class Matrix;

    template <typename U, typename T>
    friend Matrix<std::common_type_t<U, T>> operator*(const Matrix<U>& A, const Matrix<T>& B);

    template <typename T, typename U>
    friend Matrix<std::common_type_t<U, T>> operator*(const Matrix<T>&, U);

    template<typename T, typename U>
    friend Matrix<std::common_type_t<U, T>> operator/(U, Matrix<T>);

private:
    Matrix(size_t rows, size_t cols, std::vector<T> tvec = std::vector<T>()) :
        m_rows(rows), m_cols(cols), m_data(tvec)
    {}

    static void align(size_t, size_t, std::vector<T>&);

    size_t m_rows, m_cols;
    std::vector<T> m_data;
 
};













template <typename T>
Matrix<T> Matrix<T>::makeLinSpace(T begin, T end, size_t n)
{
    Matrix A = Matrix::makeMatrix(1, n);

    double delta1;
    size_t tmp;
    size_t i;
    size_t k;
    double delta2;

    delta1 = static_cast<double>(n);
    A.m_rows = 1;

    i = static_cast<size_t>(std::floor(delta1));
    A.m_cols = i;

    if (i >= 1)
    {
        tmp = i - 1;
        A.m_data[tmp] = end;
        if (A.m_cols >= 2)
        {
            A.m_data[0] = begin;
            if (A.m_cols >= 3)
            {
                if ((begin == -end) && (i > 2))
                {
                    for (k = 2; k <= tmp; k++) {
                        A.m_data[k - 1] = static_cast<T>(end * ((double)(((k << 1) - i) - 1) / ((double)i - 1.0)));
                    }

                    if ((i & 1) == 1)
                    {
                        A.m_data[i >> 1] = static_cast<T>(0.0);
                    }
                }
                else if ((begin < (T)0.0) != (end < (T)0.0))
                {
                    delta1 = begin / (static_cast<double>(A.m_cols) - 1.0);
                    delta2 = end / (static_cast<double>(A.m_cols) - 1.0);
                    tmp = A.m_cols;
                    for (k = 0; k <= tmp - 3; k++)
                    {
                        A.m_data[k + 1] = static_cast<T>((begin + delta2 * (static_cast<double>(k) + 1.0)) - delta1 *
                            (static_cast<double>(k) + 1.0));
                    }
                }
                else
                {
                    delta1 = (end - begin) / (static_cast<double>(A.m_cols) - 1.0);
                    tmp = A.m_cols;
                    for (k = 0; k <= tmp - 3; k++)
                    {
                        A.m_data[k + 1] = static_cast<T>(begin + (static_cast<double>(k) + 1.0) * delta1);
                    }
                }
            }
        }
    }
    return A;
}



template<typename T>
template<typename U>
Matrix<T>& Matrix<T>::operator*=(U scalar)
{
    *this = *this * static_cast<T>(scalar);
    return *this;
}



template<typename T>
template<typename U>
Matrix<T>& Matrix<T>::operator/=(U scalar)
{
    size_t size = m_rows * m_cols;
    for (size_t i{}; i < size; ++i)
        m_data[i] /= scalar;

    return *this;
}



template<typename T, typename U>
inline Matrix<std::common_type_t<U, T>> operator/(Matrix<T> A, U scalar)
{
    return A /= scalar;
}

template<typename T, typename U>
Matrix<std::common_type_t<U, T>> operator/(U scalar, Matrix<T> A)
{    
    auto C = Matrix<std::common_type_t<U, T>>::makeMatrix(A.m_rows, A.m_cols, 
        std::vector<std::common_type_t<U, T>>(A.m_rows * A.m_cols));

    for (size_t i{}; i < A.m_rows; ++i)
        for (size_t j{}; j < A.m_cols; ++j)
            C(i, j) = scalar / A(i, j);

    return C;
}



template<typename T>
template<typename U>
Matrix<T>& Matrix<T>::operator+=(const Matrix<U>& A)
{
    auto [Arows, Acols] = A.size();
    size_t size = m_rows * m_cols;

    if ((Arows != m_rows) || (Acols != m_cols))
        return *this; // todo: exception throw 

    for (size_t i{}; i < size; ++i)
        m_data[i] += A.m_data[i];

    return *this;
}



template<typename T>
template<typename U>
Matrix<T>& Matrix<T>::operator+=(U scalar)
{
    size_t size = m_rows * m_cols;
    for (size_t i{}; i < size; ++i)
        m_data[i] += scalar;
    return *this;
}



template<typename T, typename U>
inline Matrix<std::common_type_t<U, T>> operator+(Matrix<T> A, U scalar)
{
    return Matrix<std::common_type_t<U, T>>(A += scalar);
}


template<typename T, typename U>
inline Matrix<std::common_type_t<U, T>> operator+(U scalar, Matrix<T> A)
{
    return Matrix<std::common_type_t<U, T>>(A += scalar);
}



template <typename T, typename U>
inline Matrix<std::common_type_t<U, T>> operator+(Matrix<T> A, const Matrix<U>& B)
{
    // if A is temporary object, result is copy elision
    return Matrix<std::common_type_t<U, T>>(A += B);
}



template<typename T>
inline Matrix<T> Matrix<T>::operator+() const
{
    return *this;
}



template<typename T>
template<typename U>
Matrix<T>& Matrix<T>::operator-=(const Matrix<U>& A)
{
    return *this += -A;
}

template<typename T>
template<typename U>
Matrix<T>& Matrix<T>::operator-=(U scalar)
{

    return *this += -scalar;
}

template<typename T, typename U>
inline Matrix<std::common_type_t<U, T>> operator-(Matrix<T> A, U scalar)
{
    return Matrix<std::common_type_t<U, T>>(A -= scalar);
}


template<typename T, typename U>
inline Matrix<std::common_type_t<U, T>> operator-(U scalar, Matrix<T> A)
{
    A -= scalar;
    A = -A;
    return Matrix<std::common_type_t<U, T>>(A);
}



template <typename T, typename U>
inline Matrix<std::common_type_t<U, T>> operator-(Matrix<T> A, const Matrix<U>& B)
{
    // if A is temporary object, result is copy elision
    return Matrix<std::common_type_t<U, T>>(A -= B);
}

template<typename T>
inline Matrix<T> Matrix<T>::operator-() const
{
    std::vector<T> vec(m_rows * m_cols);
    auto C{ Matrix::makeMatrix(m_rows, m_cols, vec) };
    
    for (size_t i{}; i < m_rows; ++i)
        for (size_t j{}; j < m_cols; ++j)
            C(i, j) = -(*this)(i, j);

    return C;
}









template <typename T>
inline bool operator==(const Matrix<T>& A, const Matrix<T>& B)
{
    auto [Arows, Acols] = A.size();
    auto [Brows, Bcols] = B.size();

    if ((Arows != Brows) || (Acols != Bcols))
        return false;

    for (size_t i{}; i < Arows; ++i)
        for (size_t j{}; j < Acols; ++j)
            if (std::abs(A(i, j) - B(i, j)) > FLT_EPSILON)
                return false;
    return true;
    
}



template <typename T>
std::pair<size_t, size_t> Matrix<T>::size() const 
{
    return { m_rows, m_cols };
}



template<typename T>
template<typename U>
Matrix<T>::operator Matrix<U>() const 
{
    std::vector<U> data(m_data.begin(), m_data.end());
    return Matrix<U>(m_rows, m_cols, data);
}



template<typename T>
const T& Matrix<T>::operator()(size_t row, size_t col) const
{
    return m_data[row + m_rows * col];
}



template<typename T>
T& Matrix<T>::operator()(size_t row, size_t col)
{
    return m_data[row + m_rows * col];
}



template <typename U, typename T>
 Matrix<std::common_type_t<U, T>> operator*(const Matrix<U>& A, const Matrix<T>& B)
{
    size_t m, i, inner, n, j, coffset, boffset, b_i, k, aoffset;
    std::common_type_t<U, T> temp;
    
    if (A.m_cols != B.m_rows)
    {
        return Matrix<std::common_type_t<U, T>>::makeMatrix(0, 0);
    }

    auto C = Matrix<std::common_type_t<U, T>>::makeMatrix(A.m_rows, B.m_cols);

    m = A.m_rows;

    if (A.m_cols == 1 || B.m_rows == 1)
    {
        for (i = 0; i < m; i++)
        {
            inner = B.m_cols;
            for (n = 0; n < inner; n++)
            {
                C.m_data[i + C.m_rows * n] = static_cast<std::common_type_t<U, T>>(0.0);
                j = A.m_cols;

                for (coffset = 0; coffset < j; coffset++)
                    C.m_data[i + C.m_rows * n] += static_cast<std::common_type_t<U, T>>(A.m_data[i + A.m_rows * coffset]) *
                                                  static_cast<std::common_type_t<U, T>>(B.m_data[coffset + B.m_rows * n]);

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
                C.m_data[coffset + b_i] = static_cast < std::common_type_t<U, T>>(0.0F);

            for (k = 0; k < inner; k++)
            {
                aoffset = k * m;
                temp = static_cast<std::common_type_t<U, T>>(B.m_data[boffset + k]);
                for (b_i = 0; b_i < m; b_i++)
                {
                    i = coffset + b_i;
                    C.m_data[i] += temp * static_cast<std::common_type_t<U, T>>(A.m_data[aoffset + b_i]);
                }
            }
        }
    }
    return C;
}



template <typename T, typename U>
Matrix<std::common_type_t<U, T>> operator*(const Matrix<T>& A, U scalar)
{

    size_t size = A.m_rows * A.m_cols;
    auto C = Matrix<std::common_type_t<U, T>>::makeMatrix(A.m_rows, A.m_cols);

    for (size_t i = 0; i < size; ++i) {
        C.m_data[i] = A.m_data[i] * static_cast<std::common_type_t<U, T>>(scalar);
    }
    return C;
}



template <typename T, typename U>
inline Matrix<std::common_type_t<U, T>> operator*(U scalar, const Matrix<T>& A)
{
    return A * scalar;
}



template <typename T>
void Matrix<T>::align(size_t row, size_t col, std::vector<T>& tvec)
{
    auto localvec{ std::vector<T>(row * col) };

    localvec = tvec;

    for (size_t i{}, k{}; i < row; ++i)
        for (size_t j{}; j < col; ++j)
            tvec[i + row * j] = localvec[k++];

}


