#pragma once
#include <vector>
#include <cmath>




template <typename T>
class Matrix 
{
public:
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
        return Matrix(rows, cols, std::vector<T>(rows * cols, val));
    }

    static Matrix makeLinSpace(T begin, T end, size_t n);


    std::pair<size_t, size_t> size() const;

    


    template<typename U>
    explicit operator Matrix<U>() const;

    explicit operator bool() const;


    const T& operator()(size_t, size_t) const;
    T& operator()(size_t, size_t);

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

    Matrix& operator++();
    Matrix operator++(int);

    Matrix& operator--();
    Matrix operator--(int);


    // friends
    friend class Matrix;

    template <typename U, typename T>
    friend Matrix<std::common_type_t<U, T>> operator*(Matrix<U> A, const Matrix<T>& B);
    
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
    if (n < 1)
        return Matrix::makeMatrix(0, 0); //todo: exception

    Matrix C = Matrix::makeMatrix(1, n);

    C.m_data[n - 1] = end;
    C.m_data[0] = begin;

    if (n == 1)
        return C;

    
    if ((begin == -end) && (n > 2))
    {
        for (size_t i{2}; i <= n; ++i)
            C.m_data[i - 1] = static_cast<T>(end * ((double)(((i << 1) - n) - 1) / ((double)n - 1.0)));

        if ((n & 1) == 1)
        {
            C.m_data[n >> 1] = (T)0.0;
        }
    }
    else if ((begin < (T)0.0) != (end < (T)0.0))
    {
        size_t val = C.m_cols - 1;
        double delta1 = begin / static_cast<double>(val);
        double delta2 = end / static_cast<double>(val);

        val = C.m_cols - 3;
        for (size_t i{}; i <= val; ++i)
            C.m_data[i + 1] = static_cast<T>((begin + delta2 * (static_cast<double>(i) + 1.0)) - delta1 * (static_cast<double>(i) + 1.0));
    }
    else
    {
        size_t val = C.m_cols - 1;
        double delta1 = (end - begin) / static_cast<double>(val);
        
        val = C.m_cols - 3;
        for (size_t i{}; i <= val; ++i)
            C.m_data[i + 1] = static_cast<T>(begin + (static_cast<double>(i) + 1.0) * delta1);
    }
 
    return C;
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
    size_t size = A.m_rows * A.m_cols;
    auto C = Matrix<std::common_type_t<U, T>>::makeMatrix(A.m_rows, A.m_cols);

    for (size_t i{}; i < size; ++i)
        C.m_data[i] = scalar / A.m_data[i];

    return C;
}

template<typename T, typename U>
inline Matrix<std::common_type_t<U, T>> operator/(Matrix<T> A, Matrix<U> B)
{
    auto [Arows, Acols] = A.size();
    auto [Brows, Bcols] = B.size();

    if ((Arows != Brows) || (Acols != Bcols))
        return Matrix<std::common_type_t<U, T>>::makeMatrix(0, 0); //todo: exception


    auto C = Matrix<std::common_type_t<U, T>>::makeMatrix(Arows, Acols);

    for (int i{}; i < Arows; ++i)
        for (int j{}; j < Acols; ++j)
        {
            if(B(i, j) <= FLT_EPSILON)
                return Matrix<std::common_type_t<U, T>>::makeMatrix(0, 0);//todo exception
            C(i, j) = A(i, j) / B(i, j);
        }
            

    return C;
}






template<typename T>
Matrix<T>& Matrix<T>::operator++()
{
    *this += 1;
    return *this;
}

template<typename T>
Matrix<T> Matrix<T>::operator++(int)
{
    auto C = Matrix(m_rows, m_cols, this->m_data);
    *this += 1;
    return C;
}

template<typename T>
Matrix<T>& Matrix<T>::operator--()
{
    *this -= 1;
    return *this;
}

template<typename T>
Matrix<T> Matrix<T>::operator--(int)
{
    auto C = Matrix(m_rows, m_cols, this->m_data);
    *this -= 1;
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
    return Matrix<std::common_type_t<U, T>>(A) += B;
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
    return Matrix<std::common_type_t<U, T>>(A) -= scalar;
}


template<typename T, typename U>
inline Matrix<std::common_type_t<U, T>> operator-(U scalar, Matrix<T> A)
{
    return Matrix<std::common_type_t<U, T>>(-A + scalar);
}


template <typename T, typename U>
inline Matrix<std::common_type_t<U, T>> operator-(Matrix<T> A, const Matrix<U>& B)
{
    // if A is temporary object, result is copy elision
    return Matrix<std::common_type_t<U, T>>(A) -= B;
}

template<typename T>
inline Matrix<T> Matrix<T>::operator-() const
{
    std::vector<T> vec(m_rows * m_cols);

    for (size_t i{}; i < m_rows * m_cols; ++i)
        vec[i] = -this->m_data[i];

    auto C{ Matrix(m_rows, m_cols, std::move(vec)) };

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
    return Matrix<U>(m_rows, m_cols, std::vector<U>(m_data.begin(), m_data.end()));
}

template<typename T>
Matrix<T>::operator bool() const
{
    auto size{ m_rows * m_cols };
  
    for (size_t i{}; i < size; ++i)
        if (m_data[i] != 0)
            return true;
            
    return  false;
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












template<typename T>
template<typename U>
Matrix<T>& Matrix<T>::operator*=(U scalar)
{
    *this = *this * (T)scalar;
    return *this;
}



template <typename U, typename T> 
Matrix<std::common_type_t<U, T>> operator*(Matrix<U> A, const Matrix<T>& B)
{
    std::common_type_t<U, T> val;
    auto C = Matrix<std::common_type_t<U, T>>::makeMatrix(A.m_rows, B.m_cols);

    if (A.m_cols != B.m_rows)
        return Matrix<std::common_type_t<U, T>>::makeMatrix(0, 0); // todo: exception

    for (size_t i = 0; i < A.m_rows; i++)
        for (size_t j = 0; j < B.m_cols; j++) {
            val = (std::common_type_t<U, T>)0.0;
            for (size_t k = 0; k < A.m_cols; k++)
                val += A.m_data[i + A.m_rows * k] * B.m_data[k + B.m_rows * j];

            C.m_data[i + C.m_rows * j] = val;
        }
    return C;
}



template <typename T, typename U>
inline Matrix<std::common_type_t<U, T>> operator*(const Matrix<T>& A, U scalar)
{

    size_t size = A.m_rows * A.m_cols;
    auto C = Matrix<std::common_type_t<U, T>>::makeMatrix(A.m_rows, A.m_cols);

    for (size_t i = 0; i < size; ++i) {
        C.m_data[i] = A.m_data[i] * scalar;
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


