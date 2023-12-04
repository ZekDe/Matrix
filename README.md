### Member Functions
```cpp
Matrix makeMatrix(size_t rows, size_t cols, std::vector<T>&& tvec);

Matrix makeMatrix(size_t rows, size_t cols, const T* data);

Matrix makeMatrix(size_t rows, size_t cols, T val = 0);

Matrix makeMatrix(size_t rows, size_t cols, std::initializer_list<T>&& il);

Matrix makeRandomMatrix(size_t rows, size_t cols, T min = (T)0.0, T max = (T)1.0);

Matrix makeLinSpace(T begin, T end, size_t n);

Matrix makeLinInc(T begin, T interval, T end);


std::pair<size_t, size_t> size() const;


const T& operator()(size_t, size_t) const;

T& operator()(size_t, size_t);

std::pair<size_t, size_t> operator()(T) const;

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
```
### Friends
```cpp
friend class Matrix;

template <typename U, typename T>
friend Matrix<std::common_type_t<U, T>> operator*(const Matrix<U>& A, const Matrix<T>& B);

template <typename T, typename U>
friend Matrix<std::common_type_t<U, T>> operator*(const Matrix<T>&, U);

template<typename T, typename U>
friend Matrix<std::common_type_t<U, T>> operator/(U, const Matrix<T>&);
```
### Non-Members
```cpp
template<typename T, typename U>
Matrix<std::common_type_t<U, T>> operator+(const Matrix<T>& A, U scalar);

template<typename T, typename U>
Matrix<std::common_type_t<U, T>> operator+(U scalar, const Matrix<T>& A);

template <typename T, typename U>
Matrix<std::common_type_t<U, T>> operator+(const Matrix<T>& A, const Matrix<U>& B);

template<typename T, typename U>
Matrix<std::common_type_t<U, T>> operator-(Matrix<T> A, U scalar);

template<typename T, typename U>
Matrix<std::common_type_t<U, T>> operator-(U scalar, Matrix<T> A);

template <typename T, typename U>
Matrix<std::common_type_t<U, T>> operator-(Matrix<T> A, const Matrix<U>& B);

template <typename T>
bool operator==(const Matrix<T>& A, const Matrix<T>& B);

template <typename T, typename U>
Matrix<std::common_type_t<U, T>> operator*(U scalar, const Matrix<T>& A);

template<typename T, typename U>
inline Matrix<std::common_type_t<U, T>> operator/(const Matrix<T>& A, U scalar);

template<typename T, typename U>
Matrix<std::common_type_t<U, T>> operator/(const Matrix<T>& A, const Matrix<U>& B);

```