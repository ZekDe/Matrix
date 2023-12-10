### Member Functions

namespace MathLab

```cpp
Matrix<T> makeMatrix(size_t rows, size_t cols, std::vector<T>&& tvec);

Matrix<T> makeMatrix(size_t rows, size_t cols, T val = 0);

Matrix<T> makeMatrix(size_t rows, size_t cols, std::initializer_list<T>&& il);

Matrix<T> makeRandomMatrix(size_t rows, size_t cols, T min = (T)0.0, T max = (T)1.0);

static Matrix<T> makeEyeMatrix(size_t rows, size_t cols);
static Matrix<T> makeEyeMatrix(size_t n); 

// linearly spaced vector generate n points between begin-end
Matrix<T> makeLinSpace(T begin, T end, size_t n);

// linearly incremented vector generate a range by dt
Matrix<T> makeLinInc(T begin, T dt, T end);


std::pair<size_t, size_t> size() const;
void setSize(size_t rows, size_t cols);

const T& operator()(size_t, size_t) const;

T& operator()(size_t, size_t);

std::pair<size_t, size_t> operator()(T) const;

template <typename U>
Matrix<T>& operator*=(U);

template<typename U>
Matrix<T>& operator/=(U);

template<typename U>
Matrix<T>& operator+=(const Matrix<U>&);

template<typename U>
Matrix<T>& operator+=(U);

template<typename U>
Matrix<T>& operator-=(const Matrix<U>&);

template<typename U>
Matrix<T>& operator-=(U);

Matrix<T> operator-() const;
Matrix<T> operator+() const;

Matrix<T>& operator++();
Matrix<T> operator++(int);

Matrix<T>& operator--();
Matrix<T> operator--(int);
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

### Example - Exceptions
Same as vector class. Calls to Allocator::allocate may throw. 
```cpp
Matrix<T> makeMatrix(size_t rows, size_t cols, T val = 0);
Matrix<T> makeMatrix(size_t rows, size_t cols, std::initializer_list<T>&& il);
Matrix<T> makeRandomMatrix(size_t rows, size_t cols, T min = (T)0.0, T max = (T)1.0);
static Matrix<T> makeEyeMatrix(size_t rows, size_t cols);
static Matrix<T> makeEyeMatrix(size_t n); 

Matrix<std::common_type_t<U, T>> operator/(U scalar, const Matrix<T>& A);
Matrix<std::common_type_t<U, T>> operator/(const Matrix<T>& A, U scalar);
Matrix<std::common_type_t<U, T>> operator/(const Matrix<T>& A, const Matrix<U>& B);
Matrix<std::common_type_t<U, T>> operator+(const Matrix<T>& A, U scalar);
Matrix<std::common_type_t<U, T>> operator+(U scalar, const Matrix<T>& A);
Matrix<std::common_type_t<U, T>> operator+(const Matrix<T>& A, const Matrix<U>& B);
Matrix<T> Matrix<T>::operator+() const;
Matrix<T> operator++(int);
Matrix operator--(int);
Matrix<std::common_type_t<U, T>> operator-(Matrix<T> A, U scalar);
Matrix<T>::operator-() const;
Matrix<T>::setSize(size_t rows, size_t cols);
Matrix<T>::operator Matrix<U>() const;
Matrix<std::common_type_t<U, T>> operator*(const Matrix<U>& A, const Matrix<T>& B);
Matrix<std::common_type_t<U, T>> operator*(const Matrix<T>& A, U scalar);
Matrix<std::common_type_t<U, T>> operator*(U scalar, const Matrix<T>& A);
```

### Example - Create Matrix
```cpp

auto B = Matrix<int>::makeMatrix(4, 4, 
{ 1,2,4,-80,
-5,2,0,-10,
1,10,3,20,
1,2, 2,3 });

B.setSize(9, 8);

auto C = Matrix<float>::makeMatrix(3, 2, std::vector<float>
{8.1,5.2,
6.3,7.444,
8.5,9.6 });

#define PI 3.1415926535
auto D = Matrix<float>::makeLinSpace(-PI, +PI, 10);
auto E = Matrix<float>::makeLinInc(-PI, 0.5, PI);
auto F = Matrix<float>::makeRandomMatrix(3, 3);
auto I = Matrix<float>::makeEyeMatrix(3, 3);
auto I1 = Matrix<float>::makeEyeMatrix(4);


```
### Example - Matrix Operations
```cpp
cout << A * B;
cout << A / B;
A+=B;
cout << A;
A -=B;
cout << (A + B);
cout << 2 - A + 2;
cout << 2 / A / 2;
A++;
++A;
A--;
--A;
// ...

cout << A;
A(1, 2) = 2;
cout << A << A(1,2);

auto[row, col] = A(10); // find 10 in A
cout << "A(" << row << "," << col << ")";

auto [rows, cols] = A.size();

cout << det(A);
cout << inv(A);
cout << transpose(A);


```
