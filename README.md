# Matrix<T>

* public static member functions

Matrix makeMatrix(size_t rows, size_t cols, std::vector<T>&& tvec);

Matrix makeMatrix(size_t rows, size_t cols, const T* data);

Matrix makeMatrix(size_t rows, size_t cols, T val = 0);

Matrix makeMatrix(size_t rows, size_t cols, std::initializer_list<T>&& il);

Matrix makeRandomMatrix(size_t rows, size_t cols, T min = (T)0.0, T max = (T)1.0);

Matrix makeLinSpace(T begin, T end, size_t n);

Matrix makeLinInc(T begin, T interval, T end);

* public functions

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

* friends

friend class Matrix;

template <typename U, typename T>
friend Matrix<std::common_type_t<U, T>> operator*(const Matrix<U>& A, const Matrix<T>& B);

template <typename T, typename U>
friend Matrix<std::common_type_t<U, T>> operator*(const Matrix<T>&, U);

template<typename T, typename U>
friend Matrix<std::common_type_t<U, T>> operator/(U, const Matrix<T>&);

Example 1 - Matrix Creation


auto A = Matrix<int>::makeMatrix(4, 4, { 1,2,4,-80,
                                        -5,2,0,-10,
                                         1,10,3,20,
                                         1,2, 2,3 });
	
auto B = Matrix<float>::makeMatrix(3, 2, std::vector<float>
										{8.1,5.2,
										 6.3,7.444,
										 8.5,9.6 });
		 
float arr[] = { 1, 2, 4, -80,
               -5, 2, 0, -10,
                1, 10, 3, 20,
                1, 2, 2, 3, };

auto C = Matrix<float>::makeMatrix(4, 4, arr);		 