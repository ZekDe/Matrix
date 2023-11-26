//template <typename U, typename T>
//Matrix<std::common_type_t<U, T>> operator*(const Matrix<U> A, const Matrix<T>& B)
//{
//    size_t m, i, inner, n, j, coffset, boffset, b_i, k, aoffset;
//    std::common_type_t<U, T> temp;
//
//    if (A.m_rows != B.m_cols)
//    {
//        return Matrix<std::common_type_t<U, T>>::makeMatrix(0, 0);
//    }
//
//    auto C = Matrix<std::common_type_t<U, T>>::makeMatrix(A.m_rows, B.m_cols);
//
//    m = A.m_rows;
//
//    if (A.m_cols == 1 || B.m_rows == 1)
//    {
//        for (i = 0; i < m; i++)
//        {
//            inner = B.m_cols;
//            for (n = 0; n < inner; n++)
//            {
//                C.m_data[i + C.m_rows * n] = static_cast<std::common_type_t<U, T>>(0.0);
//                j = A.m_cols;
//
//                for (coffset = 0; coffset < j; coffset++)
//                    C.m_data[i + C.m_rows * n] += static_cast<std::common_type_t<U, T>>(A.m_data[i + A.m_rows * coffset]) *
//                    static_cast<std::common_type_t<U, T>>(B.m_data[coffset + B.m_rows * n]);
//
//            }
//        }
//    }
//    else {
//
//        inner = A.m_cols;
//        n = B.m_cols;
//
//        for (j = 0; j < n; j++)
//        {
//            coffset = j * m;
//            boffset = j * inner;
//            for (b_i = 0; b_i < m; b_i++)
//                C.m_data[coffset + b_i] = static_cast <std::common_type_t<U, T>>(0.0F);
//
//            for (k = 0; k < inner; k++)
//            {
//                aoffset = k * m;
//                temp = static_cast<std::common_type_t<U, T>>(B.m_data[boffset + k]);
//                for (b_i = 0; b_i < m; b_i++)
//                {
//                    i = coffset + b_i;
//                    C.m_data[i] += temp * static_cast<std::common_type_t<U, T>>(A.m_data[aoffset + b_i]);
//                }
//            }
//        }
//    }
//    return C;
//}