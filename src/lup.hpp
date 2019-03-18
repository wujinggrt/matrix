#ifndef LUP_H__
#define LUP_H__

#include "mat_type.hpp"

namespace wj{

template<typename T = float>
class Mat;

/*************************************************************
 * Solving linear system
 * 参考《算法导论》的矩阵操作章节
*************************************************************/

template<typename MatValueType>
std::tuple<Mat<MatValueType>, Mat<MatValueType>> LUDecomposition(const Mat<MatValueType>& mat) {
    auto a = mat.Clone();
    auto n = a.RowSize();
    Mat<MatValueType> l(a.RowSize(), a.ColSize());
    Mat<MatValueType> u(a.RowSize(), a.ColSize());
    for (int32_t i = 0; i < l.RowSize(); ++i) {
        for (int32_t j = 0; j < l.ColSize(); ++j) {
            if (i == j) {
                l[i][j] = 1;
            }
        }
    }
    for (int32_t k = 0; k < n; ++k) {
        u[k][k] = a[k][k];
        for (int32_t i = k + 1; i < n; ++i) {
            // col-vector
            l[i][k] = a[i][k] / u[k][k];
            // row-vector
            u[k][i] = a[k][i];
        }
        for (int32_t i = k + 1; i < n; ++i) {
            for (int32_t j = k + 1; j < n; ++j) {
                a[i][j] = a[i][j] - l[i][k] * u[k][j];
            }
        }
    }
    return std::make_tuple(l, u);
}

// 奇异矩阵的话会抛出exception:invalid_argument
template<typename MatValueType>
std::tuple<Mat<MatValueType>, Mat<MatValueType>> LUPDecomposition(const Mat<MatValueType>& mat) {
    auto a = mat.Clone();
    auto n = a.RowSize();
    Mat<MatValueType> pi(mat.RowSize(), 1);
    for (int32_t i = 0; i < pi.RowSize(); ++i) {
        pi[i][0] = static_cast<MatValueType>(i);
    }
    for (int32_t k = 0; k < n; ++k) {
        auto p = 0.;
        auto k2 = k;
        // find the max absolute value.
        for (int32_t i = k; i < n; ++i) {
            if (std::abs(a[i][k]) > p) {
                p = std::abs(a[i][k]);
                k2 = i;
            }
        }
        if (p == 0.) {
            throw std::invalid_argument("singular matrix");
        }
        std::swap(pi[k][0], pi[k2][0]);
        // swap rows
        for (int32_t i = 0; i < n; ++i) {
            std::swap(a(k, i), a(k2, i));
        }
        for (int32_t i = k + 1; i < n; ++i) {
            // col-vector
            a(i, k) = a(i, k) / a(k, k);
            for (int32_t j = k + 1; j < n; ++j) {
                a(i, j) = a(i, j) - a(i, k) * a(k, j);
            }
        }
    }

    return std::make_tuple(pi, a);
}

template<typename MatValueType>
Mat<MatValueType> LUPSolve(Mat<MatValueType>& l, Mat<MatValueType>& u, Mat<MatValueType>& pi, Mat<MatValueType>& b) {
    auto n = l.RowSize();
    Mat<MatValueType> x(n, 1);
    Mat<MatValueType> y(n, 1);
    for (int32_t i = 0; i < n; ++i) {
        auto sum = 0.;
        for (int32_t j = 0; j < i; ++j) {
            sum += l[i][j] * y[j][0];
        }
        y[i][0] = b[pi[i][0]][0] - sum;
    }
    for (int32_t i = n - 1; i >= 0; --i) {
        auto sum = 0.;
        for (int32_t j = i + 1; j < n; ++j) {
            sum += u[i][j] * x[j][0];
        }
        x[i][0] = (y[i][0] - sum) / u[i][i];
    }
    return x;
}

/*************************************************************
 * Solving linear system
 * end
*************************************************************/

}
#endif