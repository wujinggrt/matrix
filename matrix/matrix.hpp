#ifndef MATRIX_H__
#define MATRIX_H__

#include <cmath>

#include <vector>
#include <string>
#include <iostream>
#include <initializer_list>
#include <stdexcept>
#include <tuple>
#include <random>
#include <functional>

namespace wj {

template<typename T = float>
class Mat;

/*************************************************************
 * Operation declaration
*************************************************************/

template<typename MatValueType>
std::ostream& operator<<(std::ostream& os, const Mat<MatValueType>& mat);

template<typename MatValueType>
std::string to_string(const Mat<MatValueType>& mat);

/*************************************************************
 * Binary operator overload declaration
 * end
*************************************************************/

template<typename MatValueType>
std::tuple<Mat<MatValueType>, Mat<MatValueType>> LUDecomposition(const Mat<MatValueType>& mat);

template<typename MatValueType>
std::tuple<Mat<MatValueType>, Mat<MatValueType>> LUPDecomposition(const Mat<MatValueType>& mat);

template<typename MatValueType>
Mat<MatValueType> LUPSolve(Mat<MatValueType>& l, Mat<MatValueType>& u, Mat<MatValueType>& pi, Mat<MatValueType>& b);

/*************************************************************
 * template class: Mat
*************************************************************/
template<typename T>
class Mat {
private:
    std::vector<std::vector<T>> data_;

public:
    using iterator = decltype(data_.begin());
    using value_type = T;
    using size_type = std::size_t;
    // using iterator = typename std::vector<std::vector<T>>::iterator;

    template<typename MatValueType>
    friend std::ostream& operator<<(std::ostream& os, const Mat<MatValueType>& mat);
    template<typename MatValueType>
    friend std::string to_string(const Mat<MatValueType>& mat);

    template<typename MatValueType>
    friend std::tuple<Mat<MatValueType>, Mat<MatValueType>> LUDecomposition(const Mat<MatValueType>& mat);
    template<typename MatValueType>
    friend std::tuple<Mat<MatValueType>, Mat<MatValueType>> LUPDecomposition(const Mat<MatValueType>& mat);

public:
    Mat(std::initializer_list<std::vector<T>> ls)
        : data_{ls} {
    }

    Mat(std::size_t rows = 1, std::size_t cols = 1)
        : data_(rows, std::vector<T>(cols)) {
        if (rows == 0 || cols == 0) {
            throw std::invalid_argument("argument can not be 0!");
        }
    }

    static Mat<T> Eye(std::size_t row) {
        Mat<T> ret(row, row);
        for (std::size_t i = 0; i < row; ++i) {
            ret[i][i] = 1;
        }
        return ret;
    }

    static Mat<T> Random(T low, T high, std::size_t row = 1, std::size_t col = 1) {
        Mat<T> ret(row, col);
        std::default_random_engine engine;
        std::uniform_real_distribution<T> num_distribution(low, high);
        for (std::size_t i = 0; i < row; ++i) {
            for (std::size_t j = 0; j < col; ++j) {
                ret[i][j] = num_distribution(engine);
            }
        }
        return ret;
    }

    /*
    * output the matrix contents to stdout.
    */
    void Print(std::ostream& os = std::cout) const {
        for (auto e: data_) {
            bool first = true;
            for (auto ee: e) {
                if (first) {
                    first = false;
                } else {
                    os << " ";
                }
                os << ee;
            }
            os << '\n';
        }
        os << '\n';
    }

    // throw range
    std::vector<T>& operator[](size_type i) {
        if (i >= data_.size()) {
            throw std::out_of_range("Mat:data_:index out of range!");
        }
        return data_[i];
    }

            // throw range
    const std::vector<T>& operator[](size_type i) const {
        if (i >= data_.size()) {
            throw std::out_of_range("Mat:data_:index out of range!");
        }
        return data_[i];
    }

    iterator begin() {
        std::begin(data_);
    }

    iterator end() {
        std::end(data_);
    }

    size_type RowSize() const {
        return data_.size();
    }

    size_type ColSize() const {
        return data_[0].size();
    }
    
    // +,-,*,/操作在行、列不匹配的时候则会抛出exception:invalid_argument
    Mat<T> operator+(const Mat<T>& other) const {
        CheckSize(other);

        Mat<T> ret(RowSize(), ColSize());
        for (size_type i = 0; i < RowSize(); ++i) {
            for (int j = 0; j < ColSize(); ++j) {
                ret[i][j] = data_[i][j] + other.data_[i][j];
            }
        }
        return ret;
    }

    Mat<T> operator-(const Mat<T>& other) const {
        CheckSize(other);

        Mat<T> ret(RowSize(), ColSize());
        for (size_type i = 0; i < RowSize(); ++i) {
            for (size_type j = 0; j < ColSize(); ++j) {
                ret[i][j] = data_[i][j] - other.data_[i][j];
            }
        }
        return ret;
    }

    Mat<T> operator*(const Mat<T>& other) const {
        if (ColSize() != other.RowSize()) {
            throw std::invalid_argument(std::string("incompatible dimensions\n") + 
                "right-matrix rows, cols:" + std::to_string(RowSize()) + ", " + std::to_string(ColSize()) +
                "\nleft-matrx rows, cols:" + std::to_string(other.RowSize()) + ", " + std::to_string(other.ColSize())
                );
        }

        Mat<T> ret(data_.size(), other.ColSize());
        // this rows
        for (size_type i = 0; i < RowSize(); ++i) {
            // other cols
            for (size_type j = 0; j < other.ColSize(); ++j) {
                ret[i][j] = 0.;
                for (size_type k = 0; k < ColSize(); ++k) {
                    ret[i][j] += data_[i][k] * other.data_[k][j];
                }
            }
        }
        return ret;
    }

    Mat<T> operator/(const Mat<T>& other) const {
        CheckSize(other);

        Mat<T> ret(RowSize(), ColSize());
        for (size_type i = 0; i < RowSize(); ++i) {
            for (size_type j = 0; j < ColSize(); ++j) {
                ret[i][j] = data_[i][j] / other.data_[i][j];
            }
        }
        return ret;
    }

    T& operator()(std::size_t row_index, size_t col_index) {
        return data_[row_index][col_index];
    }

    T operator()(std::size_t row_index, size_t col_index) const {
        return data_[row_index][col_index];
    }

    Mat<T> DotProduct(const Mat<T>& other) const {
        CheckSize(other);

        Mat<T> ret(RowSize(), ColSize());
        for (size_type i = 0; i < RowSize(); ++i) {
            for (size_type j = 0; j < ColSize(); ++j) {
                ret[i][j] = data_[i][j] * other.data_[i][j];
            }
        }
        return ret;
    }

    Mat<T> Transpose() const {
        Mat<T> ret(ColSize(), RowSize());
        for (size_type i = 0; i < RowSize(); ++i) {
            for (size_type j = 0; j < ColSize(); ++j) {
                ret[j][i] = data_[i][j];
            }
        }
        return ret;
    }

    Mat<T> Clone() const {
        Mat<T> ret(RowSize(), ColSize());
        for (size_type i = 0; i < RowSize(); ++i) {
            for (size_type j = 0; j < ColSize(); ++j) {
                ret.data_[i][j] = data_[i][j];
            }
        }
        return ret;
    }

    Mat<T> Inverse() const {
        Mat<T> ret(RowSize(), RowSize());
        auto r = LUPDecomposition(*this);
        Mat<T> l(RowSize(), ColSize());
        Mat<T> u(RowSize(), ColSize());
        Mat<T> pi = std::get<0>(r);
        Mat<T> I = Mat<T>::Eye(RowSize());
        for (size_type k = 0; k < RowSize(); ++k) {
            Mat<T> b(RowSize(), 1);
            for (size_type i = 0; i < RowSize(); ++i) {
                b[i][0] = I[i][k];
            }
            for (size_type i = 0; i < RowSize(); ++i) {
                for (size_type j = 0; j < l.ColSize(); ++j) {
                    if (i == j) {
                        l[i][j] = 1.;
                        u[i][j] = std::get<1>(r)[i][j];
                    } else if (i < j) {
                        u[i][j] = std::get<1>(r)[i][j];
                    } else
                    {
                        l[i][j] = std::get<1>(r)[i][j];
                    }
                }
            }
            auto x = wj::LUPSolve(l, u, pi, b);
            for (size_type i = 0; i < RowSize(); ++i) {
                ret[i][k] = x[i][0];
            }
        }
        return ret;
    }

private:
    // if the matrix other is not equal to this,
    // throw std::invalid_argument
    void CheckSize(const Mat<T>& other) const {
        if (RowSize() != other.RowSize() || ColSize() != other.ColSize()) {
            throw std::invalid_argument(std::string("incompatible dimensions\n") +
                                        "right-matrix rows, cols:" + std::to_string(RowSize()) + ", " + std::to_string(ColSize()) +
                                        "\nleft-matrx rows, cols:" + std::to_string(other.RowSize()) + ", " + std::to_string(other.ColSize())
                                        );
        }
    }
};

/*************************************************************
 * template class: Mat
 * end
*************************************************************/


/*************************************************************
 * binary operator overload implementation
*************************************************************/

// To do
// using SFINAE 来限制num必须是arithmetic 
// 这个函数主要用来完成单个数加上这个Mat，
// 然后给Mat中的元素与他进行进行op(num, each src element)
template<typename NumType, typename MatValueType, typename BinaryOperation>
Mat<MatValueType> DoBinaryOperate(
                                  NumType num,
                                  const Mat<MatValueType>& src,
                                  BinaryOperation op,
                                  std::enable_if_t<
                                      std::is_arithmetic<NumType>::value
                                      >* = nullptr) { 
    Mat<MatValueType> ret(src.RowSize(), src.ColSize());
    MatValueType casted_num = static_cast<MatValueType>(num);
    for (std::size_t i = 0; i < src.RowSize(); ++i) {
        for (std::size_t j = 0; j < src.ColSize(); ++j) {
            ret[i][j] = op(casted_num, src(i, j));
        }
    }
    return ret;
}

template<typename NumType, typename MatValueType>
Mat<MatValueType> operator+(NumType num, const Mat<MatValueType>& mat) {
    return DoBinaryOperate(num, mat, std::plus<MatValueType>());
}

template<typename NumType, typename MatValueType>
Mat<MatValueType> operator-(NumType num, const Mat<MatValueType>& mat) {
    return DoBinaryOperate(num, mat, std::minus<MatValueType>());
}

template<typename NumType, typename MatValueType>
Mat<MatValueType> operator*(NumType num, const Mat<MatValueType>& mat) {
    return DoBinaryOperate(num, mat, std::multiplies<MatValueType>());
}

template<typename NumType, typename MatValueType>
Mat<MatValueType> operator/(NumType num, const Mat<MatValueType>& mat) {
    return DoBinaryOperate(num, mat, std::divides<MatValueType>());
}

template<typename NumType, typename MatValueType, typename BinaryOperation>
Mat<MatValueType> DoReverseBinaryOperate(
                                         NumType num,
                                         const Mat<MatValueType>& src,
                                         BinaryOperation op,
                                         std::enable_if_t<
                                              std::is_arithmetic<NumType>::value
                                              >* = nullptr) {
    Mat<MatValueType> ret(src.RowSize(), src.ColSize());
    MatValueType casted_num = static_cast<MatValueType>(num);
    for (std::size_t i = 0; i < src.RowSize(); ++i) {
        for (std::size_t j = 0; j < src.ColSize(); ++j) {
            ret[i][j] = op(src(i, j), num);
        }
    }
    return ret;
}

template<typename NumType, typename MatValueType>
Mat<MatValueType> operator+(const Mat<MatValueType>& mat, NumType num) {
    return DoReverseBinaryOperate(num, mat, std::plus<MatValueType>());
}

template<typename NumType, typename MatValueType>
Mat<MatValueType> operator-(const Mat<MatValueType>& mat, NumType num) {
    return DoReverseBinaryOperate(num, mat, std::minus<MatValueType>());
}

template<typename NumType, typename MatValueType>
Mat<MatValueType> operator*(const Mat<MatValueType>& mat, NumType num) {
    return DoReverseBinaryOperate(num, mat, std::multiplies<MatValueType>());
}

template<typename NumType, typename MatValueType>
Mat<MatValueType> operator/(const Mat<MatValueType>& mat, NumType num) {
    return DoReverseBinaryOperate(num, mat, std::divides<MatValueType>());
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const Mat<T>& mat) {
    os << "matrix:\nvalue_type: " 
       << typeid(typename Mat<T>::value_type).name()
       << "\nsize:       " 
       << mat.RowSize() 
       << " x " 
       << mat.ColSize() 
       << '\n';
    os << "[";
    for (int i = 0; i < mat.RowSize(); ++i) {
        os << (i == 0 ? "" : " ") << "[";
        for (int j = 0; j < mat.ColSize(); ++j) {
            os << mat.data_[i][j] << (j == mat.ColSize() - 1 ? "" : ", ");
        }
        os << "]" << (i == mat.RowSize() - 1 ? "" : "\n");
    }
    os << "]\n";
    return os;
}

template<typename T>
std::string to_string(const Mat<T>& mat) {
    std::string ret;
    ret.reserve(mat.RowSize() * mat.ColSize());
    ret += "[";
    for (int i = 0; i < mat.RowSize(); ++i) {
        ret += (i == 0 ? "" : " ") + std::string("[ ");
        for (int j = 0; j < mat.ColSize(); ++j) {
            ret += std::to_string(mat.data_[i][j]) + std::string(", ");
        }
        ret += "]" + std::string(i == mat.RowSize() - 1 ? "" : "\n");
    }
    ret += "]\n";
    return ret;
}
/*************************************************************
 * binary operator overload implementation
 * end
*************************************************************/

/*************************************************************
 * Solving linear system
*************************************************************/
template<typename MatValueType>
std::tuple<Mat<MatValueType>, Mat<MatValueType>> LUDecomposition(const Mat<MatValueType>& mat) {
    auto a = mat.Clone();
    auto n = a.RowSize();
    Mat<MatValueType> l(a.RowSize(), a.ColSize());
    Mat<MatValueType> u(a.RowSize(), a.ColSize());
    for (int i = 0; i < l.RowSize(); ++i) {
        for (int j = 0; j < l.ColSize(); ++j) {
            if (i == j) {
                l[i][j] = 1;
            }
        }
    }
    for (int k = 0; k < n; ++k) {
        u[k][k] = a[k][k];
        for (int i = k + 1; i < n; ++i) {
            // col-vector
            l[i][k] = a[i][k] / u[k][k];
            // row-vector
            u[k][i] = a[k][i];
        }
        for (int i = k + 1; i < n; ++i) {
            for (int j = k + 1; j < n; ++j) {
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
    for (int i = 0; i < pi.RowSize(); ++i) {
        pi[i][0] = static_cast<MatValueType>(i);
    }
    for (int k = 0; k < n; ++k) {
        auto p = 0.;
        auto k2 = k;
        // find the max absolute value.
        for (int i = k; i < n; ++i) {
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
        for (int i = 0; i < n; ++i) {
            std::swap(a.data_[k][i], a.data_[k2][i]);
        }
        for (int i = k + 1; i < n; ++i) {
            // col-vector
            a.data_[i][k] = a.data_[i][k] / a.data_[k][k];
            for (int j = k + 1; j < n; ++j) {
                a.data_[i][j] = a.data_[i][j] - a.data_[i][k] * a.data_[k][j];
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
    for (int i = 0; i < n; ++i) {
        auto sum = 0.;
        for (int j = 0; j < i; ++j) {
            sum += l[i][j] * y[j][0];
        }
        y[i][0] = b[pi[i][0]][0] - sum;
    }
    for (int i = n - 1; i >= 0; --i) {
        auto sum = 0.;
        for (int j = i + 1; j < n; ++j) {
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

using Matd = wj::Mat<double>;
using Matf = wj::Mat<float>;
using Mati = wj::Mat<int>;
}

#endif