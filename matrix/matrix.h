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
 * Binary operator overload declaration
*************************************************************/

template<typename T>
std::ostream& operator<<(std::ostream& os, const Mat<T>& mat);

template<typename T>
std::string to_string(const Mat<T>& mat);

/*************************************************************
 * Binary operator overload declaration
 * end
*************************************************************/

template<typename T>
std::tuple<Mat<T>, Mat<T>> LU_decomposition(const Mat<T>& mat);

template<typename T>
std::tuple<Mat<T>, Mat<T>> LUP_decomposition(const Mat<T>& mat);

template<typename T>
Mat<T> lup_solve(Mat<T>& l, Mat<T>& u, Mat<T>& pi, Mat<T>& b);

/*************************************************************
 * template class: Mat
*************************************************************/
template<typename T>
class Mat {
private:
    std::vector<std::vector<T>> vec_;

public:
    using iterator = decltype(vec_.begin());
    using value_type = T;
    using size_type = std::size_t;
    // using iterator = typename std::vector<std::vector<T>>::iterator;

    friend std::ostream& operator<<<>(std::ostream& os, const Mat<T>& mat);
    friend std::string to_string<>(const Mat<T>& mat);

    friend std::tuple<Mat<T>, Mat<T>> LU_decomposition<>(const Mat<T>& mat);
    friend std::tuple<Mat<T>, Mat<T>> LUP_decomposition<>(const Mat<T>& mat);

public:
    Mat(std::initializer_list<std::vector<T>> ls)
        : vec_{ls} {
    }

    Mat(size_type rows = 1, size_type cols = 1)
        : vec_{} {
        if (rows == 0 || cols == 0) {
            throw std::invalid_argument("argument can not be 0!");
        }
        vec_.reserve(rows);
        for (size_type i = 0; i < rows; ++i) {
            vec_.push_back(std::vector<T>(cols));
        }
    }

    static Mat<T> eye(size_type row) {
        Mat<T> ret(row, row);
        for (int i = 0; i < row; ++i) {
            ret[i][i] = 1;
        }
        return ret;
    }

    static Mat<T> random(T low, T high, size_type row = 1, size_type col = 1) {
        Mat<T> ret(row, col);
        std::default_random_engine engine;
        std::uniform_real_distribution<T> num_distribution(low, high);
        for (int i = 0; i < row; ++i) {
            for (int j = 0; j < col; ++j) {
                ret[i][j] = num_distribution(engine);
            }
        }
        return ret;
    }

    /*
    * output the matrix contents to stdout.
    */
    void print(std::ostream& os = std::cout) const {
        for (auto e: vec_) {
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
    decltype(vec_[0])& operator[](size_type i) {
        if (i >= vec_.size()) {
            throw std::out_of_range("Mat:vec_:index out of range!");
        }
        return vec_[i];
    }

            // throw range
    const decltype(vec_[0])& operator[](size_type i) const {
        if (i >= vec_.size()) {
            throw std::out_of_range("Mat:vec_:index out of range!");
        }
        return vec_[i];
    }

    iterator begin() {
        begin(vec_);
    }

    iterator end() {
        end(vec_);
    }

    size_type row_size() const {
        return vec_.size();
    }

    size_type col_size() const {
        return vec_[0].size();
    }
    
    // +,-,*,/操作在行、列不匹配的时候则会抛出exception:invalid_argument
    Mat<T> operator+(const Mat<T>& other) const {
        check_size(other);

        Mat<T> ret(row_size(), col_size());
        for (size_type i = 0; i < row_size(); ++i) {
            for (int j = 0; j < col_size(); ++j) {
                ret[i][j] = vec_[i][j] + other.vec_[i][j];
            }
        }
        return ret;
    }

    Mat<T> operator-(const Mat<T>& other) const {
        check_size(other);

        Mat<T> ret(row_size(), col_size());
        for (size_type i = 0; i < row_size(); ++i) {
            for (size_type j = 0; j < col_size(); ++j) {
                ret[i][j] = vec_[i][j] - other.vec_[i][j];
            }
        }
        return ret;
    }

    Mat<T> operator*(const Mat<T>& other) const {
        if (col_size() != other.row_size()) {
            throw std::invalid_argument(std::string("incompatible dimensions\n") + 
                "right-matrix rows, cols:" + std::to_string(row_size()) + ", " + std::to_string(col_size()) +
                "\nleft-matrx rows, cols:" + std::to_string(other.row_size()) + ", " + std::to_string(other.col_size())
                );
        }

        Mat<T> ret(vec_.size(), other.col_size());
        // this rows
        for (size_type i = 0; i < row_size(); ++i) {
            // other cols
            for (size_type j = 0; j < other.col_size(); ++j) {
                ret[i][j] = 0.;
                for (size_type k = 0; k < col_size(); ++k) {
                    ret[i][j] += vec_[i][k] * other.vec_[k][j];
                }
            }
        }
        return ret;
    }

    Mat<T> operator/(const Mat<T>& other) const {
        check_size(other);

        Mat<T> ret(row_size(), col_size());
        for (size_type i = 0; i < row_size(); ++i) {
            for (size_type j = 0; j < col_size(); ++j) {
                ret[i][j] = vec_[i][j] / other.vec_[i][j];
            }
        }
        return ret;
    }

    T& operator()(std::size_t row_index, size_t col_index) {
        return vec_[row_index][col_index];
    }

    T operator()(std::size_t row_index, size_t col_index) const {
        return vec_[row_index][col_index];
    }

    Mat<T> dot_product(const Mat<T>& other) const {
        check_size(other);

        Mat<T> ret(row_size(), col_size());
        for (size_type i = 0; i < row_size(); ++i) {
            for (size_type j = 0; j < col_size(); ++j) {
                ret[i][j] = vec_[i][j] * other.vec_[i][j];
            }
        }
        return ret;
    }

    void prsize_type_type() const {
        std::cout << typeid(T).name() << '\n'; 
    }

    Mat<T> trans() const {
        Mat<T> ret(col_size(), row_size());
        for (size_type i = 0; i < row_size(); ++i) {
            for (size_type j = 0; j < col_size(); ++j) {
                ret[j][i] = vec_[i][j];
            }
        }
        return ret;
    }

    Mat<T> clone() const {
        Mat<T> ret(row_size(), col_size());
        for (size_type i = 0; i < row_size(); ++i) {
            for (size_type j = 0; j < col_size(); ++j) {
                ret.vec_[i][j] = vec_[i][j];
            }
        }
        return ret;
    }

    Mat<T> inv() const {
        Mat<T> ret(row_size(), row_size());
        auto r = LUP_decomposition(*this);
        Mat<T> l(row_size(), col_size());
        Mat<T> u(row_size(), col_size());
        Mat<T> pi = std::get<0>(r);
        Mat<T> I = Mat<T>::eye(row_size());
        for (size_type k = 0; k < row_size(); ++k) {
            Mat<T> b(row_size(), 1);
            for (size_type i = 0; i < row_size(); ++i) {
                b[i][0] = I[i][k];
            }
            for (size_type i = 0; i < row_size(); ++i) {
                for (size_type j = 0; j < l.col_size(); ++j) {
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
            auto x = wj::lup_solve(l, u, pi, b);
            for (size_type i = 0; i < row_size(); ++i) {
                ret[i][k] = x[i][0];
            }
        }
        return ret;
    }

private:
    // if the matrix other is not equal to this,
    // throw std::invalid_argument
    void check_size(const Mat<T>& other) const {
        if (row_size() != other.row_size() || col_size() != other.col_size()) {
            throw std::invalid_argument(std::string("incompatible dimensions\n") +
                                        "right-matrix rows, cols:" + std::to_string(row_size()) + ", " + std::to_string(col_size()) +
                                        "\nleft-matrx rows, cols:" + std::to_string(other.row_size()) + ", " + std::to_string(other.col_size())
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
#define BINARY_OPERATOR_OVERLOAD_RETURN(n, op, mat) \
    for (std::size_t i = 0; i < mat.row_size(); ++i) { \
        for (std::size_t j = 0; j < mat.col_size(); ++j) { \
            ret[i][j] = n op mat.vec_[i][j]; \
        } \
    }

// To do
// using SFINAE 来限制num必须是arithmetic 
// 这个函数主要用来完成单个数加上这个Mat，
// 然后给Mat中的元素与他进行进行op(num, each src element)
template<typename NumType, typename MatValueType, typename BinaryOperation>
Mat<MatValueType> do_binary_operate(
                                    NumType num,
                                    const Mat<MatValueType>& src,
                                    BinaryOperation op,
                                    std::enable_if_t<
                                        std::is_arithmetic_v<NumType>
                                        >* = nullptr) { 
    Mat<MatValueType> ret(src.row_size(), src.col_size());
    MatValueType casted_num = static_cast<MatValueType>(num);
    for (std::size_t i = 0; i < src.row_size(); ++i) {
        for (std::size_t j = 0; j < src.col_size(); ++j) {
            ret[i][j] = op(casted_num, src(i, j));
        }
    }
    return ret;
}

template<typename NumType, typename MatValueType>
Mat<MatValueType> operator+(NumType num, const Mat<MatValueType>& mat) {
    return do_binary_operate(num, mat, std::plus<MatValueType>());
}

template<typename NumType, typename MatValueType>
Mat<MatValueType> operator-(NumType num, const Mat<MatValueType>& mat) {
    return do_binary_operate(num, mat, std::minus<MatValueType>());
}

template<typename NumType, typename MatValueType>
Mat<MatValueType> operator*(NumType num, const Mat<MatValueType>& mat) {
    return do_binary_operate(num, mat, std::multiplies<MatValueType>());
}

template<typename NumType, typename MatValueType>
Mat<MatValueType> operator/(NumType num, const Mat<MatValueType>& mat) {
    return do_binary_operate(num, mat, std::divides<MatValueType>());
}

template<typename NumType, typename MatValueType, typename BinaryOperation>
Mat<MatValueType> do_reverse_binary_operate(
                                           NumType num,
                                           const Mat<MatValueType>& src,
                                           BinaryOperation op,
                                           std::enable_if_t<
                                               std::is_arithmetic_v<NumType>
                                               >* = nullptr) {
    Mat<MatValueType> ret(src.row_size(), src.col_size());
    MatValueType casted_num = static_cast<MatValueType>(num);
    for (std::size_t i = 0; i < src.row_size(); ++i) {
        for (std::size_t j = 0; j < src.col_size(); ++j) {
            ret[i][j] = op(src(i, j), num);
        }
    }
    return ret;
}

template<typename NumType, typename MatValueType>
Mat<MatValueType> operator+(const Mat<MatValueType>& mat, NumType num) {
    return do_reverse_binary_operate(num, mat, std::plus<MatValueType>());
}

template<typename NumType, typename MatValueType>
Mat<MatValueType> operator-(const Mat<MatValueType>& mat, NumType num) {
    return do_reverse_binary_operate(num, mat, std::minus<MatValueType>());
}

template<typename NumType, typename MatValueType>
Mat<MatValueType> operator*(const Mat<MatValueType>& mat, NumType num) {
    return do_reverse_binary_operate(num, mat, std::multiplies<MatValueType>());
}

template<typename NumType, typename MatValueType>
Mat<MatValueType> operator/(const Mat<MatValueType>& mat, NumType num) {
    return do_reverse_binary_operate(num, mat, std::divides<MatValueType>());
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const Mat<T>& mat) {
    os << "matrix:\nvalue_type: " 
       << typeid(typename Mat<T>::value_type).name()
       << "\nsize:       " 
       << mat.row_size() 
       << " x " 
       << mat.col_size() 
       << '\n';
    os << "[";
    for (int i = 0; i < mat.row_size(); ++i) {
        os << (i == 0 ? "" : " ") << "[";
        for (int j = 0; j < mat.col_size(); ++j) {
            os << mat.vec_[i][j] << (j == mat.col_size() - 1 ? "" : ", ");
        }
        os << "]" << (i == mat.row_size() - 1 ? "" : "\n");
    }
    os << "]\n";
    return os;
}

template<typename T>
std::string to_string(const Mat<T>& mat) {
    std::string ret;
    ret.reserve(mat.row_size() * mat.col_size());
    ret += "[";
    for (int i = 0; i < mat.row_size(); ++i) {
        ret += (i == 0 ? "" : " ") + std::string("[ ");
        for (int j = 0; j < mat.col_size(); ++j) {
            ret += std::to_string(mat.vec_[i][j]) + std::string(", ");
        }
        ret += "]" + std::string(i == mat.row_size() - 1 ? "" : "\n");
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
template<typename T>
std::tuple<Mat<T>, Mat<T>> LU_decomposition(const Mat<T>& mat) {
    auto a = mat.clone();
    auto n = a.row_size();
    Mat<T> l(a.row_size(), a.col_size());
    Mat<T> u(a.row_size(), a.col_size());
    for (int i = 0; i < l.row_size(); ++i) {
        for (int j = 0; j < l.col_size(); ++j) {
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
template<typename T>
std::tuple<Mat<T>, Mat<T>> LUP_decomposition(const Mat<T>& mat) {
    auto a = mat.clone();
    auto n = a.row_size();
    Mat<T> pi(mat.row_size(), 1);
    for (int i = 0; i < pi.row_size(); ++i) {
        pi[i][0] = static_cast<T>(i);
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
            std::swap(a.vec_[k][i], a.vec_[k2][i]);
        }
        for (int i = k + 1; i < n; ++i) {
            // col-vector
            a.vec_[i][k] = a.vec_[i][k] / a.vec_[k][k];
            for (int j = k + 1; j < n; ++j) {
                a.vec_[i][j] = a.vec_[i][j] - a.vec_[i][k] * a.vec_[k][j];
            }
        }
    }

    return std::make_tuple(pi, a);
}

template<typename T>
Mat<T> lup_solve(Mat<T>& l, Mat<T>& u, Mat<T>& pi, Mat<T>& b) {
    auto n = l.row_size();
    Mat<T> x(n, 1);
    Mat<T> y(n, 1);
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