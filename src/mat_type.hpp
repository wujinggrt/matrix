#ifndef MAT_TYPE_H__
#define MAT_TYPE_H__

#include <cmath>

#include <iterator>
#include <vector>
#include <string>
#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <tuple>
#include <random>
#include <functional>

namespace wj {

template<typename T>
class Mat;

/*************************************************************
 * stream and to_string declaration
*************************************************************/

template<typename MatValueType>
std::ostream& operator<<(std::ostream& os, const Mat<MatValueType>& mat);

template<typename MatValueType>
std::string to_string(const Mat<MatValueType>& mat);

/*************************************************************
 * Sovling linear declaration
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

/**
 * 系统合成的 copy/remove constructor/assignment operator 够用。
 * 目前不用去重写。
 * */
template<typename T>
class Mat {
private:
    std::vector<std::vector<T>> data_;

public:
    using iterator = decltype(data_.begin());
    using value_type = T;
    using size_type = std::size_t;

    template<typename MatValueType>
    friend std::ostream& operator<<(std::ostream& os, const Mat<MatValueType>& mat);
    template<typename MatValueType>
    friend std::string to_string(const Mat<MatValueType>& mat);

public:
    Mat(std::initializer_list<std::vector<T>> ls)
        : data_{ls} {
    }

    Mat(std::size_t rows = 1, std::size_t cols = 1)
        : data_(rows, std::vector<T>(cols)) {
        if (rows == 0 || cols == 0) {
            throw std::invalid_argument("Occured in cunstructor:\nargument can not be 0!");
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
        if constexpr (std::is_floating_point_v<T>) {
            std::uniform_real_distribution<T> num_distribution(low, high);
            for (std::size_t i = 0; i < row; ++i) {
                for (std::size_t j = 0; j < col; ++j) {
                    ret[i][j] = num_distribution(engine);
                }
            }
        } else {
            std::uniform_int_distribution<T> num_distribution(low, high);
            for (std::size_t i = 0; i < row; ++i) {
                for (std::size_t j = 0; j < col; ++j) {
                    ret[i][j] = num_distribution(engine);
                }
            }
        }
        return ret;
    }

    /**
     * output the matrix contents to stdout.
     * */
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

    bool operator==(const Mat<T>& other) {
        return data_ == other.data_;
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
        CheckSize(other, "Occured in addition:\n");

        Mat<T> ret(RowSize(), ColSize());
        for (size_type i = 0; i < RowSize(); ++i) {
            for (int32_t j = 0; j < ColSize(); ++j) {
                ret[i][j] = data_[i][j] + other.data_[i][j];
            }
        }
        return ret;
    }

    Mat<T> operator-(const Mat<T>& other) const {
        CheckSize(other, "Occured in subtraction:\n");

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
            throw std::invalid_argument(std::string("Occured in mupltiplication:\nincompatible dimensions\n") + 
                "right-matrix rows, cols:" + std::to_string(RowSize()) + ", " + std::to_string(ColSize()) +
                "\nleft-matrx rows, cols:" + std::to_string(other.RowSize()) + ", " + std::to_string(other.ColSize())
                );
        }

        Mat<T> ret(data_.size(), other.ColSize());
        // this rows
        for (std::size_t i = 0; i < RowSize(); ++i) {
            // other cols
            for (std::size_t j = 0; j < other.ColSize(); ++j) {
                ret[i][j] = 0.;
                for (std::size_t k = 0; k < ColSize(); ++k) {
                    ret[i][j] += data_[i][k] * other.data_[k][j];
                }
            }
        }
        return ret;
    }

    Mat<T> operator/(const Mat<T>& other) const {
        CheckSize(other, "Occured in division:\n");

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
        CheckSize(other, "Occured in DotProduct:\n");

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
    void CheckSize(const Mat<T>& other, const std::string& place) const {
        if (RowSize() != other.RowSize() || ColSize() != other.ColSize()) {
            throw std::invalid_argument(
                                        place + 
                                        std::string("incompatible dimensions\n") +
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

}

#endif