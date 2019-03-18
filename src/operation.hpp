#ifndef OPERATION_H__
#define OPERATION_H__

#include "mat_type.hpp"

namespace wj{



/*************************************************************
 * binary operator overload implementation
*************************************************************/

// SFINAE 来限制num必须是arithmetic 
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
    for (int32_t i = 0; i < mat.RowSize(); ++i) {
        os << (i == 0 ? "" : " ") << "[";
        for (int32_t j = 0; j < mat.ColSize(); ++j) {
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
    for (int32_t i = 0; i < mat.RowSize(); ++i) {
        ret += (i == 0 ? "" : " ") + std::string("[ ");
        for (int32_t j = 0; j < mat.ColSize(); ++j) {
            ret += std::to_string(mat.data_[i][j]) + std::string(", ");
        }
        ret += "]" + std::string(i == mat.RowSize() - 1 ? "" : "\n");
    }
    ret += "]\n";
    return ret;
}

/**
 * To do
 * 使用动态规划来提高矩阵链式乘法的效率
 * 限制传入的类型都是一个Mat<T>的
 * */
template<typename T>
std::enable_if_t<std::is_arithmetic<T>::value, Mat<T>>
AuxOptimizedMatrixChainMultiplication(std::initializer_list<Mat<T>> il) {
    auto iter = std::begin(il);
    Mat<T> ret(iter->RowSize(), std::rbegin(il)->ColSize());
    return ret;
}

/**
 * Return Mat type
 * */
template<typename... Args>
auto OptimizedMatrixChainMultiplication(Args&&... args) {
    // 包扩展到il上，然后转发
    return AuxOptimizedMatrixChainMultiplication({std::forward<Args>(args)...});
}


/*************************************************************
 * binary operator overload implementation
 * end
*************************************************************/
}

#endif