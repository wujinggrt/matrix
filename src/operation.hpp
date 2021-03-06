#ifndef OPERATION_H__
#define OPERATION_H__

#include <stack>
#include <variant>
#include <memory>

#include "mat_type.hpp"

namespace wj{

/*************************************************************
 * binary operator overload implementation
*************************************************************/

// 使用 SFINAE 来限制 num 必须是 arithmetic, 
// 也就是能够与矩阵中元素进行四则运算的。
// 这个函数主要用来完成 num 加上这个 Mat 中的所有元素。
// 然后给Mat中的元素与他进行进行op(num, each src element)
template<typename NumType, typename MatValueType, typename BinaryOperation>
inline
Mat<MatValueType> DoBinaryOperate(
                                  NumType num,
                                  const Mat<MatValueType>& src,
                                  BinaryOperation op,
                                  std::enable_if_t<
                                      std::is_arithmetic<NumType>::value
                                      >* = nullptr) {
    Mat<MatValueType> ret(src.RowSize(), src.ColSize());
    MatValueType casted_num;
    if constexpr (std::is_same_v<NumType, MatValueType>) {
        casted_num = num;
    } else {
        casted_num = static_cast<MatValueType>(num);
    }
    for (std::size_t i = 0; i < src.RowSize(); ++i) {
        for (std::size_t j = 0; j < src.ColSize(); ++j) {
            ret(i, j) = op(casted_num, src(i, j));
        }
    }
    return ret;
}

template<typename NumType, typename MatValueType>
inline
Mat<MatValueType> operator+(NumType num, const Mat<MatValueType>& mat) {
    return DoBinaryOperate(num, mat, std::plus<MatValueType>());
}

template<typename NumType, typename MatValueType>
inline
Mat<MatValueType> operator-(NumType num, const Mat<MatValueType>& mat) {
    return DoBinaryOperate(num, mat, std::minus<MatValueType>());
}

template<typename NumType, typename MatValueType>
inline
Mat<MatValueType> operator*(NumType num, const Mat<MatValueType>& mat) {
    return DoBinaryOperate(num, mat, std::multiplies<MatValueType>());
}

template<typename NumType, typename MatValueType>
inline
Mat<MatValueType> operator/(NumType num, const Mat<MatValueType>& mat) {
    return DoBinaryOperate(num, mat, std::divides<MatValueType>());
}

template<typename NumType, typename MatValueType, typename BinaryOperation>
inline
Mat<MatValueType> DoReverseBinaryOperate(
                                         NumType num,
                                         const Mat<MatValueType>& src,
                                         BinaryOperation op,
                                         std::enable_if_t<
                                              std::is_arithmetic<NumType>::value
                                              >* = nullptr) {
    Mat<MatValueType> ret(src.RowSize(), src.ColSize());
    MatValueType casted_num;
    if constexpr (std::is_same_v<NumType, MatValueType>) {
        casted_num = num;
    } else {
        casted_num = static_cast<MatValueType>(num);
    }
    for (std::size_t i = 0; i < src.RowSize(); ++i) {
        for (std::size_t j = 0; j < src.ColSize(); ++j) {
            ret(i, j) = op(src(i, j), casted_num);
        }
    }
    return ret;
}

template<typename NumType, typename MatValueType>
inline
Mat<MatValueType> operator+(const Mat<MatValueType>& mat, NumType num) {
    return DoReverseBinaryOperate(num, mat, std::plus<MatValueType>());
}

template<typename NumType, typename MatValueType>
inline
Mat<MatValueType> operator-(const Mat<MatValueType>& mat, NumType num) {
    return DoReverseBinaryOperate(num, mat, std::minus<MatValueType>());
}

template<typename NumType, typename MatValueType>
inline
Mat<MatValueType> operator*(const Mat<MatValueType>& mat, NumType num) {
    return DoReverseBinaryOperate(num, mat, std::multiplies<MatValueType>());
}

template<typename NumType, typename MatValueType>
inline
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
            os << mat(i, j) << (j == mat.ColSize() - 1 ? "" : ", ");
        }
        os << "]" << (i == mat.RowSize() - 1 ? "" : "\n");
    }
    os << "]\n";
    return os;
}

/*************************************************************
 * binary operator overload implementation
 * end
*************************************************************/

template<typename T>
std::string to_string(const Mat<T>& mat) {
    std::string ret;
    ret.reserve(mat.RowSize() * mat.ColSize());
    ret += "[";
    for (int32_t i = 0; i < mat.RowSize(); ++i) {
        ret += (i == 0 ? "" : " ") + std::string("[ ");
        for (int32_t j = 0; j < mat.ColSize(); ++j) {
            ret += std::to_string(mat(i, j)) + std::string(", ");
        }
        ret += "]" + std::string(i == mat.RowSize() - 1 ? "" : "\n");
    }
    ret += "]\n";
    return ret;
}

/**
 * 使用动态规划来提高矩阵链式乘法的效率
 * 限制传入的类型都是一个Mat<T>的
 * */


/**
 * print the position of parentheisization for viewing the priority of multiplication.
 * 
 * Mati a1(30, 35);
 * Mati a2(35, 15);
 * Mati a3(15, 5);
 * Mati a4(5, 10);
 * Mati a5(10, 20);
 * Mati a6(20, 25);
 * 
 * for such matrix data, this function will print as follow:
 * ((A0(A1A2))((A3A4)A5))
 * as 0 start pos.
 * */ 
template<typename Iter>
auto ApplyMultiplicationInOrder(const std::vector<Iter>& matrices, 
                                const std::vector<std::vector<int32_t>>& total_multiply_operation,
                                const std::vector<std::vector<int32_t>>& parenthesis_pos) {
    // dereference 之后为引用类型，即&, 所以remove掉,
    // 也就是实际的 Mat<T>,对应类型的T
    using actual_mat_type =
        std::remove_reference_t<decltype(*(matrices.front()))>;
    actual_mat_type product;
    // product 是迭代器指向的Mat类型
    // variant中存放的是能够 dereference 为对应Mat类型
    using dereference_to_mat = 
        std::variant<Iter, 
            std::shared_ptr<actual_mat_type>>;
    // stack中保存那些能够被解引用的值，
    // 然后解析乘法顺序，弹出两个操作数相乘，然后添加结果item到stack上
    // 直到最后，得到结果。stack最后的值
    std::stack<dereference_to_mat> matrix_order;
    std::function<void(int32_t, int32_t)> 
        process_parenthesis = [&] (int32_t pos_former, int32_t pos_latter) {
            if (pos_former == pos_latter) {
                matrix_order.push(matrices[pos_former]);
            } else {
                process_parenthesis(pos_former, parenthesis_pos[pos_former][pos_latter]);
                process_parenthesis(parenthesis_pos[pos_former][pos_latter] + 1, pos_latter);
                // pop twice and operate.
                // latter 和 former 都是variant类型
                auto latter = matrix_order.top();
                matrix_order.pop();
                auto former = matrix_order.top();
                matrix_order.pop();
                actual_mat_type* ptr_former = nullptr;
                actual_mat_type* ptr_latter = nullptr;
                std::visit([&] (auto&& arg) {
                    ptr_former = &(*arg);
                }, former);
                std::visit([&] (auto&& arg) {
                    ptr_latter = &(*arg);
                }, latter);
                // 相乘
                auto sp_mat = 
                    std::make_shared<actual_mat_type>((*ptr_former) * (*ptr_latter));
                matrix_order.push(sp_mat);
            }
    };
    process_parenthesis(0, matrices.size() - 1);
    // 返回实际的矩阵，需要处理内存泄漏。
    // 不过也不用担心，因为没有使用new申请内存。
    // stack里面的内容都是迭代器，迭代器在这个函数是不会失效的，因为是initializer_list的。
    // shared_ptr也不会失效，如果当做value返回的话。
    return std::visit([] (auto&& arg) -> actual_mat_type {
        return *arg;
    }, matrix_order.top());
}

/**
 * @param matrices: 指向initializer_list的迭代器vector
 * */
template<typename Iter>
auto DetermineMatrixChainOrderAndApply(const std::vector<Iter>& matrices) {
    int32_t num_matrix = matrices.size();
    // 一个二维数组，[i][j]即是 i, j 最少需要多少次，矩阵中每一位数乘法运算
    std::vector<std::vector<int32_t>> 
        total_multiply_operation(num_matrix, 
                                 std::vector<int32_t>(num_matrix, 0));
    // 二维数组括号的位置，也就是乘法优先级
    std::vector<std::vector<int32_t>> 
        parenthesis_pos(num_matrix - 1, 
                        std::vector<int32_t>(num_matrix - 1, 0));
    // j 相对 i 的偏移量
    // 保证乘法次序，所以从 offset = 1 开始，即最少两个matrix相乘
    for (int32_t offset = 1; offset < num_matrix; ++offset) {
        // matrix_i 和 matrix_j所能够取到的值
        // 保证 j = i + offset 不会越界
        for (int32_t i = 0; i < num_matrix - offset; ++i) {
            // j = i + offset < n, j 不会越界，且能够取到最后的值
            auto j = i + offset;
            total_multiply_operation[i][j] = INT32_MAX;
            // k 从 i 处开始分割，
            // using m = total_multiply_operation for convenience
            // m[i][j] =
            //          0 if i == j
            //          for k = i to j - 1 // k can reach j - 1
            //              min{m[i][k] + [k + 1][j] + i_rowsize * k_colsize * j_colsize}
            for (auto k = i; k < j; ++k) {
                // 以i为起点，所以从第i个的行size * 第k个的列size * 第j的ColSize
                // 因为左边的ColSize == 右边的RowSize
                // 跳过k部分，一个个的试
                // 最初的m[0][0] == 0
                // m[0][0] + m[1][1] + 30 * 35 *
                auto q = total_multiply_operation[i][k] +
                         total_multiply_operation[k + 1][j] +
                         matrices[i]->RowSize() * 
                         matrices[k]->ColSize() * 
                         matrices[j]->ColSize();
                if (q < total_multiply_operation[i][j]) {
                    total_multiply_operation[i][j] = q;
                    parenthesis_pos[i][j] = k;
                }
            }
        }
    }
    return ApplyMultiplicationInOrder(matrices, total_multiply_operation, parenthesis_pos);
}

template<typename MatValueType>
std::enable_if_t<std::is_arithmetic<MatValueType>::value, Mat<MatValueType>>
AuxOptimizedChainMultiply(std::initializer_list<Mat<MatValueType>> il) {
    // 这个vector里面包含了指向各个Mat迭代器, 使用花括号初始化
    std::vector<decltype(std::cbegin(il))> matrices{};
    matrices.reserve(il.size());
    for (auto citer = std::cbegin(il); citer != std::cend(il); ++citer) {
        matrices.push_back(citer);
    }
    return DetermineMatrixChainOrderAndApply(matrices);
}

/**
 * Return Mat type
 * */
template<typename... Mats>
auto OptimizedChainMultiply(Mats&&... mats) {
    // 包扩展到il上，然后转发
    return AuxOptimizedChainMultiply({std::forward<Mats>(mats)...});
}

}

#endif