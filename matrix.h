#ifndef MATRIX_H__
#define MATRIX_H__

#include <vector>
#include <string>
#include <iostream>
#include <initializer_list>
#include <stdexcept>
#include <tuple>
#include <cmath>
#include <random>

namespace wj
{

    template<typename T = float>
    class Mat;

/*************************************************************
 * Binary operator overload declaration
*************************************************************/
    template<typename T>
    Mat<T> operator+(double n, const Mat<T> &m);

    template<typename T>
    Mat<T> operator-(double n, const Mat<T> &m);

    template<typename T>
    Mat<T> operator*(double n, const Mat<T> &m);

    template<typename T>
    Mat<T> operator/(double n, const Mat<T> &m);

    template<typename T>
    Mat<T> operator+(const Mat<T> &m, double n);

    template<typename T>
    Mat<T> operator-(const Mat<T> &m, double n);

    template<typename T>
    Mat<T> operator*(const Mat<T> &m, double n);

    template<typename T>
    Mat<T> operator/(const Mat<T> &m, double n);

    template<typename T>
    ostream& operator<<(ostream& os, const Mat<T> &m);

    template<typename T>
    std::string to_string(const Mat<T> &m);

/*************************************************************
 * Binary operator overload declaration
 * end
*************************************************************/

    template<typename T>
    std::tuple<Mat<T>, Mat<T>> LU_decomposition(const Mat<T> &m);
    
    template<typename T>
    std::tuple<Mat<T>, Mat<T>> LUP_decomposition(const Mat<T> &m);

    template<typename T>
    Mat<T> LUP_solve(Mat<T> &l, Mat<T> &u, Mat<T> &pi, Mat<T> &b);

/*************************************************************
 * template class: Mat
*************************************************************/
    template<typename T>
    class Mat
    {
    private:
        std::vector<std::vector<T>> vec_;

    public:
        using iterator = decltype(vec_.begin());
        using value_type = T;
        // using iterator = typename std::vector<std::vector<T>>::iterator;

        friend Mat<T> operator+<>(double n, const Mat<T> &m);
        friend Mat<T> operator-<>(double n, const Mat<T> &m);
        friend Mat<T> operator*<>(double n, const Mat<T> &m);
        friend Mat<T> operator/<>(double n, const Mat<T> &m);

        friend Mat<T> operator+<>(const Mat<T> &m, double n);
        friend Mat<T> operator-<>(const Mat<T> &m, double n);
        friend Mat<T> operator*<>(const Mat<T> &m, double n);
        friend Mat<T> operator/<>(const Mat<T> &m, double n);

        friend ostream& operator<<<>(ostream& os, const Mat<T> &m);
        friend std::string to_string<>(const Mat<T> &m);

        friend std::tuple<Mat<T>, Mat<T>> LU_decomposition<>(const Mat<T> &m);
        friend std::tuple<Mat<T>, Mat<T>> LUP_decomposition<>(const Mat<T> &m);

    public:
        Mat(std::initializer_list<std::vector<T>> ls)
            : vec_{ls}
        {
        }

        Mat(std::size_t rows = 1, std::size_t cols = 1)
            : vec_{}
        {
            if (!rows || !cols)
            {
                throw invalid_argument("argument can't be 0!");
            }
            vec_.reserve(rows);
            for (auto i = 0; i < rows; ++i)
            {
                vec_.push_back(std::vector<T>(cols));
            }
        }

        static Mat<T> eye(std::size_t row)
        {
            Mat<T> ret(row, row);
            for (int i = 0; i < row; ++i)
            {
                ret[i][i] = 1;
            }
            return ret;
        }

        static Mat<T> random(T low, T high, std::size_t row = 1, std::size_t col = 1)
        {
            Mat<T> ret(row, col);
            std::default_random_engine e;
            if (typeid(T) == typeid(double))
            {
                std::uniform_real_distribution<double> ud(low, high);
                for (int i = 0; i < row; ++i)
                {
                    for (int j = 0; j < col; ++j)
                    {
                        ret[i][j] = ud(e);
                    }
                }
            }
            if (typeid(T) == typeid(float))
            {
                std::uniform_real_distribution<float> uf(low, high);
                for (int i = 0; i < row; ++i)
                {
                    for (int j = 0; j < col; ++j)
                    {
                        ret[i][j] = uf(e);
                    }
                }
            }
            if (typeid(T) == typeid(int))
            {
                std::uniform_int_distribution<int> ui(low, high);
                for (int i = 0; i < row; ++i)
                {
                    for (int j = 0; j < col; ++j)
                    {
                        ret[i][j] = ui(e);
                    }
                }
            }

            return ret;
        }

        void print() const
        {
            for (auto e: vec_)
            {
                for (auto ee: e)
                {
                    std::cout << ee << ' ';
                }
                std::cout << '\n';
            }
            std::cout << '\n';
        }

        // throw range
        decltype(vec_[0])& operator[](std::size_t i)
        {
            if (i >= vec_.size())
            {
                throw out_of_range("Mat:vec_:index out of range!");
            }
            return vec_[i];
        }

        iterator begin()
        {
            vec_.begin();
        }

        iterator end()
        {
            vec_.end();
        }

        std::size_t row_size() const 
        {
            return vec_.size();
        }

        std::size_t col_size() const 
        {
            return vec_[0].size();
        }

#define THROW_EXCEPTION_FOR_INCOMPATIBLE_MATRIX(other) \
    if (row_size() != other.row_size() || col_size() != other.col_size()) \
    { \
        throw invalid_argument(std::string("incompatible dimensions\n") +  \
            "right-matrix rows, cols:" + std::to_string(row_size()) + ", " + std::to_string(col_size()) + \
            "\nleft-matrx rows, cols:" + std::to_string(other.row_size()) + ", " + std::to_string(other.col_size()) \
            ); \
    }

        // +,-,*,/操作在行、列不匹配的时候则会抛出exception:invalid_argument
        Mat<T> operator+(const Mat<T> &other) const
        {
            THROW_EXCEPTION_FOR_INCOMPATIBLE_MATRIX(other);

            Mat<T> ret(row_size(), col_size());
            for (int i = 0; i < row_size(); ++i)
            {
                for (int j = 0; j < col_size(); ++j)
                {
                    ret[i][j] = vec_[i][j] + other.vec_[i][j];
                }
            }
            return ret;
        }

        Mat<T> operator-(const Mat<T> &other) const
        {
            THROW_EXCEPTION_FOR_INCOMPATIBLE_MATRIX(other);

            Mat<T> ret(row_size(), col_size());
            for (int i = 0; i < row_size(); ++i)
            {
                for (int j = 0; j < col_size(); ++j)
                {
                    ret[i][j] = vec_[i][j] - other.vec_[i][j];
                }
            }
            return ret;
        }

        Mat<T> operator*(const Mat<T> &other) const
        {
            if (col_size() != other.row_size())
            {
                throw invalid_argument(std::string("incompatible dimensions\n") + 
                    "right-matrix rows, cols:" + std::to_string(row_size()) + ", " + std::to_string(col_size()) +
                    "\nleft-matrx rows, cols:" + std::to_string(other.row_size()) + ", " + std::to_string(other.col_size())
                    );
            }

            Mat<T> ret(vec_.size(), other.col_size());
            // this rows
            for (int i = 0; i < row_size(); ++i)
            {
                // other cols
                for (int j = 0; j < other.col_size(); ++j)
                {
                    ret[i][j] = 0.;
                    for (int k = 0; k < col_size(); ++k)
                    {
                        ret[i][j] += vec_[i][k] * other.vec_[k][j];
                    }
                }
            }
            return ret;
        }

        Mat<T> operator/(const Mat<T> &other) const
        {
            THROW_EXCEPTION_FOR_INCOMPATIBLE_MATRIX(other);

            Mat<T> ret(row_size(), col_size());
            for (int i = 0; i < row_size(); ++i)
            {
                for (int j = 0; j < col_size(); ++j)
                {
                    ret[i][j] = vec_[i][j] / other.vec_[i][j];
                }
            }
            return ret;
        }

        Mat<T> dot_product(const Mat<T> &other) const
        {
            THROW_EXCEPTION_FOR_INCOMPATIBLE_MATRIX(other);

            Mat<T> ret(row_size(), col_size());
            for (int i = 0; i < row_size(); ++i)
            {
                for (int j = 0; j < col_size(); ++j)
                {
                    ret[i][j] = vec_[i][j] * other.vec_[i][j];
                }
            }
            return ret;
        }

        void print_type() const
        {
            std::cout << typeid(T).name() << '\n'; 
        }

        Mat<T> trans() const 
        {
            Mat<T> ret(col_size(), row_size());
            for (int i = 0; i < row_size(); ++i)
            {
                for (int j = 0; j < col_size(); ++j)
                {
                    ret[j][i] = vec_[i][j];
                }
            }
            return ret;
        }

        Mat<T> clone() const
        {
            Mat<T> ret(row_size(), col_size());
            for (int i = 0; i < row_size(); ++i)
            {
                for (int j = 0; j < col_size(); ++j)
                {
                    ret.vec_[i][j] = vec_[i][j];
                }
            }

            return ret;
        }

        Mat<T> inv() const
        {
            Mat<T> ret(row_size(), row_size());
            auto r = LUP_decomposition(*this);
            Mat<T> l(row_size(), col_size());
            Mat<T> u(row_size(), col_size());
            Mat<T> pi = get<0>(r);
            Mat<T> I = Mat<T>::eye(row_size());
            for (int k = 0; k < row_size(); ++k)
            {
                Mat<T> b(row_size(), 1);
                for (int i = 0; i < row_size(); ++i)
                {
                    b[i][0] = I[i][k];
                }
                for (int i = 0; i < row_size(); ++i)
                {
                    for (int j = 0; j < l.col_size(); ++j)
                    {
                        if (i == j)
                        {
                            l[i][j] = 1.;
                            u[i][j] = get<1>(r)[i][j];
                        }
                        else if (i < j)
                        {
                            u[i][j] = get<1>(r)[i][j];
                        }
                        else
                        {
                            l[i][j] = get<1>(r)[i][j];
                        }
                    }
                }
                auto x = wj::LUP_solve(l, u, pi, b);
                for (int i = 0; i < row_size(); ++i)
                {
                    ret[i][k] = x[i][0];
                }
            }
            
            return ret;
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
    Mat<T> ret(mat.row_size(), mat.col_size()); \
    for (int i = 0; i < mat.row_size(); ++i) \
    { \
        for (int j = 0; j < mat.col_size(); ++j) \
        { \
            ret[i][j] = n op mat.vec_[i][j]; \
        } \
    } \
    return ret;

    template<typename T>
    Mat<T> operator+(double n, const Mat<T> &m)
    {
        BINARY_OPERATOR_OVERLOAD_RETURN(n, +, m)
    }
    
    template<typename T>
    Mat<T> operator-(double n, const Mat<T> &m)
    {
        BINARY_OPERATOR_OVERLOAD_RETURN(n, -, m)
    }

    template<typename T>
    Mat<T> operator*(double n, const Mat<T> &m)
    {
        BINARY_OPERATOR_OVERLOAD_RETURN(n, *, m)
    }
    
    template<typename T>
    Mat<T> operator/(double n, const Mat<T> &m)
    {
        BINARY_OPERATOR_OVERLOAD_RETURN(n, /, m)
    }
    
#define REVERSE_BINARY_OPERATOR_OVERLOAD_RETURN(mat, op, n) \
    Mat<T> ret(mat.row_size(), mat.col_size()); \
    for (int i = 0; i < mat.row_size(); ++i) \
    { \
        for (int j = 0; j < mat.col_size(); ++j) \
        { \
            ret[i][j] = mat.vec_[i][j] op n; \
        } \
    } \
    return ret;

    template<typename T>
    Mat<T> operator+(const Mat<T> &m, double n)
    {
        REVERSE_BINARY_OPERATOR_OVERLOAD_RETURN(m, +, n)
    }

    template<typename T>
    Mat<T> operator-(const Mat<T> &m, double n)
    {
        REVERSE_BINARY_OPERATOR_OVERLOAD_RETURN(m, -, n)
    }

    template<typename T>
    Mat<T> operator*(const Mat<T> &m, double n)
    {
        REVERSE_BINARY_OPERATOR_OVERLOAD_RETURN(m, *, n)
    }

    template<typename T>
    Mat<T> operator/(const Mat<T> &m, double n)
    {
        REVERSE_BINARY_OPERATOR_OVERLOAD_RETURN(m, /, n)
    }
    
    template<typename T>
    ostream& operator<<(ostream& os, const Mat<T> &m)
    {
        os << "matrix:\nvalue_type: " << (typeid(typename Mat<T>::value_type).name()) 
            << "\nsize:       " << m.row_size() << " x " << m.col_size() 
            << '\n';
        os << "[";
        for (int i = 0; i < m.row_size(); ++i)
        {
            os << (i == 0 ? "" : " ") << "[";
            for (int j = 0; j < m.col_size(); ++j)
            {
                os << m.vec_[i][j] << (j == m.col_size() - 1 ? "" : ", ");
            }
            os << "]" << (i == m.row_size() - 1 ? "" : "\n");
        }
        os << "]\n";
        return os;
    }

    template<typename T>
    std::string to_string(const Mat<T> &m)
    {
        std::string ret;
        ret += "[";
        for (int i = 0; i < m.row_size(); ++i)
        {
            ret += (i == 0 ? "" : " ") + std::string("[ ");
            for (int j = 0; j < m.col_size(); ++j)
            {
                ret += std::to_string(m.vec_[i][j]) + std::string(", ");
            }
            ret += "]" + std::string(i == m.row_size() - 1 ? "" : "\n");
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
    std::tuple<Mat<T>, Mat<T>> LU_decomposition(const Mat<T> &m)
    {
        auto a = m.clone();
        auto n = a.row_size();
        Mat<T> l(a.row_size(), a.col_size());
        Mat<T> u(a.row_size(), a.col_size());
        for (int i = 0; i < l.row_size(); ++i)
        {
            for (int j = 0; j < l.col_size(); ++j)
            {
                if (i == j)
                {
                    l[i][j] = 1;
                }
            }
        }

        for (int k = 0; k < n; ++k)
        {
            u[k][k] = a[k][k];
            for (int i = k + 1; i < n; ++i)
            {
                // col-vector
                l[i][k] = a[i][k] / u[k][k];
                // row-vector
                u[k][i] = a[k][i];
            }
            for (int i = k + 1; i < n; ++i)
            {
                for (int j = k + 1; j < n; ++j)
                {
                    a[i][j] = a[i][j] - l[i][k] * u[k][j];
                }
            }
        }
        return std::make_tuple(l, u);
    }

    // 奇异矩阵的话会抛出exception:invalid_argument
    template<typename T>
    std::tuple<Mat<T>, Mat<T>> LUP_decomposition(const Mat<T> &m)
    {
        auto a = m.clone();
        auto n = a.row_size();
        Mat<T> pi(m.row_size(), 1);
        for (int i = 0; i < pi.row_size(); ++i)
        {
            pi[i][0] = static_cast<T>(i);
        }
        for (int k = 0; k < n; ++k)
        {
            auto p = 0.;
            auto k2 = k;
            // find the max absolute value.
            for (int i = k; i < n; ++i)
            {
                if (std::abs(a[i][k]) > p)
                {
                    p = abs(a[i][k]);
                    k2 = i;
                }
            }
            if (p == 0.)
            {
                throw invalid_argument("singular matrix");
            }
            swap(pi[k][0], pi[k2][0]);
            // swap rows
            for (int i = 0; i < n; ++i)
            {
                swap(a.vec_[k][i], a.vec_[k2][i]);
            }
            for (int i = k + 1; i < n; ++i)
            {
                // col-vector
                a.vec_[i][k] = a.vec_[i][k] / a.vec_[k][k];
                for (int j = k + 1; j < n; ++j)
                {
                    a.vec_[i][j] = a.vec_[i][j] - a.vec_[i][k] * a.vec_[k][j];
                }
            }
        }

        return std::make_tuple(pi, a);
    }

    template<typename T>
    Mat<T> LUP_solve(Mat<T> &l, Mat<T> &u, Mat<T> &pi, Mat<T> &b)
    {
        auto n = l.row_size();
        Mat<T> x(n, 1);
        Mat<T> y(n, 1);
        for (int i = 0; i < n; ++i)
        {
            auto sum = 0.;
            for (int j = 0; j < i; ++j)
            {
                sum += l[i][j] * y[j][0];
            }
            y[i][0] = b[pi[i][0]][0] - sum;
        }
        for (int i = n - 1; i >= 0; --i)
        {
            auto sum = 0.;
            for (int j = i + 1; j < n; ++j)
            {
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