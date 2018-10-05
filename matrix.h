#ifndef MATRIX_H__
#define MATRIX_H__

#include <vector>
#include <string>
#include <iostream>
#include <initializer_list>
#include <stdexcept>

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

        Mat<T> inv() const;
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
            os << (i == 0 ? "" : " ") << "[ ";
            for (int j = 0; j < m.col_size(); ++j)
            {
                os << m.vec_[i][j] << ", ";
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

}

#endif