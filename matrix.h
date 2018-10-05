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

    template<typename T>
    Mat<T> operator*(double n, const Mat<T> &m);

    template<typename T>
    Mat<T> operator+(double n, const Mat<T> &m);

    template<typename T>
    class Mat
    {
    private:
        std::vector<std::vector<T>> vec_;

    public:
        using iterator = decltype(vec_.begin());
        // using iterator = typename std::vector<std::vector<T>>::iterator;

        friend Mat<T> operator*<>(double n, const Mat<T> &m);
        friend Mat<T> operator+<>(double n, const Mat<T> &m);
        
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
    if (col_size() != other.row_size()) \
    { \
        throw invalid_argument(std::string("incompatible dimensions\n") +  \
            "right-matrix cols:" + to_string(vec_[0].size()) + \
            "left-matrx rows:" + to_string(other.vec_.size()) \
            ); \
    }

        // 行、列不匹配则会抛出exception:invalid_argument
        Mat<T> operator+(const Mat<T> &other) const
        {
            THROW_EXCEPTION_FOR_INCOMPATIBLE_MATRIX(other);

            Mat<T> ret(row_size(), other.col_size());
            for (int i = 0; i < row_size(); ++i)
            {
                for (int j = 0; j < other.col_size(); ++j)
                {
                    for (int k = 0; k < col_size(); ++k)
                    {
                        ret[i][j] = vec_[i][j] + other.vec_[i][j];
                    }
                }
            }
            return ret;
        }

        // 行、列不匹配则会抛出exception:invalid_argument
        Mat<T> operator*(const Mat<T> &other) const
        {
            THROW_EXCEPTION_FOR_INCOMPATIBLE_MATRIX(other);

            Mat<T> ret(vec_.size(), other.vec_[0].size());
            // this rows
            for (int i = 0; i < vec_.size(); ++i)
            {
                // other cols
                for (int j = 0; j < other.vec_[0].size(); ++j)
                {
                    ret[i][j] = 0.;
                    for (int k = 0; k < vec_[i].size(); ++k)
                    {
                        ret[i][j] += vec_[i][k] * other.vec_[k][j];
                    }
                }
            }
            return ret;
        }

        Mat<T> operator*(const double &n) const
        {
            return n * (*this);
        }

    };

    template<typename T>
    Mat<T> operator*(double n, const Mat<T> &m)
    {
        Mat<T> ret(m.row_size(), m.col_size());
        for (int i = 0; i < m.row_size(); ++i)
        {
            for (int j = 0; j < m.col_size(); ++j)
            {
                ret[i][j] = n * m.vec_[i][j];
            }
        }
        return ret;
    }

}

#endif