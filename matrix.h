#ifndef MATRIX_H__
#define MATRIX_H__

template<typename T = float>
class Mat
{
private:
    std::vector<std::vector<T>> vec_;

public:
    using iterator = decltype(vec_.begin());
    // using iterator = typename std::vector<std::vector<T>>::iterator;
    
public:
    Mat(std::initializer_list<std::vector<T>> ls)
        : vec_{ls}
    {
    }

    void print()
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

    // throw range_error
    decltype(vec_[0])& operator[](std::size_t i)
    {
        if (i >= vec_.size())
        {
            throw range_error("Mat:vec_:index out of range!");
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

};

#endif