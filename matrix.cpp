
using namespace std;

#include "matrix.h"

template<typename T>
void print_test(initializer_list<initializer_list<T>> ls)
{
    for (auto e: ls)
    {
        for (auto ee: e)
        {
            cout << ee << ' ';
        }
        cout << '\n';
    }
    cout << '\n';
}

int main()
{
    wj::Mat<int> m{{1}, {2}, {3}};
    wj::Mat<int> mm{{1, 2, 3},
                {4, 5, 6},
                {7, 8, 9}};
    auto mmm = mm * m;
    mmm.print();
    auto a = 5.0 * mm;
    auto b = mm * 5.0;
    b.print();
    auto c = a + b;
    c.print();

    return 0;
}