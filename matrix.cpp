
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

void test()
{
        wj::Mat<int> m{{1}, {2}, {3}};
    wj::Mat<int> mm{{1, 2, 3},
                {4, 5, 6},
                {7, 8, 9}};
    auto mmm = mm * m;
    cout << "*:\n";
    mmm.print();
    auto a = 5.0 * mm;
    auto b = mm * 5.0;
    cout << "a:\n";
    a.print();
    auto c = a + b;
    cout << "+:c:\n";
    c.print();

    auto d = c - a;
    cout << "-:\n";
    d.print();

    auto e = c / a;
    cout << "/:\n";
    e.print();

    auto f = a.dot_product(e);
    cout << "dot_product:\n";
    f.print();
    
    auto g = 1000 + a;
    cout << "g + a:\n";
    g.print();

    auto h = a + 1000;
    cout << "a + 1000:\n";
    h.print();

    auto ii = h - 1000;
    cout << "h - 1000:\n";
    ii.print();

    auto kk = a / 2;
    cout << "a / 2:\n";
    kk.print();

    auto l = 2 / a;
    cout << "2 / a:\n";
    l.print();
}

void test_double()
{
    wj::Mat<double> m{{1., 2., 3.},
                {4., 5., 6.},
                {7., 8., 9.}};   
    auto a = 2 * m;
    cout << "2 * a:\n ";
    a.print();

    auto b = 2 / m;
    cout << "2 / m:\n";
    b.print();
    cout << m << '\n';
}

int main()
{
    test_double();
    
    return 0;
}