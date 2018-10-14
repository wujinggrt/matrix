

#include "matrix.h"
using namespace std;

using namespace wj;

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
    auto c = a + b + a;
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
    cout << m.trans() << endl;

    string s = to_string(m);
    cout << s << endl;
}

void test_lu()
{
    wj::Mat<double> m{{2, 3, 1, 5},
                    {6, 13, 5, 19},
                    {2, 19, 10, 23},
                    {4, 10, 11, 31}};
    auto rr = wj::LUP_decomposition(m);
    cout << get<0>(rr) << get<1>(rr);
    
    auto r = wj::LU_decomposition(m);
    cout << get<0>(r);
    cout << get<1>(r);


}

void test_lup()
{
    
    wj::Mat<double> mm{{2, 0, 2, 0.6},
                        {3, 3, 4, -2},
                        {5, 5, 4, 2},
                        {-1, -2, 3.4, -1}};
    auto rr = wj::LUP_decomposition(mm);
    cout << get<0>(rr);
    cout << get<1>(rr);
}

void test_sln()
{
    wj::Mat<double> a{{1, 2, 0},
                      {3, 4, 4},
                      {5, 6, 3}};
    wj::Mat<double> b{{3}, {7}, {8}};
    auto r = wj::LUP_decomposition(a);
    cout << get<0>(r) << get<1>(r);
    Matd l(a.row_size(), a.col_size());
    Matd u(a.row_size(), a.col_size());
    Matd pi = get<0>(r);
    for (int i = 0; i < a.row_size(); ++i)
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
    cout << l << u << x;
/*
matrix:
value_type: d
size:       3 x 1
[[2]
 [0]
 [1]]
matrix:
value_type: d
size:       3 x 3
[[5, 6, 3]
 [0.2, 0.8, -0.6]
 [0.6, 0.5, 2.5]]
matrix:
value_type: d
size:       3 x 3
[[1, 0, 0]
 [0.2, 1, 0]
 [0.6, 0.5, 1]]
matrix:
value_type: d
size:       3 x 3
[[5, 6, 3]
 [0, 0.8, -0.6]
 [0, 0, 2.5]]
matrix:
value_type: d
size:       3 x 1
[[-1.4]
 [2.2]
 [0.6]]
 */
}

void test_inv()
{
    cout << wj::Mat<int>::eye(1);
    Matd m{{1, 2, 0},
          {3, 4, 4},
          {5, 6, 3}};
    cout << m.inv();
}

void test_readme()
{
    Matd a{{1, 2, 3},
            {4, 5, 6},
            {7, 8, 9}};
    cout << a;
    wj::Mat<int> b(3, 4);
    auto c = Matd::eye(3);

    Matd d = a + a;
    cout << d;

    Matd e = 2 + a;
    cout << e;
    Matd f = c * a;
    cout << f;

    Matd g = a.dot_product(d);
    cout << g;
}

void test_31()
{
    Matd L = {{1, 1}};
    L = L.trans();
    Matd mu_Y{{0}, 
              {0}};
    Matd D_Y{{2, 0}, 
             {0, 2}};
    Matd D_delta{{2, 0}, 
                 {0, 2}};
    Matd D_Y_delta{{0, -1}, 
                   {0, 0}};
    Matd B{{-1, -1},
           {-1, 0}};
    // D_YY * B.trans() + D_Y_delta
    auto tmp1 = D_Y * B.trans() + D_Y_delta;
    auto tmp2 = D_delta + B * D_Y_delta + D_Y_delta.trans() * B.trans() + B * D_Y * B.trans();
    auto tmp3 = tmp2.inv();
    auto Y_hat = mu_Y + tmp1 * tmp3 * (L - B * mu_Y);
    cout << "tmp3(inverse):\n" << tmp3;
    cout << "Y_hat :\n" << Y_hat;
/*
tmp3(inverse):
matrix:
value_type: d
size:       2 x 2
[[0.222222, -0.111111]
 [-0.111111, 0.222222]]
Y_hat :
matrix:
value_type: d
size:       2 x 1
[[-0.555556]
 [-0.222222]]
*/
}

int main()
{
    test_31();
    
    return 0;
}