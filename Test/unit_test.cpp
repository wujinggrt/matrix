#include "unit_test.h"
#include "../matrix/matrix.hpp"

void TestBasicArithmetic() {
    wj::Mat<int> vector_3by1{{1}, {2}, {3}};
    wj::Mat<int> matrix_3by3{
                             {1, 2, 3},
                             {4, 5, 6},
                             {7, 8, 9}};
    printf("我是a[0][0]:%i\n 我是a(1, 1):%i\n", matrix_3by3[0][0], matrix_3by3(1, 1));
    wj::Mat<int> result = matrix_3by3 * vector_3by1;
    std::cout << "乘法*:\n";
    result.Print();
    wj::Mat<int> a = 5.0 * matrix_3by3;
    wj::Mat<int> ab = 5.0 * matrix_3by3;
    wj::Mat<int> b = matrix_3by3 * 5.0;
    std::cout << "a:\n";
    a.Print();
    wj::Mat<int> c = a + b + a;
    std::cout << "加法+:\n";
    c.Print();

    auto d = c - a;
    std::cout << "减法-:\n";
    d.Print();

    auto e = c / a;
    std::cout << "除法/:\n";
    e.Print();

    auto f = a.DotProduct(e);
    std::cout << "DotProduct:\n";
    f.Print();
    
    auto g = 1000 + a;
    std::cout << "g + a:\n";
    g.Print();

    auto h = a + 1000;
    std::cout << "a + 1000:\n";
    h.Print();

    auto ii = h - 1000;
    std::cout << "h - 1000:\n";
    ii.Print();

    auto kk = a / 2;
    std::cout << "a / 2:\n";
    kk.Print();

    auto l = 2 / a;
    std::cout << "2 / a:\n";
    l.Print();
}

TEST_CASE(TestBasicArithmetic) {
    TestBasicArithmetic();
}

void TestDouble() {
    wj::Mat<double> m{
                      {1., 2., 3.},
                      {4., 5., 6.},
                      {7., 8., 9.}};   
    auto a = 2 * m;
    std::cout << "2 * a:\n";
    a.Print();

    auto b = 2 / m;
    std::cout << "2 / m:\n";
    b.Print();
    std::cout << m << '\n';
    std::cout << m.Transpose() << std::endl;

    std::string s = to_string(m);
    std::cout << s << std::endl;
}

TEST_CASE(TestDouble) {
    TestDouble();
}

void TestLU() {
    wj::Mat<double> m{
                      {2, 3, 1, 5},
                      {6, 13, 5, 19},
                      {2, 19, 10, 23},
                      {4, 10, 11, 31}};
    auto rr = wj::LUPDecomposition(m);
    std::cout << std::get<0>(rr) << std::get<1>(rr);
    
    auto r = wj::LUDecomposition(m);
    std::cout << std::get<0>(r);
    std::cout << std::get<1>(r);
}

TEST_CASE(TestLU) {
    TestLU();
}

void TestLUP() {
    wj::Mat<double> mm{
                       {2, 0, 2, 0.6},
                       {3, 3, 4, -2},
                       {5, 5, 4, 2},
                       {-1, -2, 3.4, -1}};
    auto rr = wj::LUPDecomposition(mm);
    std::cout << std::get<0>(rr);
    std::cout << std::get<1>(rr);
}

TEST_CASE(TestLUP) {
    TestLUP();
}


// Terminal will output below:
// matrix:
// value_type: d
// size:       3 x 1
// [[2]
//  [0]
//  [1]]
// matrix:
// value_type: d
// size:       3 x 3
// [[5, 6, 3]
//  [0.2, 0.8, -0.6]
//  [0.6, 0.5, 2.5]]
// matrix:
// value_type: d
// size:       3 x 3
// [[1, 0, 0]
//  [0.2, 1, 0]
//  [0.6, 0.5, 1]]
// matrix:
// value_type: d
// size:       3 x 3
// [[5, 6, 3]
//  [0, 0.8, -0.6]
//  [0, 0, 2.5]]
// matrix:
// value_type: d
// size:       3 x 1
// [[-1.4]
//  [2.2]
//  [0.6]]
//
void TestLUPSolve() {
    wj::Mat<double> a{
                      {1, 2, 0},
                      {3, 4, 4},
                      {5, 6, 3}};
    wj::Mat<double> b{{3}, {7}, {8}};
    auto r = wj::LUPDecomposition(a);
    std::cout << std::get<0>(r) << std::get<1>(r);
    wj::Matd l(a.RowSize(), a.ColSize());
    wj::Matd u(a.RowSize(), a.ColSize());
    wj::Matd pi = std::get<0>(r);
    for (int i = 0; i < a.RowSize(); ++i) {
        for (int j = 0; j < l.ColSize(); ++j) {
            if (i == j) {
                l[i][j] = 1.;
                u[i][j] = std::get<1>(r)[i][j];
            } else if (i < j) {
                u[i][j] = std::get<1>(r)[i][j];
            } else {
                l[i][j] = std::get<1>(r)[i][j];
            }
        }
    }
    auto x = wj::LUPSolve(l, u, pi, b);
    std::cout << l << u << x;
}

TEST_CASE(TestLUPSolve) {
    TestLUPSolve();
}

void TestInverse() {
    std::cout << wj::Mat<int>::Eye(1);
    wj::Matd m{
               {1, 2, 0},
               {3, 4, 4},
               {5, 6, 3}};
    std::cout << m.Inverse();
}

TEST_CASE(TestInverse) {
    TestInverse();
}

void TestDotProduct() {
    wj::Matd a{
               {1, 2, 3},
               {4, 5, 6},
               {7, 8, 9}};
    std::cout << a;
    wj::Mat<int> b(3, 4);
    auto c = wj::Matd::Eye(3);

    wj::Matd d = a + a;
    std::cout << d;

    wj::Matd e = 2 + a;
    std::cout << e;
    wj::Matd f = c * a;
    std::cout << f;

    wj::Matd g = a.DotProduct(d);
    std::cout << g;
}

TEST_CASE(TestDotProduct) {
    TestDotProduct();
}

// Terminal will output below:
// tmp3(inverse):
// matrix:
// value_type: d
// size:       2 x 2
// [[0.222222, -0.111111]
//  [-0.111111, 0.222222]]
// Y_hat :
// matrix:
// value_type: d
// size:       2 x 1
// [[-0.555556]
//  [-0.222222]]
void TestFilterInGeometic() {
    wj::Matd L = {{1, 1}};
    L = L.Transpose();
    wj::Matd mu_Y{{0}, {0}};
    wj::Matd D_Y{
                 {2, 0}, 
                 {0, 2}};
    wj::Matd D_delta{
                     {2, 0}, 
                     {0, 2}};
    wj::Matd D_Y_delta{
                       {0, -1}, 
                       {0, 0}};
    wj::Matd B{
               {-1, -1},
               {-1, 0}};
    // D_YY * B.Transpose() + D_Y_delta
    auto tmp1 = D_Y * B.Transpose() + D_Y_delta;
    auto tmp2 = D_delta + B * D_Y_delta + D_Y_delta.Transpose() * B.Transpose() + B * D_Y * B.Transpose();
    auto tmp3 = tmp2.Inverse();
    auto Y_hat = mu_Y + tmp1 * tmp3 * (L - B * mu_Y);
    std::cout << "tmp3(inverse):\n" << tmp3;
    std::cout << "Y_hat :\n" << Y_hat;

}

TEST_CASE(TestFilterInGeometic) {
    TestFilterInGeometic();
}