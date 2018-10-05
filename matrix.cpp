#include <vector>
#include <string>
#include <iostream>
#include <initializer_list>

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
    Mat<int> m{{1}, {2}, {3}};
    for (auto e: m)
    {
        for (auto ee: e)
        {
            cout << ee << '\n';
        }
        cout << '\n';
    }
    return 0;
}