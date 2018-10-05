#include <iostream>
#include <vector>
#include <initializer_list>
using namespace std;

class A
{
public:
    A(initializer_list<int> ls)
        :_vec{ls}
    {}
    vector<int> _vec;
};

int main()
{
    A a{1, 2, 3, 4};
    cout << (a._vec)[2] << endl;
    
    return 0;
}