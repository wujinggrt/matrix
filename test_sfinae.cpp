#include <iostream>
#include <type_traits>
using namespace std;

struct MyInt {
    int32_t i = 0;
    void increase() {
        ++i;
    }
};

template<typename T>
void foo(T&& t, std::decay_t<decltype(++t)>* = nullptr) {
    ++t;
}

template<typename T>
void foo(T&& t, std::enable_if_t<std::is_same_v<MyInt, remove_reference_t<T>>>* = nullptr) {
    t.increase();
}

int main() {
    int i = 32;
    foo(i);
    cout << i << "\n";
    MyInt my_int;
    foo(my_int);
    cout << my_int.i << "\n";
}