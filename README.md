# 矩阵库
这是一个使用了C++标准库写出来的矩阵相关操作的库，能够完成矩阵求逆、转置、基本的算数运算。并且重载了大部分操作符，所以使用起来比较方便。

缺点是在矩阵乘法使用的是原生的算法，还没有优化，时间复杂度应该是O(n^3)。

## 使用