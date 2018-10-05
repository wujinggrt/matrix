# 矩阵库
这是一个使用了C++标准库写出来的矩阵相关操作的库，能够完成矩阵求逆、转置、基本的算数运算。并且重载了大部分操作符，所以使用起来比较方便。

缺点是在矩阵乘法使用的是原生的算法，还没有优化，时间复杂度应该是O(n^3)。

## 使用方法

在使用的CPP文件目录下，仿佛matrix.h文件，并且在使用的文件```#include "matrix.h"```

### 1 基本类型

矩阵的基础类型为```wj::Mat<T>```, ```T```为```double, float, int```，分别对于类型```Matd, Matf, Mati```。

```C++
using Matd = wj::Mat<double>;
using Matf = wj::Mat<float>;
using Mati = wj::Mat<int>;
```

#### 1.1 构造函数

```C++
// 一个3x3的矩阵
Matd a{{1, 2, 3},
       {4, 5, 6},
       {7, 8, 9}};
// 一个3x4的矩阵，默认为元素全为0的矩阵
wj::Mat<int> b(3, 4);
```

#### 1.2 单位矩阵

```C++
// c 是一个3x3的单位矩阵
auto c = Matd::eye(3);
```

#### 1.3 使用流进行(cout)输出

```C++
cout << a;
/* 输出：
matrix:
value_type: d
size:       3 x 3
[[1, 2, 3]
 [4, 5, 6]
 [7, 8, 9]]
 */
```

#### 1.4 clone()

```C++
Matd cl = a.clone()
```

#### 1.5 行列大小

```C++
std::size_t rows = a.row_size();
std::size_t cols = a.col_size();
```

### 2 加法

如果矩阵之间的操作，行列不匹配则会抛出异常(invalid_argument)。

#### 2.1 矩阵之间加法

```C++
Matd d = a + a;
/*
matrix:
value_type: d
size:       3 x 3
[[2, 4, 6]
 [8, 10, 12]
 [14, 16, 18]]
 */
```

#### 2.2 矩阵和常数的加减法

```C++
Matd e = 2 + a;
// e = a + 2; 相同的结果
/*
matrix:
value_type: d
size:       3 x 3
[[3, 4, 5]
 [6, 7, 8]
 [9, 10, 11]]
 */
```

### 3 减法、除法、乘法同样

### 4 矩阵乘法

#### 4.1 矩阵乘法

```C++
Matd f = c * a;
/*
matrix:
value_type: d
size:       3 x 3
[[1, 2, 3]
 [4, 5, 6]
 [7, 8, 9]]
 */
```

#### 4.2 点乘

```C++
// a点乘d
Matd g = a.dot_product(d);
/*
matrix:
value_type: d
size:       3 x 3
[[2, 8, 18]
 [32, 50, 72]
 [98, 128, 162]]
 */
```

### 5 转置和求逆

转置和求逆，以及各种操作都是产生新的矩阵(Mat对象)。

```C++
Matd m{{1, 2, 0},
        {3, 4, 4},
        {5, 6, 3}};

cout << m.trans();
/*
matrix:
value_type: d
size:       3 x 3
[[1, 3, 5]
 [2, 4, 6]
 [0, 4, 3]]
*/

cout << m.inv();
/*
matrix:
value_type: d
size:       3 x 3
[[-1.2, -0.6, 0.8]
 [1.1, 0.3, -0.4]
 [-0.2, 0.4, -0.2]]
 */
 ```