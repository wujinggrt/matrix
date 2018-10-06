# 矩阵库
这是一个使用了C++标准库写出来的矩阵相关操作的库，能够完成矩阵求逆、转置、基本的算数运算。并且重载了大部分操作符，所以使用起来比较方便。

缺点是在矩阵乘法使用的是原生的算法，还没有优化，时间复杂度应该是O(n^3)。

## 已完成的功能

>v1.0

|算数运算|跳转|
|---|---|
|基本类型|[1. 基本类型](#1.)|
|加法|[2. 加法](#2.)|
|减法、除法、(与常数的)乘法|[3. 减法、除法、(与常数的)乘法同样](#3.)|
|矩阵乘法|[4. 矩阵乘法](#4.)
|转置和求逆|[5. 转置和求逆](#5.)|
|随机数矩阵|[6. 随机数矩阵](#6.)|
|例子|[7. 滤波章节例子3.1](#7.)|

## 使用方法

在使用的CPP文件目录下， 拷贝matrix.h，在使用的cpp```#include "matrix.h"```

<h3 id="1.">1. 基本类型</h3>

矩阵的基础类型为```wj::Mat<T>```, ```T```为```double, float, int```，分别对于类型```Matd, Matf, Mati```。

```C++
using Matd = wj::Mat<double>;
using Matf = wj::Mat<float>;
using Mati = wj::Mat<int>;
```

#### 1.1构造函数

```C++
// 一个3x3的矩阵
wj::Matd a{{1, 2, 3},
       {4, 5, 6},
       {7, 8, 9}};
// 一个3x4的矩阵，默认为元素全为0的矩阵
wj::Mat<int> b(3, 4);
```

如果不想在每一个对象/函数前用```wj::```，可以在```#include "matrix.h"```之后使用```using namespace wj;```。（以下例子默认使用using...）

#### 1.2 单位矩阵

```C++
// c 是一个3x3的单位矩阵
Matd c = Matd::eye(3);
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

<h3 id="2.">2. 加法</h3>

如果矩阵之间的操作，行列不匹配则会抛出异常(invalid_argument)。

#### 下标索引

```C++
// d 是第一行第二个元素
double d = a[0][1];
```

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

<h3 id="3.">3. 减法、除法、(与常数的)乘法同样</h3>

<h3 id="4.">4. 矩阵乘法</h3>

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

<h3 id="5.">5. 转置和求逆</h3>

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

 <h3 id="6.">6. 随机数矩阵</h3>

随机数矩阵使用```Mat<T>```的静态函数```Mat<T>::random(low, high, row, col)```来生成。是随机数均匀分布，生成的类型是```T```。

例子：</br>
```C++
Matd m = wj::Matd::random(-100., 100., 5, 5);
cout << m;
/*
matrix:
value_type: d
size:       5 x 5
[[-73.6924, -8.26997, -56.2082, 35.7729, 86.9386]
 [3.88327, -93.0856, 5.94004, -98.4604, -86.6316]
 [37.3545, 86.0873, 5.38576, 30.7838, 40.2381]
 [52.4396, -90.5071, -34.3532, 51.2821, -26.9323]
 [96.5101, 50.6712, -85.4628, 76.9414, -12.7177]]
 */

Matd mm = wj::Mat<int>::random(-100., 100., 5, 5);
cout << mm;
/*
matrix:
value_type: i
size:       5 x 5
[[-100, -74, 51, -8, 7]
 [-56, -91, 36, 36, 87]
 [-23, 4, 67, -94, -90]
 [6, 34, -99, -23, -87]
 [-17, 38, 18, 87, 70]]
 */
```

<h3 id="7.">7. 滤波章节例子3.1</h3>

```C++
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
Matd tmp1 = D_Y * B.trans() + D_Y_delta;
Matd tmp2 = D_delta + B * D_Y_delta + D_Y_delta.trans() * B.trans() + B * D_Y * B.trans();
Matd tmp3 = tmp2.inv();
Matd Y_hat = mu_Y + tmp1 * tmp3 * (L - B * mu_Y);
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
```