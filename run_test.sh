#!/bin/bash

g++ -std=c++17 -O2 ./test/unit_test.cpp ./test/main.cpp -o a.out
if [ -f ./a.out ];then
    ./a.out > output
fi
