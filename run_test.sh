#!/bin/bash

g++ -std=c++17 ./TEST/unit_test.cpp ./TEST/main.cpp
if [ -f ./a.out ];then
    ./a.out > output
fi