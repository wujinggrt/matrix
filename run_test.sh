#!/bin/bash

clang++ -std=c++17 ./test/test.cpp ./test/main.cpp
./a.out > ./output.txt