#!/bin/bash

clang++ -std=c++14 ./test/unit_test.cpp ./test/main.cpp
./a.out > output
