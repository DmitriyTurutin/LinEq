#!/bin/bash


clang-format -style='{IndentWidth: 4}' -i *.c
clang -Wall -o cg cg.c 
