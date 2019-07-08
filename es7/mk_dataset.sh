#!/bin/bash


c++ dataset.cpp prng.cpp xorshiro.cpp `root-config --glibs --libs --cflags` -o dataset
