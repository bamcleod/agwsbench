#!/bin/sh

echo "DFS ramp processing for 3000 frames"
time ./dfs_preprocess
echo ""

echo "DFS fft sum"
./dfs_timing 64 
echo ""

echo "SH48 48x48 subaps 8x8 pix"
./shcentroids 48 8 | grep FFT
echo ""

echo "SH24 24x24 subaps; 12x12 pix"
./shcentroids 24 12 | grep FFT

