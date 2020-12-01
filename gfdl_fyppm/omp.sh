#!/bin/bash

thread=2

export OMP_NUM_THREADS=$thread

srun -A hackathon -p batch -N 1 -c $thread --gres=gpu:1 -t 0:30:00  ./exec.omp/fyppm.x
