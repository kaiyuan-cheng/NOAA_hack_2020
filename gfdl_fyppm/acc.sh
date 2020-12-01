#!/bin/bash

export PGI_ACC_TIME=1

thread=2

srun -N 1 --ntasks=$thread -p batch -t 00:05:00 --gres=gpu:$thread ./exec.acc/fyppm.x
