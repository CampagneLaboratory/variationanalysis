#!/usr/bin/env bash
GPU_SPEC=$1
shift
COMMAND="$@"

sem -j 1 --id gpu_${GPU_SPEC} "export CUDA_VISIBLE_DEVICES=${GPU_SPEC} && $COMMAND"