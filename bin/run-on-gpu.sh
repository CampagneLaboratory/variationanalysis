#!/usr/bin/env bash
GPU_SPEC=$1
shift
COMMAND="$@"

sem --fg -j 1 --id "gpu_"${USER}"_"${GPU_SPEC} "export CUDA_VISIBLE_DEVICES=${GPU_SPEC} && $COMMAND"