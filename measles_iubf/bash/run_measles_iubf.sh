#!/bin/bash
# nohup bash run_measles_iubf.sh run_level &
export LC_ALL=en_US.UTF-8
Rscript ../analyses/iubf_run.R $1
