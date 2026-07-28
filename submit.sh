#!/bin/bash

time=$1
mem_per_cpu=$2
cpu_per_task=$3
shift 3
other_args="$@"

convert_minutes_to_hms() {
    mins=$1
    hrs=$(( mins / 60 ))
    mins=$(( mins % 60 ))
    printf "%02d:%02d:00\n" "$hrs" "$mins"
}


total_mem=$(( mem_per_cpu * cpu_per_task ))

qsub -cwd \
     -l h_rt=$(convert_minutes_to_hms ${time}) \
     -l mem_free=${total_mem} \
     -pe thread ${cpu_per_task} \
     ${other_args}
