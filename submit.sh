#!/bin/bash
# SGE submit wrapper used by config/profiles/Migale (cluster-generic executor).
# Called as: submit.sh <runtime_min> <mem_mb> <cpus> [qsub args...]
#
# - <mem_mb> is the TOTAL memory of the job in MB (snakemake resources.mem_mb),
#   so it is passed as is with the "M" unit (SGE reads a value without unit as bytes).
# - "-q a.q,b.q": SGE sorts queues by load, so the order of the list gives no
#   priority. The list is kept as hard request and its first queue is added as
#   soft request, so that it is preferred and the others are only used when it is full.

time=$1
mem_mb=$2
cpu_per_task=$3
shift 3

convert_minutes_to_hms() {
    mins=$1
    hrs=$(( mins / 60 ))
    mins=$(( mins % 60 ))
    printf "%02d:%02d:00\n" "$hrs" "$mins"
}

queue_args=()
other_args=()
while [[ $# -gt 0 ]]; do
    if [[ "$1" == "-q" && $# -ge 2 ]]; then
        queue_args=(-hard -q "$2")
        if [[ "$2" == *,* ]]; then
            queue_args+=(-soft -q "${2%%,*}" -hard)
        fi
        shift 2
    else
        other_args+=("$1")
        shift
    fi
done

qsub -cwd \
     -r n \
     -l h_rt=$(convert_minutes_to_hms ${time}) \
     -l mem_free=${mem_mb}M \
     -pe thread ${cpu_per_task} \
     "${queue_args[@]}" \
     "${other_args[@]}"
