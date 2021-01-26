#!/bin/bash

# Need to source .profile before running so you are using bash 4.4

# ls -p appends a slash to any directories
# grep -v / tells grep to return only lines not containing a slash
FILES=$(ls -p | grep -v /) # Excludes @eaDir directory
coverage=$1
num_jobs=$2
reference=$3
bed=$4
job_ids=($(seq -w $num_jobs))

# Create array to hold all commands to execute
neat_cmds=()

# Go through all files, generate all python commands and store in array
for file in ${FILES[@]}; do
  sample_name=$(basename $file .vcf) # i.e. CF_1, Hemo_1, P_2
  dirname="${sample_name}_${coverage}x"
  mkdir $dirname
  curr_dir=$(pwd)
  output="${curr_dir}/$dirname/$dirname"
  for id in ${job_ids[@]}; do
    build_job_cmd=("python /home/NBS/Read_Simulation/neat-genreads/genReads.py \
                -r $reference \
                -R 150 \
                -o $output \
                --bam \
                --vcf \
                --pe 300 30 \
                -t $bed \
                -to 0 \
                -c $coverage \
                --gz \
                -v $file \
                -M 0 \
                --job $id ${#job_ids[@]}")
    neat_cmds+=("$build_job_cmd")
  done
done

max_jobs=40; curr_jobs=0
for c in "${neat_cmds[@]}"; do
  # If true (number of jobs is greater than max_jobs),
  # wait until the next background job finishes to continue
  (($curr_jobs >= $max_jobs)) && wait -n
  # Print out command you are going to execute
  echo $c
  # Execute command and increment curr_jobs
  $c & ((++curr_jobs))
done
# Wait until all jobs are complete
# In other words, don't exit out of this script
wait
