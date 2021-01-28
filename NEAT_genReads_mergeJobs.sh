#!/bin/bash

# Need to source .profile before running so you are using bash 4.4
source ~/.profile

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
mergeJobs_cmds=()

# Go through all files, generate all python commands and store in array
for file in ${FILES[@]}; do
  sample_name=$(basename $file .vcf) # i.e. CF_1, Hemo_1, P_2
  dirname="${sample_name}_${coverage}x"
  mkdir $dirname
  curr_dir=$(pwd)
  output="${curr_dir}/$dirname/$dirname"
  mergeJobs_output="${curr_dir}/${dirname}_merged"
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
  build_merge_cmd=("python /home/NBS/Read_Simulation/neat-genreads/mergeJobs.py \
                      -i $output \
                      -o $mergeJobs_output \
                      -s /usr/local/bin/samtools")
  mergeJobs_cmds+=("$build_merge_cmd")
done

max_genReads_jobs=40; curr_genReads_jobs=0
# Generate all read sets
for c in "${neat_cmds[@]}"; do
  # If true (number of jobs is greater than max_jobs),
  # wait until the next background job finishes to continue
  (($curr_genReads_jobs >= $max_genReads_jobs)) && wait -n
  # Print out command you are going to execute
  echo $c
  # Execute command and increment curr_jobs
  $c & ((++curr_genReads_jobs))
done

# For a single mergeJobs command, it creates 4 commands:
# cat all R1 fastqs
# cat all R2 fastqs
# cat all bams
# cat all vcfs
max_mergeReads_jobs=10; curr_mergeReads_jobs=0
# Concatenate reads, vcf, bams for each sample
for c in "${mergeJobs_cmds[@]}"; do
  (($curr_mergeReads_jobs >= $max_mergeReads_jobs)) && wait -n
  echo $c
  $c & ((++curr_mergeReads_jobs))
done
# Wait until all jobs are complete
# In other words, don't exit out of this script
wait
