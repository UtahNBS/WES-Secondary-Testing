#!/bin/bash

# Need to source .profile before running so you are using bash 4.4
source ~/.profile

# List all directories for the current directory, excluding @eaDir
DIRS=$(ls -d */ | grep -v '@eaDir')

# Create array to hold all commands to execute
mergeJobs_cmds=()
curr_dir=$(pwd)

# Go through all files, generate all python commands and store in array
for d in ${DIRS[@]}; do
  input="${curr_dir}/${d}/${d}"
  output="${curr_dir}/${d}/${d}_merged"
  build_merge_cmd=("python /home/NBS/Read_Simulation/neat-genreads/mergeJobs.py \
                      -i $input \
                      -o $output \
                      -s /usr/local/bin/samtools")
  mergeJobs_cmds+=("$build_merge_cmd")
done

# For a single mergeJobs command, it creates 4 commands:
# cat all R1 fastqs
# cat all R2 fastqs
# cat all bams
# cat all vcfs
max_jobs=10; curr_jobs=0
for c in "${mergeJobs_cmds[@]}"; do
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
