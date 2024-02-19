#!/bin/bash
#SBATCH --job-name=api-job       # create a short name for your job
#SBATCH --partition=gpu          # use gpu partition
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=1G         # memory per cpu-core
#SBATCH --time=00:05:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all          # send email for all job steps
#SBATCH --mail-user=username@domain.com

# File that we submit for the job scheduler and we give the argument

# /usr/bin/python3 /home/username/job.py "$1"

