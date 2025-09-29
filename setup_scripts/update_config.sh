#!/usr/bin/env bash

# This script update the SLURM configuration based on the current hardware

# Reconfigure SLURM
cpu=`nproc`  # Get number of CPUs
available_mem=`free -m | awk 'NR==2{print $2-4096}'`  # Get total RAM in MB; we remove 4GB for safety
sockets=`lscpu | grep "Socket(s):" | awk '{print $2}'`  # Get number of sockets
corespersocket=`lscpu | grep "Core(s) per socket:" | awk '{print $NF}'`  # Get number of cores per socket
threads_core=`lscpu | grep "Thread(s) per core:" | awk '{print $NF}'`  # Get number of threads per core
config="NodeName=localhost CPUs=${cpu} RealMemory=${available_mem} Sockets=${sockets} CoresPerSocket=${corespersocket} ThreadsPerCore=${threads_core} State=UNKNOWN"  # Build replacing line
sudo sed -i -e "s/NodeName.*/$config/" /etc/slurm/slurm.conf  # Replace config; only works on a server with a single node

# Restart SLURM to load new config
sudo systemctl restart slurmd
sudo systemctl restart slurmctld
