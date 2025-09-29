#!/usr/bin/env bash

# This script installs all the necessary packages for the Snakemake course

# Don't forget to clone the repo before running this script
# git clone https://github.com/sib-swiss/containers-snakemake-training.git
# cd containers-snakemake-training

# Update repo list
echo "Installing base packages" > server_setup.log
sudo apt update
sudo apt install -y software-properties-common

# Add Apptainer PPA and install it
echo "Installing Apptainer" >> server_setup.log
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt update
sudo apt install -y apptainer

# Install SLURM
echo "Installing SLURM" >> server_setup.log
sudo apt install -y slurmd slurmctld

# Configure SLURM
echo "Configuring SLURM" >> server_setup.log
cpu=`nproc`  # Get number of CPUs
sed -i "s/NB_CPU/$cpu/" setup_scripts/slurm.conf  # Replace actual number of CPUs
available_mem=`free -m | awk 'NR==2{print $2-4096}'`  # Get total RAM in MB; we remove 4GB for safety
sed -i "s/AVAILABLE_RAM/$available_mem/" setup_scripts/slurm.conf  # Replace actual available RAM
sockets=`lscpu | grep "Socket(s):" | awk '{print $2}'`  # Get number of sockets
sed -i "s/NB_SOCKETS/$sockets/" setup_scripts/slurm.conf  # Replace actual number of sockets
corespersocket=`lscpu | grep "Core(s) per socket:" | awk '{print $NF}'`  # Get number of cores per socket
sed -i "s/CORES_SOCKET/$corespersocket/" setup_scripts/slurm.conf  # Replace actual number of cores per socket
threads_core=`lscpu | grep "Thread(s) per core:" | awk '{print $NF}'`  # Get number of threads per core
sed -i "s/THREADS_CORE/$threads_core/" setup_scripts/slurm.conf  # Replace actual number of threads per core
sudo cp setup_scripts/slurm.conf /etc/slurm/slurm.conf

# Add SLURM alias for all users
echo "Adding SLURM alias" >> server_setup.log
echo "" | sudo tee -a /etc/bash.bashrc > /dev/null
echo "alias Squeue='squeue --user=\$USER --format=\"%11u %12i %10j %4Q %8q %3t %20R %20V %20S %11M %11l %6D %6C %7m\"'" | sudo tee -a /etc/bash.bashrc > /dev/null

# Start SLURM service and make node ready to receive jobs
echo "Starting SLURM" >> server_setup.log
sudo systemctl start slurmd
sudo scontrol update nodename=localhost state=idle

# Install miniforge in /opt/miniforge3
echo "Installing Miniforge" >> server_setup.log
sudo curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
sudo bash Miniforge3-$(uname)-$(uname -m).sh -b -p /opt/miniforge3
sudo rm Miniforge3-$(uname)-$(uname -m).sh

# Make conda available for all users
echo "Initialising conda for all users" >> server_setup.log
sudo /opt/miniforge3/bin/conda init --system

# Install Snakemake course environment
echo "Creating Snakemake environment"
sudo /opt/miniforge3/bin/conda env create -y -f conda/snakemake_course.yaml
echo "Setup complete!" >> server_setup.log

