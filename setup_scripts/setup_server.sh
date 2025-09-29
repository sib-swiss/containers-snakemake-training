#!/usr/bin/env bash

# This script installs all the necessary packages for the Snakemake course

# Don't forget to clone the repo before running this script
# git clone https://github.com/sib-swiss/containers-snakemake-training.git
# cd containers-snakemake-training

# Update repo list
sudo apt update
sudo apt install -y software-properties-common

# Add Apptainer PPA and install it
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt update
sudo apt install -y apptainer

# Install SLURM
sudo apt install -y slurmd slurmctld

# Configure SLURM
cpu=`nproc`  # Get number of CPUs
sed -i "s/NB_CPU/$cpu/" setup_scripts/slurm.conf  # Replace actual number of CPUs
available_mem=`free -m | awk 'NR==2{print $2-2048}'`  # Get total RAM in MB; we remove 2GB for safety
sed -i "s/AVAILABLE_RAM/$available_mem/" setup_scripts/slurm.conf  # Replace actual available RAM
sudo cp setup_scripts/slurm.conf /etc/slurm/slurm.conf

# Add SLURM alias for all users
sudo echo -e "\nalias Squeue = 'squeue -u $USER'" >> /etc/bash.bashrc

# Start SLURM service and make node ready to receive jobs
sudo systemctl start slurmd
sudo scontrol update nodename=localhost state=idle

# Install miniforge in /opt/miniforge3
sudo curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
sudo bash Miniforge3-$(uname)-$(uname -m).sh -b -p /opt/miniforge3
sudo rm Miniforge3-$(uname)-$(uname -m).sh

# Install Snakemake course environment
sudo mamba env create -f conda/snakemake_course.yaml

