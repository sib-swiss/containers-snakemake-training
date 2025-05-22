#!/usr/bin/env bash

# This script installs the necessary packages for the Snakemake course

# Update repo list
sudo apt update
sudo apt install -y software-properties-common

# Add Apptainer PPA and install it
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt update
sudo apt install -y apptainer

# Install miniforge in /opt/miniforge3
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh

# Install snake_course environment
git clone https://github.com/sib-swiss/containers-snakemake-training.git
mamba env create -f conda/snakemake_course.yaml

