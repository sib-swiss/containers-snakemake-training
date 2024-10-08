## Workshop goal

Over the course of the workshop, you will implement and improve a workflow to trim bulk RNAseq reads, align them to a genome, perform some quality checks (QC), count mapped reads, and identify Differentially Expressed Genes (DEG). The goal of the workshop is that after the last series of exercises, you will have implemented a simple workflow with commonly used Snakemake features. You will be able to use this workflow as a reference to implement your own workflows in the future.

## Software

All the software needed in this workflow is either:

* Already installed in the `snake_course` conda environment
* Already installed in a Docker container
* Will be installed via a conda environment during today's exercises

All information of this course is based on the [official documentation](https://snakemake.readthedocs.io/en/v7.32.3/) for Snakemake version `7.32.3`. 

## Exercises

Each series of exercises is divided into multiple questions. We first provide a general explanation on the context behind each question; we then explicitly describe the task and provide details when they are required. We also provide hints that should help you with the most challenging parts of some questions. You should first try to solve the problems without using these hints! Do not hesitate to modify and overwrite your code from previous questions when specified in an exercise, as the solutions for each series of exercises are provided. If something is not clear at any point, please call us and we will do our best to answer your questions. You can also check the [official Snakemake documentation](https://snakemake.readthedocs.io/en/stable/index.html) for more information.
