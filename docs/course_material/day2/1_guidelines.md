## Course goal

Throughout the course, you will implement and improve a workflow to trim bulk RNAseq reads, align them to a genome, perform some quality checks (QC), count mapped reads, and identify Differentially Expressed Genes (DEG). After the last series of exercises, you will have implemented a simple workflow with commonly used Snakemake features. You will be able to use this workflow as a reference to implement your own workflows in the future.

## Software

All the software needed in this workflow is either:

* Already installed in the `snake_course` conda environment
* Already installed in a Docker container
* Will be installed via a conda environment during today's exercises

All information of this course is based on the [official documentation](https://snakemake.readthedocs.io/en/v8.20.3/) for Snakemake version `8.20.3`.

## Website colour code explanation

We tried to use a colour code throughout the website to make the different pieces of information easily distinguishable. Here's a quick summary about the colour blocks you will encounter:

!!! info "This is a supplementary piece of information"

!!! tip "This is a tip to help you solve an exercise"

!!! success "This is the answer to an exercise"

!!! warning "This is a warning about a potential problem"

!!! bug "This is an explanation about a common bug/error"

## Exercises

Each series of exercises is divided into multiple questions. Each question provides a background explanation, a description of the task at hand and additional details when required.

!!! tip "Hints for challenging questions"
    For the most challenging questions, hints will be provided. However, you should first try to solve the problems without tips!

## Answers

Do not hesitate to modify and overwrite your code from previous answers as the difficulty is incremental. The questions are designed to incite you to build your answers upon the previous ones.

!!! info "Restarting from a clean Snakefile"
    * If you feel that you drifted too far apart from the solution, you can always restart from the files provided in the [solutions folder](https://github.com/sib-swiss/containers-snakemake-training/tree/main/docs/solutions_day2/) of the course repository
    * At the start of sessions 3 and 4, you will also find a short note showing the command to download the Snakefile that should have been written at the end of the previous session

If something is not clear at any point, please call us and we will do our best to answer your questions! You can also check the [official Snakemake documentation](https://snakemake.readthedocs.io/en/stable/index.html) for more information.

## Computing environment

!!! warning "Development and computation"
    You can develop and write your scripts in a distant folder (using an `ssh` connection via VS code) or locally (if you do so, you will need to copy them on the server (with `scp`) before running them), but **remember that all the computations should be performed on the server, so don't forget to log in with `ssh -i key_username.pem username@18.195.170.182`**.

!!! bug "`Error: Command not found`"
    If you try to run a command and get an error such as `Command 'snakemake' not found`, you are probably in the wrong conda environment:

    * To list available conda environments, use `conda env list`
    * To activate an environment, use `conda activate <env_name>`
    * To deactivate an environment, use `conda deactivate`
    * To list packages installed in an environment, activate it and use `conda list`. **The computing environment on the server is called `snake_course`.**
