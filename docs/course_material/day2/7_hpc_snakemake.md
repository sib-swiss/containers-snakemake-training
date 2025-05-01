## TODOs

- Benchmark, resources
- Create and add slides
- Sync solutions

## Learning outcomes

**After having completed this chapter you will be able to:**

- Use the `benchmark` directive to understand the resource usage of each rule.
- Optimise resource usage in your workflow using the `resources` directive.
- Run your Snakemake workflow in an HPC environment using SLURM.
- Set up workflow-specific configuration profiles.
- Define rule-specific resource settings in the Snakefile.

## Material

[:fontawesome-solid-file-pdf: Download the presentation](../../assets/pdf/day2/7_snakemake_in_hpc.pdf){: .md-button }

## Snakefile from previous session

You can get the full RNA-seq workflow we have been working by running the following commands:

```sh
git clone https://github.com/sib-swiss/containers-snakemake-training.git  # Clone repository
cd containers-snakemake-training/  # Go in repository
cp -r docs/solutions_day2/session4 destination_path  # Copy folder where needed; adapt destination_path to required path
```

## Exercises

In the following exercises you will learn how to get a better idea of the resources used by your workflow, optimising the workflow resources usage, and the necessary command-line arguments, files and rule directives to be able to run your workflow on an HPC environment using SLURM. You will learn how to set up configuration profiles and further develop the Snakefile to define rule-specific settings. 

### Workflow benchmarking

Knowing how many resources are being used by each job can be very useful to optimise your workflow both for local and remote execution. The `benchmark` directive allows you to do exactly that.

**Exercise:** add the and `benchmark` directive to the rules. You can follow the same logic as you did with the `log` directive.

??? tip "Benchmarks"
    `benchmark` directives must contain the same `wildcards` as the `output` directive, here `sample`

??? success "Answer"
    ```python linenums="1" hl_lines="13-14"
    rule reads_quantification_genes:
        """
        This rule quantifies the reads of a bam file mapping on genes and produces
        a count table for all genes of the assembly.
        """
        input:
            bam_once_sorted = rules.sam_to_bam.output.bam_sorted,
        output:
            gene_level = 'results/{sample}/{sample}_genes_read_quantification.tsv',
            gene_summary = 'results/{sample}/{sample}_genes_read_quantification.summary'
        log:  
            'logs/{sample}/{sample}_genes_read_quantification.log'  # Path of log file
        benchmark:  # benchmark directive
            'benchmarks/{sample}/{sample}_genes_read_quantification.txt'  # Path of benchmark file
        container:
            'https://depot.galaxyproject.org/singularity/subread%3A2.0.6--he4a0461_2'
        shell:
            '''
            featureCounts -t exon -g gene_id -s 2 -p --countReadPairs \
            -B -C --largestOverlap --verbose -F GTF \
            -a resources/Scerevisiae.gtf -o {output.gene_level} {input.bam_once_sorted} &>> {log}
            mv {output.gene_level}.summary {output.gene_summary}
            '''
    ```

### Controlling memory usage and runtime

Another way to optimise resource usage in a workflow is to control resources such as the amount of memory or runtime of each job. This ensures your instance (computer, cluster...) won't run out of memory during computation (which could interrupt jobs or even crash the instance) and that your jobs will run in a reasonable amount of time (a job taking more time to run than usual might be a sign that something is going on).

??? info "Resource usage and schedulers"
    Optimising resource usage is especially important when submitting jobs to a scheduler (for instance on a cluster), as it allows a better and more precise definition of your job priority: jobs with low threads/memory/runtime requirements often have a higher priority than heavy jobs, which means they will often start first.

Controlling memory usage and runtime in Snakemake is very simple: you only need to need to use the `resources` directive with the `memory` or `runtime` keywords and in most software, you don't even need to specify the amount of memory available to the software via a parameter. Determining how much memory and runtime to use is also easier... **after the first run** of your workflow. The benchmark files created in the previous section are especially useful for that:

|  **s**  | **h: m: s** | **max_rss** | **max_vms** | **max_uss** | **max_pss** | **io_in** | **io_out** | **mean_load** | **cpu_time** |
|:-------:|:---------:|:-----------:|:-----------:|:-----------:|:-----------:|:---------:|:----------:|:-------------:|:------------:|
| 31.3048 |  0:00:31  |    763.04   |    904.29   |    757.89   |    758.37   |    1.81   |   230.18   |     37.09     |     11.78    |

* `s` and `h:m:s` give the job wall clock time (in seconds and hours-minutes-seconds, respectively), which is the actual time taken from the start of a software to its end. You can use these results to infer a safe value for the `runtime` keyword
* Likewise, you can use `max_rss` (shown in megabytes) to figure out how much memory was used by the job and use this value in the `memory` keyword

??? info "What are the other columns?"
    In case you are wondering about the other columns of the table, the [official documentation](https://snakemake.readthedocs.io/en/v8.20.5/snakefiles/rules.html#benchmark-rules) has detailed explanations about their content.

Here are some suggested values for the current workflow:

* `fastq_trim`: 500 MB
* `read_mapping`: 2 GB
* `sam_to_bam`: 250 MB
* `reads_quantification_genes`: 500 MB

**Exercise:** Implement memory usage limit in the workflow rules.

??? tip "Two ways to declare memory values "
    There are two ways to declare memory values in `resources`:

    1. `mem_<unit> = n`
    1. `mem = 'n<unit>'`: in this case, you must pass a **string**, so you have to enclose `n<unit>` with quotes `''`

    `<unit>` is a unit in [B, KB, MB, GB, TB, PB, KiB, MiB, GiB, TiB, PiB] and `n` is a **float**.

??? success "Answer"
    We implemented memory usage control in all the rules so that you can check everything. Rules `fastq_trim` and `read_mapping` have the first format while rules `sam_to_bam` and `reads_quantification_genes` have the second one. Feel free to copy this in your Snakefile:
    ```python linenums="1" hl_lines="18 19 47 48 76 77 107 108"
    rule fastq_trim:
        '''
        This rule trims paired-end reads to improve their quality. Specifically, it removes:
        - Low quality bases
        - A stretches longer than 20 bases
        - N stretches
        '''
        input:
            reads1 = 'data/{sample}_1.fastq',
            reads2 = 'data/{sample}_2.fastq',
        output:
            trim1 = 'results/{sample}/{sample}_atropos_trimmed_1.fastq',
            trim2 = 'results/{sample}/{sample}_atropos_trimmed_2.fastq'
        log:
            'logs/{sample}/{sample}_atropos_trimming.log'
        benchmark:
            'benchmarks/{sample}/{sample}_atropos_trimming.txt'
        resources:  # Add directive
            mem_mb = 500  # Add keyword and value with format 1
        threads: 2
        container:
            'https://depot.galaxyproject.org/singularity/atropos%3A1.1.32--py312hf67a6ed_2'
        shell:
            '''
            echo "Trimming reads in <{input.reads1}> and <{input.reads2}>" > {log}
            atropos trim -q 20,20 --minimum-length 25 --trim-n --preserve-order --max-n 10 \
            --no-cache-adapters -a "A{{20}}" -A "A{{20}}" --threads {threads} \
            -pe1 {input.reads1} -pe2 {input.reads2} -o {output.trim1} -p {output.trim2} &>> {log}
            echo "Trimmed files saved in <{output.trim1}> and <{output.trim2}> respectively" >> {log}
            echo "Trimming report saved in <{log}>" >> {log}
            '''

    rule read_mapping:
        '''
        This rule maps trimmed reads of a fastq onto a reference assembly.
        '''
        input:
            trim1 = rules.fastq_trim.output.trim1,
            trim2 = rules.fastq_trim.output.trim2
        output:
            sam = 'results/{sample}/{sample}_mapped_reads.sam',
            report = 'results/{sample}/{sample}_mapping_report.txt'
        log:
            'logs/{sample}/{sample}_mapping.log'
        benchmark:
            'benchmarks/{sample}/{sample}_mapping.txt'
        resources:  # Add directive
            mem_gb = 2  # Add keyword and value with format 1
        threads: 4
        container:
            'https://depot.galaxyproject.org/singularity/hisat2%3A2.2.1--hdbdd923_6'
        shell:
            '''
            echo "Mapping the reads" > {log}
            hisat2 --dta --fr --no-mixed --no-discordant --time --new-summary --no-unal \
            -x resources/genome_indices/Scerevisiae_index --threads {threads} \
            -1 {input.trim1} -2 {input.trim2} -S {output.sam} --summary-file {output.report} 2>> {log}
            echo "Mapped reads saved in <{output.sam}>" >> {log}
            echo "Mapping report saved in <{output.report}>" >> {log}
            '''

    rule sam_to_bam:
        '''
        This rule converts a sam file to bam format, sorts it and indexes it.
        '''
        input:
            sam = rules.read_mapping.output.sam
        output:
            bam = 'results/{sample}/{sample}_mapped_reads.bam',
            bam_sorted = 'results/{sample}/{sample}_mapped_reads_sorted.bam',
            index = 'results/{sample}/{sample}_mapped_reads_sorted.bam.bai'
        log:
            'logs/{sample}/{sample}_mapping_sam_to_bam.log'
        benchmark:
            'benchmarks/{sample}/{sample}_mapping_sam_to_bam.txt'
        resources:  # Add directive
            mem = '250MB'  # Add keyword and value with format 2
        threads: 2
        container:
            'https://depot.galaxyproject.org/singularity/samtools%3A1.21--h50ea8bc_0'
        shell:
            '''
            echo "Converting <{input.sam}> to BAM format" > {log}
            samtools view {input.sam} --threads {threads} -b -o {output.bam} 2>> {log}
            echo "Sorting .bam file" >> {log}
            samtools sort {output.bam} --threads {threads} -O bam -o {output.bam_sorted} 2>> {log}
            echo "Indexing sorted .bam file" >> {log}
            samtools index -b {output.bam_sorted} -o {output.index} 2>> {log}
            echo "Sorted file saved in <{output.bam_sorted}>" >> {log}
            '''

    rule reads_quantification_genes:
        '''
        This rule quantifies the reads of a bam file mapping on genes and produces
        a count table for all genes of the assembly. The strandedness parameter
        is determined by get_strandedness().
        '''
        input:
            bam_once_sorted = rules.sam_to_bam.output.bam_sorted,
        output:
            gene_level = 'results/{sample}/{sample}_genes_read_quantification.tsv',
            gene_summary = 'results/{sample}/{sample}_genes_read_quantification.summary'
        log:
            'logs/{sample}/{sample}_genes_read_quantification.log'
        benchmark:
            'benchmarks/{sample}/{sample}_genes_read_quantification.txt'
        resources:  # Add directive
            mem = '500MB'  # Add keyword and value with format 2
        threads: 2
        container:
            'https://depot.galaxyproject.org/singularity/subread%3A2.0.6--he4a0461_2'
        shell:
            '''
            echo "Counting reads mapping on genes in <{input.bam_once_sorted}>" > {log}
            featureCounts -t exon -g gene_id -s 2 -p --countReadPairs \
            -B -C --largestOverlap --verbose -F GTF \
            -a resources/Scerevisiae.gtf -T {threads} -o {output.gene_level} {input.bam_once_sorted} &>> {log}
            echo "Renaming output files" >> {log}
            mv {output.gene_level}.summary {output.gene_summary}
            echo "Results saved in <{output.gene_level}>" >> {log}
            echo "Report saved in <{output.gene_summary}>" >> {log}
            '''
    ```

### Using configuration profiles 

Configuration profiles allow you to specify default and rule-specific SLURM settings. Currently, Snakemake supports two types of profiles: global and workflow-specific. In this tutorial we will focus on workflow-specific profiles. 

You can create a workflow profile named `slurm_profile` by:

1. Creating a directory named `slurm_profile` the project directory.
2. Creating a `config.yaml` file inside `slurm_profile`. 

Below you can find an example file structure:

```hl_lines="3 4"
│── config
│   └── config.yaml
│── slurm_profile
│   └── config.yaml
└── workflow
    │── Snakefile
    │── envs
    │   │── tool1.yaml
    │   └── tool2.yaml
    │── rules
    │   │── module1.smk
    │   └── module2.smk
    └── scripts
        │── script1.py
        └── script2.R
```

The `config.yaml` file inside the `slurm_profile` directory can contain multiple settings related to the execution of the workflow. Below you can find an example of a configuration file:

```yaml
executor: cluster-generic
cluster-generic-submit-cmd: 'sbatch --job-name={rule} --cpus-per-task={threads}'
jobs: 2
```

**What does each of the settings in the `slurm_profile/config.yaml` file do?**

??? success "Answer"
    These settings tell Snakemake how to communicate with the job scheduler. If you check out the Snakemake help with `snakemake --help` you will see that all the settings in the configuration profile are arguments you can pass to the `snakemake` command:
    
    * `executor` indicates which executor plugin to use. The `cluster-generic` pluggin allows us to send jobs by specifying the submission command.
    * `cluster-generic-submit-cmd` indicates the command to use to submit the jobs. Since we are using SLURM, the command is `sbatch`. This string can contain any valid `sbatch` arguments and use values from the `Snakefile` (i.e. using `threads` specified in the Snakefile as a value for the `--cpus-per-task` argument). 
    * `jobs` indicates how many jobs should be sent at a time.

**Exercise:** add the required SLURM argument to log the run of each job to a file in `slurm_logs/{rule}/{rule}_{wildcards}.log`. **Hint:** you can find the SLURM argument to use by running `sbatch -h` or by checking this [SLURM cheatsheet](https://slurm.schedmd.com/pdfs/summary.pdf).

??? warning "Important!"
    In older SLURM versions, the directory `logs` needs to exist before running the workflow! In order to be able to save the logs into a directory, you will need to have created it before running Snakemake. Two ways to go about this are:

    * Creating it before running the workflow.
    * Adding the directory creation command in `cluster-generic-submit-cmd`. This option ensures that no problems arise if you are running the workflow for the very first time.

??? success "Answer"

    ```yaml
    executor: cluster-generic
    cluster-generic-submit-cmd: 
        mkdir -p slurm_logs/{rule} &&
        sbatch 
            --job-name={rule}
            --cpus-per-task={threads}
            --output=slurm_logs/{rule}/{rule}_{wildcards}.log
    jobs: 2
    ```

### Passing rule-specific resources to the job scheduler

As we have seen, it is possible to specify rule-specific resources to have better control over how the workflow uses the available resources. This is especially important when running your workflow in an HPC to avoid problems such as RAM-intensive jobs crashing, or unnecessarily allocating a lot of resources for small jobs. As we saw before, you can define rule-specific resources using the `resources` directive as shown below:

```python linenums="1" hl_lines="6 7"
rule myrule:
    input:
        ...
    output:
        ...
    resources:
        mem_mb=100
    shell:
        """
        cp {input} {output};
        """
```

This will ensure that jobs from rule `myrule` never use more than 100 MB of RAM. Note that the jobs will crash if they surpass this limit, so try to account for some wiggle room when you set resources. 

**Exercise:** update the configuration profile to define the amount of memory per job using the values specified in each rule through the `resources` directive. Note that some rules use `mem_mb` while some others use `mem_gb`, so change them all first to `mem_mb`. 

??? success "Answer"

    ```yaml
    executor: cluster-generic
    cluster-generic-submit-cmd: 
        mkdir -p logs/{rule} &&
        sbatch
            --cpus-per-task={threads}
            --mem={resources.mem_mb}
            --job-name={rule}-{wildcards}
            --output=slurm_logs/{rule}/{rule}-{wildcards}.log
            --export=ALL
            --chdir=$PWD
    cluster-generic-cancel-cmd: scancel
    jobs: 2
    software-deployment-method: 
        - conda
        - apptainer
    printshellcmds: True
    ```

### Adapting the workflow to the available resources

Before launching our worklfow, it is very important to understand what are the available resources in the machine where we are going to run it. This information is typically provided by the HPC user guide. In this case, you can use the following commands to know what is the amount of available cores and RAM:

```sh
nproc --all     # Number of cores
free -g         # Amount of memory in GB under the column named "total"
```

**Exercise:** based on the available resources, figure out if you need to update any of the rule resources.

??? warning "Important!"
    A very important role of schedulers is to ensure the currently running jobs do not surpass the total amount of available resources. When not enough resources are available, Snakemake will display the message `Waiting for more resources.` to indicate that a job is waiting until enough resources are available for it to run. However, if a job is asking for more resources than a machine has (i.e. the job asks for 50 CPUs, but the machine only has 20), Snakemake will get stuck and display the same message but will never manage to launch the job.

### Final exercise

To conclude, we will run our workflow in the HPC environment with the following command:

```sh
snakemake  --profile slurm_profile
```

This is a great time to get yourself a coffee/tea and celebrate! :coffee: :tea:

Congratulations, you made it to the end! You are now able to create a Snakemake workflow, make it reproducible thanks to Conda and Docker/Apptainer and even run it in an HPC! To make things even better, have a look at [Snakemake's best practices](https://snakemake.readthedocs.io/en/v8.20.5/snakefiles/best_practices.html)!