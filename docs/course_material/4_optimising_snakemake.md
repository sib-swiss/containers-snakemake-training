## Learning outcomes

**After having completed this chapter you will be able to:**

* Use non-file parameters and config files in rules
* Modularise a workflow
* Make a workflow process list of inputs instead of one input at a time
* Aggregate outputs in a target rule
* (Optimise CPU usage)

## Material

[:fontawesome-solid-file-pdf: Download the presentation](../assets/pdf/4_optimising_snakemake.pdf){: .md-button }

## Snakefile from previous session

If you didn't finish the previous part or didn't do the optional exercises, you can restart from a fully commented Snakefile, with log messages implemented in all rules. You can download it [here](https://raw.githubusercontent.com/sib-swiss/containers-snakemake-training/main/docs/solutions/session2/workflow/Snakefile) or download it in your current directory with:
```sh
wget https://raw.githubusercontent.com/sib-swiss/containers-snakemake-training/main/docs/solutions/session2/workflow/Snakefile
```

## Exercises

This series of exercises focuses on how to improve the workflow that you developed in the previous session. As a result, you will add only one rule to your workflow. But, fear not, it's a crucial one!

??? tip "Development and back-up"
    During this session, you will modify your Snakefile quite heavily, so it may be a good idea to make back-ups from time to time (with `cp` or a simple copy/paste) or use a versioning system. As a general rule, if you have a doubt on the code you are developing, do not hesitate to make a back-up beforehand.

### Using non-file parameters and config files

#### Non-file parameters

As you have seen, Snakemake execution revolves around input and output files. However, a lot of software also use **non-file parameters** to run. In the previous presentation and series of exercises, we advocated against using hard-coded file paths. Yet, if you look back at previous rules, you will find two occurrences of this behaviour in `shell` directives:

* In rule `read_mapping`, the index parameter:
```sh
-x resources/genome_indices/Scerevisiae_index
```
* In rule `reads_quantification_genes`, the annotation parameter:
```sh
-a resources/Scerevisiae.gtf
```

This reduces readability and makes it very hard to change the values of these parameters, because this requires to change the `shell` directive code.

The `params` directive was (partly) designed to solve this problem: it contains parameters and variables that can be accessed in the `shell` directive. It allows to specify additional non-file parameters instead of hard-coding them into shell commands or using them as inputs/outputs.

??? info "Main properties of parameters from the `params` directive"
    * Their values can be of any type (integer, string, list...)
    * Their values can depend on wildcard values and use input functions (explained [here](5_reproducibility_snakemake.md#gathering-the-input-files)). This means that parameters can be changed conditionally, for example depending on the value of a wildcard
        * In contrast to the `input` directive, the `params` directive can take more arguments than only `wildcards`, namely `input`, `output`, `threads`, and `resources`
    * Similarly to `{input}` and `{output}` placeholders, they can be accessed from within the `shell` directive with the `{params}` placeholder
    * Multiple parameters can be defined in a rule (do not forget the comma between each entry!) and they can also be named. While it isn't mandatory, un-named parameters are not explicit at all, so you should **always name your parameters**

Here is an example of `params` utilisation:
```python linenums="1"
rule get_header:
    input:
        'data/example.txt'
    output:
        'results/header.txt'
    params:
        lines = 5
    shell:
        'head -n {params.lines} {input} > {output}'
```

**Exercise:** Pick one of the two hard-coded paths mentioned earlier and replace it using `params`.

??? success "Answer"
    You need to add a `params` directive to the rule, name the parameter and replace the path by the placeholder in the `shell` directive. We did this for both rules so that you can check everything. Feel free to copy this in your Snakefile. For clarity, only lines that changed are shown below:

    * `read_mapping`:
    ```python linenums="1" hl_lines="5"
    params:
        index = 'resources/genome_indices/Scerevisiae_index'
    shell:
        'hisat2 --dta --fr --no-mixed --no-discordant --time --new-summary --no-unal \
        -x {params.index} \  # Parameter was replaced here
        -1 {input.trim1} -2 {input.trim2} -S {output.sam} --summary-file {output.report} 2>> {log}'
    ```

    * `reads_quantification_genes`:
    ```python linenums="1" hl_lines="6"
    params:
        annotations = 'resources/Scerevisiae.gtf'
    shell:
        'featureCounts -t exon -g gene_id -s 2 -p --countReadPairs \
        -B -C --largestOverlap --verbose -F GTF \
        -a {params.annotations} -o {output.gene_level} {input.bam_once_sorted} &>> {log}'  # Parameter was replaced here
    ```

But doing this only shifted the problem: now, hard-coded paths are in `params` instead of `shell`. This is better, but not by much! Luckily, there is an even better way to handle parameters: instead of hard-coding parameter values in the Snakefile, Snakemake can use parameters (and values) defined in config files.

#### Config files

Config files are stored in the `config` subfolder and written in [JSON](https://en.wikipedia.org/wiki/JSON) or [YAML](https://en.wikipedia.org/wiki/YAML) format. You will use the latter for this course as it is more user-friendly. In .yaml files:

* Parameters are defined with the syntax `<name>: <value>`
* Values can be strings, integers, booleans...
    * For a complete overview of available value types, see [this list](https://learnxinyminutes.com/docs/yaml/)
* A parameter can have multiple values, each value being on an indented line starting with "**-**"
    * These values will be stored in a Python list when Snakemake parses the config file
* Parameters can be named and nested to have a hierarchical structure, each sub-parameter and its value being on an indented lines
    * These parameters will be stored as a dictionary when Snakemake parses the config file

Config files will be parsed by Snakemake when executing the workflow, and parameters and their values will be stored in a [Python dictionary](https://docs.python.org/3/tutorial/datastructures.html#dictionaries) named `config`. The config file path can be specified in the Snakefile with `configfile: <path/to/file.yaml>` at the top of the file, or at runtime with the execution parameter `--configfile <path/to/file.yaml>`.

The example below shows a parameter with a single value (`lines_number`), a parameter with multiple values (`samples`), and nested parameters (`resources`):
```yaml linenums="1"
lines_number: 5  # Parameter with single value (string, int, float, bool ...)
samples:  # Parameter with multiple values
    - sample1
    - sample2
resources:  # Nested parameters
    threads: 4
    memory: 4G
```

Then, each parameter can be accessed in Snakemake with:
```python linenums="1"
config['lines_number']  # --> 5
config['samples']  # --> ['sample1', 'sample2']  # A list of parameters becomes a list
config['resources']  # --> {'threads': 4, 'memory': '4G'}  # A list of named parameters becomes a dictionary
config['resources']['threads']  # --> 4
```

??? warning "Accessing config values in `shell`"
    You cannot use values from the `config` dictionary directly in a `shell` directive. If you need to access parameter value in `shell`, first define it in `params` and assign its value from the dictionary, then use `params.<name>` in `shell`.

**Exercise:** Create a config file in YAML format and fill it with variables and values to replace one of the two hard-coded parameters mentioned before. Then replace the hard-coded parameter values by variables from the config file. Finally, add the config file import on top of your Snakefile.

??? success "Answer"
    First, create an empty config file:
    ```sh
    touch config/config.yaml  # Create empty config file
    ```

    Then, fill it with the desired values:
    ```yaml linenums="1"
    # Configuration options of RNAseq-analysis workflow
    # Location of genome indices
    index: 'resources/genome_indices/Scerevisiae_index'
    # Location of annotation file
    annotations: 'resources/Scerevisiae.gtf'
    ```

    Then, replace the `params` values in the Snakefile. We did this for both rules so that you can check everything. Feel free to copy this in your Snakefile. For simplicity, only lines that changed are shown below:

    * `read_mapping`:
    ```python linenums="1"
    params:
        index = config['index']
    ```

    * `reads_quantification_genes`:
    ```python linenums="1"
    params:
        annotations = config['annotations']
    ```

    Finally, add the file path on top of the Snakefile:
    `configfile: 'config/config.yaml'`

From now on, if you need to change these values, you can easily do it in the config file instead of modifying the code!

### Modularising a workflow

If you develop a large workflow, you are bound to encounter some cluttering problems. Have a look at your current Snakefile: with only four rules, it is already almost 150 lines long. Imagine what happens when your workflow has dozens of rules? The Snakefile may (will?) become messy and harder to maintain and edit. This is why it quickly becomes crucial to modularise your workflow. This approach also makes it easier to re-use pieces of one workflow into another. Snakemake can be modularised on three different levels:

1. The most fine-grained level is wrappers

    ??? info "More information on wrappers"
        Wrappers allow to quickly use popular tools and libraries in Snakemake workflows, thanks to the `wrapper` directive. Wrappers are automatically downloaded and deploy a conda environment when running the workflow, which increases reproducibility. However their implementation can be 'rigid' and sometimes it may be better to write your own rule. See the [official documentation](https://snakemake.readthedocs.io/en/v9.11.6/snakefiles/modularization.html#wrappers) for more explanations

1. For larger parts belonging to the same workflow, it is recommended to split the main Snakefile into smaller snakefiles, each containing rules with a common topic. Smaller snakefiles are then integrated into the main Snakefile with the `include` statement. In this case, all rules share a common config file. See the [official documentation](https://snakemake.readthedocs.io/en/v9.11.6/snakefiles/modularization.html#includes) for more explanations

    ??? info "Rules organisation"
        There is no official guideline on how to regroup rules, but a simple and logic approach is to create "thematic" snakefiles, _i.e._ place rules related to the same topic in the same file. Modularisation is a common practice in programming in general: it is often easier to group all variables, functions, classes... related to a common theme into a single script, package, software...

1. The final level of modularisation is modules

    ??? info "More on modules"
        It enables combination and re-use of rules in the same workflow and between workflows. This is done with the `module` statement, similarly to Python `import`. See the [official documentation](https://snakemake.readthedocs.io/en/v9.11.6/snakefiles/modularization.html#modules) for more explanations

In this course, you will only use the second level of modularisation. Briefly, the idea is to write a main Snakefile in `workflow/Snakefile`, to place the other snakefile containing rules in the subfolder `workflow/rules` (these 'sub-Snakefile' should end with `.smk`, the recommended file extension of Snakemake) and to tell Snakemake to import the modular snakefile in the main Snakefile with the `include: <path/to/snakefile.smk>` syntax.


**Exercise:** Move your current Snakefile into the subfolder `workflow/rules` and rename it to `read_mapping.smk`. Then create a new Snakefile in `workflow/` and import `read_mapping.smk` in it using the `include` syntax. You should also move the importation of the config file from the modular Snakefile to the main one.

??? success "Answer"
    First, move and rename the main Snakefile:
    ```sh linenums="1"
    mv workflow/Snakefile workflow/rules/read_mapping.smk  # Move and rename main Snakefile to modular snakefile
    touch workflow/Snakefile  # Recreate main Snakefile
    ```

    Then, add `include` and `configfile` statements to the new Snakefile. It should resemble this:
    ```python linenums="1"
    '''
    Main Snakefile of RNAseq analysis workflow. This workflow can clean and
    map reads, and perform Differential Expression Analyses.
    '''

    # Config file path
    configfile: 'config/config.yaml'

    # Rules to execute workflow
    include: 'rules/read_mapping.smk'
    ```

    Finally, remove the config file import (`configfile: 'config/config.yaml'`) from the modular snakefile (`workflow/rules/read_mapping.smk`).

??? info "Relative paths"
    * Include statements are relative to the directory of the Snakefile in which they occur. For example, if the Snakefile is in `workflow`, then Snakemake will search for included snakefiles in `workflow/path/to/other/snakefile`, regardless of the working directory
    * You can place snakefiles in a sub-directory without changing input and output paths, because these paths are relative to the working directory
    * **However, you will need to edit paths to external scripts and conda environments, because these paths are relative to the snakefile from which they are called** (this will be discussed in the [last series of exercises](5_reproducibility_snakemake.md#deciding-how-to-run-the-script))

If you have trouble visualising what an `include` statement does, you can imagine that the entire content of the included file gets copied into the Snakefile. As a consequence, syntaxes like `rules.<rule_name>.output.<output_name>` can still be used in modular snakefiles, even if the rule `<rule_name>` is defined in another snakefile. However, you have to make sure that **the snakefile in which `<rule_name>` is defined is included before the snakefile that uses `rules.<rule_name>.output`**. This is also true for input functions and checkpoints.

### Using a target rule instead of a target file

Modularisation also offers a great opportunity to facilitate workflows execution. By default, if no target is given in the command line, Snakemake executes the first rule in the Snakefile. So far, you have always executed the workflow with a target file to avoid this behaviour. But we can actually use this property to make execution easier by writing a pseudo-rule (also called target-rule and usually named rule `all`) which contains all the desired outputs files as inputs in the Snakefile. This rule will look like this:
```python linenums="1"
rule all:
    input:
        'path/to/ouput1',
        'path/to/ouput2',
        '...'
```

**Exercise:** Implement a rule `all` in your Snakefile to generate the final outputs by default when running `snakemake` without specifying a target. Then, test your workflow with a dry-run and the `-F` parameter. How many rules does Snakemake run?

??? tip "Content of a rule `all`"
    * A rule is not required to have an output nor a `shell` directive
    * Inputs of rule `all` should be the final outputs that you want to generate, here those of rule `reads_quantification_genes`

??? success "Answer"
    `reads_quantification_genes` is currently creating the last workflow outputs (with the same example as before, `results/highCO2_sample1/highCO2_sample1_genes_read_quantification.tsv`). We need to use these files as inputs of rule `all`:
    ```python linenums="1"
    # Master rule used to launch workflow
    rule all:
        '''
        Dummy rule to automatically generate required outputs.
        '''
        input:
            'results/highCO2_sample1/highCO2_sample1_genes_read_quantification.tsv'
    ```

    You can launch a dry-run with:
    ```sh
    snakemake -c 4 -F -p -n
    ```

    You should see all the rules appearing, including rule `all`, and the job stats:
    ```sh
    rule all:
        input: results/highCO2_sample1/highCO2_sample1_genes_read_quantification.tsv
        jobid: 0
        reason: Forced execution
        resources: tmpdir=/tmp

    Job stats:
    job                           count
    --------------------------  -------
    all                               1
    fastq_trim                        1
    read_mapping                      1
    reads_quantification_genes        1
    sam_to_bam                        1
    total                             5
    ```
    Snakemake runs 5 rules in total: the 4 of the previous session and the rule `all`.

After several (dry-)runs, you may have noticed that the rule order is not always the same: apart from Snakemake considering the first rule of the workflow as a default target, the order of rules in Snakefile/snakefiles is arbitrary and does not influence the DAG of jobs.

### Aggregating outputs to process lists of files

Using a target rule like the one presented in the previous paragraph gives another opportunity to make things easier. In the previous rule `all`, inputs are still hard-coded... and you know that this is not an optimal solution, especially if there are many samples to process. The [`expand()` function](https://snakemake.readthedocs.io/en/v9.11.6/snakefiles/rules.html#the-expand-function) will solve both problems.

`expand()` is used to generate a list of output files by automatically expanding a wildcard expression to several values. In other words, it will replace a wildcard in an expression by all the values of a list, successively. For instance, `expand('{sample}.tsv', sample=['A', 'B', 'C'])` will create the list of files `['A.tsv', 'B.tsv', 'C.tsv']`.

**Exercise:** Use an expand syntax to transform rule `all` to generate a list of final outputs **with all the samples**. Then, test your workflow with a dry-run and the `-F` parameter.

??? tip "Two things are required for an expand syntax"
    1. A Python list of values that will replace a wildcard; here, a sample list
    1. An output path with a wildcard that can be turned into an `expand()` function to create all the required outputs

??? success "Answer"
    First, we need to add a sample list in the Snakefile, **before rule `all`**. This list contains all the values that the wildcard will be replaced with:
    ```python
    SAMPLES = ['highCO2_sample1', 'highCO2_sample2', 'highCO2_sample3', 'lowCO2_sample1', 'lowCO2_sample2', 'lowCO2_sample3']
    ```

    Then, we need to transform the rule `all` inputs to use the `expand` function:
    ```python
    expand('results/{sample}/{sample}_genes_read_quantification.tsv', sample=SAMPLES)
    ```

    The Snakefile should like this:
    ```python linenums="1"
    # Sample list
    SAMPLES = ['highCO2_sample1', 'highCO2_sample2', 'highCO2_sample3', 'lowCO2_sample1', 'lowCO2_sample2', 'lowCO2_sample3']

    # Master rule used to launch workflow
    rule all:
        '''
        Dummy rule to automatically generate required outputs.
        '''
        input:
            expand('results/{sample}/{sample}_genes_read_quantification.tsv', sample=SAMPLES)
    ```

    You can launch a dry-run with the same command as before:
    ```sh
    snakemake -c 4 -F -p -n
    ```

    You should see rule `all` requiring 6 inputs and all the other rules appearing 6 times (1 for each sample):
    ```sh
    rule all:
        input: results/highCO2_sample1/highCO2_sample1_genes_read_quantification.tsv, results/highCO2_sample2/highCO2_sample2_genes_read_quantification.tsv, results/highCO2_sample3/highCO2_sample3_genes_read_quantification.tsv, results/lowCO2_sample1/lowCO2_sample1_genes_read_quantification.tsv, results/lowCO2_sample2/lowCO2_sample2_genes_read_quantification.tsv, results/lowCO2_sample3/lowCO2_sample3_genes_read_quantification.tsv
        jobid: 0
        reason: Forced execution
        resources: tmpdir=/tmp

    Job stats:
    job                           count
    --------------------------  -------
    all                               1
    fastq_trim                        6
    read_mapping                      6
    reads_quantification_genes        6
    sam_to_bam                        6
    total                            25
    ```

But there is an even better solution! At the moment, samples are defined as a list in the Snakefile. Processing other samples still means having to locate and change some chunks of code. To further improve workflow usability, samples can instead be defined in config files, so they can easily be added, removed, or modified by users without actually modifying code.

**Exercise:** Implement a parameter in the config file to specify sample names and modify the rule `all` to use this parameter instead of the `SAMPLES` variable in the `expand()` syntax.

??? success "Answer"
    First, we need to add the sample names to the config file:
    ```yaml linenums="1"
    # Configuration options of RNAseq-analysis workflow

    # Location of genome indices
    index: 'resources/genome_indices/Scerevisiae_index'

    # Location of annotation file
    annotations: 'resources/Scerevisiae.gtf'

    # Sample names
    samples:
      - highCO2_sample1
      - highCO2_sample2
      - highCO2_sample3
      - lowCO2_sample1
      - lowCO2_sample2
      - lowCO2_sample3
    ```

    Then, we need to remove `SAMPLES` from the Snakefile and use the config file in `expand()` instead:
    ```python linenums="1"
    # Master rule used to launch workflow
    rule all:
        '''
        Dummy rule to automatically generate required outputs.
        '''
        input:
            expand('results/{sample}/{sample}_genes_read_quantification.tsv', sample=config['samples'])
    ```
    Here, `config['samples']` is a Python list containing strings, each string being a sample name. This is because a list of parameters become a list when the config file is parsed. You can launch a dry-run with the same command as before (`snakemake -c 4 -F -p -n`) and you should see the same jobs.

??? info "An even more Snakemake-idiomatic solution"
    There is an even better and more Snakemake-idiomatic version of the `expand()` syntax in rule `all`:
    ```python
    expand(rules.reads_quantification_genes.output.gene_level, sample=config['samples'])
    ```
    This entirely removes the need to write output paths, even though it might be less easy to understand at first sight.

**Exercise:** Run the workflow on the other samples and generate the workflow DAG and filegraph. If you implemented parallelisation and multithreading in all the rules, the execution should take less than 10 min in total to process all the samples, otherwise it will be a few minutes longer.

??? success "Answer"
    You can run the workflow by removing `-F` and `-n` from the dry-run command, which makes a very simple command:
    ```sh
    snakemake -c 4 -p --sdm apptainer
    ```

    To generate the DAG, you can use:
    ```sh
    snakemake -c 1 -F -p --dag | dot -Tpng > images/all_samples_dag.png
    ```
    <p align="center">
      <img src="../../../assets/images/all_samples_dag.png"/>
    </p>

    If needed, open the picture in a new tab to zoom in. Then, you can generate the filegraph with:
    ```sh
    snakemake -c 1 -F -p --filegraph | dot -Tpng > images/all_samples_filegraph.png
    ```
    <p align="center">
      <img src="../../../assets/images/all_samples_filegraph.png" width="50%"/>
    </p>
    You probably noticed that these two figures have an extra rule, `fastq_qc_sol4`. It is the rule implemented in the supplementary exercise below.

### Optimising resource usage in a workflow by multithreading

!!! warning "This part is about CPU usage in Snakemake. It is quite long, so do it only if you finished all the other exercises."

When working with real, larger, datasets, some processes can take a long time to run. Fortunately, computation time can be decreased by running jobs in parallel and using several [threads](https://en.wikipedia.org/wiki/Thread_(computing)) or [cores](https://en.wikipedia.org/wiki/Multi-core_processor) for a single job.

**Exercise:** What are the two things you need to add to a rule to enable multithreading?

??? success "Answer"
    You need to add:

    1. The `threads` directive to tell Snakemake that it needs to allocate several threads to this rule
    1. The software-specific parameters in the `shell` directive to tell a software that it can use several threads

    If you add only the first element, the software will not be aware of the number of threads allocated to it and will use its default number of threads (usually one). If you add only the second element, Snakemake will allocate only one thread to the software, which means that it will run slowly or crash (the software expects multiple threads but gets one).

Usually, you need to read documentation to identify which software can make use of multithreading and which parameters control multithreading. We did it for you to save some time:

* `atropos trim`, `hisat2`, `samtools view`, and `samtools sort` can parallelise with the `--threads <nb_thread>` parameter
* `featureCounts` can parallelise with the `-T <nb_thread>` parameter
* `samtools index` can't parallelise. Remember that multithreading only applies to software that were developed to this end, Snakemake itself cannot parallelise a software!

Unfortunately, there is no easy way to find the optimal number of threads for a job. It depends on tasks, datasets, software, resources you have at your disposal... It often takes of few rounds of trial and error to see what works best. We already decided the number of threads to use for each software:

* 4 threads for `hisat2`
* 2 for all the other software

**Exercise:** Implement multithreading in a rule of your choice (it's usually best to start by multithreading the longest job, here `read_mapping`, but the example dataset is small, so it doesn't really matter).

??? tip "`threads` placeholder"
    `threads` can also be replaced by a `{threads}` placeholder in the `shell` directive.

??? success "Answer"
    We implemented multithreading in all the rules so that you can check everything. Feel free to copy this in your Snakefile:
    ```python linenums="1" hl_lines="16 23 43 50 68 74 76 95 105"
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
        threads: 2  # Add directive
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
        params:
            index = config['index']
        threads: 4  # Add directive
        container:
            'https://depot.galaxyproject.org/singularity/hisat2%3A2.2.1--hdbdd923_6'
        shell:
            '''
            echo "Mapping the reads" > {log}
            hisat2 --dta --fr --no-mixed --no-discordant --time --new-summary --no-unal \
            -x {params.index} --threads {threads} \
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
        threads: 2  # Add directive
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
        threads: 2  # Add directive
        params:
            annotations = config['annotations']
        container:
            'https://depot.galaxyproject.org/singularity/subread%3A2.0.6--he4a0461_2'
        shell:
            '''
            echo "Counting reads mapping on genes in <{input.bam_once_sorted}>" > {log}
            featureCounts -t exon -g gene_id -s 2 -p --countReadPairs \
            -B -C --largestOverlap --verbose -F GTF \
            -a {params.annotations} -T {threads} -o {output.gene_level} {input.bam_once_sorted} &>> {log} 
            echo "Renaming output files" >> {log}
            mv {output.gene_level}.summary {output.gene_summary}
            echo "Results saved in <{output.gene_level}>" >> {log}
            echo "Report saved in <{output.gene_summary}>" >> {log}
            '''
    ```

**Exercise:** Do you need to change anything in the `snakemake` command to run the workflow with 4 cores?

??? success "Answer"
    You need to provide additional cores to Snakemake with the parameter `-c 4`. Using the same sample as before (`highCO2_sample1`), the workflow can be run with:
    ```sh
    snakemake -c 4 -F -p results/highCO2_sample1/highCO2_sample1_genes_read_quantification.tsv --sdm apptainer
    ```

    The number of threads allocated to all jobs running at a given time cannot exceed the value specified with `--cores`, so if you use `-c 1`, Snakemake will not be able to use multiple threads. Conversely, if you ask for more threads in a rule than what was provided with `--cores`, Snakemake will cap rule threads at `--cores` to avoid requesting too many. Another benefit of increasing `--cores` is to allow Snakemake to run multiple jobs in parallel (for example, here, running two jobs using two threads each).

If you run the workflow from scratch with multithreading in all rules, it should take ~6 min, compared to ~10 min before (_i.e._ a 40% decrease!). This gives you an idea of how powerful multithreading is when datasets and computing power get bigger!

??? warning "Things to keep in mind when using parallel execution"
    * On-screen output from parallel jobs will be mixed, so save any output to log files instead
    * Parallel jobs use more RAM. If you run out then either your OS will swap data to disk (which slows data access), or a process will die (which can crash Snakemake)
    * Parallelising is not without consequences and has a cost. This is a topic too wide for this course, but just know that using too many cores on a dataset that is too small can slow down computation, as explained [here](https://stackoverflow.com/questions/45256953/why-is-multiprocess-pool-slower-than-a-for-loop)
