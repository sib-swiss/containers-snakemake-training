## Learning outcomes

**After having completed this chapter you will be able to:**

* Use an input function to work with an unknown number of files
* Run scripts from other languages (Python and R)
* Deploy a rule-specific conda environment
* Deploy a rule-specific Docker/Singularity container

## Material

[:fontawesome-solid-file-pdf: Download the presentation](../../assets/pdf/day2/5_reproducibility_snakemake.pdf){: .md-button }

## Exercises

In this series of exercises, you will create the last two rules of the workflow. Each rule will execute a script (one in Python and one in R; don't worry, this is not a programming course, so we wrote the scripts for you!), and both rules will have dedicated environments that you will need to take into account in the snakefiles.

??? tip "Snakefile from previous session"
    If you want to restart from a fully commented workflow (including a supplementary rule for .fastq files quality check), with multithreading and resource usage control implemented in all rules, you can manually get it [here](https://github.com/sib-swiss/containers-snakemake-training/blob/2024_update/docs/solutions_day2/session3) or download it in your current directory with:
    ```
    wget https://raw.githubusercontent.com/sib-swiss/containers-snakemake-training/refs/heads/2024_update/docs/solutions_day2/session2/workflow/Snakefile
    ```

??? tip "Development and back-up"
    During this session, you will modify your Snakefile quite heavily, so it may be a good idea to make back-ups from time to time (with `cp` or a simple copy/paste) or use a versioning system. As a general rule, if you have a doubt on the code you are developing, do not hesitate to make a back-up beforehand.

### Creating a rule to gather read counts

To perform a Differential Expression Analysis (DEA), it is easier to have a single file gathering all the read counts from the different samples. The next rule that you will create will both find the required files and merge them, thanks to an input function and a Python script.

#### Building the overall rule structure

We already wrote the common elements of the rule so that you can focus on the most interesting parts (the missing `input` and the missing elements at the end):

```python
rule count_table:
    '''
    This rule merges gene count tables of an assembly into one table.
    '''
    input:
        ?
    output:
        count_table = 'results/total_count_table.tsv'
    log:
        'logs/total_count_table.log'
    benchmark:
        'benchmarks/total_count_table.txt'
    resources:
        mem_mb = 500
    threads: 1
    ?:
        ?
    ?:
        ?
```

??? tip "Explicit is better than implicit"
    Even if a software cannot multithread, it is useful to add `threads: 1` to keep the syntax consistent between rules and clearly state that the software works with a single thread.

**Exercise:** Given that this rule and the next one will be quite different from the previous ones, it is a good idea to implement them in a new snakefile. Create a new file `workflow/rules/analysis.smk` and copy the previous rule structure in it. Do you need to change anything else to use this rule in your workflow?

??? success "Answer"
    To actually use this rule, Snakemake needs to be aware that it exists: this is done with the `include` statement. We need to add the following line to the main Snakefile:
    ```
    include: 'rules/analyses.smk'
    ```
    The Snakefile should now look like this:
    ```python
    '''
    Main Snakefile of RNAseq analysis workflow. This workflow can clean and
    map reads, and perform Differential Expression Analyses.
    '''

    # Config file path
    configfile: 'config/config.yaml'

    # Rules to execute workflow
    include: 'rules/read_mapping.smk'
    include: 'rules/analyses.smk'

    # Master rule used to launch workflow
    rule all:
        '''
        Dummy rule to automatically generate required outputs.
        '''
        input:
            expand(rules.reads_quantification_genes.output.gene_level, sample=config['samples'])
    ```

Now, let's start filling those missing elements!

#### Gathering input files

This task is quite complex: we need a way to identify all the rule inputs and gather them in a Python list. Here, there are only six samples, so in theory, you could list them directly... However, it isn't good practice and it quickly becomes unmanageable when the number of sample increases. Fortunately, there is a much more elegant and convenient solution: an input function, which provides the added benefit of scaling up very well.

We wrote one for you:
```python
# Input function used in rule count_table
def get_gene_counts(wildcards):
    '''
    This function lists count tables from every sample in the config file
    '''
    return [f"results/{sample}/{sample}_genes_read_quantification.tsv"
            for sample in config['samples']]
```
??? warning "Snakemake wildcard vs Python f-string"
    This input function is pure Python code: in the return statement, {sample} isn't a wildcard, it is an f-string variable! This shows that you can natively use Python constructions in a workflow: Snakemake will still be able compute them. This is because Snakemake was built on top of Python.

This function will loop over the list of samples in the config file and replace `{sample}` with the current sample name of the iteration to create a string which is the output path from the rule `reads_quantification_genes` of said sample. Then, it will aggregate all the paths in a list and return this list.

??? info "To go further on input functions"
    * Input functions take the `wildcards` global object as **single argument**
    * You can access wildcard values inside an input function with the syntax `{wildcards.wildcards_name}`
    * Input functions can return a list of files, which will then be automatically handled like multiple inputs by Snakemake
        * Input functions can also return a dictionary; in this case, the function should be called with the syntax:
        ```
        input: unpack(<function_name>)
        ```
        The dictionary keys will be interpreted as input names and the dictionary values will be interpreted as input values, providing a list of named inputs

??? warning "Input functions and output directory"
    Input functions are evaluated **before** the workflow is executed. As a consequence, they cannot be used to list the content of an output directory, since the directory does not exist before the workflow is executed!

**Exercise:** Insert the function `get_gene_counts()` in `workflow/rules/analysis.smk` and adapt the input value of `count_table` accordingly. Do you need to insert the function in a specific location?

??? success "Answer"
    First step is to add the input function to the file. However, it needs to appear **before the rule `count_table`**, otherwise you will see the error `name 'get_gene_counts' is not defined`. In other words, the function needs to be defined **before** Snakemake looks for it when it parses the input. Then, you need to set the function name as the rule input value.

    `workflow/rules/analysis.smk` should now resemble this:
    ```python
    # Input function used in rule count_table
    def get_gene_counts(wildcards):
        '''
        This function lists count tables from every sample in the config file
        '''
        return [f"results/{sample}/{sample}_genes_read_quantification.tsv"
                for sample in config['samples']]

    rule count_table:
        '''
        This rule merges gene count tables of an assembly into one table.
        '''
        input:
            get_gene_counts
        output:
            count_table = 'results/total_count_table.tsv'
        log:
            'logs/total_count_table.log'
        benchmark:
            'benchmarks/total_count_table.txt'
        conda:
            '../envs/py.yaml'
        resources:
            mem_mb = 500
        threads: 1
        ?:
            ?
        ?:
            ?
    ```
    An important thing to remember is that you don't need to use parentheses or specify any argument when you call an input function in the `input` directive.


<!--

AT. 0 exercises

Check usage of we/you
Check usage of the
Remove 'note' usage

The checked up to ### Running the workflow on the other samples

Add this somewhere:
Note that we do not need to use `wildcards` in the output name: only one file will be created, and its name will not change depending on the sample name, because we use **all** the samples to create it.

 -->




### TMP


The counts will be concatenated thanks to a Python script available




**Exercise:** Implement a rule to list and merge read count files (coming from rule `reads_quantification_genes`) into a single file using the Python script provided [here](https://raw.githubusercontent.com/sib-swiss/containers-snakemake-training/main/scripts/solutions/day2/session4/workflow/scripts/count_table.py).

??? info "Information on the script to compute the table"
    * You can download the script with `wget https://raw.githubusercontent.com/sib-swiss/containers-snakemake-training/main/docs/solutions_day2/session4/workflow/scripts/count_table.py`
    * It goes in a special directory in your workflow
    * It is written in Python
    * It takes a list of files as input, each file being the read count output of featureCounts
    * It produces one output, a tab-separated table containing all the read counts of the different sample gathered by gene

??? tip
    While the goal of this rule is quite easy to grasp, setting it up requires using several advanced notions of Snakemake, so here is a little outline of the steps you should take:

    1. Think about the directive you want to use to run the script
        * Looking at the script length with `wc -l <path/to/script>` could help you decide
    1. Think about the location/path of the script
    1. Check the beginning of the script to see if you need any special Python packages. You can do that with `head <path/to/script>`. If you see lines containing `import <package_name>`, it means that the script is using external Python packages
        * If the script is using external packages, think on how can you provide them

Now, let's solve these problems one by one!



#### Getting the Python script and running it

**Exercise:** Download the script and place it the proper folder. Remember that per the [official documentation](https://snakemake.readthedocs.io/en/v8.20.5/snakefiles/deployment.html), scripts should be stored in a subfolder `workflow/scripts`.

??? success "Answer"
    ```sh
    wget https://raw.githubusercontent.com/sib-swiss/containers-snakemake-training/main/scripts/solutions/day2/session4/workflow/scripts/count_table.py  # Download script
    mkdir -p workflow/scripts  # Create appropriate folder
    mv count_table.py workflow/scripts  # Move script in the newly created folder
    ```

**Exercise:** Check how long the script is.

??? success "Answer"
    If you run `wc -l workflow/scripts/count_table.py`, you will see that the script is 67 lines long. This is too much to use a `run` directive, so we will use the `script` directive instead. This means that we need to add the following to our rule `count_table`:
    ```python
    script:
        '../scripts/count_table.py'
    ```

??? info "Script path"
    If you placed included files in subfolders (like `rules/analysis.smk`), you need to change relative paths for external script files, hence the `../` in the script path.

??? tip
    In many cases, it would be nice to have a script that can be called by Snakemake but also work with standard Python, so that the code can be reused in other projects. There are several ways to do that:

    * You could implement most of the functionalities in a module and use this module in a simple script called by Snakemake
    * You could test for the existence of a `snakemake` object and handle parameter values differently (_e.g._ command-line arguments) if the object does not exist

    Inside the script, you have access to an object `snakemake` that provides access to the same objects that are available in the `run` and `shell` directives (`input`, `output`, `params`, `wildcards`, `log`, `threads`, `resources`, config), _e.g._ you can use `snakemake.input[0]` to access the first input file of a rule, or `snakemake.input.input_name` to access a specific named input.

**Exercise:** Check the script content to see whether it requires specific packages.

??? success "Answer"
    If you run `head workflow/scripts/count_table.py`, you will see several `import` commands at the start of the script, including `import pandas as pd`. [`pandas`](https://pandas.pydata.org/docs/index.html) is a great package, but it is not part of the default packages natively shipped with Python. This means that we need to find a solution to provide it to the rule. The easiest way to do that is to create a conda environment dedicated to the rule. Conda environments should be stored in a subfolder `workflow/envs`.

    Create the appropriate folder: `mkdir -p workflow/envs`. Then, write the following configuration in the environment file, `workflow/envs/py.yaml`:
    ```yaml
    # Environment file to perform data processing with python
    name : python
    channels :
        - conda-forge
        - bioconda
    dependencies :
        - python >= 3.10
        - pandas == 1.4.3
    ```

    This means that you need to add the following to your rule `count_table`:
    ```python
    conda:
        '../envs/py.yaml'
    ```

??? info "Environment file path"
    If you placed included files in subfolders (like `rules/analysis.smk`), you need to change relative paths for conda environments files as well, hence the `../` in the environment file path.

Using conda environments improves reproducibility for many reasons, including version control and the fact that users do not need to manually manage software dependencies. **Note: the first execution of the workflow after adding Conda environments will take some time, because Conda will have to download and install all the software**.
<!--
#### Identifying and listing the input files

The before-last step is the most complex one: identifying all the inputs of the rule and gathering them in a list. Here, there are only six samples, so in theory, you could list them directly... However, by now you should that is not a good solution. Fortunately, there is a much more elegant and convenient way to do this: an input function, which provides the added benefit of scaling up very well if the number of samples increase.

We already wrote the input function for you:

```python
# Input function used in rule count_table
def get_gene_counts(wildcards):
    '''
    This function lists count tables from every sample in the config file
    '''
    # Note that here {sample} is not a wildcard, it is an f-string variable!
    return [f"results/{sample}/{sample}_genes_read_quantification.tsv"
            for sample in config['samples']]
```

This function will loop over the list of samples in the config file, replace `{sample}` with the current sample name of the iteration to create a string which is the output path from the rule `reads_quantification_genes` of said sample. Then, it will aggregate all the paths in a list and return this list.

??? info "Details on input functions"
    * Input functions take the `wildcards` global object as **single argument**
    * You can access wildcard values inside an input function with the syntax `{wildcards.wildcards_name}`
    * Input and output functions can return a list of files, which will then be automatically handled like multiple inputs or outputs by Snakemake. These functions can also return a dictionary; in this case, the function should be called with the syntax `input: unpack(<function_name>)`. The dictionary's keys will be interpreted as input/output names and the dictionary's values will be interpreted as input/output values
    * Functions are evaluated **before** the workflow is executed. As a consequence, they cannot be used to list the content of an output directory, since the directory does not exist before the workflow is executed!

**Exercise:** Insert the function in the proper Snakefile and adapt the input value of the rule accordingly.

??? success "Answer"
    There are two things to do:

    1. Insert the input function in `workflow/rules/analysis.smk`, before the rule, otherwise you will get a `name 'get_gene_counts' is not defined` error (the function needs to be defined **before** Snakemake looks for it when it parses the rule input)
    1. Use the function name as value for the input directive

    Your rule and function should resemble this:
    ```yaml
    # Input function used in rule count_table
    def get_gene_counts(wildcards):
        '''
        This function lists count tables from every sample in the config file
        '''
        return [f"results/{sample}/{sample}_genes_read_quantification.tsv"
                for sample in config['samples']]
    ```

    ```python
    rule count_table:
        '''
        This rule merges gene count tables of an assembly into one table.
        '''
        input:
            get_gene_counts
        output:
            count_table = 'results/total_count_table.tsv'
        log:
            'logs/total_count_table.log'
        benchmark:
            'benchmarks/total_count_table.txt'
        conda:
            '../envs/py.yaml'
        resources:
            mem_mb = 500
        threads: 1
        script:
            '../scripts/count_table.py'
    ```

??? info "Input functions and parentheses"
    You don't need to use parentheses or specify any argument when you call an input function in the `input` directive.
 -->
#### Adapting the Snakefile and running the rule

Now, all that is left is to run the rule to create the table.

**Exercise:** Which command should you use to create the output? Is there anything else to do beforehand?

??? success "Answer"
    It turns out that we cannnot launch the workflow directly: we need to include the new rule file in the Snakefile and adapt the output of the rule `all`! Your Snakefile should now resemble this:

    ```python
    '''
    Main Snakefile of the RNAseq analysis workflow. This workflow can clean and
    map reads, and perform Differential Expression Analyses.
    '''

    # Config file path
    configfile: 'config/config.yaml'

    # Rules to execute workflow
    include: 'rules/read_mapping.smk'
    include: 'rules/analyses.smk'

    # Master rule used to launch workflow
    rule all:
        '''
        Dummy rule to automatically generate required outputs.
        '''
        input:
            'results/total_count_table.tsv'
    ```

    Finally, run the workflow with `snakemake -c 4 -p --use-conda`

??? info "--use-conda"
    Do not forget `--use-conda` otherwise Snakemake will not use the environments!!!

### Creating a rule to detect Differentially Expressed Genes

It is now time to write the final rule of the workflow. This rule will perform the DEA using the global count table you previously created.

**Exercise:** Implement a rule to perform DEA using the R script provided [here](https://raw.githubusercontent.com/sib-swiss/containers-snakemake-training/main/scripts/solutions/day2/session4/workflow/scripts/DESeq2.R).

??? info "Information on the script to compute the table"
    * You can download the script with `wget https://raw.githubusercontent.com/sib-swiss/containers-snakemake-training/main/scripts/solutions/day2/session4/workflow/scripts/DESeq2.R`
    * It goes in a special directory in your workflow
    * It is written in R
    * It takes a global read counts table as input
    * It produces two outputs:
        1. A tab-separated table containing the DEG and the associated statistical results
        1. A pdf file containing control plots of the analysis

??? tip
    While not being trivial, this rule is much easier than the previous one and some things work similarly. Still, here is a little outline of the steps you should take:

    1. Build the basic structure of your rule: name, input, outputs, log, benchmark
        * Memory should be set at 1 GB
        * Threads should be set at 2
    1. Think about the directive you want to use to run the script
        * Remember that there is only one way to run an R script
    1. Think about the location/path of the script
    1. This script is using a lot of external packages. Fortunately, all these packages are available in a _certain Docker image you worked with yesterday_

Now, let's write this last rule!

#### Building the rule structure

You have done that a few times already, so it should not be too difficult.

**Exercise:** Set the input, outputs, log, benchmark, resources and thread values.

??? success "Answer"
    ```python
    rule differential_expression:
        '''
        This rule detects DEGs and plots associated visual control graphs (PCA,
        heatmaps...).
        '''
        input:
            table = rules.count_table.output.table
        output:
            deg = 'results/deg_list.tsv',
            pdf = 'results/deg_plots.pdf'
        log:
            'logs/differential_expression.log'
        benchmark:
            'benchmarks/differential_expression.txt'
        ?:
            ?
        resources:
            mem_gb = 1
        threads: 2
        ?:
            ?

    ```

    Note that we do not need to use `wildcards` in this rule, because all the files are precisely defined.

#### Getting the R script and running it

**Exercise:** Download the script and place it the proper folder. Remember that per the [official documentation](https://snakemake.readthedocs.io/en/v8.20.5/snakefiles/deployment.html), scripts should be stored in a subfolder `workflow/scripts`.

??? success "Answer"
    ```sh
    wget https://raw.githubusercontent.com/sib-swiss/containers-snakemake-training/main/scripts/solutions/day2/session4/workflow/scripts/DESeq2.R  # Download the script
    mv DESeq2.R workflow/scripts  # Move the script in the newly created folder
    ```

**Exercise:** Find a way to run the script.

??? success "Answer"
    There is only one way to run an R script: use the `script` directive. This means that we need to add the following to our rule `differential_expression`:
    ```python
    script:
        '../scripts/DESeq2.R'
    ```

??? info "Script path"
    If you placed included files in subfolders (like `rules/analysis.smk`), you need to change relative paths for external script files, hence the `../` in the script path.

??? tip
    Inside the script, an S4 object named `snakemake` analogous to the Python case available and allows access to input and output files and other parameters. Here the syntax follows that of S4 classes with attributes that are R lists, _e.g._ you can access the first input file with `snakemake@input[[1]]` (note that the first file does not have index 0 here, because R starts counting from 1). Named input and output files can be accessed in the same way, by just providing the name instead of an index, _e.g._ `snakemake@input[["myfile"]]`.

**Exercise:** Find an efficient way to create a computing environment for the rule.

??? tip
    Remember what you did during Day 1, session 3 "Working with Dockerfiles"!

??? success "Answer"
    During Day 1, you built your own docker image, called deseq2. This image actually contains everything we need to run DEA, so let's use it again, but with Snakemake this time! This means that you need to add the following to your rule `differential_expression`:
    ```python
    container:
        'docker://geertvangeest/deseq2:v1'
    ```

??? info "Your own Docker image"
    First try with your own image. If it doesn't work, then you can use Geert's image: `geertvangeest/deseq2:v1`.

After all these modifcations, this is what your final rule should look like:

??? success "Answer"
    ```python
    rule differential_expression:
        '''
        This rule detects DEGs and plots associated visual control graphs (PCA,
        heatmaps...).
        '''
        input:
            table = rules.count_table.output.table
        output:
            deg = 'results/deg_list.tsv',
            pdf = 'results/deg_plots.pdf'
        log:
            'logs/differential_expression.log'
        benchmark:
            'benchmarks/differential_expression.txt'
        container:
            'docker://geertvangeest/deseq2:v1'
        resources:
            mem_gb = 1
        threads: 2
        script:
            '../scripts/DESeq2.R'
    ```

#### Adapting the Snakefile and running the rule

Now, all that is left is to run the rule to create the DEG list.

**Exercise:** Which command should you use to create the output? Is there anything else to do beforehand?

??? success "Answer"
    It turns out that we cannnot launch the workflow directly: we need to include the new rule file in the Snakefile and adapt the output of the rule `all`! Your Snakefile should now resemble this:

    ```python
    '''
    Main Snakefile of the RNAseq analysis workflow. This workflow can clean and
    map reads, and perform Differential Expression Analyses.
    '''

    # Config file path
    configfile: 'config/config.yaml'

    # Rules to execute workflow
    include: 'rules/read_mapping.smk'
    include: 'rules/analyses.smk'

    # Master rule used to launch workflow
    rule all:
        '''
        Dummy rule to automatically generate required outputs.
        '''
        input:
            'results/deg_list.tsv'
    ```

    Finally, run the workflow with `snakemake -c 4 -p --use-singularity`. You should see new Snakemake information messages:
    ```
    Rscript --vanilla /path/to/snakemake_rnaseq/.snakemake/scripts/tmpge97d_lz.DESeq2.R
    Activating singularity image /path/to/snakemake_rnaseq/.snakemake/singularity/8bfdbe93244feb95887ab5d33a705017.simg
    INFO:    squashfuse not found, will not be able to mount SIF
    INFO:    fuse2fs not found, will not be able to mount EXT3 filesystems
    INFO:    gocryptfs not found, will not be able to use gocryptfs
    INFO:    Converting SIF file to temporary sandbox...
    INFO:    Cleaning up image...
    ```

??? info "--use-singularity"
    Do not forget `--use-singularity` otherwise Snakemake will not pull the image!!!

??? tip
    If you want to see how a Snakemake-generated Dockerfile looks like, use: `snakemake -c 1 --containerize > Dockerfile`.

**Exercise:** How many DEG were detected?

??? success "Answer"
    Have a look at the list that was just created: `cat results/deg_list.tsv`. 8 genes are differentially expressed!

**Exercise:** If you had to re-run the entire workflow from scratch, what command would you use?

??? success "Answer"
    You would need to execute `snakemake -c 4 -p --use-conda --use-singularity -F`.
    * `-F` is to force the execution of the entire workflow
    * Don't forget `--use-conda --use-singularity`! Otherwise, you will lack some software and packages and the workflow will crash!

**Exercise:** Visualise the DAG of the entire workflow.

??? success "Answer"
    You should now be used to this. `snakemake -c 1 -p -F --dag | dot -T png > images/total_dag.png`

This is the DAG you should see:

<figure align="center">
  <img src="../../../assets/images/total_dag.png" width="100%"/>
</figure>

Congratulations, you are now able to create a Snakemake workflow and make it reproducible thanks to conda/mamba and Docker/Singularity! To make things even better, have a look at [Snakemake's best practices](https://snakemake.readthedocs.io/en/v8.20.5/snakefiles/best_practices.html)!
