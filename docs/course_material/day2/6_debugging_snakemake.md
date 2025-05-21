## Learning outcomes

**After having completed this chapter you will be able to:**

* Understand the inner workings and order of operations in Snakemake
* Efficiently Design and debug a Snakemake workflow
* Use flags to delete or protect outputs
* Use flags to consider a directory as output

## Material

[:fontawesome-solid-file-pdf: Download the presentation](../../assets/pdf/day2/6_additional_concepts.pdf){: .md-button }

## Designing a workflow

There are many ways to design a new workflow, but these few pieces of advice will be useful in most cases:

* Start with a pen and paper: try to find out how many rules you will need and how they depend on each other. In other terms, start by sketching the DAG of your workflow!
    * Remember that Snakemake has a bottom-up approach (it goes from the final outputs to the first input), so it may be easier for you to work in that order as well and write your last rule first
    * Determine which rules (if any) aggregate or split inputs and create input functions accordingly (the topic is tackled in [this](5_reproducibility_snakemake.md) series of exercises)
* Make sure your input and output directives are right before worrying about anything else, especially the shell sections. There is no point in executing commands with the wrong inputs/outputs!
    * Remember that Snakemake builds the DAG before running the shell commands, so you can use the `-n/--dry-run/--dryrun` parameter to test the workflow before running it. You can even do that without writing all the shell commands!
* List any parameters or settings that might need to be adjusted later
* Choose meaningful and easy-to-understand names for your rules, inputs, outputs, parameters, `wildcards`... to make your Snakefile as readable as possible. This is true for every script, piece of code, variable... and Snakemake is no exception! Have a look at [The Zen of Python](https://peps.python.org/pep-0020/) for more information

## Debugging a workflow

It is very likely you will see bugs and errors the first time you try to run a new Snakefile: don’t be discouraged, this is normal!

### Order of operations in Snakemake

The topic was approached when we discussed [DAGs](2_introduction_snakemake.md#chaining-rules), but to efficiently debug a workflow, it is worth taking a deeper look at what Snakemake does when you execute the command `snakemake -c 1 <target>`. There are 3 main phases:

1. Prepare to run:
    1. Read all the rule definitions from the Snakefile
1. Resolve the DAG (happens when Snakemake says 'Building DAG of jobs'):
    1. Check what output(s) are required
    1. Look for matching input(s) by looking at the outputs of all the rules
    1. Fill in the `wildcards` to determine the exact input(s) of the matching rule(s)
    1. Check whether this(these) input(s) is(are) available; if not, repeat Step 2 until everything is resolved
1. Run:
    1. If needed, create the output(s) folder path
    1. If needed, remove the outdated output(s)
    1. Fill in the placeholders and run the shell commands
    1. Check that the commands ran without errors and produced the expected output(s)

### Debugging advice

Sometimes, Snakemake will give you a precise error report, but other times... less so. Try to identify which phase of execution failed (see previous paragraph on [order of operations](6_debugging_snakemake.md#order-of-operations-in-snakemake)) and double-check the most common error causes for that phase:

1. Parsing phase failures (phase 1):
    * Syntax errors, among which (but not limited to):
        * Missing commas/colons/semicolons
        * Unbalanced quotes/brackets/parenthesis
        * Wrong indentation
        * _These errors can be easily solved using a text editor with Python/Snakemake text colouring_
    * Failure to evaluate expressions
        * Problems in functions (`expand()`, input functions...) in input/output directives
        * Python logic added outside of rules
    * Other problems with rule definition
        * Invalid rule names/directives
        * Invalid wildcard names
        * Mismatched `wildcards`
1. DAG building failures (phase 2, before Snakemake tries to run any job):
    * Failure to determine the target
    * Ambiguous rules making the same output(s)
    * On the contrary, no rule making the required output(s)
    * Circular dependency (violating the **'Acyclic'** property of a D**A**G).
    * Write-protected output(s)
1. DAG running failures (phase 3, `--dry-run` works and builds the DAG, but the jobs execution fails):
    * Shell command returning non-zero status
    * Missing output file(s) after the commands have run
    * Reference to a `$shell_variable` before it was set
    * Use of a wrong/unknown placeholder inside `{ }`
    * _When a job fails, Snakemake reports an error, deletes all output file(s) for that job (because of potential corruption), and stops_

## Using non-conventional outputs

_This part is an extra-exercise about non-conventional Snakemake outputs. It is quite long, so do it only if you finished the rest of the course beforehand._

Snakemake has several built-in utilities to assign properties to outputs that are deemed 'special'. These properties are listed in the table below:

| Property  |              Syntax              |                                                    Function                                                   |
| --------- | -------------------------------- | ------------------------------------------------------------------------------------------------------------- |
| Temporary | `temp('file.txt')`       | File is deleted as soon as it is not required by a future jobs                                              |
| Protected | `protected('file.txt')`  |File cannot be overwritten after job ends (useful to prevent erasing a file by mistake, for example files requiring heavy computation)                   |
| Ancient   | `ancient('file.txt')`    | Ignore file timestamp and assume file is older than any outputs: file will not be re-created when re-running workflow, except when `--force` options are used       |
| Directory | `directory('directory')` | Output is a directory instead of a file (use `touch()` instead if possible)                                     |
| Touch     | `touch('file.txt')`      | Create an empty flag file `file.txt` regardless of shell commands (only if commands finished without errors) |

The next paragraphs will show how to use some of these properties.

### Removing and safeguarding outputs

This part shows a few examples using `temp()` and `protected()` flags.

**Exercise:** Can you think of a convenient use of the `temp()` flag?

??? success "Answer"
    `temp()` is extremely useful to automatically remove intermediary outputs that are no longer needed.

**Exercise:** In your workflow, identify outputs that are intermediary and mark them as temporary with `temp()`.

??? success "Answer"
    Unsorted .bam and .sam outputs seem like great candidates to be marked as temporary. One could argue that trimmed .fastq files could also be marked as temporary, but we won't bother with them here. For clarity, only lines that changed are shown below:

    * rule `read_mapping`:
    ```python linenums="1"
    output:
        sam = temp('results/{sample}/{sample}_mapped_reads.sam'),
    ```

    * rule `sam_to_bam`:
    ```python linenums="1"
    output:
        bam = temp('results/{sample}/{sample}_mapped_reads.bam'),
    ```

??? warning "Advantages and drawbacks of `temp()`"
    On one hand, removing temporary outputs is a great way to save storage space. If you look at the size of your current `results/` folder (`du -bchd0 results/`), you will notice that it drastically increased. Just removing these two files would allow to save ~1 GB. While it may not seem like much, remember that you usually have much bigger files and many more samples!

    On the other hand, using temporary outputs might force you to re-run more jobs than necessary if an input changes, so think carefully before using `temp()`.

**Exercise:** On the contrary, is there a file of your workflow that you would like to protect? If so, mark it with `protected()`.

??? success "Answer"
    Sorted .bam files from rule `sam_to_bam` seem like good candidates for protection:
    ```python linenums="1"
    output:
        bam_sorted = protected('results/{sample}/{sample}_mapped_reads_sorted.bam'),
    ```
    If you set this output as protected, be careful when you want to re-run your workflow to recreate the file!

### Using an output directory: the FastQC example

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a program designed to spot potential problems in high-througput sequencing datasets. It is a very popular tool, notably because it is fast and does not require a lot of configuration. It runs a set of analyses on one or more sequence files in FASTQ or BAM format and produces a quality report with plots that summarise the results. It highlights potential problems that might require a closer look in your dataset.

As such, it would be interesting to run FastQC on the original and trimmed .fastq files to check whether trimming actually improved read quality. FastQC can be run interactively or in batch mode, during which it saves results as an HTML file and a ZIP file. You will later see that running FastQC in batch mode is a bit harder than it looks.

??? info "Data types and FastQC"
    FastQC does not differentiate between sequencing techniques and as such can be used to look at libraries coming from a large number of experiments (Genomic Sequencing, ChIP-Seq, RNAseq, BS-Seq etc...).

If you run `fastqc -h`, you will notice something a bit surprising, but not unusual in bioinformatics:
```sh
    -o --outdir     Create all output files in the specified output directory.
                    Please note that this directory must exist as the program
                    will not create it. If this option is not set then the
                    output file for each sequence file is created in the same
                    directory as the sequence file which was processed.

   -d --dir         Selects a directory to be used for temporary files written when
                    generating report images. Defaults to system temp directory if
                    not specified.
```

Two files are produced for each .fastq file and these files appear in the same directory as the input file: FastQC does not allow to specify the output file names! However, you can set an alternative output directory, even though **it needs to be manually created before FastQC is run**.

There are different solutions to this problem:

1. Work with the default file names produced by FastQC and leave the reports in the same directory as the input files
1. Create the outputs in a new directory and leave the reports with their default name
1. Create the outputs in a new directory and tell Snakemake that the directory itself is the output
1. Force a naming convention by manually renaming the FastQC output files within the rule

It would be too long to test all four solutions, so you will work on the 3<sup>rd</sup> **or** the 4<sup>th</sup> solution. Here is a brief summary of solutions 1 and 2:

1. This could work, but it's better not to put reports in the same directory as input sequences. As a general principle, when writing Snakemake rules, we prefer to be in charge of the output names and to have all the files linked to a sample in the same directory
1. This involves manually constructing the output directory path to use with the `-o` parameter, which works but isn't very convenient

The first part of the FastQC command is:
```sh
fastqc --format fastq --threads {threads} --outdir {folder_path} --dir {folder_path} <input_fastq1> <input_fastq2>
```

??? info "Explanation of FastQC parameters"
    * `-t/--threads`: specify how many files can be processed simultaneously. Here, it will be 2 because inputs are paired-end files
    * `-o/--outdir`: create output files in specified output directory
    * `-d/--dir`: select a directory to be used for temporary files

**Choose between solution 3 or 4 and implement it. You will find more information and help below.**

=== "Solution 3"
    This option is equivalent to tell Snakemake not to worry about individual files at all and consider an entire directory as the rule output.

    **Exercise:** Implement a single rule to run FastQC on both the original and trimmed .fastq files (four files in total) using directories as ouputs with the `directory()` flag.

    * The container image can be found at `https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0`

    ??? success "Answer"
        This makes the rule definition quite 'simple' compared to solution 4:
        ```python linenums="1" hl_lines="25 31"
        rule fastq_qc_sol3:
            '''
            This rule performs a QC on paired-end .fastq files before and after trimming.
            '''
            input:
                reads1 = rules.fastq_trim.input.reads1,
                reads2 = rules.fastq_trim.input.reads2,
                trim1 = rules.fastq_trim.output.trim1,
                trim2 = rules.fastq_trim.output.trim2
            output:
                before_trim = directory('results/{sample}/fastqc_reports/before_trim/'),
                after_trim = directory('results/{sample}/fastqc_reports/after_trim/')
            log:
                'logs/{sample}/{sample}_fastqc.log'
            threads: 2
            container:
                'https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0'
            shell:
                '''
                echo "Creating output directory <{output.before_trim}>" > {log}
                mkdir -p {output.before_trim} 2>> {log}  # FastQC doesn't create output directories so we have to do it manually
                echo "Performing QC of reads before trimming in <{input.reads1}> and <{input.reads2}>" >> {log}
                fastqc --format fastq --threads {threads} --outdir {output.before_trim} \
                --dir {output.before_trim} {input.reads1} {input.reads2} &>> {log}
                echo "Results saved in <{output.before_trim}>" >> {log}
                echo "Creating output directory <{output.after_trim}>" >> {log}
                mkdir -p {output.after_trim} 2>> {log}  # FastQC doesn't create output directories so we have to do it manually
                echo "Performing QC of reads after trimming in <{input.trim1}> and <{input.trim2}>" >> {log}
                fastqc --format fastq --threads {threads} --outdir {output.after_trim} \
                --dir {output.after_trim} {input.trim1} {input.trim2} &>> {log}
                echo "Results saved in <{output.after_trim}>" >> {log}
                '''
        ```

    ??? info " `.snakemake_timestamp`"
        When `directory()` is used, Snakemake creates an empty file called `.snakemake_timestamp` in the output directory. This is the marker file it uses to know whether it needs to re-run the rule producing the directory.

    Overall, this rule works quite well and allows for an easy rule definition. However, in this case, individual files are not explicitly named as outputs and this may cause problems to chain rules later. Also, remember that some software won’t give you any control at all over the outputs, which is why you need a back-up plan, _i.e._ solution 4: the most powerful solution is still to use shell commands to move and/or rename files to names you want. Also, the Snakemake developers advise to use `directory()` only as a last resort and to use [`touch()`](https://snakemake.readthedocs.io/en/v8.20.5/snakefiles/rules.html#flag-files) instead.

=== "Solution 4"
    This option amounts to let FastQC follows its default behaviour but manually rename the files afterwards to obtain the exact outputs we require.

    **Exercise:** Implement a single rule to run FastQC on both the original and trimmed .fastq files (four files in total) and rename the files created by FastQC to the desired output names using the `mv <old_name> <new_name>` command.

    * The container image can be found at `https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0`

    ??? success "Answer"
        This makes the rule definition (much) more complicated than solution 3:
        ```python linenums="1" hl_lines="45 49 57"
        rule fastq_qc_sol4:
            '''
            This rule performs a QC on paired-end .fastq files before and after trimming.
            '''
            input:
                reads1 = rules.fastq_trim.input.reads1,
                reads2 = rules.fastq_trim.input.reads2,
                trim1 = rules.fastq_trim.output.trim1,
                trim2 = rules.fastq_trim.output.trim2
            output:
                # QC before trimming
                html1_before = 'results/{sample}/fastqc_reports/{sample}_before_trim_1.html',  # Forward-read report in HTML format before trimming
                zipfile1_before = 'results/{sample}/fastqc_reports/{sample}_before_trim_1.zip',  # Forward-read report in ZIP format before trimming
                html2_before = 'results/{sample}/fastqc_reports/{sample}_before_trim_2.html',  # Reverse-read report in HTML format before trimming
                zipfile2_before = 'results/{sample}/fastqc_reports/{sample}_before_trim_2.zip',  # Reverse-read report in ZIP format before trimming
                # QC after trimming
                html1_after = 'results/{sample}/fastqc_reports/{sample}_after_trim_1.html',  # Forward-read report in HTML format after trimming
                zipfile1_after = 'results/{sample}/fastqc_reports/{sample}_after_trim_1.zip',  # Forward-read report in ZIP format after trimming
                html2_after = 'results/{sample}/fastqc_reports/{sample}_after_trim_2.html',  # Forward-read report in HTML format after trimming
                zipfile2_after = 'results/{sample}/fastqc_reports/{sample}_after_trim_2.zip'  # Forward-read report in ZIP format after trimming
            params:
                wd = 'results/{sample}/fastqc_reports/',  # Temporary directory to store files before renaming
                # QC before trimming
                html1_before = 'results/{sample}/fastqc_reports/{sample}_1_fastqc.html',  # Default FastQC output name for forward-read report in HTML format before trimming
                zipfile1_before = 'results/{sample}/fastqc_reports/{sample}_1_fastqc.zip',  # Default FastQC output name for forward-read report in ZIP format before trimming
                html2_before = 'results/{sample}/fastqc_reports/{sample}_2_fastqc.html',  # Default FastQC output name for reverse-read report in HTML format before trimming
                zipfile2_before = 'results/{sample}/fastqc_reports/{sample}_2_fastqc.zip',  # Default FastQC output name for reverse-read report in ZIP format before trimming
                # QC after trimming
                html1_after = 'results/{sample}/fastqc_reports/{sample}_atropos_trimmed_1_fastqc.html',# Default FastQC output name for forward-read report in HTML format after trimming
                zipfile1_after = 'results/{sample}/fastqc_reports/{sample}_atropos_trimmed_1_fastqc.zip',# Default FastQC output name for forward-read report in ZIP format after trimming
                html2_after = 'results/{sample}/fastqc_reports/{sample}_atropos_trimmed_2_fastqc.html',# Default FastQC output name for reverse-read report in HTML format after trimming
                zipfile2_after = 'results/{sample}/fastqc_reports/{sample}_atropos_trimmed_2_fastqc.zip'# Default FastQC output name for reverse-read report in ZIP format after trimming
            log:
                'logs/{sample}/{sample}_fastqc.log'
            threads: 2
            container:
                'https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0'
            shell:
                '''
                echo "Creating output directory <{params.wd}>" > {log}
                mkdir -p {params.wd} 2>> {log}  # FastQC doesn't create output directories so we have to do it manually
                echo "Performing QC of reads before trimming in <{input.reads1}> and <{input.reads2}>" >> {log}
                fastqc --format fastq --threads {threads} --outdir {params.wd} \
                --dir {params.wd} {input.reads1} {input.reads2} &>> {log}
                echo "Renaming results from original fastq analysis" >> {log}  # Renames files because we can't choose fastqc output
                mv {params.html1_before} {output.html1_before} 2>> {log}
                mv {params.zipfile1_before} {output.zipfile1_before} 2>> {log}
                mv {params.html2_before} {output.html2_before} 2>> {log}
                mv {params.zipfile2_before} {output.zipfile2_before} 2>> {log}
                echo "Performing QC of reads after trimming in <{input.trim1}> and <{input.trim2}>" >> {log}
                fastqc --format fastq --threads {threads} --outdir {params.wd} \
                --dir {params.wd} {input.trim1} {input.trim2} &>> {log}
                echo "Renaming results from trimmed fastq analysis" >> {log}  # Renames files because we can't choose fastqc output
                mv {params.html1_after} {output.html1_after} 2>> {log}
                mv {params.zipfile1_after} {output.zipfile1_after} 2>> {log}
                mv {params.html2_after} {output.html2_after} 2>> {log}
                mv {params.zipfile2_after} {output.zipfile2_after} 2>> {log}
                echo "Results saved in <results/{wildcards.sample}/fastqc_reports/>" >> {log}
                '''
        ```

    This solution is very long and much more complicated than the other one. However, it makes up for its complexity by allowing a total control on what is happening: with this method, we can choose where temporary files are saved as well as output names. It could have been shortened by using `-o .` to tell FastQC to create files in the current working directory instead of a specific one, but this would have created another problem: if we run multiple jobs in parallel, then Snakemake may try to produce files from different jobs at the same temporary destination. In this case, the different Snakemake instances would be trying to write to the same temporary files at the same time, overwriting each other and corrupting the output files.

Three interesting things are happening in both versions of this rule:

* Similarly to outputs, it is possible to refer to inputs of a rule directly in another rule with the syntax `rules.<rule_name>.input.<input_name>`
* FastQC doesn't create output directories by itself (other programs might insist that output directories **do not** already exist), so you have to manually create it with `mkdir` in the `shell` directive before running FastQC
    * Overall, most software will create the required directories for you; FastQC is an exception
* When using output **files**, Snakemake creates the missing folders by itself. However, when using a `directory(`) output, Snakemake **will not** create the directory. Otherwise, this would guarantee the rule to succeed every time, independently of the `shell` directive status: as soon as the directory is created, the rule succeeds, even if the `shell` directive fails. The current behaviour prevents Snakemake from continuing the workflow when a command may have failed. This also shows that solution 3, albeit simple, can be risky

??? tip "Controlling execution flow"
    If you want to make sure that a certain rule is executed before another, you can write the outputs of the first rule as inputs of the second one, even if the rule doesn't use them. For example, you could force the execution of FastQC before mapping reads with one more line in rule `read_mapping`:
    ```python linenums="1"
    rule read_mapping:
        '''
        This rule maps trimmed reads of a fastq onto a reference assembly.
        '''
        input:
            trim1 = rules.fastq_trim.output.trim1,
            trim2 = rules.fastq_trim.output.trim2,
            fastqc = rules.fastq_qc_sol4.output.html1_before  # This single line will force the execution of FastQC before read mapping
        output:
            sam = 'results/{sample}/{sample}_mapped_reads.sam',
            report = 'results/{sample}/{sample}_mapping_report.txt'
        params:
            index = 'resources/genome_indices/Scerevisiae_index'
        log:
            'logs/{sample}/{sample}_mapping.log'
        threads: 4
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
    ```
