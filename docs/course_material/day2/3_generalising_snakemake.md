## Learning outcomes

**After having completed this chapter you will be able to:**

* Create rules with multiple inputs and outputs
* Make the code shorter and more general by using placeholders and wildcards
* Check the workflow behaviour
* Visualise a workflow DAG

## Material

[:fontawesome-solid-file-pdf: Download the presentation](../../assets/pdf/day2/3_generalising_snakemake.pdf){: .md-button }

## Advice and reminders

In each rule, you should try (as much as possible) to:

* Choose meaningful rule names
* Use rules dependency, with the syntax `rules.<rule_name>.output`
    * If you use named outputs (recommended), the syntax becomes `rules.<rule_name>.output.<output_name>`
    * If you use numbered outputs, the syntax becomes `rules.<rule_name>.output[n]` (with `n` starting at 0)
* Use multiple inputs/outputs when needed/possible
* Use placeholders
* Use wildcards
    * Choose meaningful wildcard names
    * You can use the same wildcard names in multiple rules for consistency and readability, but Snakemake will treat them as independent wildcards and their values will not be shared: rules are self-contained and wildcards are local to each rule (see this [very nice summary](http://ivory.idyll.org/blog/2023-snakemake-slithering-wildcards.html) on wildcards)
* Create a log file with the `log` directive and a benchmark file with the `benchmark` directive
    * The `output`, `log`, and `benchmark` directives must contain the same wildcards!

## Testing the workflow's logic

* If you have a doubt, do not hesitate to test your workflow logic with a **dry-run** (the `-n/--dry-run/--dryrun` parameter): `snakemake -c 1 -n <target>`. Snakemake will then display all the jobs required to generate the target as well as a **reason** field explaining why each job is required
* To visualise the exact commands executed by each job (with the placeholders and wildcards replaced by their values), run snakemake with the `-p/--printshellcmds` parameter: `snakemake -c 1 -p <target>`
* These two parameters are often used together to check an entire workflow: `snakemake -c 1 -n -p <target>`

## Data origin

The data we will use during the exercises was produced [in this work](https://pubmed.ncbi.nlm.nih.gov/31654410/). Briefly, the team studied the transcriptional response of a strain of baker's yeast, [_Saccharomyces cerevisiae_](https://en.wikipedia.org/wiki/Saccharomyces_cerevisiae), facing environments with different concentrations of CO<sub>2</sub>. To this end, they performed **150 bp paired-end** sequencing of mRNA-enriched samples. Detailed information on all the samples are available [here](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA550078), but just know that for the purpose of the course, we selected **6 samples** (**3 replicates per condition**, **low** and **high CO<sub>2</sub>**) and **down-sampled** them to **1 million read-pairs each** to reduce computation time.

## Exercises

One of the aims of today's course is to develop a simple, yet efficient, workflow to analyse bulk RNAseq data. This workflow takes reads coming from RNA sequencing as inputs and produces a list of genes that are differentially expressed between the two conditions. The files containing the reads are in [FASTQ format](https://en.wikipedia.org/wiki/FASTQ_format) and the output will be a tab-separated file containing a list of genes with expression changes, results of statistical tests...

In this series of exercises, we will create the 'backbone' of the workflow, _i.e._ the rules that are the most computationally expensive, namely:

* A rule to trim poor-quality reads
* A rule to map the trimmed reads on a reference genome
* A rule to convert and sort files from the SAM format to the BAM format
* A rule to count the reads mapping on each gene

??? tip "Designing and debugging a workflow"
    If you have problems designing your Snakemake workflow or debugging it, you can find some help [here](6_debugging_snakemake.md#designing-a-workflow).

At the end of this series of exercises, your workflow DAG should look like this:
<figure markdown align="center">
  ![backbone_rulegraph](../../assets/images/backbone_rulegraph.png)
  <figcaption>Workflow rulegraph at <br>the end of the session</figcaption>
</figure>

### Downloading the data and setting up the folder structure

In this part, we will download the data and start building the directory structure of our workflow according to the [official recommendations](https://snakemake.readthedocs.io/en/v8.20.5/snakefiles/deployment.html). We already starting doing so in the previous series of exercises and at the end of the course, it should resemble this:

```
│── .gitignore
│── README.md
│── LICENSE.md
│── benchmarks
│   │── sample1.txt
│   └── sample2.txt
│── config
│   │── config.yaml
│   └── some-sheet.tsv
│── data
│   │── sample1.fastq
│   └── sample2.fastq
│── images
│   │── dag.png
│   │── filegraph.png
│   └── rulegraph.png
│── logs
│   │── sample1.log
│   └── sample2.log
│── results
│   │── DEG_list.tsv
│   │── sample1
│   │   └── sample1.bam
│   └── sample2
│       └── sample2.bam
│── resources
│   │── Scerevisiae.fasta
│   │── Scerevisiae.gtf
│   └── genome_indices
|       │── Scerevisiae_index.1.ht2
|       │── Scerevisiae_index.2.ht2
|       │── Scerevisiae_index.3.ht2
|       │── Scerevisiae_index.4.ht2
|       │── Scerevisiae_index.5.ht2
|       │── Scerevisiae_index.6.ht2
|       │── Scerevisiae_index.7.ht2
|       └── Scerevisiae_index.8.ht2
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

For now, the main thing to remember is that the **code** goes into the **`workflow` subfolder**  and the **rest** is mostly **input/output files**, except for the **`config` subfolder**, which will be **explained later**. All **output files** generated in the workflow should be stored under **`results/`**.

Now, let's download the data, uncompress it and build the first part of the directory structure.

```sh linenums="1"
ssh -i ~/.ssh/key_username.pem username@18.195.170.182  # Connect to the server; don't forget to change key_username and username
wget https://containers-snakemake-training.s3.eu-central-1.amazonaws.com/snakemake_rnaseq.tar.gz  # Download data
tar -xvf snakemake_rnaseq.tar.gz  # Uncompress archive
rm snakemake_rnaseq.tar.gz  # Delete archive
cd snakemake_rnaseq/  # Start developing in a new folder
```

In the ` snakemake_rnaseq/` folder, you should now see 2 subfolders:

* `data/`, which contains the data to analyse
* `resources/`, which contains retrieved resources, here the assembly, the genome indices and the annotation file of _S. cerevisiae_. It may also contain small resources delivered along with the workflow via git

Let's create the other missing subfolders and the Snakefile:

```sh linenums="1"
mkdir -p config/ images/ workflow/envs workflow/rules workflow/scripts  # Create subfolder structure
touch workflow/Snakefile  # Create empty Snakefile
```

??? info "What does `-p` do?"
    The `-p` parameter of `mkdir` make parent directories as needed and does not return an error if the directory already exists.

The **Snakefile** marks the **entrypoint** of the workflow. It will be automatically discovered when running Snakemake from the root of the structure, here `snakemake_rnaseq/`.

??? info "Using a Snakefile from a non-default location"
    Snakemake can use a Snakefile from a non-default location thanks the `-s/--snakefile` parameter: `snakemake -c 1 -s <Snakefile_path> <target>`. However, it is **highly discouraged** as it **hampers reproducibility**.

??? warning "Relative paths in Snakemake"
    All the paths in a Snakefile are relative to the working directory in which the `snakemake` command is executed:

    * If you execute Snakemake in `snakemake_rnaseq/`, the relative path to the input files in the rule is `data/<sample>.fastq`
    * If you execute Snakemake in `snakemake_rnaseq/workflow/`, the relative path to the input files in the rule is `../data/<sample>.fastq`

If you followed the [advice](#advice-and-reminders) at the top of this page, Snakemake should create all the other missing folders by itself, so it is now time to create the rules mentioned [earlier](#exercises). If needed, you can check [here](6_debugging_snakemake.md#designing-a-workflow) for a few pieces of advice on workflow design.

??? info "'bottom-up' or 'top-down' development?"
    Even if it is often easier to start from the final outputs and work backwards to the first inputs, the next exercises are presented in the opposite direction (first inputs to last outputs) to make the session easier to understand. That being said, feel free to work and develop your code in the order you prefer!

#### Important: do not process all the samples!

**Do not try to process all the samples yet, even if we asked you to use wildards. For now, choose only one sample (which means two .fastq files because reads are paired-end). We will see an efficient way to process a list of files in the next series of exercises.**

### Creating a rule to trim reads

Usually, the first step when dealing with sequencing data is to improve read quality by removing low quality bases, stretches of As and Ns and reads that are too short.

??? info "Trimming sequencing adapters"
    In theory, trimming should also remove sequencing adapters, but we will not do it here to keep computation time low and avoid parsing other files to extract the adapter sequences.

We will use [atropos](https://peerj.com/articles/3720/) to trim the reads. The first part of the trimming command is:
```
atropos trim -q 20,20 --minimum-length 25 --trim-n --preserve-order --max-n 10 --no-cache-adapters -a "A{{20}}" -A "A{{20}}"
```

??? info "Explanation of atropos parameters"
    * `-q 20,20`: trim low-quality bases from 5' and 3' ends of each read before adapter removal
    * `--minimum-length 25`: discard trimmed reads that are shorter than 25 bp
    * `--trim-n`: trim Ns at the ends of reads
    * `--preserve-order`: preserve order of reads in input files
    * `--max-n 10`: discard reads with more than 10 Ns
    * `--no-cache-adapters`: do not cache adapters list as '.adapters' in the working directory
    * `-a "A{{20}}" -A "A{{20}}"`: remove series of 20 As in the adapter sequence (`-a` for the first read of the pair, `-A` for the second one)
        * The usual command-line syntax is `-a "A{20}"`. Here, brackets were doubled to prevent Snakemake from interpreting `{20}` as a wildcard

**Exercise:** Complete the atropos command given above with parameters to specify inputs (files to trim) and outputs (trimmed files). You can find information on how to use `atropos` and its parameters with `atropos trim -h`. Then, implement a rule containing your command to trim the reads contained in the .fastq files.

??? tip "atropos inputs and outputs"
    * The fastq files to trim are located in `data/`
    * The paths of the files to trim (_i.e._ input files, in FASTQ format) are specified with the parameters `-pe1` (first read) and `-pe2` (second read)
    * The paths of the trimmed files (_i.e._ output files, also in FASTQ format) are specified with the parameters `-o` (first read) and `-p` (second read)

??? success "Answer"
    This is one way of writing this rule, but definitely not the only way (and this is true for all the rules presented in these exercises).

    ```python linenums="1"
    rule fastq_trim:
        """
        This rule trims paired-end reads to improve their quality. Specifically, it removes:
        - Low quality bases
        - A stretches longer than 20 bases
        - N stretches
        """
        input:
            reads1 = 'data/{sample}_1.fastq',
            reads2 = 'data/{sample}_2.fastq',
        output:
            trim1 = 'results/{sample}/{sample}_atropos_trimmed_1.fastq',
            trim2 = 'results/{sample}/{sample}_atropos_trimmed_2.fastq'
        shell:
            '''
            atropos trim -q 20,20 --minimum-length 25 --trim-n --preserve-order --max-n 10 \
            --no-cache-adapters -a "A{{20}}" -A "A{{20}}" \
            -pe1 {input.reads1} -pe2 {input.reads2} -o {output.trim1} -p {output.trim2}
            '''
    ```

    There are three interesting things happening here:

    1. We added a comment between double quotes under the rule name. In Python, it is called a docstring. While it is not mandatory, it is good practice (and very recommended) to write docstrings to explain what a rule does, what are its inputs, outputs, parameters...
    1. We used the `{sample}` wildcards twice in the output paths (L12-13). This is because we prefer to have all the files from a single sample in the same directory
    1. We used a backslash `\` at the end of L16-17 to split a very long command in smaller lines. This is purely 'cosmetic' but it avoids very long lines that are painful to read, understand and debug...

**Exercise:** If you had to run the workflow by specifying only one output, what command would you use?

??? success "Answer"
    To process the sample `highCO2_sample1`, for example, you would use the command:

    ```snakemake -c 1 -p results/highCO2_sample1/highCO2_sample1_atropos_trimmed_1.fastq```

    You don't need to ask for the two outputs of the rule: asking only for one will still trigger the execution, but the workflow will complete without errors if and only if **both** outputs are present. Like with [intermediary files](2_introduction_snakemake.md#chaining-rules), this property also helps reducing the number of targets to write in the snakemake command used to execute the workflow!

### Creating a rule to map trimmed reads onto a reference genome

Once the reads are trimmed, the next step is to map those reads onto the species genome (_S. cerevisiae_ strain S288C) to eventually obtain read counts. The reference assembly used in this exercise is [RefSeq GCF_000146045.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000146045.2/) and was retrieved via the NCBI website. We will use [HISAT2](https://www.nature.com/articles/s41587-019-0201-4) to map the reads.

??? info "HISAT2 genome index"
    To align reads to a genome, HISAT2 relies on a graph-based index. To save some time, we built the genome index for you, using the command:
    ```
    hisat2-build -p 24 -f Scerevisiae.fasta resources/genome_indices/Scerevisiae_index
    ```
    The parameters are explained here:

    * `-p`: number of threads to use
    * `-f`: genomic sequence in [FASTA format](https://en.wikipedia.org/wiki/FASTA_format)
    * `Scerevisiae_index`: global name shared by all the index files

    The index files can be found in  `resources/genome_indices/`.

The first part of the mapping command is:
```
hisat2 --dta --fr --no-mixed --no-discordant --time --new-summary --no-unal
```

??? info "Explanation of HISAT2 parameters"
    * `--dta`: report alignments tailored for transcript assemblers
    * `--fr`: set alignment of -1, -2 mates to forward/reverse (position of reads in a pair relatively to each other)
    * `--no-mixed`: remove unpaired alignments for paired reads
    * `--no-discordant`: remove discordant alignments for paired reads
    * `--time`: print wall-clock time taken by search phases
    * `--new-summary`: print alignment summary in a new style
    * `--no-unal`: suppress SAM records for reads that failed to align

**Exercise:** Complete the HISAT2 command given above with parameters to specify inputs and outputs. Briefly, you will need 2 inputs, 2 outputs (in 2 different formats) and the genome indices mentioned above (should they be considered as inputs?). You can find more information on how to use HISAT2 and its parameters with `hisat2 -h`. Then, implement a rule containing your command to map the trimmed reads contained in the .fastq files.

??? tip "HISAT2 inputs and outputs"
    * The paths of the trimmed files (_i.e._ input files) are specified with the parameters `-1` (first read) and `-2` (second read)
    * The basename of the genome indices (binary format) is specified with the parameter `-x`. The files have a shared name of `resources/genome_indices/Scerevisiae_index`, which is the value you need to use for `-x`
    * The path of the mapped reads file (_i.e._ output file, in SAM format) is specified with the parameter `-S` (do not forget the .sam extension at the end of the file)
    * The path of the mapping report (_i.e._ output file, in text format) is specified with the parameter `--summary-file`

??? success "Answer"
    ```python linenums="1"
    rule read_mapping:
        """
        This rule maps trimmed reads of a fastq on a reference assembly.
        """
        input:
            trim1 = rules.fastq_trim.output.trim1,
            trim2 = rules.fastq_trim.output.trim2
        output:
            sam = 'results/{sample}/{sample}_mapped_reads.sam',
            report = 'results/{sample}/{sample}_mapping_report.txt'
        shell:
            '''
            hisat2 --dta --fr --no-mixed --no-discordant --time --new-summary --no-unal \
            -x resources/genome_indices/Scerevisiae_index \
            -1 {input.trim1} -2 {input.trim2} -S {output.sam} --summary-file {output.report}
            '''
    ```

**Exercise:** What do you think about the value of the `-x` parameter?

??? success "Answer"
    Something very interesting is happening here: to run, HISAT2 requires several genome index files. As such, they could (should?) be considered as inputs... and this is a problem:

    * On one hand, the `-x` parameter only accepts strings containing the files path and common name. This means that if you were to manually add all the index files as inputs, HISAT2 would not recognize them... and crash!
    * On the other hand, if you were add the value of `-x` as input, Snakemake would look for a file called `resources/genome_indices/Scerevisiae_index`... and crash because such a file doesn't exist!

    This highlights the fact that **inputs must be files** and this is why we directly added the value of `-x` in the `shell` command. However, this is not very convenient: we will see later a better way to deal with this problem.

Using the same sample than before (`highCO2_sample1`), the workflow can be run with:

```snakemake -c 1 -p results/highCO2_sample1/highCO2_sample1_mapped_reads.sam```

That being said, we recommend **not to run it now**, because this step is the longest of the workflow (with the current settings, it will take ~6 min to complete). If you want to run it now, you can, but you should launch it and start working on the next rules while it finishes.

### Creating a rule to convert and sort SAM files to BAM

HISAT2 only outputs mapped reads in the [SAM format](https://en.wikipedia.org/wiki/SAM_(file_format)). However, most downstream analysis tools use the [BAM format](https://en.wikipedia.org/wiki/Binary_Alignment_Map), which is the compressed binary version of the SAM format and, as such, is much smaller, easier to manipulate and transfer and allows a faster data retrieval. Additionally, many analyses require BAM files to be sorted by genomic coordinates and indexed because sorted BAM files can be processed much more easily and quickly than unsorted ones. Operations on SAM and BAM files are usually performed with [Samtools](https://doi.org/10.1093/gigascience/giab008).

??? info "What are alignment formats?"
    More information on alignment formats can be found on [samtools github repository](https://github.com/samtools/hts-specs).

We wrote a **single rule** performing all these operations for you:

```python linenums="1"
rule sam_to_bam:
    """
    This rule converts a sam file to bam format, sorts it and indexes it.
    """
    input:
        sam = rules.read_mapping.output.sam
    output:
        bam = 'results/{sample}/{sample}_mapped_reads.bam',
        bam_sorted = 'results/{sample}/{sample}_mapped_reads_sorted.bam',
        index = 'results/{sample}/{sample}_mapped_reads_sorted.bam.bai'
    shell:
        '''
        samtools view {input.sam} -b -o {output.bam}
        samtools sort {output.bam} -O bam -o {output.bam_sorted}
        samtools index -b {output.bam_sorted} -o {output.index}
        '''
```

??? info "Explanation of Samtools parameters"
    * You can find information on how to use Samtools and its parameters with `samtools --help`
    * `samtools view`:
        * `-b`: create an output in BAM format
        * `-o`: path of the output file
    * `samtools sort`:
        * `-O bam`: create an output in BAM format
        * `-o`: path of the output file
    * `samtools index`:
        * `-b`: create an index in BAI format

**Exercise:** Copy this rule in your Snakefile. Don't forget to update the output values to match the ones you used in your previous rules. Do you notice anything surprising in this rule?

??? success "Answer"
    Let's start with a quick breakdown of what the `shell` directive does:

    * L13: `samtools view` converts a file in SAM format to BAM format
    * L14: `samtools sort` sorts a BAM file by genomic coordinates
    * L15: `samtools index` indexes a sorted BAM file. The index must have the exact same basename than its associated BAM file; the only difference is that it finishes with the extension `.bam.bai` instead of `.bam`

    The interesting this is that so far, all the rules had only one command in the `shell` directive. In this rule, there are 3 commands grouped together, each with their own inputs and outputs. This means 2 things:

    1. We could have split this rule into 3 separate rules with dependencies. There is no official guideline on whether to split rules like this, but a good rule of thumb is: does it make sense to run these commands together? Is it not too computationally expensive (time, memory, cpu) to run these rules together?
    1. `{output.bam}` and `{output.bam_sorted}` are outputs of  `samtools view` and  `samtools sort`... But they are also inputs of  `samtools sort` and  `samtools index`! This means that files that are created by a command can instantly be re-used in subsequent commands **within the same rule**!

**Exercise:** Do you see any drawbacks to using a rule like this?

??? success "Answer"
    When Snakemake has a problem and crashes, it removes all the rule outputs to avoid further computation with corrupted files. This means that if the first two commands complete but the last command fails (here, `samtools index`) and doesn't produce the expected output, Snakemake will remove **all the outputs**, including those that were created without error (here, `{output.bam}` and `{output.bam_sorted}`), causing a waste of time and money


Using the same sample than before (`highCO2_sample1`), the workflow can be run with
```
snakemake -c 1 -p results/highCO2_sample1/highCO2_sample1_mapped_reads_sorted.bam
```
However, you will soon run the entire workflow, so it might be worth waiting!

### Creating a rule to count mapped reads

Most of the analyses happening downstream the alignment step, including Differential Expression Analyses, are starting off read counts, by exon or by gene. However, we are still missing those counts!

??? info "Genome annotations"

    * To count reads mapping on genomic features, we first need a definition of those features, called annotations. Here, we chose one of the best-known model organism, _S. cerevisiae_, which has been annotated for a long time. Its annotations are easily findable on the [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000146045.2/) or the [Saccharomyces Genome Database](https://www.yeastgenome.org/). If an organism has not been annotated yet, there are ways to work around this problem, but this is an entirely different field that we won't discuss here!
    * You should also know that there are two main annotations format: [GTF](http://mblab.wustl.edu/GTF2.html) and [GFF](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md). The former is lighter and easier to work with, so that is the one we will use

??? warning "Chromosome names must match between files"
    If you are working with genome sequences and annotations from different sources, remember that they **must** contain matching chromosome names, otherwise counting will not work, as the counting software will not be able to read counts to match exon/gene locations.

We already wrote a rule to count the reads mapping on each gene of the _S. cerevisiae_ genome using [featureCounts](https://academic.oup.com/bioinformatics/article/30/7/923/232889):

```python linenums="1"
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
        'logs/{sample}/{sample}_genes_read_quantification.log'
    benchmark:
        'benchmarks/{sample}/{sample}_genes_read_quantification.txt'
    shell:
        '''
        featureCounts -t exon -g gene_id -s 2 -p --countReadPairs \
        -B -C --largestOverlap --verbose -F GTF \
        -a resources/Scerevisiae.gtf -o {output.gene_level} {input.bam_once_sorted}
        mv {output.gene_level}.summary {output.gene_summary}
        '''
```

??? info "Explanation of featureCounts command"
    * You can find information on how to use featureCounts and its parameters with `featureCounts -h`
    * `-t`: specify on which feature type to count the reads
    * `-g`: specify if and how to gather feature counts. Here, reads are counted by features (exon) (`-t`) and the exon counts are gathered by 'meta-features’ (genes) (`-g`)
    * `-s`: perform strand-specific read counting
        * Strandedness is determined by looking at the mRNA library preparation kit. It can also be determined _a posteriori_ with scripts such as [infer_experiment.py](https://rseqc.sourceforge.net/#infer-experiment-py) from the [RSeQC package](http://doi.org/10.1093/bioinformatics/bts356)
    * `-p`: specify that input data contain paired-end reads
    * `--countReadPairs`: count read pairs instead of reads
    * `-B`: only count read pairs that have both ends aligned
    * `-C`: do not count read pairs that have their two ends mapping to different chromosomes or mapping on the same chromosome but on different strands
    * `--largestOverlap`: assign reads to the meta-feature/feature that has the largest number of overlapping bases
    * `--verbose`: output verbose information, such as unmatched chromosome/contig names
    * `-F`: specify format of the provided annotation file
    * `-a`: specify path of file containing the annotations (_i.e._ input files, in GTF format)
    * `-o`: specify path of file containing the count results (_i.e._ output file, in tsv format)
    * Paths of the sorted BAM file(s) (_i.e._ input file(s)) are not specified with an parameter, they are simply added at the end of the command

**Exercise:** What does L19 do? Why did we add it?

??? success "Answer"
    The `mv` command can be used to move or rename a file. Here, it does the latter. featureCounts outputs a second, separate file (in tsv format) containing summary statistics about the read counting, with the name `<output_name>.summary`. For example, if the output is `test.tsv`, the summary will be printed in `test.tsv.summary`. However, there is no parameter available to choose the filename, so if we need this file as an output, we have to manually rename it.

It would be interesting to know what is happening when featureCounts runs. This is where the `log` and `benchmark` directives play a role!

**Exercise:** Copy this rule in your Snakefile. Then, add the `log` and `benchmark` directives to the rule. Don't forget to update the directive values to match the ones you used in your previous rules.

??? tip "Logs and benchmarks"
    * The `log` and `benchmark` directives must contain the same wildcards than the `output` directive, here `sample`
    * Logs need to be handled manually, so you need to redirect what is produced by featureCounts to the log file. You can redirect both `stdout` and `stderr` streams with `&>` at the end of the command

??? success "Answer"
    ```python linenums="1"
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
            'logs/{sample}/{sample}_genes_read_quantification.log'
        benchmark:
            'benchmarks/{sample}/{sample}_genes_read_quantification.txt'
        shell:
            '''
            echo "Counting reads mapping on genes in <{input.bam_once_sorted}>" > {log}
            featureCounts -t exon -g gene_id -s 2 -p --countReadPairs \
            -B -C --largestOverlap --verbose -F GTF \
            -a resources/Scerevisiae.gtf -o {output.gene_level} {input.bam_once_sorted} &>> {log}
            echo "Renaming output files" >> {log}
            mv {output.gene_level}.summary {output.gene_summary}
            echo "Results saved in <{output.gene_level}>" >> {log}
            echo "Report saved in <{output.gene_summary}>" >> {log}
            '''
    ```

### Running the whole workflow

**Exercise:** If you have not done it after each step, it is now time to run the entire workflow on your sample of choice. What command will you use to run it? Once Snakemake has finished, check the log of the rule `reads_quantification_genes`. Does it contain useful information?

??? success "Answer"
    Because all the rules are chained together, you only need to specify one of the final outputs to trigger the execution of all the previous rules. Using the same sample than before (`highCO2_sample1`). You can add the `-F` parameter to force an entire re-run, which should take ~10 min.:
    ```
    snakemake -c 1 -F -p results/highCO2_sample1/highCO2_sample1_genes_read_quantification.tsv
    ```
    You can check the logs with the command:
    ```
    cat logs/highCO2_sample1/highCO2_sample1_genes_read_quantification.log
    ```

**Extra:** If you have time, check Snakemake's log in `.snakemake/log/`. Is everything as you expected, especially the wildcard values, input and output names etc...?

??? success "Answer"
    You can check the logs with the command:
    ```
    cat .snakemake/log/<latest_log>
    ```
    This general log says exactly what was run by Snakemake, after the placeholders, wildcards... were replaced by their actual values.

### Visualising the DAG of the workflow

We have now implemented and run the main steps of our workflow. It is always a good idea to visualise the whole process to check for errors and inconsistencies. Interestingly, Snakemake's has a built-in workflow visualisation feature to do this.



**Exercise:** Visualise the entire workflow’s Directed Acyclic Graph using the `--dag` parameter. Remember that Snakemake prints a DAG in text format, so you need to pipe its results (using `|`) into the `dot` command to transform it into a picture. Do you need to specify a target?

??? tip "Creating a DAG"
    * Try to follow the official recommendations on workflow structure, which states that images are supposed to go in the `images/` subfolder
    * The `images/` subfolder is not automatically created by Snakemake because it isn't handled as part of an actual run, so you need to create it beforehand. You did this when setting up the [workflow structure](#downloading-the-data-and-setting-up-the-folder-structure)
    * If you already computed all the outputs of the workflow, steps in the DAG will have dotted lines. To visualise the DAG before running the workflow, add `-F/--forceall` to the snakemake command to force the execution of all jobs.

??? tip "The `dot` command"
    * [`dot`](https://graphviz.org/docs/layouts/dot/) is a part of the [graphviz package](https://graphviz.org/) and is used to draw hierarchical or layered drawings of directed graphs, _i.e._ graphs in which edges (arrows) have a direction
    * `-T`: control the image format. Available formats are listed [here](https://graphviz.org/docs/outputs/)


??? success "Answer"
    To run the command without target, you can use:
    ```
    snakemake -c 1 --dag -F | dot -Tpng > images/dag.png
    ```
    You will get a `Target rules may not contain wildcards.` error. This makes sense, because if you don't give a target to Snakemake, it can't compute the wildcard values. To run the command using the same sample than before (`highCO2_sample1`), you can target one of the final outputs of the workflow:
    ```
    snakemake -c 1 --dag -F results/highCO2_sample1/highCO2_sample1_genes_read_quantification.tsv | dot -Tpng > images/dag.png
    ```
    We added the `-F` parameter to force Snakemake to compute the entire DAG and ensure all the jobs are shown. You can also use `-f <target>` to show fewer jobs.

    An important thing to remember is that the `--dag` parameter implicitly activates the `--dry-run/--dryrun/-n` parameter, which means that no jobs are executed during the plot creation.

Snakemake can also create other graphs. In total, there are 3 types:

* A DAG, created with the `--dag` parameter: dependency graph of all the jobs (rules appear once per wildcard value; wildcard value are displayed)
* A rulegraph, created with the `--rulegraph` parameter: dependency graph of all the rules (rules appear only once)
* A filegraph, created with the `--filegraph` parameter: dependency graph of all the rules with inputs and outputs (rule appears once, wildcards are shown but not replaced)

Here are the graphs of the workflow:

<figure>
  <img src="../../../assets/images/backbone_dag.png" width="30%" height="450" align="center"/>
  <img src="../../../assets/images/backbone_rulegraph.png" width="30%" height="450" align="center"/>
  <img src="../../../assets/images/backbone_filegraph.png" width="30%" height="450" align="center"/>
  <figcaption>DAG, rulegraph and filegraph (respectively) of the workflow<br>at the end of the session</figcaption>
</figure>

And the commands to create them:

* Rulegraph:
```
snakemake -c 1 --rulegraph -F results/highCO2_sample1/highCO2_sample1_genes_read_quantification.tsv | dot -Tpdf > images/rulegraph.pdf
```
* Filegraph:
```
snakemake -c 1 --filegraph -F results/highCO2_sample1/highCO2_sample1_genes_read_quantification.tsv | dot -Tjpg > images/filegraph.jpg
```
