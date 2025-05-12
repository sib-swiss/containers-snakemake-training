'''
This Snakefile contains rules to trim reads, map them to a reference genome,
convert and sort the resulting bam files and count the mapped reads.
'''

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

# rule fastq_qc_sol3:
#     '''
#     This rule performs a QC on paired-end .fastq files before and after trimming.
#     '''
#     input:
#         reads1 = rules.fastq_trim.input.reads1,
#         reads2 = rules.fastq_trim.input.reads2,
#         trim1 = rules.fastq_trim.output.trim1,
#         trim2 = rules.fastq_trim.output.trim2
#     output:
#         before_trim = directory('results/{sample}/fastqc_reports/before_trim/'),
#         after_trim = directory('results/{sample}/fastqc_reports/after_trim/')
#     log:
#         'logs/{sample}/{sample}_fastqc.log'
#     threads: 2
#     container:
#         'https://depot.galaxyproject.org/singularity/fastqc%3A0.12.1--hdfd78af_0'
#     shell:
#         '''
#         echo "Creating output directory <{output.before_trim}>" > {log}
#         mkdir -p {output.before_trim} 2>> {log}
#         echo "Performing QC of reads before trimming in <{input.reads1}> and <{input.reads2}>" >> {log}
#         fastqc --format fastq --threads {threads} --outdir {output.before_trim} --dir {output.before_trim} {input.reads1} {input.reads2} &>> {log}
#         echo "Results saved in <{output.before_trim}>" >> {log}
#         echo "Creating output directory <{output.after_trim}>" >> {log}
#         mkdir -p {output.after_trim} 2>> {log}
#         echo "Performing QC of reads after trimming in <{input.trim1}> and <{input.trim2}>" >> {log}
#         fastqc --format fastq --threads {threads} --outdir {output.after_trim} --dir {output.after_trim} {input.trim1} {input.trim2} &>> {log}
#         echo "Results saved in <{output.after_trim}>" >> {log}
#         '''

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


rule read_mapping:
    '''
    This rule maps trimmed reads of a fastq onto a reference assembly.
    '''
    input:
        trim1 = rules.fastq_trim.output.trim1,
        trim2 = rules.fastq_trim.output.trim2,
        fastqc = rules.fastq_qc_sol4.output.html1_before
    output:
        sam = 'results/{sample}/{sample}_mapped_reads.sam',
        report = 'results/{sample}/{sample}_mapping_report.txt'
    params:
        index = config['index']
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
    a count table for all genes of the assembly
    '''
    input:
        bam_once_sorted = rules.sam_to_bam.output.bam_sorted,
    output:
        gene_level = 'results/{sample}/{sample}_genes_read_quantification.tsv',
        gene_summary = 'results/{sample}/{sample}_genes_read_quantification.summary'
    params:
        annotations = config['annotations']
    log:
        'logs/{sample}/{sample}_genes_read_quantification.log'
    threads: 2
    container:
        'https://depot.galaxyproject.org/singularity/subread%3A2.0.6--he4a0461_2'
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
