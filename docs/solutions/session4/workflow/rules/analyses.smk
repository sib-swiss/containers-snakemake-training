'''
This Snakefile contains rules used to gather files containing read counts into
a single table and perform Differential Expression Analyses.
'''


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
        table = 'results/total_count_table.tsv'
    log:
        'logs/total_count_table.log'
    threads: 1
    conda:
        '../envs/py.yaml'
    script:
        '../scripts/count_table.py'


rule differential_expression:
    '''
    This rule detects DEGs and plots control graphs (PCA, heatmaps...).
    '''
    input:
        table = rules.count_table.output.table
    output:
        deg = 'results/deg_list.tsv',
        pdf = 'results/deg_plots.pdf'
    log:
        'logs/differential_expression.log'
    threads: 2
    container:
        'docker://athiebaut/deseq2:v3'
    script:
        '../scripts/DESeq2.R'
