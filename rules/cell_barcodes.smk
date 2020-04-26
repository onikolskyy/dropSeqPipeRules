whitelist_opts = pd.read_table("whitelist_opts.csv", sep=",").set_index("index", drop=False)

localrules:
    get_cell_whitelist,
#    extend_barcode_whitelist,
#    extend_barcode_top

rule merge_and_repair:
    input:
        # define input paths
        R1='{results_dir}/samples/{sample}/trimmmed_repaired_R1.fastq.gz',
        R2='{results_dir}/samples/{sample}/Aligned.out.bam'
    params:
        cell_barcode_length=(config['FILTER']['cell-barcode']['end'] - config['FILTER']['cell-barcode']['start'] + 1),
        umi_barcode_length=(config['FILTER']['UMI-barcode']['end'] - config['FILTER']['UMI-barcode']['start'] + 1),
        num_cells=lambda wildcards: round(int(samples.loc[wildcards.sample,'expected_cells'])*1.2),
    output:
        repaired_bam='{results_dir}/samples/{sample}/Aligned.repaired.bam',
        mapping='{results_dir}/samples/{sample}/top_barcodes.csv'
    conda: 'envs/merge_and_repair.yaml'
    benchmark : '{results_dir}/sample/MERGE_AND_REPAIR.{sample}.txt'
    script: 'scripts/merge_and_repair.py'

# for plotting
rule get_cell_whitelist:
    input:
        '{results_dir}/samples/{sample}/top_barcodes.csv'
    output:
        '{results_dir}/samples/{sample}/barcodes.csv'
    benchmark: '{results_dir}/benchmarks/get_cell_whitelist.{sample}.txt'
    shell:
        """cat {input} | cut -f 1 > {output}"""


