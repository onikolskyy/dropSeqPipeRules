"""Get fastqc reports"""

#Which rules will be run on the host computer and not sent to nodes
localrules:
    multiqc_fastqc_reads,
    multiqc_fastqc_barcodes

rule fastqc_barcodes:
    """Create fastqc report"""
    input: 
        get_R1_files,
    output:
        html='{results_dir}/logs/fastqc/{sample}_R1_fastqc.html',
        zip='{results_dir}/logs/fastqc/{sample}_R1_fastqc.zip'
    params: '--extract'
    benchmark: '{results_dir}/benchmarks/fastqc_barcodes.benchmark.{sample}.txt'
    wrapper:
        '0.27.1/bio/fastqc'

rule fastqc_reads:
    """Create fastqc report"""
    input: 
        get_R2_files,
    output:
        html='{results_dir}/logs/fastqc/{sample}_R2_fastqc.html',
        zip='{results_dir}/logs/fastqc/{sample}_R2_fastqc.zip'
    benchmark: '{results_dir}/benchmarks/fastqc_reads.benchmark.{sample}.txt'
    params: '--extract'
    wrapper:
        '0.27.1/bio/fastqc'


rule multiqc_fastqc_barcodes:
    input:
        expand('{results_dir}/logs/fastqc/{sample}_R1_fastqc.html', sample=samples.index, results_dir=results_dir)
    output:
        html='{results_dir}/reports/fastqc_barcodes.html'
    params: '-m fastqc --ignore *_R2*'
    benchmark: '{results_dir}/benchmarks/multiqc_fastqc_barcodes.txt'
    wrapper:
        '0.27.1/bio/multiqc'

rule multiqc_fastqc_reads:
    input: 
        expand('{results_dir}/logs/fastqc/{sample}_R2_fastqc.html', sample=samples.index, results_dir=results_dir)
    output:
        html='{results_dir}/reports/fastqc_reads.html'
    benchmark: '{results_dir}/benchmarks/multiqc_fastqc_reads.txt'
    params: '-m fastqc --ignore *_R1*'
    wrapper:
        '0.27.1/bio/multiqc'
