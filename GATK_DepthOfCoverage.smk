'''
GATK DepthOfCoverage
Not available in GATK4 so GATK3 needs to be used to run this tool
-geneList option produces *.sample_gene_summary file
'''
import glob

(DIRECTORIES, ANALYSES, FILES, SAMPLES) = glob_wildcards("{directory}/GATK{analysis}SeqPipeline_{file}/bqsr_readgroups_markdups_sorted_{sample}.bam")

rule all:
  input:
    expand("{directory}/GATK{analysis}SeqPipeline_{file}/GATK_DepthOfCoverage/{sample}", zip, directory=DIRECTORIES, analysis=ANALYSES, file=FILES, sample=SAMPLES)

rule GATK_DepthOfCoverage:
  input:
    fa = "/home/NBS/WES_Reference_Data/human_g1k_v37.fasta",
    bam="{directory}/GATK{analysis}SeqPipeline_{file}/bqsr_readgroups_markdups_sorted_{sample}.bam",
    interval_list = "/home/NBS/WES_Reference_Data/Illumine_Exome_CEX_TargetedRegions_v1.2.interval_list",
    gene_list = "/home/NBS/WES_Reference_Data/geneTrack_DepthOfCoverage_modified_refseq_sorted.txt",
    GATK3_jar_file="/home/nruiz/Desktop/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar"
  output:
    "{directory}/GATK{analysis}SeqPipeline_{file}/GATK_DepthOfCoverage/{sample}"
  params:
    "DepthOfCoverage"
  shell:
    "java -Djava.io.tmpdir=/home/NBS/new_tmp -jar {input.GATK3_jar_file} "
    "-T {params} -R {input.fa} -o {output} -I {input.bam} "
    "-geneList {input.gene_list} -L {input.interval_list}"
