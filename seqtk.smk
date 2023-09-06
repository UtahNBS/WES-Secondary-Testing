# Write it first to calculate 60x coverage

import glob

(SAMPLES, FILES) = glob_wildcards("{sample}/{sample}/{file}.fastq.gz")

rule all:
  input:
    expand("{sample}/{sample}_picard_hsmetrics.txt", sample=SAMPLES)

rule Picard_CollectHsMetrics:
  input:
    ref = "/home/NBS/WES_Reference_Data/human_g1k_v37.fasta.gz",
    bam="{sample}/bqsr_readgroups_markdups_sorted_{sample}.bam",
    interval_list = "/home/NBS/WES_Reference_Data/Illumine_Exome_CEX_TargetedRegions_v1.2.interval_list"
  output:
    "{sample}/{sample}_picard_hsmetrics.txt"
  log:
    '{sample}/{sample}_picard_hsmetrics.log'
  shell:
    'java -jar $picardJAR CollectHsMetrics TMP_DIR=/home/NBS/new_tmp I={input.bam} O={output} R={input.ref} BAIT_INTERVALS={input.interval_list} TARGET_INTERVALS={input.interval_list}'
