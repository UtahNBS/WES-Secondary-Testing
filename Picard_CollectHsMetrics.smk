# You can use a SAM or BAM file as input
import glob

(DIRECTORIES, SAMPLES, FILES) = glob_wildcards("{directory}/GATKExomeSeqPipeline_{sample}/bqsr_readgroups_markdups_sorted_{file}.bam")

rule all:
  input:
    expand("{directory}/GATKExomeSeqPipeline_{sample}/{sample}_picard_hsmetrics.txt", zip, directory=DIRECTORIES, sample=SAMPLES)

rule Picard_CollectHsMetrics:
  input:
    ref = "/home/NBS/WES_Reference_Data/human_g1k_v37.fasta.gz",
    bam="{directory}/GATKExomeSeqPipeline_{sample}/bqsr_readgroups_markdups_sorted_{sample}.bam",
    interval_list = "/home/NBS/WES_Reference_Data/Illumine_Exome_CEX_TargetedRegions_v1.2.interval_list"
  output:
    "{directory}/GATKExomeSeqPipeline_{sample}/{sample}_picard_hsmetrics.txt"
  log:
    "{directory}/GATKExomeSeqPipeline_{sample}/{sample}_picard_hsmetrics.log"
  shell:
    "java -jar $picardJAR CollectHsMetrics TMP_DIR=/home/UT_NBS/tmp "
    "I={input.bam} O={output} R={input.ref} "
    "BAIT_INTERVALS={input.interval_list} TARGET_INTERVALS={input.interval_list}"
