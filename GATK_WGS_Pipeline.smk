'''

GATK_WGS_Pipeline.smk
This is for the analysis of whole-genome sequence data.
The only modifications I made was to rename output directory folder and
remove exome targets file and padding for intervals.

Usage:
snakemake -ps GATK_WGS_Pipeline.smk --directory /dir --cores 48

Tools:
- BWA MEM version 0.7.17-r1194-dirty
- Picard 2.18.14-SNAPSHOT
- samtools version 1.8
- GATK version 4.0
- bcftools version 1.8
- SnpEff version 4.3
- VEP version 96
- vt normalize version 0.5

'''
import glob

(SAMPLES, FILES) = glob_wildcards("{sample}/{file}.fastq.gz")

FASTA = "/home/NBS/WES_Reference_Data/human_g1k_v37.fasta.gz"
dbSNP = "/home/NBS/WES_Reference_Data/dbsnp_138.b37.vcf.gz"

def get_R1(wildcards):
  R1 = glob.glob(wildcards.sample + '/' + wildcards.sample + '_R1.fastq.gz')
  return(R1)

def get_R2(wildcards):
  R2 = glob.glob(wildcards.sample + '/' + wildcards.sample + '_R2.fastq.gz')
  return(R2)

rule all:
  input:
    expand('{sample}/GATK_WGS_Pipeline/{sample}.sam.gz', sample=SAMPLES),
    expand('{sample}/GATK_WGS_Pipeline/sorted_{sample}.bam', sample=SAMPLES),
    expand('{sample}/GATK_WGS_Pipeline/markdups_sorted_{sample}.bam', sample=SAMPLES),
    expand('{sample}/GATK_WGS_Pipeline/logs/{sample}_metrics_markdups_sorted.txt', sample=SAMPLES),
    expand('{sample}/GATK_WGS_Pipeline/readgroups_markdups_sorted_{sample}.bam', sample=SAMPLES),
    expand('{sample}/GATK_WGS_Pipeline/readgroups_markdups_sorted_{sample}.bam.bai', sample=SAMPLES),
    expand('{sample}/GATK_WGS_Pipeline/recal_data_{sample}.table', sample=SAMPLES),
    expand('{sample}/GATK_WGS_Pipeline/bqsr_readgroups_markdups_sorted_{sample}.bam', sample=SAMPLES),
    expand('{sample}/GATK_WGS_Pipeline/{sample}.vcf.gz', sample=SAMPLES),
    expand('{sample}/GATK_WGS_Pipeline/{sample}_HC.bam', sample=SAMPLES),
    expand('{sample}/GATK_WGS_Pipeline/{sample}_raw_snps.vcf', sample=SAMPLES),
    expand('{sample}/GATK_WGS_Pipeline/{sample}_filtered_snps.vcf', sample=SAMPLES),
    expand('{sample}/GATK_WGS_Pipeline/{sample}_raw_indels.vcf', sample=SAMPLES),
    expand('{sample}/GATK_WGS_Pipeline/{sample}_filtered_indels.vcf', sample=SAMPLES),
    expand('{sample}/GATK_WGS_Pipeline/{sample}_snps_indels.vcf', sample=SAMPLES),
    expand('{sample}/GATK_WGS_Pipeline/{sample}_snpeff.vcf', sample=SAMPLES),
    expand('{sample}/GATK_WGS_Pipeline/{sample}_vep.vcf', sample=SAMPLES),
    expand('{sample}/GATK_WGS_Pipeline/{sample}_snpeff_norm.vcf', sample=SAMPLES),
    expand('{sample}/GATK_WGS_Pipeline/{sample}_vep_norm.vcf', sample=SAMPLES),
    expand('{sample}/GATK_WGS_Pipeline/{sample}_snpeff_norm_clinvar.vcf', sample=SAMPLES),
    expand('{sample}/GATK_WGS_Pipeline/{sample}_vep_norm_clinvar.vcf', sample=SAMPLES)

rule bwa:
  input:
    ref=FASTA,
    R1=get_R1,
    R2=get_R2
  output:
    '{sample}/GATK_WGS_Pipeline/{sample}.sam.gz'
  log:
    "{sample}/GATK_WGS_Pipeline/logs/{sample}_bwa.log"
  threads: 16
  shell:
    "bwa mem -t {threads} {input.ref} {input.R1} {input.R2} | gzip -3 > {output} 2> {log} "
    "|| true; touch {output}"

rule picard_sort:
  input:
    "{sample}/GATK_WGS_Pipeline/{sample}.sam.gz"
  output:
    "{sample}/GATK_WGS_Pipeline/sorted_{sample}.bam"
  log:
    "{sample}/GATK_WGS_Pipeline/logs/{sample}_picard_SortSam.log"
  shell:
    "java -jar $picardJAR SortSam TMP_DIR=/home/NBS/new_tmp I={input} O={output} "
    "SORT_ORDER=coordinate 2> {log} || true; touch {output}"

rule picard_markdups:
  input:
    "{sample}/GATK_WGS_Pipeline/sorted_{sample}.bam"
  output:
    bam="{sample}/GATK_WGS_Pipeline/markdups_sorted_{sample}.bam",
    metrics="{sample}/GATK_WGS_Pipeline/logs/{sample}_metrics_markdups_sorted.txt"
  params:
    tmp_dir="/home/NBS/new_tmp"
  log:
    "{sample}/GATK_WGS_Pipeline/logs/{sample}_picard_MarkDuplicates.log"
  shell:
    "java -jar $picardJAR MarkDuplicates I={input} O={output.bam} "
    "M={output.metrics} TMP_DIR={params.tmp_dir} 2> {log} || true; touch {output}"

rule picard_add_groups:
  input:
    "{sample}/GATK_WGS_Pipeline/markdups_sorted_{sample}.bam"
  output:
    "{sample}/GATK_WGS_Pipeline/readgroups_markdups_sorted_{sample}.bam"
  log:
    "{sample}/GATK_WGS_Pipeline/logs/{sample}_picard_AddOrReplaceReadGroups.log"
  shell:
    # TODO: Change the default parameters?
    "java -jar $picardJAR AddOrReplaceReadGroups I={input} O={output} RGLB=LaneX "
    "RGPL=illumina RGPU=NONE RGSM=Sample1 2> {log} || true; touch {output}"

rule indexBam:
  input:
    "{sample}/GATK_WGS_Pipeline/readgroups_markdups_sorted_{sample}.bam"
  output:
    "{sample}/GATK_WGS_Pipeline/readgroups_markdups_sorted_{sample}.bam.bai"
  shell:
    "samtools index {input} {output}"

rule BaseRecalibrator:
  input:
    index="{sample}/GATK_WGS_Pipeline/readgroups_markdups_sorted_{sample}.bam.bai",
    fa=FASTA,
    bam="{sample}/GATK_WGS_Pipeline/readgroups_markdups_sorted_{sample}.bam",
    dbsnp=dbSNP,
    indels_1000G='/home/NBS/WES_Reference_Data/1000G_phase1.indels.b37.vcf.gz',
    mills='/home/NBS/WES_Reference_Data/Mills_and_1000G_gold_standard.indels.b37.vcf.gz',
  output:
    "{sample}/GATK_WGS_Pipeline/recal_data_{sample}.table"
  params:
    tmp_dir="/home/NBS/new_tmp"
  log:
    "{sample}/GATK_WGS_Pipeline/logs/{sample}_GATK_BaseRecalibrator.log"
  shell:
    "gatk BaseRecalibrator -R {input.fa} -I {input.bam} -O {output} "
    "--known-sites {input.dbsnp} --known-sites {input.indels_1000G} --known-sites "
    "{input.mills} --tmp-dir {params.tmp_dir} 2> {log} || true; touch {output}"

rule ApplyBQSR:
  input:
    fa=FASTA,
    bam="{sample}/GATK_WGS_Pipeline/readgroups_markdups_sorted_{sample}.bam",
    table="{sample}/GATK_WGS_Pipeline/recal_data_{sample}.table"
  output:
    "{sample}/GATK_WGS_Pipeline/bqsr_readgroups_markdups_sorted_{sample}.bam"
  log:
    "{sample}/GATK_WGS_Pipeline/logs/{sample}_GATK_ApplyBQSR.log"
  shell:
    "gatk ApplyBQSR -R {input.fa} -I {input.bam} --bqsr-recal-file "
    "{input.table} -O {output} 2> {log} || true; touch {output}"

rule HaplotypeCaller:
  input:
    fa=FASTA,
    bam="{sample}/GATK_WGS_Pipeline/bqsr_readgroups_markdups_sorted_{sample}.bam",
  output:
    vcf="{sample}/GATK_WGS_Pipeline/{sample}.vcf.gz",
    bam="{sample}/GATK_WGS_Pipeline/{sample}_HC.bam"
  params:
    tmp_dir="/home/NBS/new_tmp"
  log:
    "{sample}/GATK_WGS_Pipeline/logs/{sample}_GATK_HaplotypeCaller.log"
  shell:
    "gatk HaplotypeCaller -R {input.fa} -I {input.bam} -O {output.vcf} "
    "-bamout {output.bam} --tmp-dir {params.tmp_dir} 2> {log} || true; touch {output}"

rule SelectVariantsSNP:
  input:
    fa=FASTA,
    vcf="{sample}/GATK_WGS_Pipeline/{sample}.vcf.gz"
  output:
    "{sample}/GATK_WGS_Pipeline/{sample}_raw_snps.vcf"
  log:
    "{sample}/GATK_WGS_Pipeline/logs/{sample}_GATK_SelectVariants_SNPs.log"
  shell:
    "gatk SelectVariants -R {input.fa} -V {input.vcf} -select-type SNP -O {output} 2> {log}"

rule SNPFilter:
  input:
    fa=FASTA,
    vcf="{sample}/GATK_WGS_Pipeline/{sample}_raw_snps.vcf"
  output:
    "{sample}/GATK_WGS_Pipeline/{sample}_filtered_snps.vcf"
  log:
    "{sample}/GATK_WGS_Pipeline/logs/{sample}_GATK_VariantFiltration_SNPs.log"
  params:
    QD="QD < 2.0",
    FS="FS > 60.0",
    MQ="MQ < 40.0",
    MQRankSum="MQRankSum < -12.5",
    ReadPosRankSum="ReadPosRankSum < -8.0",
    filtername="snp_filter"
  shell:
    "gatk VariantFiltration -R {input.fa} -V {input.vcf} --filter-expression "
    "'{params.QD} || {params.FS} || {params.MQ} || {params.MQRankSum} || "
    "{params.ReadPosRankSum}' --filter-name '{params.filtername}' -O {output} "
    "2> {log} || true; touch {output}"

rule SelectVariantsIndels:
  input:
    fa=FASTA,
    vcf="{sample}/GATK_WGS_Pipeline/{sample}.vcf.gz"
  output:
    "{sample}/GATK_WGS_Pipeline/{sample}_raw_indels.vcf"
  log:
    "{sample}/GATK_WGS_Pipeline/logs/{sample}_GATK_SelectVariants_Indels.log"
  shell:
    "gatk SelectVariants -R {input.fa} -V {input.vcf} -select-type INDEL "
    "-O {output} 2> {log} || true; touch {output}"

rule IndelFilter:
  input:
    fa=FASTA,
    vcf="{sample}/GATK_WGS_Pipeline/{sample}_raw_indels.vcf"
  output:
    "{sample}/GATK_WGS_Pipeline/{sample}_filtered_indels.vcf"
  log:
    "{sample}/GATK_WGS_Pipeline/logs/{sample}_GATK_VariantFiltration_Indels.log"
  params:
    QD="QD < 2.0",
    FS="FS > 200.0",
    ReadPosRankSum="ReadPosRankSum < -20.0",
    filtername="indel_filter"
  shell:
    "gatk VariantFiltration -R {input.fa} -V {input.vcf} --filter-expression "
    "'{params.QD} || {params.FS} || {params.ReadPosRankSum}' --filter-name "
    "'{params.filtername}' -O {output} 2> {log} || true; touch {output}"

rule CombineVCFs:
  input:
    snps="{sample}/GATK_WGS_Pipeline/{sample}_filtered_snps.vcf",
    indels="{sample}/GATK_WGS_Pipeline/{sample}_filtered_indels.vcf"
  output:
    "{sample}/GATK_WGS_Pipeline/{sample}_snps_indels.vcf"
  shell:
    "perl /home/NBS/scripts/combine_vcfs.pl {input} {output} || true; touch {output}"

rule bcftools_sort:
  input:
    "{sample}/GATK_WGS_Pipeline/{sample}_snps_indels.vcf"
  output:
    "{sample}/GATK_WGS_Pipeline/{sample}_snps_indels.vcf"
  log:
    "{sample}/GATK_WGS_Pipeline/logs/{sample}_bcftools_sort.log"
  shell:
    "bcftools sort {input} > tmp && mv tmp {output} 2> {log} || true; touch {output}"

rule variant_annotation:
  input:
    vcf="{sample}/GATK_WGS_Pipeline/{sample}_snps_indels.vcf",
    fa=FASTA
  output:
    snpeff="{sample}/GATK_WGS_Pipeline/{sample}_snpeff.vcf",
    vep="{sample}/GATK_WGS_Pipeline/{sample}_vep.vcf",
    vep_stats="{sample}/GATK_WGS_Pipeline/{sample}_vep_stats.html"
  log:
    snpeff="{sample}/GATK_WGS_Pipeline/logs/{sample}_snpeff.log",
    vep="{sample}/GATK_WGS_Pipeline/logs/{sample}_vep.log"
  params:
    "GRCh37.75"
  run:
    shell("java -Xmx8g -jar $snpeffJAR -v {params} -noStats {input.vcf} > {output.snpeff} 2> {log.snpeff}")
    shell("vep --cache --offline --refseq --dir_cache /home/NBS/.vep --fasta {input.fa} --format vcf "
          "--fork 4 --vcf --hgvs --hgvsg -i {input.vcf} -o {output.vep} "
          "-sf {output.vep_stats} 2> {log.vep}")

rule vt_normalize:
  '''
  Normalize variant representation in VCF file using vt normalize
  '''
  input:
    snpeff="{sample}/GATK_WGS_Pipeline/{sample}_snpeff.vcf",
    vep="{sample}/GATK_WGS_Pipeline/{sample}_vep.vcf",
    fa=FASTA
  output:
    snpeff_norm="{sample}/GATK_WGS_Pipeline/{sample}_snpeff_norm.vcf",
    vep_norm="{sample}/GATK_WGS_Pipeline/{sample}_vep_norm.vcf"
  log:
    snpeff_log="{sample}/GATK_WGS_Pipeline/logs/{sample}_snpeff_vt_normalize.log",
    vep_log="{sample}/GATK_WGS_Pipeline/logs/{sample}_vep_vt_normalize.log"
  run:
    shell("vt normalize {input.snpeff} -o {output.snpeff_norm} -r {input.fa} &> {log.snpeff_log} || true; touch {output.snpeff_norm}")
    shell("vt normalize {input.vep} -o {output.vep_norm} -r {input.fa} &> {log.vep_log} || true; touch {output.vep_norm}")

rule clinvar_annotation:
  input:
    snpeff="{sample}/GATK_WGS_Pipeline/{sample}_snpeff_norm.vcf",
    vep="{sample}/GATK_WGS_Pipeline/{sample}_vep_norm.vcf"
  output:
    snpeff_clinvar="{sample}/GATK_WGS_Pipeline/{sample}_snpeff_norm_clinvar.vcf",
    vep_clinvar="{sample}/GATK_WGS_Pipeline/{sample}_vep_norm_clinvar.vcf"
  run:
    shell("perl /home/NBS/scripts/filter_gene_specific_variants_clinvar.pl {input.snpeff} {output.snpeff_clinvar}")
    shell("perl /home/NBS/scripts/filter_gene_specific_variants_clinvar.pl {input.vep} {output.vep_clinvar}")
