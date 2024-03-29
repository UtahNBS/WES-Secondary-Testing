'''

GATKExomeSeqPipeline.smk

Usage:
snakemake -ps GATKExomeSeqPipeline.smk --directory \
/home/NBS/Sequencing_Runs/Raw/NBS_Exome-89857768/WES_analysis --cores 48

TOOLS:
- BWA MEM (version 0.7.17-r1194-dirty)
- Picard
- Samtools
- GATK (version 4.0)
- SnpEff (version 4.3)
- VEP (version 96)
- vt normalize (version 0.5)

'''
import glob

(SAMPLES, FILES, TOOLS, READS) = glob_wildcards("{sample}/{file}_{tool}_{read}.fastq.gz")


# TODO: Add in command line option to specify reference set (RefSeq or Ensembl)

FASTA = "/home/NBS/WES_Reference_Data/human_g1k_v37.fasta.gz"
dbSNP = "/home/NBS/WES_Reference_Data/dbsnp_138.b37.vcf.gz"
ExomeTargets = "/home/NBS/WES_Reference_Data/Illumine_Exome_CEX_TargetedRegions_v1.2.bed"
Indels_1000G_phase1 = "/home/NBS/WES_Reference_Data/1000G_phase1.indels.b37.vcf.gz"
Mills_1000G_gold_standard_indels = "/home/NBS/WES_Reference_Data/Mills_and_1000G_gold_standard.indels.b37.vcf.gz"

def get_R1(wildcards):
  R1 = glob.glob(wildcards.sample + '/' + wildcards.sample + '_' + wildcards.tool + '_R1.fastq.gz')
  return(R1)

def get_R2(wildcards):
  R2 = glob.glob(wildcards.sample + '/' + wildcards.sample + '_' + wildcards.tool + '_R2.fastq.gz')
  return(R2)

rule all:
  input:
    #expand('{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}.sam.gz', zip, sample=SAMPLES, tool=TOOLS),
    #expand('{sample}/GATKExomeSeqPipeline_{sample}_{tool}/sorted_{sample}_{tool}.bam', zip, sample=SAMPLES, tool=TOOLS),
    #expand('{sample}/GATKExomeSeqPipeline_{sample}_{tool}/markdups_sorted_{sample}_{tool}.bam', zip, sample=SAMPLES, tool=TOOLS),
    expand('{sample}/GATKExomeSeqPipeline_{sample}_{tool}/logs/{sample}_{tool}_metrics_markdups_sorted.txt', zip, sample=SAMPLES, tool=TOOLS),
    expand('{sample}/GATKExomeSeqPipeline_{sample}_{tool}/readgroups_markdups_sorted_{sample}_{tool}.bam', zip, sample=SAMPLES, tool=TOOLS),
    expand('{sample}/GATKExomeSeqPipeline_{sample}_{tool}/readgroups_markdups_sorted_{sample}_{tool}.bam.bai', zip, sample=SAMPLES, tool=TOOLS),
    #expand('{sample}/GATKExomeSeqPipeline_{sample}_{tool}/recal_data_{sample}_{tool}.table', zip, sample=SAMPLES, tool=TOOLS),
    expand('{sample}/GATKExomeSeqPipeline_{sample}_{tool}/bqsr_readgroups_markdups_sorted_{sample}_{tool}.bam', zip, sample=SAMPLES, tool=TOOLS),
    expand('{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}.vcf.gz', zip, sample=SAMPLES, tool=TOOLS),
    expand('{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}_HC.bam', zip, sample=SAMPLES, tool=TOOLS),
    #expand('{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}_vep.vcf', zip, sample=SAMPLES, tool=TOOLS),
    #expand('{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}_vep_snpeff.vcf', zip, sample=SAMPLES, tool=TOOLS),
    expand('{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}_vep_snpeff_norm.vcf', zip, sample=SAMPLES, tool=TOOLS),
    expand('{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}_vep_snpeff_norm_clinvar.vcf', zip, sample=SAMPLES, tool=TOOLS)
  #shell:
    #"multiqc -f ."

rule bwa:
  input:
    ref=FASTA,
    R1=get_R1,
    R2=get_R2
  output:
    temp("{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}.sam.gz")
  log:
    "{sample}/GATKExomeSeqPipeline_{sample}_{tool}/logs/{sample}_{tool}_bwa.log"
  threads: 6
  shell:
    "bwa mem -t {threads} {input.ref} {input.R1} {input.R2} | gzip -3 > {output} 2> {log} "
    "|| true; touch {output}"

rule picard_sort_sam:
  input:
    "{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}.sam.gz"
  output:
    temp("{sample}/GATKExomeSeqPipeline_{sample}_{tool}/sorted_{sample}_{tool}.bam")
  params:
    tmp_dir="/home/UT_NBS/tmp"
  log:
    "{sample}/GATKExomeSeqPipeline_{sample}_{tool}/logs/{sample}_{tool}_picard_SortSam.log"
  shell:
    "java -jar $picardJAR SortSam TMP_DIR={params.tmp_dir} I={input} O={output} "
    "SORT_ORDER=coordinate 2> {log} || true; touch {output}"

rule picard_markdups:
  input:
    "{sample}/GATKExomeSeqPipeline_{sample}_{tool}/sorted_{sample}_{tool}.bam"
  output:
    bam=temp("{sample}/GATKExomeSeqPipeline_{sample}_{tool}/markdups_sorted_{sample}_{tool}.bam"),
    metrics="{sample}/GATKExomeSeqPipeline_{sample}_{tool}/logs/{sample}_{tool}_metrics_markdups_sorted.txt"
  params:
    tmp_dir="/home/UT_NBS/tmp"
  log:
    "{sample}/GATKExomeSeqPipeline_{sample}_{tool}/logs/{sample}_{tool}_picard_MarkDuplicates.log"
  shell:
    "java -jar $picardJAR MarkDuplicates I={input} O={output.bam} "
    "M={output.metrics} TMP_DIR={params.tmp_dir} 2> {log} || true; touch {output}"

rule picard_add_groups:
  input:
    "{sample}/GATKExomeSeqPipeline_{sample}_{tool}/markdups_sorted_{sample}_{tool}.bam"
  output:
    "{sample}/GATKExomeSeqPipeline_{sample}_{tool}/readgroups_markdups_sorted_{sample}_{tool}.bam"
  log:
    "{sample}/GATKExomeSeqPipeline_{sample}_{tool}/logs/{sample}_{tool}_picard_AddOrReplaceReadGroups.log"
  shell:
    # TODO: Change the default parameters?
    "java -jar $picardJAR AddOrReplaceReadGroups I={input} O={output} RGLB=LaneX "
    "RGPL=illumina RGPU=NONE RGSM=Sample1 2> {log} || true; touch {output}"

# BAM file has to be indexed before using BaseRecalibrator since I am specifying target intervals with the BED file
rule indexBam:
  input:
    "{sample}/GATKExomeSeqPipeline_{sample}_{tool}/readgroups_markdups_sorted_{sample}_{tool}.bam"
  output:
    "{sample}/GATKExomeSeqPipeline_{sample}_{tool}/readgroups_markdups_sorted_{sample}_{tool}.bam.bai"
  shell:
    "samtools index {input} {output}"

rule BaseRecalibrator:
  input:
    index="{sample}/GATKExomeSeqPipeline_{sample}_{tool}/readgroups_markdups_sorted_{sample}_{tool}.bam.bai",
    fa=FASTA,
    bam="{sample}/GATKExomeSeqPipeline_{sample}_{tool}/readgroups_markdups_sorted_{sample}_{tool}.bam",
    dbsnp=dbSNP,
    indels_1000G=Indels_1000G_phase1,
    mills=Mills_1000G_gold_standard_indels,
    targets=ExomeTargets
  output:
    temp("{sample}/GATKExomeSeqPipeline_{sample}_{tool}/recal_data_{sample}_{tool}.table")
  params:
    tmp_dir="/home/UT_NBS/tmp",
    padding="100"
  log:
    "{sample}/GATKExomeSeqPipeline_{sample}_{tool}/logs/{sample}_{tool}_GATK_BaseRecalibrator.log"
  shell:
    "gatk BaseRecalibrator -R {input.fa} -I {input.bam} -O {output} "
    "--known-sites {input.dbsnp} --known-sites {input.indels_1000G} --known-sites "
    "{input.mills} -L {input.targets} -ip {params.padding} --tmp-dir "
    "{params.tmp_dir} 2> {log} || true; touch {output}"

rule ApplyBQSR:
  input:
    fa=FASTA,
    bam="{sample}/GATKExomeSeqPipeline_{sample}_{tool}/readgroups_markdups_sorted_{sample}_{tool}.bam",
    table="{sample}/GATKExomeSeqPipeline_{sample}_{tool}/recal_data_{sample}_{tool}.table"
  output:
    "{sample}/GATKExomeSeqPipeline_{sample}_{tool}/bqsr_readgroups_markdups_sorted_{sample}_{tool}.bam"
  log:
    "{sample}/GATKExomeSeqPipeline_{sample}_{tool}/logs/{sample}_{tool}_GATK_ApplyBQSR.log"
  shell:
    "gatk ApplyBQSR -R {input.fa} -I {input.bam} --bqsr-recal-file "
    "{input.table} -O {output} 2> {log} || true; touch {output}"

rule HaplotypeCaller:
  input:
    fa=FASTA,
    bam="{sample}/GATKExomeSeqPipeline_{sample}_{tool}/bqsr_readgroups_markdups_sorted_{sample}_{tool}.bam",
    targets=ExomeTargets
  output:
    vcf="{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}.vcf.gz",
    bam="{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}_HC.bam"
  params:
    tmp_dir="/home/UT_NBS/tmp",
    padding="100"
  log:
    "{sample}/GATKExomeSeqPipeline_{sample}_{tool}/logs/{sample}_{tool}_GATK_HaplotypeCaller.log"
  shell:
    "gatk HaplotypeCaller -R {input.fa} -I {input.bam} -O {output.vcf} "
    "-bamout {output.bam} -L {input.targets} -ip {params.padding} --tmp-dir {params.tmp_dir} "
    "2> {log} || true; touch {output}"

rule SelectVariantsSNP:
  input:
    fa=FASTA,
    vcf="{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}.vcf.gz"
  output:
    temp("{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}_raw_snps.vcf")
  log:
    "{sample}/GATKExomeSeqPipeline_{sample}_{tool}/logs/{sample}_{tool}_GATK_SelectVariants_SNPs.log"
  shell:
    "gatk SelectVariants -R {input.fa} -V {input.vcf} -select-type SNP -O {output} 2> {log}"

rule SNPFilter:
  input:
    fa=FASTA,
    vcf="{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}_raw_snps.vcf"
  output:
    temp("{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}_filtered_snps.vcf")
  log:
    "{sample}/GATKExomeSeqPipeline_{sample}_{tool}/logs/{sample}_{tool}_GATK_VariantFiltration_SNPs.log"
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
    vcf="{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}.vcf.gz"
  output:
    temp("{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}_raw_indels.vcf")
  log:
    "{sample}/GATKExomeSeqPipeline_{sample}_{tool}/logs/{sample}_{tool}_GATK_SelectVariants_Indels.log"
  shell:
    "gatk SelectVariants -R {input.fa} -V {input.vcf} -select-type INDEL "
    "-O {output} 2> {log} || true; touch {output}"

rule IndelFilter:
  input:
    fa=FASTA,
    vcf="{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}_raw_indels.vcf"
  output:
    temp("{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}_filtered_indels.vcf")
  log:
    "{sample}/GATKExomeSeqPipeline_{sample}_{tool}/logs/{sample}_{tool}_GATK_VariantFiltration_Indels.log"
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
    snps="{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}_filtered_snps.vcf",
    indels="{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}_filtered_indels.vcf"
  output:
    temp("{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}_snps_indels.vcf")
  shell:
    "perl /home/NBS/scripts/combine_vcfs.pl {input} {output} || true; touch {output}"

rule picard_sort:
  input:
    "{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}_snps_indels.vcf"
  output:
    temp("{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}_snps_indels_sorted.vcf")
  shell:
    "java -jar $picardJAR SortVcf I={input} O={output} || true; touch {output}"

rule VEP:
  input:
    vcf="{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}_snps_indels_sorted.vcf",
    fa=FASTA
  output:
    temp("{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}_vep.vcf")
  shell:
    "vep --cache --offline --refseq --dir_cache /home/NBS/.vep "
    "--fasta {input.fa} --format vcf --fork 4 --vcf --hgvsg -e "
    "-i {input.vcf} -o {output} --no_stats || true; touch {output}"

rule SnpEff:
  input:
    "{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}_vep.vcf"
  output:
    temp("{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}_vep_snpeff.vcf")
  shell:
    "java -Xmx8g -jar $snpeffJAR -v hg19 -noStats {input} > {output} "
    "|| true; touch {output}"

rule vt_normalize:
  '''
  Normalize variant representation in VCF file using vt normalize
  '''
  input:
    vcf="{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}_vep_snpeff.vcf",
    fa=FASTA
  output:
    "{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}_vep_snpeff_norm.vcf"
  log:
    "{sample}/GATKExomeSeqPipeline_{sample}_{tool}/logs/{sample}_{tool}_vep_snpeff_vt_normalize.log"
  run:
    shell("vt normalize {input.vcf} -o {output} -r {input.fa} &> {log} || true; touch {output}")

rule clinvar_annotation:
  input:
    "{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}_vep_snpeff_norm.vcf"
  output:
    "{sample}/GATKExomeSeqPipeline_{sample}_{tool}/{sample}_{tool}_vep_snpeff_norm_clinvar.vcf"
  run:
    shell("perl /home/NBS/scripts/filter_gene_specific_variants_clinvar.pl {input} {output}")
