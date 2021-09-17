'''

CF_Pipeline_unmasked.smk - GATK-based targeted analysis exome sequencing pipeline with CFTR polyT/polyTG mask removed.

Pipeline aligns entire exome (GRCh37) and then restricts the rest of the analysis
to gene(s) of interest.

Usage:
snakemake -ps CF_Pipeline_unmasked.smk --directory /path/to/dir --cores 48 > pipeline.log &

Tools:
- BWA MEM (version 0.7.17-r1194-dirty)
- Picard
- Samtools
- GATK (version 4.0)
- SnpEff (version 4.3)
- VEP (version 96)
- vt normalize (version 0.5)

'''
import glob
import vcf

(SAMPLES, FILES, TOOLS, READS) = glob_wildcards("{sample}/{file}_{tool}_{read}.fastq.gz")

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
    expand('{sample}/CF_Pipeline_unmasked_{sample}_{tool}/logs/{sample}_metrics_markdups_sorted.txt', zip, sample=SAMPLES, tool=TOOLS),
    expand('{sample}/CF_Pipeline_unmasked_{sample}_{tool}/bqsr_readgroups_markdups_sorted_{sample}.bam', zip, sample=SAMPLES, tool=TOOLS),
    expand('{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}_HC.bam', zip, sample=SAMPLES, tool=TOOLS),
    expand('{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}_vep_snpeff_filtered_norm.vcf', zip, sample=SAMPLES, tool=TOOLS),
    expand('{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}_vep_snpeff_filtered_norm_clinvar.vcf', zip, sample=SAMPLES, tool=TOOLS)

rule bwa:
  input:
    ref = FASTA,
    R1 = get_R1,
    R2 = get_R2
  output:
    temp('{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}.sam.gz')
  log:
    "{sample}/CF_Pipeline_unmasked_{sample}_{tool}/logs/{sample}_bwa.log"
  threads: 6
  shell:
    "bwa mem -t {threads} {input.ref} {input.R1} {input.R2} | gzip -3 > {output} "
    "2> {log} || true; touch {output}"

rule picard_sort_sam:
  input:
    "{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}.sam.gz"
  output:
    temp("{sample}/CF_Pipeline_unmasked_{sample}_{tool}/sorted_{sample}.bam")
  params:
    tmp_dir="/home/UT_NBS/tmp"
  log:
    "{sample}/CF_Pipeline_unmasked_{sample}_{tool}/logs/{sample}_picard_SortSam.log"
  shell:
    "java -jar $picardJAR SortSam TMP_DIR={params.tmp_dir} I={input} O={output} "
    "SORT_ORDER=coordinate 2> {log} || true; touch {output}"

rule index_bam_select_gene:
  input:
    "{sample}/CF_Pipeline_unmasked_{sample}_{tool}/sorted_{sample}.bam"
  output:
    temp("{sample}/CF_Pipeline_unmasked_{sample}_{tool}/sorted_{sample}.bam.bai")
  shell:
    "samtools index {input} {output}"

rule select_target_reads:
  input:
    index="{sample}/CF_Pipeline_unmasked_{sample}_{tool}/sorted_{sample}.bam.bai",
    bam="{sample}/CF_Pipeline_unmasked_{sample}_{tool}/sorted_{sample}.bam"
  output:
    temp("{sample}/CF_Pipeline_unmasked_{sample}_{tool}/sorted_target_reads_{sample}.bam")
  log:
    "{sample}/CF_Pipeline_unmasked_{sample}_{tool}/logs/{sample}_select_target_reads.log"
  shell:
    "samtools view -b -h {input.bam} '7:117120017-117308719' > {output} 2> {log} "
    "|| true; touch {output}"

rule picard_markdups:
  input:
    "{sample}/CF_Pipeline_unmasked_{sample}_{tool}/sorted_target_reads_{sample}.bam"
  output:
    bam=temp("{sample}/CF_Pipeline_unmasked_{sample}_{tool}/markdups_sorted_{sample}.bam"),
    metrics="{sample}/CF_Pipeline_unmasked_{sample}_{tool}/logs/{sample}_metrics_markdups_sorted.txt"
  params:
    tmp_dir="/home/UT_NBS/tmp"
  log:
    "{sample}/CF_Pipeline_unmasked_{sample}_{tool}/logs/{sample}_picard_MarkDuplicates.log"
  shell:
    "java -jar $picardJAR MarkDuplicates I={input} O={output.bam} "
    "M={output.metrics} TMP_DIR={params.tmp_dir} 2> {log} || true; touch {output}"

rule picard_add_groups:
  input:
    "{sample}/CF_Pipeline_unmasked_{sample}_{tool}/markdups_sorted_{sample}.bam"
  output:
    temp("{sample}/CF_Pipeline_unmasked_{sample}_{tool}/readgroups_markdups_sorted_{sample}.bam")
  log:
    "{sample}/CF_Pipeline_unmasked_{sample}_{tool}/logs/{sample}_picard_AddOrReplaceReadGroups.log"
  shell:
    "java -jar $picardJAR AddOrReplaceReadGroups I={input} O={output} "
    "RGLB=LaneX RGPL=illumina RGPU=NONE RGSM=Sample1 2> {log} || true; touch {output}"

rule index_bam_for_BaseRecalibrator:
  input:
    "{sample}/CF_Pipeline_unmasked_{sample}_{tool}/readgroups_markdups_sorted_{sample}.bam"
  output:
    temp("{sample}/CF_Pipeline_unmasked_{sample}_{tool}/readgroups_markdups_sorted_{sample}.bam.bai")
  shell:
    "samtools index {input} {output}"

rule BaseRecalibrator:
  input:
    index="{sample}/CF_Pipeline_unmasked_{sample}_{tool}/readgroups_markdups_sorted_{sample}.bam.bai",
    fa=FASTA,
    bam="{sample}/CF_Pipeline_unmasked_{sample}_{tool}/readgroups_markdups_sorted_{sample}.bam",
    dbsnp=dbSNP,
    indels_1000G=Indels_1000G_phase1,
    mills=Mills_1000G_gold_standard_indels,
    targets=ExomeTargets
  output:
    temp("{sample}/CF_Pipeline_unmasked_{sample}_{tool}/recal_data_{sample}.table")
  params:
    tmp_dir="/home/UT_NBS/tmp",
    padding="100"
  log:
    "{sample}/CF_Pipeline_unmasked_{sample}_{tool}/logs/{sample}_GATK_BaseRecalibrator.log"
  shell:
    "gatk BaseRecalibrator -R {input.fa} -I {input.bam} -O {output} "
    "--known-sites {input.dbsnp} --known-sites {input.indels_1000G} --known-sites "
    "{input.mills} -L {input.targets} -ip {params.padding} --tmp-dir "
    "{params.tmp_dir} 2> {log} || true; touch {output}"

rule ApplyBQSR:
  input:
    fa=FASTA,
    bam="{sample}/CF_Pipeline_unmasked_{sample}_{tool}/readgroups_markdups_sorted_{sample}.bam",
    table="{sample}/CF_Pipeline_unmasked_{sample}_{tool}/recal_data_{sample}.table"
  output:
    "{sample}/CF_Pipeline_unmasked_{sample}_{tool}/bqsr_readgroups_markdups_sorted_{sample}.bam"
  log:
    "{sample}/CF_Pipeline_unmasked_{sample}_{tool}/logs/{sample}_GATK_ApplyBQSR.log"
  shell:
    "gatk ApplyBQSR -R {input.fa} -I {input.bam} --bqsr-recal-file {input.table} "
    "-O {output} 2> {log} || true; touch {output}"

rule HaplotypeCaller:
  input:
    fa=FASTA,
    bam="{sample}/CF_Pipeline_unmasked_{sample}_{tool}/bqsr_readgroups_markdups_sorted_{sample}.bam",
    targets=ExomeTargets
  output:
    vcf=temp("{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}.vcf.gz"),
    bam="{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}_HC.bam"
  params:
    tmp_dir="/home/UT_NBS/tmp"
  log:
    "{sample}/CF_Pipeline_unmasked_{sample}_{tool}/logs/{sample}_GATK_HaplotypeCaller.log"
  shell:
    "gatk HaplotypeCaller -R {input.fa} -I {input.bam} -O {output.vcf} "
    "-bamout {output.bam} -L {input.targets} --tmp-dir {params.tmp_dir} "
    "2> {log} || true; touch {output}"

rule SelectVariantsSNP:
  input:
    fa=FASTA,
    vcf="{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}.vcf.gz"
  output:
    temp("{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}_raw_snps.vcf")
  log:
    "{sample}/CF_Pipeline_unmasked_{sample}_{tool}/logs/{sample}_GATK_SelectVariants_SNPs.log"
  shell:
    "gatk SelectVariants -R {input.fa} -V {input.vcf} -select-type SNP "
    "-O {output} 2> {log} || true; touch {output}"

rule SNPFilter:
  input:
    fa=FASTA,
    vcf="{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}_raw_snps.vcf"
  output:
    temp("{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}_filtered_snps.vcf")
  log:
    "{sample}/CF_Pipeline_unmasked_{sample}_{tool}/logs/{sample}_GATK_VariantFiltration_SNPs.log"
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
    vcf="{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}.vcf.gz"
  output:
    temp("{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}_raw_indels.vcf")
  log:
    "{sample}/CF_Pipeline_unmasked_{sample}_{tool}/logs/{sample}_GATK_SelectVariants_Indels.log"
  shell:
    "gatk SelectVariants -R {input.fa} -V {input.vcf} -select-type INDEL "
    "-O {output} 2> {log} || true; touch {output}"

rule IndelFilter:
  input:
    fa=FASTA,
    vcf="{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}_raw_indels.vcf"
  output:
    temp("{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}_filtered_indels.vcf")
  log:
    "{sample}/CF_Pipeline_unmasked_{sample}_{tool}/logs/{sample}_GATK_VariantFiltration_Indels.log"
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
    snps="{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}_filtered_snps.vcf",
    indels="{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}_filtered_indels.vcf"
  output:
    temp("{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}_snps_indels.vcf")
  shell:
    "perl /home/NBS/scripts/combine_vcfs.pl {input} {output} || true; touch {output}"

rule picard_sort:
  input:
    "{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}_snps_indels.vcf"
  output:
    temp("{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}_snps_indels_sorted.vcf")
  shell:
    "java -jar $picardJAR SortVcf I={input} O={output} || true; touch {output}"

rule VEP:
  input:
    vcf="{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}_snps_indels_sorted.vcf",
    fa=FASTA
  output:
    temp("{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}_vep.vcf")
  shell:
    "vep --cache --offline --refseq --dir_cache /home/NBS/.vep --fasta "
    "{input.fa} --format vcf --fork 4 --vcf --hgvsg -e -i {input.vcf} "
    "-o {output} --no_stats || true; touch {output}"

rule SnpEff:
  input:
    "{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}_vep.vcf"
  output:
    temp("{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}_vep_snpeff.vcf")
  shell:
    "java -Xmx8g -jar $snpeffJAR -v hg19 -noStats {input} > {output} "
    "|| true; touch {output}"

rule filter_overlapping_genes:
  '''
  Uses script to filter out genes that are not associated with NBS disorders but
  overlap genes that are associated with NBS disorders.
  This is done post-variant annotation because the annotations will be made for
  all genes that the variant affects.
  '''
  input:
    vcf="{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}_vep_snpeff.vcf",
    bed="{sample}/{sample}.bed"
  output:
    temp("{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}_vep_snpeff_filtered.vcf")
  shell:
    "perl /home/NBS/scripts/filter_overlapping_genes.pl {input.vcf} {input.bed} {output}"

rule vt_normalize:
  '''
  Normalize variant representation in VCF file using vt normalize
  '''
  input:
    vcf="{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}_vep_snpeff_filtered.vcf",
    fa=FASTA
  output:
    "{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}_vep_snpeff_filtered_norm.vcf"
  log:
    "{sample}/CF_Pipeline_unmasked_{sample}_{tool}/logs/{sample}_vt_normalize.log"
  shell:
    "vt normalize {input.vcf} -o {output} -r {input.fa} &> {log} || true; touch {output}"

rule clinvar_annotation:
  '''
  Use script to annotate variants using ClinVar VCF file. This will soon be
  replaced by querying the variant database.
  '''
  input:
    "{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}_vep_snpeff_filtered_norm.vcf"
  output:
    "{sample}/CF_Pipeline_unmasked_{sample}_{tool}/{sample}_vep_snpeff_filtered_norm_clinvar.vcf"
  shell:
    "perl /home/NBS/scripts/filter_gene_specific_variants_clinvar.pl {input} {output}"