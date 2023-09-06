import glob

(SAMPLES, FILES) = glob_wildcards("{sample}/{file}.fastq.gz")

def get_R1(wildcards):
  R1 = glob.glob(wildcards.sample + '/' + wildcards.sample + '_R1.fastq.gz')
  return(R1)

def get_R2(wildcards):
  R2 = glob.glob(wildcards.sample + '/' + wildcards.sample + '_R2.fastq.gz')
  return(R2)

rule all:
  input:
    expand('{sample}/seqyclean/{sample}_seqyclean', sample=SAMPLES)

rule seqyclean:
  input:
    R1 = get_R1,
    R2 = get_R2
  output:
    "{sample}/seqyclean/{sample}_seqyclean"
  log:
    "{sample}/seqyclean/seqyclean.error"
  shell:
    "seqyclean -1 {input.R1} -2 {input.R2} -qual -o {output} -c /home/NBS/WES_Reference_Data/phix.fasta > {log}"
