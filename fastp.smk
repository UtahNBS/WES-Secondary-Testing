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
    expand("{sample}/GATKTargetedSeqPipeline_{sample}_fastp", sample=SAMPLES),
    expand("{sample}/{sample}_fastp_R1.fastq.gz", sample=SAMPLES),
    expand("{sample}/{sample}_fastp_R2.fastq.gz", sample=SAMPLES),
    expand("{sample}/{sample}_fastp.html", sample=SAMPLES),
    expand("{sample}/{sample}_fastp.json", sample=SAMPLES)

rule fastp:
  input:
    R1 = get_R1,
    R2 = get_R2,
  output:
    dir = directory("{sample}/GATKTargetedSeqPipeline_{sample}_fastp"),
    Out1 = "{sample}/{sample}_fastp_R1.fastq.gz",
    Out2 = "{sample}/{sample}_fastp_R2.fastq.gz",
    html = "{sample}/{sample}_fastp.html",
    json = "{sample}/{sample}_fastp.json"
  log:
    "{sample}/{sample}_fastp.log"
  run:
    shell("mkdir {output.dir}") # This should be ok because I'm always just running for Targeted Analysis. I don't need to rerun PosCon or NTC.
    shell("fastp -i {input.R1} -I {input.R2} -o {output.Out1} -O {output.Out2} -h {output.html} -j {output.json} 2> {log}")
