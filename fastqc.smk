import glob

(SAMPLES, FILES) = glob_wildcards("{sample}/{file}.fastq.gz")

rule all:
  input:
    expand("{sample}/fastqc/{file}_fastqc/", zip, sample=SAMPLES, file=FILES)

rule make_fastqc_dirs:
  run:
    commands = []
    for f in FILES:
      commands.append('mkdir' + ' ' + f[:-3]+'/'+f+'_fastqc')
    for c in commands:
      shell(c)

rule fastqc:
  input:
    lambda wildcards: [os.path.join(SAMPLES[i], x + '.fastq.gz') for i,x in enumerate(FILES) if x == wildcards.file]
  output:
    directory('{sample}/fastqc/{file}_fastqc/')
  log:
    '{sample}/fastqc/logs/{file}_fastqc.log'
  params:
    "--threads 8"
  shell:
    'fastqc {params} {input} -o {output} 2> {log}'
