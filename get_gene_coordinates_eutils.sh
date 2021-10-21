#!/bin/bash

# TODO: Genes PSAP and G6PD have multiple coordinates listed for the same version and chr
# This is because the gene symbols are aliases for other genes so those coordinates are included in the output

# Store parsed gene coordinate data retrieved from NCBI
updated_coordinates=()

# Map of chromosome accession numbers to chromosome number
declare -A chr_accession_map
chr_accession_map[NC_000001.10]=1
chr_accession_map[NC_000002.11]=2
chr_accession_map[NC_000003.11]=3
chr_accession_map[NC_000004.11]=4
chr_accession_map[NC_000005.9]=5
chr_accession_map[NC_000006.11]=6
chr_accession_map[NC_000007.13]=7
chr_accession_map[NC_000008.10]=8
chr_accession_map[NC_000009.11]=9
chr_accession_map[NC_000010.10]=10
chr_accession_map[NC_000011.9]=11
chr_accession_map[NC_000012.11]=12
chr_accession_map[NC_000013.10]=13
chr_accession_map[NC_000014.8]=14
chr_accession_map[NC_000015.9]=15
chr_accession_map[NC_000016.9]=16
chr_accession_map[NC_000017.10]=17
chr_accession_map[NC_000018.9]=18
chr_accession_map[NC_000019.9]=19
chr_accession_map[NC_000020.10]=20
chr_accession_map[NC_000021.8]=21
chr_accession_map[NC_000022.10]=22
chr_accession_map[NC_000023.10]=X
chr_accession_map[NC_000024.9]=Y
chr_accession_map[NC_012920.1]=MT

gene_list="MASTER_NBS_gene_disorder_list.bed"
new_bed_file="updated_NBS_gene_disorder_list.bed"
while IFS=$'\t' read -u 9 -r chr start end gene disorder
do
  output=$(esearch -db gene -query "$gene [GENE] AND human [ORGN]" | \
  efetch -format docsum | \
  xtract -pattern LocationHistType -element AnnotationRelease ChrAccVer -inc ChrStart ChrStop)
  readarray -t gene_info <<< "$output"
  for line in "${gene_info[@]}"
  do
    read version accession start end <<< ${line}
    if [[ $version == "105.20201022" && $accession == *"NC_"* ]]; then
      if [ ${chr_accession_map[$accession]} == $chr ]; then
        # Check if start is greater than end (reverse strand genes).
        # I want to make sure the smaller value is start and the greater value is end
        if (( start > end )); then
          temp=$start
          start=$end
          end=$temp
        fi
        updated_coordinates+=(${chr_accession_map[$accession]}$'\t'$start$'\t'$end$'\t'$gene$'\t'$disorder)
      fi
    fi
  done
done 9< "$gene_list"

for i in "${updated_coordinates[@]}"
do
  printf "$i\n" >> $new_bed_file
done
