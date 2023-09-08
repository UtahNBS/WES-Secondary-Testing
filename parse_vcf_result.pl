#!/usr/bin/perl
use strict;
use warnings;

# Parse VCFs and combine output into txt file
# Usage: perl /home/NBS/scripts/parse_vcf_result.pl annotated.vcf out.txt

my $vcf = $ARGV[0];
my $out = $ARGV[1];
open (my $fh1, "<", $vcf) or die $!;
open (my $fh2, ">", $out) or die $!;

my %zygosity_map = ('0/1' => 'Heterozygous',
                    '1/0' => 'Heterozygous',
                    '1/1' => 'Homozygous',
                    '0/0' => 'Homozygous',
                    '1/2' => 'Heterozygous');
my %variants;

# Parse VCF
while (<$fh1>) {
  chomp $_;
  next if $_ =~ /^#/;
  # Catches files with no variants
  if ($_) {

    my ($chr, $pos, $id, $ref, $alt, $qual, 
        $filter, $info, $gt_format, $gt) = split("\t", $_);
    my $key = $chr.":".$pos.":".$ref.":".$alt;
    my ($zygosity, $ad, $dp, $gq, $pl) = split(":", $gt);
    my ($ref_depth, $alt_depth) = split(",", $ad);
    $variants{$key}{'chr'} = $chr;
    $variants{$key}{'pos'} = $pos;
    $variants{$key}{'ref'} = $ref;
    $variants{$key}{'alt'} = $alt; # Using alt from alt col instead of annotation. 
    $variants{$key}{'zygosity'} = $zygosity;
    $variants{$key}{'mapped_zygosity'} = $zygosity_map{$zygosity};
    $variants{$key}{'read_dp'} = $dp;
    $variants{$key}{'ref_dp'} = $ref_depth;
    $variants{$key}{'alt_dp'} = $alt_depth;
    # ClinVar annotations if available
    $variants{$key}{'clinvar_CLNSIG'} = "";
    $variants{$key}{'clinvar_ALLELEID'} = "";
    $variants{$key}{'clinvar_CLNHGVS'} = "";
    $variants{$key}{'clinvar_CLNDN'} = "";

    my @info = split(";", $info);

    foreach my $i (@info) {

      my ($k, $v) = split("=", $i);

      # Parse VEP annotations
      if ($k eq 'CSQ') {
        my @annos = split(",", $v);
        #foreach my $a (@annos) {
          #my @ann = split(/\|/, $a);
          # VEP can output multiple annotations for a single variant because it will annotate 
          # multiple transcripts. In this script, I am only looking at the first annotation.
          # Generally we are only interested in the first annotation however this might not 
          # always be the case.
          my @ann = split(/\|/, $annos[0]);
          $variants{$key}{'gene'} = $ann[3]; # only keeping last gene name
	      # Variant effect
          $variants{$key}{'vep_effect'} = $ann[1];
          # Impact
          $variants{$key}{'impact'} = $ann[2];
          # Exon
          if ($ann[8]) {
            $variants{$key}{'exon'} = $ann[8];
          }
          else {
            $variants{$key}{'exon'} = "";
          }
	      # coding HGVS
          if ($ann[10]) {
            $variants{$key}{'hgvs_c'} = $ann[10];
          }
          else {
            $variants{$key}{'hgvs_c'} = "";
          }
	      # protein HGVS
          if ($ann[11]) {
            $variants{$key}{'hgvs_p'} = $ann[11];
          }
          else {
            $variants{$key}{'hgvs_p'} = "";
          }
	      # dbSNP ID
	      if ($ann[17]) {
	        my @ids = split(/&/, $ann[17]);
	        foreach (@ids) {
	          if ($_ =~ /^rs/) {
                $variants{$key}{'dbsnp_id'} = $_;
	          }
 	          else {
		        $variants{$key}{'dbsnp_id'} = "";
	          }
            }
          }
          else {
	        $variants{$key}{'dbsnp_id'} = "";
	      }
	      # variant type
	      if ($ann[21]) {
	        $variants{$key}{'type'} = $ann[21];
	      }
          else {
            $variants{$key}{'type'} = "";
	      }
	      # SIFT
	      if ($ann[37]) {
	        $variants{$key}{'sift'} = $ann[37];
	      }
	      else {
	        $variants{$key}{'sift'} = "";
	      }	
	      # PolyPhen2
	      if ($ann[38]) {
	        $variants{$key}{'polyphen2'} = $ann[38];
	      }
	      else {
            $variants{$key}{'polyphen2'} = "";
          }	
          # 1000 Genomes allele frequency
          if ($ann[43]) {
	        $variants{$key}{'1000G_af'} = $ann[43];
	      }
	      else {
	        $variants{$key}{'1000G_af'} = "";
          }
          # gnomAD allele frequency
          if ($ann[51]) {
            $variants{$key}{'gnomAD_af'} = $ann[51];
          }
          else {
            $variants{$key}{'gnomAD_af'} = "";
          }
        #}
      }

      # Parse SnpEff annotation
      if ($k eq 'ANN') {
        my @anno = split(/\|/, $v);
        # variant effect
        $variants{$key}{'snpeff_effect'} = $anno[1];
      }

      # Parse ClinVar annotations if present
      # ClinVar interpretation
      if ($k eq 'CLNSIG') {
        if (exists $variants{$key}) {
          if ($v) {
            $variants{$key}{'clinvar_CLNSIG'} = $v;
          }
          else {
            $variants{$key}{'clinvar_CLNSIG'} = "";
          }
        }
        else {
          $variants{$key}{'clinvar_CLNSIG'} = "VARIANT NOT IN ORIGINAL VCF FILE\n";
        }
      }
      # ClinVar allele ID
      if ($k eq 'ALLELEID') {
	    if (exists $variants{$key}) {
          if ($v) {
            $variants{$key}{'clinvar_ALLELEID'} = $v;
	      }
	    else {
            $variants{$key}{'clinvar_ALLELEID'} = "";
	    }	
    	}
        else {
	      $variants{$key}{'clinvar_ALLELEID'} = "VARIANT NOT IN ORIGINAL VCF FILE\n";
        } 
      }
      # ClinVar genomic HGVS annotation
      if ($k eq 'CLNHGVS') {
        if (exists $variants{$key}) {
          if ($v) {
            $variants{$key}{'clinvar_CLNHGVS'} = $v;
          }
          else {
            $variants{$key}{'clinvar_CLNHGVS'} = "";
          }
        } 
        else {
          $variants{$key}{'clinvar_CLNHGVS'} = "VARIANT NOT IN ORIGINAL VCF FILE\n";
        }
      }
      # ClinVar associated disorders
      if ($k eq 'CLNDN') {
	    if (exists $variants{$key}) {
	      if ($v) {
	        $variants{$key}{'clinvar_CLNDN'} = $v;
	      }
	      else {
	        $variants{$key}{'clinvar_CLNDN'} = "";
	      }
	    }
        else {
	      $variants{$key}{'clinvar_CLNDN'} = "VARIANT NOT IN ORIGINAL VCF FILE\n";
        }
      }
    }
  }

  else {
    # Prints message to screen that no variants are in file
    print "No variants in file\n";
  }

}

print $fh2 "chromosome\tposition\treference_allele\talternate_allele\tgene\texon\timpact\t".
           "vep_effect\tsnpeff_effect\tinterpretation\ttype\trsid\t".
           "clinvar_allele_id\thgvs_g\thgvs_c\thgvs_p\tzygosity\treference_allele_coverage_depth\t".
           "alternate_allele_coverage_depth\ttotal_variant_coverage_depth\talternate_allele_fraction\t".
           "sift\tpolyphen2\tgnomad_af\t1000G_af\tassociated_disorders(CLNDN)\n";

foreach my $key (keys %variants) {
  # Calculate alternate allelic fraction
  my $alt_allelic_fraction = ($variants{$key}{'alt_dp'})/($variants{$key}{'read_dp'}) * 100;

  print $fh2 "$variants{$key}{'chr'}\t".
             "$variants{$key}{'pos'}\t".
             "$variants{$key}{'ref'}\t".
             "$variants{$key}{'alt'}\t".
             "$variants{$key}{'gene'}\t".
             "$variants{$key}{'exon'}\t".
             "$variants{$key}{'impact'}\t".
             "$variants{$key}{'vep_effect'}\t".
             "$variants{$key}{'snpeff_effect'}\t".
             "$variants{$key}{'clinvar_CLNSIG'}\t".
             "$variants{$key}{'type'}\t".
             "$variants{$key}{'dbsnp_id'}\t".
             "$variants{$key}{'clinvar_ALLELEID'}\t".
             "$variants{$key}{'clinvar_CLNHGVS'}\t".
             "$variants{$key}{'hgvs_c'}\t".
             "$variants{$key}{'hgvs_p'}\t".
             "$variants{$key}{'mapped_zygosity'}\t".
             "$variants{$key}{'ref_dp'}\t".
             "$variants{$key}{'alt_dp'}\t".
             "$variants{$key}{'read_dp'}\t".
             "$alt_allelic_fraction%\t".
             "$variants{$key}{'sift'}\t".
             "$variants{$key}{'polyphen2'}\t".
             "$variants{$key}{'gnomAD_af'}\t".
             "$variants{$key}{'1000G_af'}\t".
             "$variants{$key}{'clinvar_CLNDN'}\n";
}


close $fh1;
close $fh2;