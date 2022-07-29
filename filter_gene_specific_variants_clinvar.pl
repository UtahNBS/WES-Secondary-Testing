#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

# This script annotates variants with Clinvar
# Usage: perl /home/NBS/scripts/filter_gene_specific_variants_clinvar.pl sample_vep_snpeff.vcf sample_vep_snpeff_clinvar.vcf

my $clinvar = '/home/NBS/WES_Reference_Data/clinvar.vcf.gz';
my $vcf = $ARGV[0];
my $out = $ARGV[1];
#my $chr_interest = $ARGV[2];
#my $gene_interest = $ARGV[3];
open (my $fh1, '-|', 'gzip -dc /home/NBS/WES_Reference_Data/clinvar.vcf.gz') or die $!;
open (my $fh2, '<', $vcf) or die $!;
open (my $fh3, '>', $out) or die $!;

my %clinvar_gene;

# Build hash using Chr:VCFStart:Ref:Alt as key
# If Chr is not included, ClinVar variants will be annotated to the wrong variants in the vcf file.
while (<$fh1>) {
  chomp $_;
  next if $_ =~ /^#/;
  my ($chr, $start, $variation_id, $ref, $alt, $qual, $filter, $info) = split("\t", $_);
  #next if $chr ne $ARGV[2];
  my @info = split(';', $info);
  my $hash_key = $chr.":".$start.":".$ref.":".$alt;
  foreach (@info) {
    my ($key, $value) = split('=', $_);
    if ($key eq 'ALLELEID') {
      $clinvar_gene{$hash_key}{'ALLELEID'} = $value;
    }
    if ($key eq 'CLNSIG') {
      $clinvar_gene{$hash_key}{'CLNSIG'} = $value;
    }
    if ($key eq 'CLNDN') {
      $clinvar_gene{$hash_key}{'CLNDN'} = $value;
    }
    if ($key eq 'CLNDISDB') {
      $clinvar_gene{$hash_key}{'CLNDISDB'} = $value;
    }
    if ($key eq 'CLNHGVS') {
      $clinvar_gene{$hash_key}{'CLNHGVS'} = $value;
    }
    if ($key eq 'CLNVC') {
      $clinvar_gene{$hash_key}{'CLNVC'} = $value;
    }
    if ($key eq 'CLNVCSO') {
      $clinvar_gene{$hash_key}{'CLNVCSO'} = $value;
    }
    if ($key eq 'RS') {
      $clinvar_gene{$hash_key}{'RS'} = $value;
    }
    if ($key eq 'GENEINFO') {
      $clinvar_gene{$hash_key}{'GENEINFO'} = $value;
    }
  }
}

# Filter annotated sample file for variants of interest
while (<$fh2>) {
  chomp $_;
  if ($_ =~ '^#') {
    print $fh3 "$_\n";
  }
  else {
    my @cols = split("\t", $_);
    my $lookup = $cols[0].":".$cols[1].":".$cols[3].":".$cols[4];
    if (exists $clinvar_gene{$lookup}) {
      my @clinvar_info;
      foreach (keys $clinvar_gene{$lookup}) {
        push (@clinvar_info, $_."=".$clinvar_gene{$lookup}{$_});
      }
      my $new_info = join(';', $cols[7], @clinvar_info);
      # Add edited info column back to original VCF line
      splice @cols, 7, 1, $new_info;
      print $fh3 join("\t", @cols);
      print $fh3 "\n";
    }
  }
}

close $fh1;
close $fh2;
close $fh3;
