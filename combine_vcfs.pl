#!/usr/bin/perl
use strict;
use warnings;

# Combine SNP and Indel VCFs after filtering from HaplotypeCaller vcf output
my $snps = $ARGV[0];
my $indels = $ARGV[1];
my $combo = $ARGV[2];
open (my $fh1, "<", $snps) or die "Can't open $snps $!\n";
open (my $fh2, "<", $indels) or die "Can't open $indels $!\n";
open (my $fh3, ">", $combo) or die "Can't open $combo $!\n";

# added defined() to suppress
# warning Value of <HANDLE> construct can be "0"; test with defined()
# this is a known bug
while (defined((my $line1 = <$fh1>)) and defined((my $line2 = <$fh2>))) {
  chomp $line1;
  chomp $line2;
  # Check if header line
  # These are the ones I want to compare
  if (($line1 =~ /^#/) && ($line2 =~ /^#/)) {
    if ($line1 eq $line2) {
      print $fh3 "$line1\n";
    }
    else {
      print $fh3 "$line1\n";
      print $fh3 "$line2\n";
    }
  }
}

close $fh1;
close $fh2;

open (my $fh4, "<", $snps) or die "Can't open $snps $!\n";
open (my $fh5, "<", $indels) or die "Can't open $indels $!\n";

while (<$fh4>) {
  chomp $_;
  next if $_ =~ /^#/;
  print $fh3 "$_\n";
}

while (<$fh5>) {
  chomp $_;
  next if $_ =~ /^#/;
  print $fh3 "$_\n";
}

close $fh4;
close $fh5;
close $fh3;

# Removed this from script and added as another rule in snakefile
# Sort VCF file using bcftools
#system("bcftools sort filtered_snps_indels.vcf > tmp && mv tmp filtered_snps_indels.vcf");
