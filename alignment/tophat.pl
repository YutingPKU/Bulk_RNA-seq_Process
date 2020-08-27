#!/usr/bin/perl

use strict;
use warnings;
use Cwd;

my $dir = getcwd;
mkdir "tophat";
mkdir "tophatBatch";
my $fat4 = "cn-short";
my $f4 = "cns";
my $cores = 20;
# my $fat4 = "fat4way";
# my $f4 = "f4w";
open OUT2,">","tophatBatch/tophatBatch.sh" or die $!;
print OUT2 "#!/bin/bash\n\n";
open IN,"<","sampleInfo.csv" or die $!;
while(<IN>){
	chomp;
	# my @in = split/\,/;
	# my $id = $in[0];
	my $id = $_;
	my $out = "$id"."_thout";
	my $fq1 = "$id"."_1.clean.fq.gz";
	my $fq2 = "$id"."_2.clean.fq.gz";
	# if($fat4 eq "fat4long"){
	# 	$fat4 = "fat4way";
	# 	$f4 = "f4w";
	# }else{
	# 	$fat4 = "fat4long";
	# 	$f4 = "f4l";
	# }
	open OUT,">","tophatBatch/$id.sh" or die "$!";
	print OUT "#!/bin/bash\n";
	print OUT "#SBATCH -J $id\n";
	print OUT "#SBATCH -o $id.%j.out\n";
	print OUT "#SBATCH -e $id.%j.err\n";
	print OUT "#SBATCH -p $fat4\n";
	print OUT "#SBATCH -N 1\n";
	print OUT "#SBATCH --ntasks-per-node=$cores\n";
	print OUT "#SBATCH --no-requeue\n";
	print OUT "#SBATCH -A lch3000_g1\n";
	print OUT "#SBATCH --qos=lch3000$f4\n\n";
	print OUT "export PATH=/lustre1/lch3000_pkuhpc/liuyf/python/Python-2.7.11:\$PATH\n";
#	print OUT "export PYTHONPATH=/lustre1/lch3000_pkuhpc/liuyf/python/Python-2.7.11\n";
#	print OUT "export PATH=\$HOME/.local/bin:\$PATH\n";
#	print OUT "export PYTHON_EGG_CACHE=\n";
#	print OUT "export PATH=/lustre1/lch3000_pkuhpc/liuyf/app/bowtie2-2.2.9:\$PATH\n";
	print OUT "export PATH=/lustre1/lch3000_pkuhpc/liuyt/software/tophat-2.1.1.Linux_x86_64:\$PATH\n\n";
	print OUT "tophat -G /lustre1/lch3000_pkuhpc/liuyf/ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf \\\n";
	print OUT "       -p $cores \\\n";
	print OUT "       -o $dir/tophat/$out \\\n";
	print OUT "       /lustre1/lch3000_pkuhpc/liuyt/genome/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome \\\n";
	print OUT "       $dir/cleanData/$fq1 \\\n";
	print OUT "       $dir/cleanData/$fq2\n";
	close OUT;

	print OUT2 "pkubatch $id.sh\n";
}
close IN;
close OUT2;

