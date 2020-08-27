#!/usr/bin/perl

use strict;
use warnings;
use Cwd;

my $dir = getcwd;
mkdir "cuffquant";
mkdir "cuffquantBatch";
open OUT2,">","cuffquantBatch/cuffquantBatch.sh" or die $!;
print OUT2 "#!/bin/bash\n\n";
open IN,"<","sampleInfo.csv" or die $!;
my $fat4 = "cn-long";
my $f4 = "cnl";
# my $fat4 = "fat4way";
# my $f4 = "f4w";
my $cores = 20;
while(<IN>){
	chomp;
	# my @in = split/\,/;
	# my $id = $in[0];
	my $id = $_;
	my $in = "$dir/tophat/"."$id"."_thout/accepted_hits.bam";
	my $out = "$dir/cuffquant/"."$id"."_cqout";
	# if($fat4 eq "fat4long"){
	# 	$fat4 = "fat4way";
	# 	$f4 = "f4w";
	# }else{
	# 	$fat4 = "fat4long";
	# 	$f4 = "f4l";
	# }
	open OUT,">","cuffquantBatch/$id.sh" or die "$!";
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
	print OUT "export PATH=/lustre1/lch3000_pkuhpc/liuyt/liuyf/python/Python-2.7.11:\$PATH\n";
#	print OUT "export PYTHONPATH=/lustre1/lch3000_pkuhpc/liuyt/liuyf/python/Python-2.7.11\n";
#	print OUT "export PATH=\$HOME/.local/bin:\$PATH\n";
#	print OUT "export PYTHON_EGG_CACHE=\n";
	print OUT "export PATH=/lustre1/lch3000_pkuhpc/liuyt/liuyf/app/cufflinks-2.2.1.Linux_x86_64:\$PATH\n\n";
	print OUT "cuffquant -no-update-check -p $cores \\\n";
	print OUT "          -o $out \\\n";
#	print OUT "          /lustre1/lch3000_pkuhpc/liuyt/genome/Homo_sapiens/UCSC/hg19/Annotation/Archives/archive-2014-06-02-13-47-56/Genes/genes.gtf \\\n";
	print OUT "	    /lustre1/lch3000_pkuhpc/liuyf/ref/Homo_sapiens/UCSC/hg19/Annotation/Genes/genes.gtf  \\\n";
	print OUT "          $in\n";
	close OUT;

	print OUT2 "pkubatch $id.sh\n";
}
close IN;
close OUT2;

