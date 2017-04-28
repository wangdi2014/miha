#!/usr/bin/perl
use threads; 

#Program searches for SNP associated with hematopoetic processes
#Input - file with human SNPs (NCBI) and list of mutations arranged by tissues (proteinatlas)
#Output - vcf list of hematopoetic SNPs resulting in amino acid change;
#Fasta sequences of peptides with aa substitutions is also given
#pbzip2(optional) and snpEff (required) are needed

#Written by Andrey Shelenkov, fallandar@gmail.com
#Last update: April 28, 2017

#source files needed:
#ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b149_GRCh37p13/VCF/All_20161121.vcf.gz
#http://www.proteinatlas.org/download/normal_tissue.csv.zip
#http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
#snpEff program is needed (http://snpeff.sourceforge.net/download.html)

#convert to bz2 before running like
#system ("gunzip All_20161121.vcf.gz && pbzip2 All_20161121.vcf");
#system ("unzip normal_tissue.csv.zip && pbzip2 normal_tissue.csv");

#output files - All_20161121_common_hemato_grch37_aachange.vcf.bz2 and All_20161121_common_hemato_grch37_aachange.fasta.bz2

sub open_file($)
{
    my ($filename) = @_;
    if ($filename=~/\.gz$/)
   {open ($fd,"zcat $filename |") || die "Can not open $filename: $!\n";}
    elsif	($filename=~/\.bz2$/)
    {open ($fd,"pbzip2 -c -d $filename |") || die "Can not open $filename: $!\n";}
    else
   {	open ($fd,$filename) || die "Can not open $filename: $!\n";}

    return $fd;
}

sub read_vcf_by_chrom
{
	($name,$chrom)=@_;
	$fh = open_file($name);
	open (OUT,">$chrom.out") || die "Can not open $chrom.out for writing";
  while (<$fh>)
  {
  	#next if (index($_,"#")==0);
  	next if (index($_,$chrom)!=0);
  	chomp;
  	@arr=split /\t/;
  	$rec=$arr[7];
  	foreach $key (keys (%genes))
  	{
  	if ($rec=~/GENEINFO=$key\:/ || $rec=~/\|$key\:/)
  		{print OUT "$_\n"; last;}
  	}
  	
  }
  close $fh;
  close OUT;
}

#get common mutations
system "pbzip2 -cd All_20161121.vcf.bz2 | grep \"COMMON=1\" > All_20161121_common.vcf" unless (-e "All_20161121_common.vcf");
#filter by hemato localization and approved status
system "pbzip2 -cd normal_tissue.csv.bz2 | grep -P \"hemato.*High.*(Approved|Supported)\"  > normal_tissue_hemato.csv" unless (-e "normal_tissue_hemato.csv");

#read gene names from csv
$fh = open_file("normal_tissue_hemato.csv");
%genes=();
$i=0;
while (<$fh>)
{
	unless($i)
	{$i++;next;}
	@arr=split /,/;
	$arr[1]=~s/\"//g;
	$genes{$arr[1]}++;

}
close $fh;

#get vcf records for hemato genes

unless ($Config{useithreads})
#threads are not supported
{
	for ($i=1;$i<=22;$i++)
	{read_vcf_by_chrom("All_20161121_common.vcf",$i);}
	read_vcf_by_chrom("All_20161121_common.vcf","Y");
	read_vcf_by_chrom("All_20161121_common.vcf","X");
} 
else 
{
 @threads=();
  for ($i=1;$i<=22;$i++)
  {$threads[$i]=threads->create(\&read_vcf_by_chrom,"All_20161121_common.vcf",$i);}
  $threads[23]=threads->create(\&read_vcf_by_chrom,"All_20161121_common.vcf","Y");
  $threads[24]=threads->create(\&read_vcf_by_chrom,"All_20161121_common.vcf","X");
  for ($i=1;$i<=24;$i++)
  {$threads[$i]->join;} 
}

#join 'out' files to new vcf and add header from initital vcf
system "cat *out > All_20161121_common_hemato.vcf";
system "cat All_20161121_common_hemato.vcf | sort -V > All_20161121_common_hemato.vcf1";
system "mv All_20161121_common_hemato.vcf1 All_20161121_common_hemato.vcf";
system "pbzip2 -cd All_20161121.vcf.bz2 | head -n 56 > 1.vcf";
system "cat 1.vcf All_20161121_common_hemato.vcf > All_20161121_common_hemato1.vcf";
system "mv All_20161121_common_hemato1.vcf All_20161121_common_hemato.vcf";

#run snpEff to get amino acid changes for SNPs; GRCh37.75 should be downloaded
system "java -Xmx36g -jar snpEff.jar -v GRCh37.75 -classic -no-intergenic ./All_20161121_common_hemato.vcf > All_20161121_common_hemato_grch37.vcf";

#filter by amino acid changes (skip synonymous)
system "grep \"EFF=NON_SYNONYMOUS_CODING\" All_20161121_common_hemato_grch37.vcf > All_20161121_common_hemato_grch37_aachange.vcf";

#add header to vcf
system "cat 1.vcf All_20161121_common_hemato_grch37_aachange.vcf > All_20161121_common_hemato_grch37_aachangex.vcf";
system "mv All_20161121_common_hemato_grch37_aachangex.vcf All_20161121_common_hemato_grch37_aachange.vcf";

if (-e "Homo_sapiens.GRCh38.pep.all.fa.gz")
{system "gunzip Homo_sapiens.GRCh38.pep.all.fa.gz";}

#print amino acid sequences for variants to separate file
system "python ./snpeff_to_peptides.py -i All_20161121_common_hemato_grch37_aachange.vcf -p Homo_sapiens.GRCh38.pep.all.fa > All_20161121_common_hemato_grch37_aachange.fasta 2>/dev/null";

#clean and archive
system "rm -f *out";
system "pbzip2 All_20161121_common.vcf All_20161121_common_hemato.vcf All_20161121_common_hemato_grch37.vcf All_20161121_common_hemato_grch37_aachange.vcf All_20161121_common_hemato_grch37_aachange.fasta";

