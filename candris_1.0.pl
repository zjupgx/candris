#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Statistics::Distributions;

=head1 Description
	
	CANDRIS: Predicting Cancer Driver Sites ; V1.0; 2018.01.08

=head1 Contact
	
	Gu Xun		xgu@iastate.edu
	Zhou Zhan	zhanzhou@zju.edu.cn
	Wu Jingcheng	21619014@zju.edu.cn
	Zhao Wenyi	21719052@zju.edu.cn
	Su Zhixi	zxsu@fudan.edu.cn
	Zou Yangyun	yyzou@fudan.edu.cn

=head1 Analysis steps
	
	driver sites prediction steps as follows: 
	
		1	choose specific version of Ensembl reference  #ICGC(75_GRCh37.p13)/TCGA(80_GRCh38.p2)/preprocessed data
		2	collect mutation data from ICGC(.tsv)/TCGA(.maf)/preprocessed data
		3	data preprocessing for CANDRIS input
		4	perform CANDRIS
		5	output files in your path 
	

=head1 usage

	Options:
		-r input ensembl reference file 
		-t threshold of posterior probability
		-i input mutation file 
		-s whether to calculate the value of omega0 and omega1 of driver gene ("T" or "F"). Omega0 is the selection pressure of passenger sites, omega1 is the selection pressure of driver sites.
		-f input file including Cn, Cs value for each gene, if -w is chosen 'True'
		-o1 output file 1, mutation gene list 
		-o2 output file 2, mutation site list
		-h help

=head1 sample

	perl candris.pl -r reference -t threshold -i ICGC(.tsv)/TCGA(.maf)/preprocessed data -s T -f Cn,Cs_value_for_each_gene -o1 gene_with_driver_mutation_sites -o2 driver_mutation_sites 

=head1 done 

=cut

my($ref,$in,$out1,$out2,$version,$threshold,$select,$cncs);
GetOptions(
	"r:s"=>\$ref,
	"t:s"=>\$threshold,
	"i:s"=>\$in,
	"s:s"=>\$select,
	"f:s"=>\$cncs,
	"o1:s"=>\$out1,
	"o2:s"=>\$out2,
);
die `pod2text $0` unless($in&&$ref);
###########################################################################################

my (%nonsynonymous_site,%nonsynonymous_site_with_mutation,%cn,%cs,%cncs,%pvalue,%qvalue);
if ($select eq "T"){
	open CNCS, "<$cncs" or die "Can not open file $cncs:$!";
	my $head_cncs=<CNCS>;
	while (<CNCS>){
		chomp;
		my @line=split/\t/,$_;
		my $enst=(split/\|/,$line[0])[1];
		$nonsynonymous_site{$enst}=$line[1];
		$nonsynonymous_site_with_mutation{$enst}=$line[5];
		$cn{$enst}=$line[6];
		$cs{$enst}=$line[7];
		$cncs{$enst}=$line[8];
		$pvalue{$enst}=$line[9];
		$qvalue{$enst}=$line[10];
	}
}

## Read gene annotation
open REF, "<$ref" or die "Can not open file $ref:$!";
my $head_ref=<REF>;
my %gene_id;
my %aa_length;

while(<REF>){
	chomp;
	my @record=split /\t/, $_;
	$aa_length{$record[1]}=$record[6];
	$gene_id{$record[1]}=$record[0]."\t".$record[1]."\t".$record[3]."\t".$record[4];
}
close REF;


open IN, "<$in" or die "Can not open file $in:$!";
open OUT1, ">$out1" or die "Can not open file $out1:$!";
open OUT2, ">$out2" or die "Can not open file $out2:$!";
if ($select eq "T"){
	print OUT1 "Gene_id\tTranscript_id\tChr\tGene_name\tMax_hit\tCount\tProtein_length\tDistribution\tmut_site\tf0\tmean\tvariance\tm0\tm1\teta\tLRT_p-value\tQ(z)\tzc\tnumber_of_predited\tCnCs0\tCnCs1\tCnCs\tchisq_p-value\tq-value\n";
}
else{
	print OUT1 "Gene_id\tTranscript_id\tChr\tGene_name\tMax_hit\tCount\tProtein_length\tDistribution\tmut_site\tf0\tmean\tvariance\tm0\tm1\teta\tLRT_p-value\tQ(z)\tzc\tnumber_of_predited\n";
}
print OUT2 "Gene_id\tTranscript_id\tChr\tGene_name\tProtein_position\tz\tQ(z)\tProtein_mutation\n";

my %aa_mut_count;
my $head=<IN>;

if ($head=~/^icgc_mutation_id/){
	my @head=split/\t/,$head;
	my ($mut_type,$aa_mut,$ensg,$enst);
	foreach my $i (0 .. $#head){
		if ($head[$i]=~/\bconsequence_type\b/){
			$mut_type=$i;
		}elsif($head[$i]=~/\baa_mutation\b/){
			$aa_mut=$i;
		}elsif($head[$i]=~/\bgene_affected\b/){
			$ensg=$i;
		}elsif($head[$i]=~/\btranscript_affected\b/){
			$enst=$i;
		}
	}
	while(<IN>){
		chomp;
		my @lines=split /\t/,$_;
		if ($lines[$mut_type] =~ /missense_variant/){
			my $id=$lines[$ensg]."\t".$lines[$enst]."\t".$lines[$aa_mut];
			$aa_mut_count{$id}+=1;
		}
	}
}elsif($head=~/^#version gdc/){
	my ($mut_type,$aa_mut,$ensg,$enst);
	while(<IN>){
		chomp;
		my $line=$_;
		my @lines=split /\t/,$_;
		if($line=~/^#/){
			next;
		}elsif($line=~/^Hugo_Symbol/){
			my @title=split/\t/,$line;
			foreach my $i (0 .. $#title){
				if ($title[$i]=~/\bConsequence\b/){
					$mut_type=$i;
				}elsif($title[$i]=~/\bHGVSp_Short\b/){
					$aa_mut=$i;
				}elsif($title[$i]=~/\bGene\b/){
					$ensg=$i;
				}elsif($title[$i]=~/\bFeature\b/){
					$enst=$i;
				}
			}
		}else{
			if ($lines[$mut_type] =~ /missense_variant/){
				my $pro_mut=(split/\./,$lines[$aa_mut])[1];
				my $id=$lines[$ensg]."\t".$lines[$enst]."\t".$pro_mut;
				$aa_mut_count{$id}+=1;
			}
		}
	}
}elsif($head=~/^gene_id/){
	while(<IN>){
		chomp;
		my @lines =split /\t/,$_;
		my $id=$lines[0]."\t".$lines[1]."\t".$lines[3];
		$aa_mut_count{$id}=$lines[4];
	}
}
close IN;

my %count;
my @site_mut;
my %site_mut_list;

#Read mutation records
foreach my $key(keys %aa_mut_count){
	my @mut = split/\t/,$key;
	my $gene=$mut[1];
	my $aa_mut=$mut[2];
	my $pos;
	if (exists $gene_id{$gene}){	
		if($aa_mut=~/[a-zA-Z]+(\d+)[a-zA-Z]+/){
			$pos=$1;
		}
		my $count=$aa_mut_count{$key};
		my $site_mut=$gene.":".$pos; if(!$pos){print $key;}
		$site_mut_list{$site_mut}.=$aa_mut.":".$count.";";

		if(!exists $count{$site_mut}){
			$count{$site_mut}=$count;
			push @site_mut, $site_mut;
		}else{
			$count{$site_mut}+=$count;
		}
	}
}

## sort site mutation
foreach my $site_mut(@site_mut){
	my %mut_count;
	my $aa_mut;
	my @mutation=split /;/, $site_mut_list{$site_mut};
	foreach my $mutation(@mutation){
		$aa_mut=(split /:/, $mutation)[0];
		if(!exists $mut_count{$aa_mut}){
			$mut_count{$aa_mut}=(split /:/, $mutation)[1];
		}else{
			$mut_count{$aa_mut}+=(split /:/, $mutation)[1];
		}
	}
	$site_mut_list{$site_mut}="";
	foreach $aa_mut (sort {$mut_count{$b}<=>$mut_count{$a}} keys %mut_count){
		$site_mut_list{$site_mut}.=$aa_mut.":".$mut_count{$aa_mut}.";";
	}
}

my %max;
my %hit;
my @gene_list;
## mutation profile for each gene
foreach my $site_mut(@site_mut){
	my $gene=(split /:/, $site_mut)[0];
	my $pos=(split /:/, $site_mut)[1];
	if(!exists $count{$gene}){
		$count{$gene}=$count{$site_mut};
		$max{$gene}=$count{$site_mut};
		$hit{$gene}=$count{$site_mut}.":".$pos."\;";
		push @gene_list, $gene;
	}else{
		$count{$gene}+=$count{$site_mut};
		$hit{$gene}.=$count{$site_mut}.":".$pos."\;";
		$max{$gene}=$count{$site_mut}, if $count{$site_mut}>$max{$gene};
	}
}

##sort mutation distribution
foreach my $gene(@gene_list){
	my %site_hit;
	my $pos;
	my @site_hit=split /;/, $hit{$gene};
	foreach my $site_hit(@site_hit){
		$pos=(split /:/, $site_hit)[1];
		$site_hit{$pos}=(split /:/, $site_hit)[0];
	}
	$hit{$gene}="";
	foreach $pos (sort {$site_hit{$b}<=>$site_hit{$a} or $a<=>$b} keys %site_hit){
		$hit{$gene}.=$site_hit{$pos}.":".$pos.";";
	}
}

## calculate paramates
foreach my $gene(@gene_list){
	my $gene_id=$gene_id{$gene};
	my $count=$count{$gene};
	my $count_trim=$count;
	my $aa=$aa_length{$gene};
	my $maxhit=$max{$gene};
	my $dis=$hit{$gene};
	my @dis=split /;/, $dis;
	my $aa_trim=$aa;
	for my $i(0 .. $#dis){
		my $dis1=(split /:/,$dis[$i])[0];
		if($dis1>100){
			$count_trim-=$dis1;
			$aa_trim-=1;
		}
	}
	my $hit;
	my @q;
	if($maxhit>2 && $aa=~/\d+/){
		my $mean=$count_trim/$aa_trim;
		my $num=$#dis+1;
		my $nohit=$aa-$num;
		my $f0=$nohit/$aa;
		my $m0=-log($f0);
		my $var=$nohit*$mean*$mean;
		for my $i(0 .. $#dis){
			my $dis1=(split /:/,$dis[$i])[0];
			$var+=($dis1-$mean)*($dis1-$mean), if $dis1<100;
		}
		$var=$var/($aa_trim-1);
		my $m1=$mean+($var-$mean)/($mean-$m0);
		my $eta=($mean-$m0)/($m1-$m0);
		if($m1<0||$eta<0||$eta>1){
			print OUT1 "$gene_id\t$maxhit\t$count\t$aa\t$dis\t$num\t$f0\t$mean\t$var\t$m0\t$m1\t$eta\n";
			next;
		}
		my $zc=(($m1-$m0)-log($eta/(1-$eta))+log($threshold/(1-$threshold)))/log($m1/$m0);
		my $qz;
		my $num_pred=0;
		my $LRT= &LRT($maxhit,$eta,$m1,$m0,$mean);
		my $p=Statistics::Distributions::chisqrprob(1,$LRT);
		foreach my $i(0 .. $#dis){
			my $dis1=(split /:/,$dis[$i])[0];
			my $dis0=(split /:/,$dis[$i-1])[0];
			my $pos=(split /:/,$dis[$i])[1];
			my $ratio=log($eta/(1-$eta))+$dis1*log($m1/$m0)-($m1-$m0);
			$q[$i]=1-1/(1+exp($ratio));
			if($q[$i]>$threshold){
				$num_pred++;
				my $site_mut=$gene.":".$pos;
				print OUT2 "$gene_id\t$pos\t$dis1\t$q[$i]\t$site_mut_list{$site_mut}\n" if ($p<0.05);
			}
			if($#dis==0){  ##the distribution of driver sites in each gene
				$hit=1;
				$qz.="$q[$i]:$dis1:$hit;";
			}elsif($i==0){
				$hit=1;
			}elsif($i<$#dis && $dis1==$dis0){
				$hit++;
			}elsif($i<$#dis && $dis1!=$dis0){
				$qz.="$q[$i-1]:$dis0:$hit;";
				$hit=1;
			}elsif($i==$#dis && $dis1==$dis0){
				$hit++;
				$qz.="$q[$i]:$dis1:$hit;";
			}elsif($i==$#dis && $dis1!=$dis0){
				$qz.="$q[$i-1]:$dis0:$hit;";
				$hit=1;
				$qz.="$q[$i]:$dis1:$hit;";
			}
		}
		my $q0=1-1/(1+($eta/(1-$eta))*exp(-($m1-$m0)));
		$qz.=$q0.":0:".$nohit.";";
		
		if ($select eq "T"){
			my $cn=$cn{$gene};
			my $cn0=-log(1-$nonsynonymous_site_with_mutation{$gene}/$nonsynonymous_site{$gene});
			my $cs=$cs{$gene};
			my $cn1=($cn-(1-$eta)*$cn0)/$eta;
			my $cncs0=$cn0/$cs;
			my $cncs1=$cn1/$cs;
			my $cncs=$cncs{$gene};
			my $pvalue=$pvalue{$gene};
			my $qvalue=$qvalue{$gene};
			print OUT1 "$gene_id\t$maxhit\t$count\t$aa\t$dis\t$num\t$f0\t$mean\t$var\t$m0\t$m1\t$eta\t$p\t$qz\t$zc\t$num_pred\t$cncs0\t$cncs1\t$cncs\t$pvalue\t$qvalue\n" if ($p<0.05);
			print OUT1 "$gene_id\t$maxhit\t$count\t$aa\t$dis\t$num\t$f0\t$mean\t$var\t$m0\t$m1\t$eta\t$p\n" if ($p>=0.05);	
		}
		else{
			print OUT1 "$gene_id\t$maxhit\t$count\t$aa\t$dis\t$num\t$f0\t$mean\t$var\t$m0\t$m1\t$eta\t$p\t$qz\t$zc\t$num_pred\n";	
		}
	}
}
close OUT1;
close OUT2;
exit;

#calculate likelihood ratio test value
sub LRT{
	my ($max_hit,$eta,$m1,$m0,$mean)=@_;
	$max_hit=50 if $max_hit>50;
	my $fun_single=1;
	my $fun_two=1;
	foreach my $i (0 .. $max_hit){
		$fun_single+=log(($mean**$i)/(exp($mean)*fac($i)));
		$fun_two+=log((1-$eta)*($m0**$i)/(exp($m0)*fac($i))+$eta*($m1**$i)/(exp($m1)*fac($i)));
	}
	return 2*($fun_two-$fun_single);
}

sub fac{
	my ($n)=@_;
	return 1 if $n==0;
	return fac($n-1)*$n;
}