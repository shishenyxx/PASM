#!/usr/bin/perl
#developed by LiJR
#version=0.9 modified by YangXX, corrected by ZhengXN
#date=2016-07-14

use Statistics::R;
my $R = Statistics::R->new();
#	$R->run(q`.libPaths("/rd/home/yangxx/R/x86_64-pc-linux-gnu-library/3.1/")`);
	$R->run(q`suppressPackageStartupMessages(library(yyxMosaicHunter))`);
	$R->run(q`suppressPackageStartupMessages(library(pryr))`);
	$R->run(q`suppressPackageStartupMessages(library(stats))`);
 @base=("A","a","T","t","C","c","G","g","N","n");
 
# print "chr\tpos\tdepth\tref\talt\taf_ref\taf_alt\tins\tdel\tp_hom_ref\tp_het\tp_hom_alt\tp_mosaic\n";
 
while (<>){
	chomp;
	@items=split("\t",$_);
	$chr=$items[0];
	$pos=$items[1];
	$ref=uc($items[2]);
	$ref_lc=lc($ref);
	$depth=$items[3];
	$id=$items[6];
if($depth>10){
	@match=split("",$items[4]);
	@quality=split("",$items[5]);
	%pos_n=();
	$m=0;
	#$A=0;$T=0;$C=0;$G=0;$N=0;
	for $b(@base) {$count{$b}=0;}
	$m_l=0;
	#$X=0;$Y=0;$Z=0;
	$I=0;$D=0;
	$o=0;$a=0;
	$star=0;
	$in_base=0;$del_base=0;$mis_base=0;

	$l=@match;
	$i=0;
        $pos_q=0;
	%end_reads=();
	$star=0;
#	@base=qw("A","a","T","t","C","c","G","g");
#	print "$l\n";
	while($i<$l){
#	if($i==1000){print "$depth\t$match[$i]\t$a\t$m\t$A\t$T\t$C\t$G\t$I\t$D\t$o\n"};
		if($match[$i] eq "^" ){
			$i+=2;
		}elsif($match[$i] eq "\$"){
			$i++;
		}elsif($match[$i] eq "."){
			$count{$ref} ++;
			push @{$pos_n{$ref}},$i;
                        $pos_q++;
			$i++;
		}elsif($match[$i] eq ','){
			$count{$ref_lc} ++;
			push @{$pos_n{$ref_lc}},$i;
                        $pos_q++;
			$i++;
		}elsif($match[$i] ~~ @base){
			$count{$match[$i]} ++;
			push @{$pos_n{$match[$i]}},$i;
                        $pos_q++;
#			print $match[$i]."\t".$count{$match[$i]}."\n";
			$i++;
		}elsif($match[$i] eq "+"){
			$I++;
			if ($match[($i+2)] =~/\d/){
				$c=10*$match[($i+1)]+$match[($i+2)];
				$i+=$c+3;
			}else{
				$c=$match[($i+1)];
                                $i+=$c+2;
			}
			 $in_base+=$c;
		}elsif($match[$i] eq "-"){
			$D++;
			if ($match[($i+2)] =~/\d/){
				$c=10*$match[($i+1)]+$match[($i+2)];
				$i+=$c+3;
			}else{
                                $c=$match[($i+1)];
				$i+=$c+2;
			}
			$del_base+=$c;
		}elsif($match[$i] eq "*"){
		        $D++;
			$i++;
		}else{
			print "ERROR: $i\t$match[$i]\n";
			exit(0);
		}
	}
#	$a=$m+$A+$T+$C+$G+$N+$star;
#	$mis_base+=$A+$T+$C+$G;
	
	for $n("A","T","C","G"){
		$count{$n}=0 unless defined $count{$n};
		$count{lc($n)}=0 unless defined $count{lc($n)};
		$all{$n}=$count{$n}+$count{lc($n)};
	}


	@key=sort {$all{$b}<=>$all{$a}} keys %all;
#	if($depth>0){
#	$ref=$key[0];$alt1=$key[1];
		if($key[0] eq $ref){
			$alt1=$key[1];
		}else{
			$alt1=$key[0];
		}
		$ref_lc=lc($ref);
		$alt1_lc=lc($alt1);
		
#		$af_ref=$all{$ref}/$depth;
	#	$af_alt=$all{$alt}/$depth;
		$quality_ref="";$quality_alt="";
		@quality_score_ref=();@quality_score_alt=();

		for $p(@{$pos_n{$ref}},@{$pos_n{$ref_lc}}){
			$quality_ref.=$quality[$p];
			$quality_score=ord($quality[$p])-33;
			push @quality_score_ref,$quality_score;
		}
		for $p(@{$pos_n{$alt1}},@{$pos_n{$alt1_lc}}){
			$quality_alt.=$quality[$p];
			$quality_score=ord($quality[$p])-33;
			push @quality_score_alt,$quality_score;
		}
		
		if ($quality_alt eq ""){
		$quality_alt=".";
		}
	$quality_ref =~ s/!//g;
                $quality_alt =~ s/!//g;
                $quality_ref =~ s/\\//g;
                $quality_alt =~ s/\\//g;

#	print $chr."_"."$pos\t$quality_ref\t$quality_alt\n";	
#	print "$chr:$pos\t$quality_ref\t$quality_alt\n";



#Re-output the interval based on depth
	$mle_calculate_only=$all{$alt1}/($all{$ref}+$all{$alt1});
	$mle_interval=1.5/sqrt($depth);
	$mle_lw=$mle_calculate_only-$mle_interval;
	$mle_up=$mle_calculate_only+$mle_interval;
	if ($mle_lw<0){
	$mle_lw=0;
	}
	if ($mle_up>1){
	$mle_up=1;
	}
#	print "$mle_calculate_only\t$mle_interval\t$mle_lw\t$mle_up\n";

    


#	if($all{$ref} >3  and $all{$alt} >3 ){
		$R->set('quality_ref_in',$quality_ref);
		$R->set('quality_alt_in',$quality_alt);
		$R->run(q`quality_ref <- as.character(quality_ref_in)`);
		$R->run(q`quality_alt <- as.character(quality_alt_in)`);
		
#Transfer the MLE interval to R		
		$R->set('int_region_mle_lw',$mle_lw);
		$R->set('int_region_mle_up',$mle_up);
		
#Transfer ref+、ref-、alt+、alt- to R
		$R->set('ref_forward',$count{$ref});
		$R->set('ref_reverse',$count{lc($ref)});
		$R->set('alt_forward',$count{$alt1});
		$R->set('alt_reverse',$count{lc($alt1)});
#Estimate AF		
		$R->run(q`mh_result <- yyx_wrapped_mosaic_hunter_for_one_site(quality_ref, quality_alt)`);
#		$R->run(q`curve(mh_result$likelihood_fun(x), 0, 1, n=1001L)`);
		$R->run(q`cred_int <- yyx_get_credible_interval(mh_result$likelihood_fun, c(int_region_mle_lw,int_region_mle_up), 0.95)`);
#Fisher's exact test in R		
		$R->run(q`p_fisher<-fisher.test(matrix(c(ref_forward,ref_reverse,alt_forward,alt_reverse),nrow=2))$p.value`);
#		$MH=$R->get('mh_result$ref_het_alt_mosaic_posterior');
#		print join ("\t",@{$MH});
	

#Transfer back and output		
#		($p_hom_ref,$p_het,$p_hom_alt,$p_mosaic)=@{$mh_result};
		print "$chr\t$pos\t$depth\t$ref\t$alt1\t$I\t$D\t$count{A}\t$count{C}\t$count{G}\t$count{T}\t$count{a}\t$count{c}\t$count{g}\t$count{t}\t";
		$MLE=$R->get('cred_int$MLE');
		print "${MLE}\t";
		$CI=$R->get('cred_int$CI');
		print join ("\t",@{$CI});
		$P_FISHER=$R->get('p_fisher');
		print "\t${P_FISHER}";
		print "\n";
#	}
	}
}

	$R->stop;

