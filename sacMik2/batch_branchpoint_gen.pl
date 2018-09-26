#!/usr/bin/perl -w 

@files = qw(
mikatae_rapa00_ls077_sacMik2
mikatae_rapa60_ls078_sacMik2
);

#$jtail = ".ih1.rd1.xpe.star.jcn.jmn.trk.bed";
#$jtail = ".ih1.xpe.star.jcn.jmc.trk.bed";
#$jtail = ".ih1.rd1.xpe.star.jcn.jmc.trk.bed";
$jtail = ".ih1.rd1.xpe.star.jcn.all.trk.bed";


$usage = <<EOF;

USAGE:

./batch_branchpoint_gen.pl annotation.bed genome.2bit outputfile.tsv conservation

inputintron.bed contains bed format from data that has introns
annotation.bed contains the annotated genes in the genome
genome.2bit is the genome sequence in 2bit format
outputfile is a tab seperated value text file sorted by total bed score
conservationfile (bigwig format , optional)

You will need twoBitToFa in your path

(and hgBigWigSummary if you are using conservation)

EOF

$numargs = $#ARGV +1 ;
if ($numargs < 4){
print $usage;
exit;}

open(IF,"sc3intronssacMik2.bed");
while($a=<IF>){
chomp($a);$intronhash{$a}++;
}



$safe_filename_characters = "a-zA-Z0-9_.-/";
$annotationbed=$ARGV[0];
$genome = $ARGV[1];
$outfile = $ARGV[2];
if($ARGV[3]){$conswig = $ARGV[3];}
else{$conswig="";}

if($annotationbed =~ /[^$safe_filename_characters]/){print "unsafe sequence characters in annotationbed \n"; print "please use simple alpha numeric characters\n";exit;}
if($genome =~ /[^$safe_filename_characters]/){print "unsafe sequence characters in genome \n"; print "please use simple alpha numeric characters\n";exit;}
if($outfile =~ /[^$safe_filename_characters]/){print "unsafe sequence characters in outfile \n"; print "please use simple alpha numeric characters\n";exit;}
if($conswig =~ /[^$safe_filename_characters]/){print "unsafe sequence characters in conservation wig file \n"; print "please use simple alpha numeric characters\n";exit;}


open(IF,"$annotationbed");
while($a=<IF>){chomp($a);$annhash{$a}++;}
close(IF);

%cnthash = (); %jutotalhash = ();
open(OF,">$outfile");
print OF "Intron_name\tlink_to_browser\tchrom\tstart\tend\tstrand\tlength\t5pss_to_bp\tbp_to3pss\ttype\tus\tbp_seq\tds\texpanded 5pss_seq\t5pss_conservation\texpanded_branch_point\tbp conservation\texpanded 3pss\t3pss conservation\tgene\tannotated_intron_overlap\ttotal_score";

foreach $f (@files){print OF "\tscore_$f";}
print OF "\tintron_sequence\n";


foreach $f (@files){
$shortf = $f;
$shortf =~ s/bayanus_rapa//;
$shortf =~ s/_sacBay2//;
$shortf =~ s/mikatae_rapa//;
$shortf =~ s/_sacMik2//;
#print "$shortf\n";exit;
$t = $f;
$t =~ s/_sacBay2//;
$t =~ s/by_//;
open(IF,"$f$jtail");
$a=<IF>;
while($a=<IF>){
($chrom,$start,$end,$name,$score,$strand,$tstart,$tend,$colors,$exons,$sizelist,$startlist)=split(/\t/,$a);

if($chrom =~ /chr/){
if($exons>1){
@sizes = split(/,/,$sizelist);
@starts = split(/,/,$startlist);
for($i=1;$i<$exons;$i++){

$justart = $start + $starts[$i-1] + $sizes[$i-1];
$juend = $start  + $starts[$i];

$jupos = "$chrom\t$justart\t$juend\t$strand";
$jutotalhash{$jupos}+=$score;
$jueachhash{$f}{$jupos}+=$score;
$cnthash{$jupos}++;
$juname{$jupos} .= "$name\_$shortf,";
}
}
}
}
}

foreach $jupos (sort {$jutotalhash{$b} <=> $jutotalhash{$a}} keys %jutotalhash){
$gene="";
$notes="";

($chrom,$start,$end,$strand) = split(/\t/,$jupos);
$length = $end-$start;

foreach $intronline (keys %intronhash){
@b=split(/\t/,$intronline);
if($b[0] eq $chrom && $b[1]<$end && $b[2]>$start){

if($b[1]==$start && $b[2]==$end && $strand eq $b[5] && $notes !~ /exact_match_/ ){ $notes .= "exact_match_$b[3],";}
elsif($b[1]==$start &&  $strand eq $b[5] && $notes !~ /exact_match/ ){if($strand eq '+'){ $notes .= "exact_matchus_$b[3],";}else{ $notes .= "exact_matchds_$b[3],";}}
elsif($b[2]==$end &&  $strand eq $b[5] && $notes !~ /exact_match/ ){if($strand eq '+'){ $notes .= "exact_matchds_$b[3],";}else{ $notes .= "exact_matchus_$b[3],";}}
elsif($b[1]<$end && $b[2]>$start && $strand eq $b[5] && $notes !~ /exact_match/ && $notes !~ /overlap_intron_$b[3],/ ){ $notes .= "overlap_intron_$b[3],";$ihit++;}
elsif($b[1]<$end && $b[2]>$start && $strand ne $b[5] && $notes !~ /exact_match/ && $notes !~ /opposing_intron_$b[3],/ ){ $notes .= "opposing_intron_$b[3],";$ihit++;}
}
}




foreach $geneline (keys %annhash){
@b=split(/\t/,$geneline);
if($b[0] eq $chrom && $b[1]<$end && $b[2]>$start){
$tgene=$b[3];
if($gene!~ /$b[3],/ && $strand eq $b[5]){$gene .= "$tgene,";}
elsif($gene !~ /$b[3],/){$gene .= "opposing_$tgene";}

if($b[9]>1){
@sizes = split(/,/,$b[10]);
@starts = split(/,/,$b[11]);
$ihit=0;
#print "@sizes\n";
#print "@starts\n";
for($k=1;$k<$b[9];$k++){
$istart = $b[1] + $starts[$k-1] + $sizes[$k-1];
$iend = $b[1]  + $starts[$k];
#print "$juname{$jupos} $jupos \n"; print "if($istart==$start && $iend==$end && $strand eq $b[5] && $notes !~ /exact_match_/ ){ $notes .= exact_match_$tgene,;$ihit++;}\n";
#print "$chrom $istart $iend\n";
if($istart==$start && $iend==$end && $strand eq $b[5] && $notes !~ /exact_match_/ ){ $notes .= "exact_match_$tgene,";$ihit++;}
elsif($istart==$start &&  $strand eq $b[5] && $notes !~ /exact_match/ ){if($strand eq '+'){ $notes .= "exact_matchus_$tgene,";}else{ $notes .= "exact_matchds_$tgene,";}$ihit++;}
elsif($iend==$end &&  $strand eq $b[5] && $notes !~ /exact_match/ ){if($strand eq '+'){ $notes .= "exact_matchds_$tgene,";}else{ $notes .= "exact_matchus_$tgene,";}$ihit++;}
elsif($istart<$end && $iend>$start && $strand eq $b[5] && $notes !~ /exact_match/ && $notes !~ /overlap_intron_$tgene,/ ){ $notes .= "overlap_intron_$tgene,";$ihit++;}
elsif($istart<$end && $iend>$start && $strand ne $b[5] && $notes !~ /opposing_intron_$tgene,/ ){ $notes .= "opposing_intron_$tgene,";$ihit++;}
}
if($ihit==0){
if($strand eq $b[5] && $notes !~ /exact/ && $notes !~ /overlap_gene_$tgene,/ ){$notes .= "overlap_gene_$tgene,";}
elsif($strand ne $b[5] && $notes !~ /opposing_gene_$tgene,/) {$notes .= "opposing_gene_$tgene,";}
}
}
else{
#if($strand eq $b[5] && $notes !~ /overlap_gene_$tgene,/ ){$notes .= "overlap_gene_$tgene,";}
#elsif($notes !~ /opposing_gene_$tgene,/) {$notes .= "opposing_gene_$tgene,";}
}
}
}


if($strand eq '+'){

$l=$start;$r=$end;
system("twoBitToFa $genome:$chrom:$l-$r tmp.fa");
open(FA1 , "tmp.fa");
$fl=<FA1>;
$full = "";while($fl=<FA1>){chomp($fl);$full .= "$fl" ;}$full =~ tr/a-z/A-Z/;
close(FA1);


$l=$start;$r=$l+6;
system("twoBitToFa $genome:$chrom:$l-$r tmp.fa");
open(FA1 , "tmp.fa");
$fl=<FA1>;
$us = "";while($fl=<FA1>){chomp($fl);$us .= "$fl" ;}$us =~ tr/a-z/A-Z/;
close(FA1);

$l=$start-5;$r=$l+6+10;
system("twoBitToFa $genome:$chrom:$l-$r tmp.fa");
open(FA1 , "tmp.fa");
$fl=<FA1>;
$exus = "";while($fl=<FA1>){chomp($fl);$exus .= "$fl" ;}$exus =~ tr/a-z/A-Z/;
close(FA1);
$usl=$l; $usr=$r;


$r=$end;$l=$r-3;
system("twoBitToFa $genome:$chrom:$l-$r tmp.fa");
open(FA1 , "tmp.fa");
$fl=<FA1>;
$ds = "";while($fl=<FA1>){chomp($fl);$ds .= "$fl" ;}$ds =~ tr/a-z/A-Z/;
close(FA1);

$r=$end+5;$l=$r-3-10;
system("twoBitToFa $genome:$chrom:$l-$r tmp.fa");
open(FA1 , "tmp.fa");
$fl=<FA1>;
$exds = "";while($fl=<FA1>){chomp($fl);$exds .= "$fl" ;}$exds =~ tr/a-z/A-Z/;
close(FA1);
$dsl=$l; $dsr=$r;

if($length>35){ $l=$start+35; $r=$end;
system("twoBitToFa $genome:$chrom:$l-$r tmp.fa");
open(FA1 , "tmp.fa");
$fl=<FA1>;
$branch = "";while($fl=<FA1>){chomp($fl);$branch .= "$fl" ;}$branch =~ tr/a-z/A-Z/;
close(FA1);
}
else{$branch="short";}

}
else{

$l=$start;$r=$end;
system("twoBitToFa $genome:$chrom:$l-$r tmp.fa");
open(FA1 , "tmp.fa");
$fl=<FA1>;
$full = "";while($fl=<FA1>){chomp($fl);$full .= "$fl" ;}$full =~ tr/a-z/A-Z/;
close(FA1);
$tmp=reverse($full);$tmp=~ tr/ACTG/TGAC/;$full=$tmp;

$l=$start;$r=$l+3;
system("twoBitToFa $genome:$chrom:$l-$r tmp.fa");
open(FA1 , "tmp.fa");
$fl=<FA1>;
$ds = "";while($fl=<FA1>){chomp($fl);$ds .= "$fl" ;}$ds =~ tr/a-z/A-Z/;
close(FA1);
$tmp=reverse($ds);$tmp=~ tr/ACTG/TGAC/;$ds=$tmp;

$l=$start-5;$r=$l+3+10;
system("twoBitToFa $genome:$chrom:$l-$r tmp.fa");
open(FA1 , "tmp.fa");
$fl=<FA1>;
$exds = "";while($fl=<FA1>){chomp($fl);$exds .= "$fl" ;}$exds =~ tr/a-z/A-Z/;
close(FA1);
$tmp=reverse($exds);$tmp=~ tr/ACTG/TGAC/;$exds=$tmp;
$dsl=$l; $dsr=$r;

$r=$end;$l=$r-6;
system("twoBitToFa $genome:$chrom:$l-$r tmp.fa");
open(FA1 , "tmp.fa");
$fl=<FA1>;
$us = "";while($fl=<FA1>){chomp($fl);$us .= "$fl" ;}$us =~ tr/a-z/A-Z/;
close(FA1);
$tmp=reverse($us);$tmp=~ tr/ACTG/TGAC/;$us=$tmp;

$r=$end+5;$l=$r-6-10;
system("twoBitToFa $genome:$chrom:$l-$r tmp.fa");
open(FA1 , "tmp.fa");
$fl=<FA1>;
$exus = "";while($fl=<FA1>){chomp($fl);$exus .= "$fl" ;}$exus =~ tr/a-z/A-Z/;
close(FA1);
$tmp=reverse($exus);$tmp=~ tr/ACTG/TGAC/;$exus=$tmp;
$usl=$l; $usr=$r;

if($length>35){ $l=$start; $r=$end-35;
system("twoBitToFa $genome:$chrom:$l-$r tmp.fa");
open(FA1 , "tmp.fa");
$fl=<FA1>;
$branch = "";while($fl=<FA1>){chomp($fl);$branch .= "$fl" ;}$branch =~ tr/a-z/A-Z/;
close(FA1);
}
else{$branch="short";}
$tmp=reverse($branch);$tmp=~ tr/ACTG/TGAC/;$branch=$tmp;

} 

$motif = substr($us,0,2) . "-". substr($ds,-2);

if($conswig ne ''){

system("hgWiggle -db=sacMik2 -position=$chrom:$usl-$usr -doStats $conswig > tmp.stats ");
open(ST,"tmp.stats");
$aa=<ST>;$aa=<ST>;
$aa=<ST>;$aa=<ST>;
if($aa){@bb= split(/\t/,$aa);
$meanus = $bb[9];
}else{$meanus=0;}
#system("hgWiggle -position=$chrom:$usl-$usr -doStats $conswig > tmp.stats ");
#system("bigWigSummary $conswig $chrom $usl $usr 1 > tmp.stats ");open(ST,"tmp.stats");$aa=<ST>;chomp($aa);if($aa){$meanus = $aa;}else{$meanus=0;}


system("hgWiggle -db=sacMik2 -position=$chrom:$dsl-$dsr -doStats $conswig > tmp.stats");
open(ST,"tmp.stats");
$aa=<ST>;$aa=<ST>;
$aa=<ST>;$aa=<ST>;
if($aa){@bb= split(/\t/,$aa);
$meands = $bb[9];
}else{$meands=0;}
#system("hgWiggle -position=$chrom:$dsl-$dsr -doStats $conswig > tmp.stats");
#system("bigWigSummary $conswig $chrom $dsl $dsr 1 > tmp.stats ");open(ST,"tmp.stats");$aa=<ST>;chomp($aa);if($aa){$meands = $aa;}else{$meands=0;}
}
else{$meanus=0;$meands=0;}

#$lbranch = length($branch);

#y = [tc]
#p = [ga]
#b = [agt]

#bp = NACTAAN or NPYTPAY

$actaapos=0; $opos=0; 
$lactaac=-1;$ractaac=0;
while($branch =~ m/[ACTG]ACTAA[ACTG]/g){$lactaac = $-[0] ; $ractaac = $+[1]; }

$lmotif1=-1;$rmotif1=0;
while($branch =~ m/[ACTG][AG][TC]T[AG]A[TC][AG]/g){$lmotif1 = $-[0] ; $rmotif1 = $+[1]; }

$lmotif2=-1;$rmotif2=0;
while($branch =~ m/[ACTG][ACTG][TC]T[AG]A[TC][ACTG]/g){$lmotif2 = $-[0] ; $rmotif2 = $+[1]; }

#motif1 = NPYTPAYP
#motif1 = NPYTPAYB

#P = A T ?
#8) ACTAAC or PYTPAY in intron at least 45bp from the 5'ss; distance from 3'ss in the field  --right now you have yes/no for ACTAAC
# 45 + $lactaac
# $length - (45 + $lactaac);
$actualbranch=""; 
$usethis=-100;
if($lactaac>=0){ $usethis = $lactaac;}
elsif($lmotif1>=0){$usethis = $lmotif1;}
elsif($lmotif2>=0){$usethis = $lmotif2;}
if($usethis != -100){
$actualbranch = substr($branch,$usethis,8);
$actaapos = $length - (45 + $usethis) + 4 ; $opos = $length - $actaapos; 

if($strand eq '+'){$l=$start+$opos-10; $r = $l + 20;}
else{$l=$end-$opos-10;$r = $l+20;}
system("twoBitToFa $genome:$chrom:$l-$r tmp.fa");
open(FA1 , "tmp.fa");
$fl=<FA1>;
$exbranch = "";while($fl=<FA1>){chomp($fl);$exbranch .= "$fl" ;}$exbranch =~ tr/a-z/A-Z/;
close(FA1);
if($strand eq '-'){$tmp=reverse($exbranch);$tmp=~ tr/ACTG/TGAC/;$exbranch=$tmp;}


if($conswig ne ''){
system("hgWiggle -db=sacMik2 -position=$chrom:$l-$r -doStats $conswig > tmp.stats");
open(ST,"tmp.stats");
$aa=<ST>;$aa=<ST>;
$aa=<ST>;$aa=<ST>;
if($aa){@bb= split(/\t/,$aa);
$meanbranch = $bb[9];
}else{$meanbranch=0;}

#system("hgWiggle -position=$chrom:$l-$r -doStats $conswig > tmp.stats");
#system("bigWigSummary $conswig $chrom $l $r 1 > tmp.stats ");
#open(ST,"tmp.stats");$aa=<ST>;chomp($aa);if($aa){$meanbranch = $aa;}else{$meanbranch=0;}

}
else{$meanbranch=0;}
}
else{$actualbranch="";$exbranch="";$meanbranch=0;}

$l=$start-4;$r=$end+4;
$gene =~ s/,$//;
$notes =~ s/,$//;
print  OF "$juname{$jupos}\t<a href=\"http://intron.ucsc.edu/cgi-bin/hgTracks?db=sacMik2&position=$chrom:$l-$r\">link</a>\t$chrom\t$start\t$end\t$strand\t$length\t$opos\t$actaapos\t$motif\t$us\t$actualbranch\t$ds\t$exus\t$meanus\t$exbranch\t$meanbranch\t$exds\t$meands\t$gene\t$notes\t$jutotalhash{$jupos}";
foreach $f (@files){if(! exists $jueachhash{$f}{$jupos}){$jueachhash{$f}{$jupos}=0;} print OF "\t$jueachhash{$f}{$jupos}";}
print OF "\t$full\n";

$meanus=0;$meands=0;$meanbranch=0;
$exbranch="";$exus="";$exds="";
#print OF "$full\n";
#678910
#12 - end

}


