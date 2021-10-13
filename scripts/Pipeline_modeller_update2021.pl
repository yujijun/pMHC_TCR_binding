#!/usr/bin/perl
use Getopt::Std;
use vars qw($opt_h $opt_c $opt_a $opt_b $opt_p $opt_A $opt_B $opt_u $opt_n);
getopts('c:a:b:p:A:B:u:n:h');

sub help{
	print "Usage: Pipeline_modeller_update2021.pl -c <mhc class>\n";
	print "                                -a <mhc alpha chain>\n";
	print "                                -b <mhc beta chain>\n";
	print "                                -p <peptide>\n";
	print "                                [-A] <tcr alpha chain>\n";
	print "                                [-B] <tcr beta chain>\n";
	print "                                -u <userDir>\n";
	print "                                -n <model number>\n";
	print "DESCRIPTION：\n";
	print "      -c  the input MHC class type must be \"MHC-I\" or \"MHC-II\"!\n";
	print "      -a  input the sequence file (.fasta) or the allele name (eg:HLA-A*01:01, HLA-DPA1*01:03) of mhc alpha chain.\n";
	print "      -b  input the sequence file (.fasta) or the allele name (eg:b2m, HLA-DPB1*01:01) of mhc beta chain.\n";
	print "      -p   input the sequence file (.fasta) of peptide.\n";
	print "      [-A]  input the sequence file (.fasta) of tcr alpha chain.\n";
	print "      [-B]  input the sequence file (.fasta) of tcr beta chain.\n";
	print "      -u   input the output direction.\n";
	print "      -n   input the output model number.\n";
}

# 读取并返回文件中的第一条链
sub readFASTA{
	open fr_seqs, "<$_[0]";
	local $/ = '>';
	while(<fr_seqs>){
		chomp $_;
		my $check = $_;
		$check =~ s/\s+//g;
		if($check eq ''){
			next;
		}

		my @data = split /\n/,$_;
		return @data[1..$#data];
	}
	close(fr_seqs);
}

my $templates = '';
my $mhc_a_chain = '';
my $mhc_b_chain = '';
my $pep_chain = '';
my $tcr_a_chain = '';
my $tcr_b_chain = '';
my $userDirName = '';
my $model_num = '';
my $modeller_prog = "./program";

if($opt_h){
	help();
	exit;
}elsif(!$opt_c || !$opt_a || !$opt_b || !$opt_p || !$opt_u || !$opt_n){
	help();
	exit;
}elsif(($opt_c ne "MHC-I") && ($opt_c ne "MHC-II")){
	print "Warning: the input MHC class type must be \"MHC-I\" or \"MHC-II\"!\n";
	help();
	exit;
}elsif($opt_c eq "MHC-I"){
	
	$templates = "./Template_MHCI";
	if($opt_a =~ m/\.fasta$/){
		$mhc_a_chain = $opt_a;
	}else{
		open fr, "<$modeller_prog/getTopAlleles/Top_MHCI_alleles.csv";
		my %allele_list;
		while(<fr>){
			chomp $_;
			my @line_split = split /,/, $_;
			$allele_list{$line_split[0]} = $line_split[2];
		}
		close(fr);
		if(exists $allele_list{$opt_a}){
			#print $allele_list{$opt_a}."\n";
			$mhc_a_chain = $allele_list{$opt_a};
		}else{
			print "Warning: Please check the input mhc alpha chain, whether the input is a fasta sequence file (the suffix is \".fasta\"); or if the input is an allele of HLA class I, please check whether the input allele name is correct!\n";
			help();
			exit;
		}
	}

	if($opt_b =~ m/\.fasta$/){
		$mhc_b_chain = $opt_b;
	}elsif($opt_b eq 'b2m'){
		open fr, "<$modeller_prog/getTopAlleles/b2m.csv";
		my %allele_list;
		while(<fr>){
			chomp $_;
			my @line_split = split /,/, $_;
			$allele_list{$line_split[0]} = $line_split[2];
		}
		close(fr);
		#print $allele_list{$opt_b}."\n";
		$mhc_b_chain = $allele_list{$opt_b};
	}else{
		print "Warning: Please check the input mhc beta chain, whether the input is a fasta sequence file (the suffix is \".fasta\"); or you need to input \"b2m\" in -b option!\n";
		help();
		exit;
	}
}elsif($opt_c eq "MHC-II"){
	$templates = "./Template_MHCII";
	if($opt_a =~ m/\.fasta$/){
		$mhc_a_chain = $opt_a;
	}else{
		open fr, "<$modeller_prog/getTopAlleles/Top_MHCII_alleles.csv";
		my %allele_list;
		while(<fr>){
			chomp $_;
			my @line_split = split /,/, $_;
			$allele_list{$line_split[0]} = $line_split[2];
		}
		close(fr);
		if(exists $allele_list{$opt_a}){
			#print $allele_list{$opt_a}."\n";
			$mhc_a_chain = $allele_list{$opt_a};
		}else{
			print "Warning: Please check the input mhc alpha chain, whether the input is a fasta sequence file (the suffix is \".fasta\"); or if the input is an allele of HLA class II, please check whether the input allele name is correct!\n";
			help();
			exit;
		}
	}

	if($opt_b =~ m/\.fasta$/){
		$mhc_b_chain = $opt_b;
	}else{
		open fr, "<$modeller_prog/getTopAlleles/Top_MHCII_alleles.csv";
		my %allele_list;
		while(<fr>){
			chomp $_;
			my @line_split = split /,/, $_;
			$allele_list{$line_split[0]} = $line_split[2];
		}
		close(fr);
		if(exists $allele_list{$opt_b}){
			#print $allele_list{$opt_b}."\n";
			$mhc_b_chain = $allele_list{$opt_b};
		}else{
			print "Warning: Please check the input mhc beta chain, whether the input is a fasta sequence file (the suffix is \".fasta\"); or if the input is is an allele of HLA class II, please check whether the input allele name is correct!\n";
			help();
			exit;
		}
	}
}

if($opt_p =~ m/\.fasta$/){
	$pep_chain = $opt_p;
}else{
	print "Warning: Please check the input peptide chain, the input must be a fasta sequence file (the suffix is \".fasta\")!\n";
	help();
	exit;
}

if($opt_A){
	if($opt_A =~ m/\.fasta$/){
		 $tcr_a_chain = $opt_A;
	}else{
		print "Warning: Please check the input TCR alpha chain, the input must be a fasta sequence file (the suffix is \".fasta\")!\n";
		help();
		exit;
	}
}

if($opt_B){
	if($opt_B =~ m/\.fasta$/){
		 $tcr_b_chain = $opt_B;
	}else{
		print "Warning: Please check the input TCR beta chain, the input must be a fasta sequence file (the suffix is \".fasta\")!\n";
		help();
		exit;
	}
}

$userDirName = $opt_u;
$model_num = $opt_n;

mkdir "./model_build";
my $userDirPath = "./model_build/$userDirName";
mkdir $userDirPath;
my $query_seqs_file = "$userDirPath/query_seqs.fasta";
system "touch $query_seqs_file";

open fw, "> $query_seqs_file";
if($mhc_a_chain =~ m/\.fasta$/){
	my $mhc_a_chain_seq = readFASTA($mhc_a_chain);
	print fw ">mhc_alpha\n";
	print fw "$mhc_a_chain_seq\n";
}else{
	$opt_a =~ s/\*|\://g;
	print fw ">$opt_a\n";
	print fw "$mhc_a_chain\n";
}

if($mhc_b_chain =~ m/\.fasta$/){
	my $mhc_b_chain_seq = readFASTA($mhc_b_chain);
	print fw ">mhc_beta\n";
	print fw "$mhc_b_chain_seq\n";
}else{
	$opt_b =~ s/\*|\://g;
	print fw ">$opt_b\n";
	print fw "$mhc_b_chain\n";
}

my $pep_chain_seq = readFASTA($pep_chain);
print fw ">peptide\n";
print fw "$pep_chain_seq\n";

if($opt_A){
	my $tcr_a_chain_seq = readFASTA($tcr_a_chain);
	print fw ">tcr_alpha\n";
	print fw "$tcr_a_chain_seq\n";
}
if($opt_B){
	my $tcr_b_chain_seq = readFASTA($tcr_b_chain);
	print fw ">tcr_beta\n";
	print fw "$tcr_b_chain_seq\n";
}
close(fw);

my $sec; my $min; my $hour;
my $day; my $mon; my $year;
my $wday; my $yday; my $isdst;
   ($sec,$min,$hour,$day,$mon,$year,$wday,$yday,$isdst)=localtime(time());
   $year+=1900;
   $mon+=1;

my @tmp = split /\//,$userDirPath;
my $last_dir = $tmp[-1];
my $new_path = join("/",@tmp[0..$#tmp-1]);
system ("touch $new_path\/status.txt");
open (f, ">>$new_path\/status.txt");
print f "$year-$mon-$day $hour:$min:$sec\n";
close f;

open fw_config, ">$userDirPath\/config.txt";

# 将输入的序列转换为pir格式
open fr_seqs, "<$query_seqs_file";
local $/ = '>';
my @query_seqs;
while(<fr_seqs>){
	chomp $_;
	my $check = $_;
	$check =~ s/\s+//g;
	if($check eq ''){
		next;
	}
	
	my @data = split /\n/,$_;
	my $query_seq_name = $data[0];
	$query_seq_name =~ s/>//g;
	#$query_seq_name =~ s/^\s+|\s+$//g;
	$query_seq_name =~ s/\s+//g;
	
	open fw, ">$userDirPath\/$query_seq_name\.ali";
	print fw ">P1;$query_seq_name\n";
	print fw "HLA\:$query_seq_name\:\:\:\:\:\:\:0.00\:0.00\n";
	print fw @data[1..$#data]."*";
	close fw;

	open fw_fasta, ">$userDirPath\/$query_seq_name\.fasta";
	print fw_fasta ">$query_seq_name\n";
	print fw_fasta @data[1..$#data];
	close fw_fasta;

	push @query_seqs, "$userDirPath\/$query_seq_name\.fasta";
}

local $/ = "\n";
# 寻找最相似的模板,并匹配链的对应关系
my @template_set = glob("$templates/sequence/*");
my %template_indent;
foreach(@template_set)
{
	my @tmp = split /\//, $_;
	my $template_id = $tmp[-1];
	my @template_seqs = glob("$templates/sequence/$template_id/*");

	if($#template_seqs < $#query_seqs){
		next;
	}

	my @ident_matrix;
	open fw_matrix, ">$userDirPath\/matrix.txt";
	for($i=0; $i<=$#template_seqs; $i++){
		my @ident_row;
		for($j=0; $j<=$#query_seqs; $j++){
			my $ident = `$modeller_prog/NWalign/src/NWalign $query_seqs[$j] $template_seqs[$i]`;
			#system "echo $template_seqs[$i]\--$query_seqs[$j]--$ident >>logfile 2>&1";
			#print $template_seqs[$i]."--".$query_seqs[$j]."--".$ident."\n";
			push @ident_row, $ident;
		}
		push @ident_matrix, [@ident_row];
	}
	for(my $i=0; $i<=$#template_seqs; $i++){
		for(my $j=0; $j<=$#query_seqs; $j++){
			print fw_matrix -$ident_matrix[$i][$j]."  ";
		}
		print fw_matrix "\n";
	}
	
	system "$modeller_prog/Hungarian_nonstandard/Hungarian $userDirPath\/matrix.txt $userDirPath\/matrix_assign.txt";
	open fr_matrix_assign, "<$userDirPath\/matrix_assign.txt";
	my $ident_aver;
	my $num = 0;
	my $T_Q_list;
	while(<fr_matrix_assign>)
	{
		if($_ =~ /^SUM/)
		{
			$num ++;
			if( $num == 2){
				last;
			}
			my @tmp = split /\s+/, $_;
			$ident_aver = (-$tmp[-1])/($#query_seqs+1);
		}else{
			chomp $_;
			$T_Q_list .= "-".$_;
		}
	}
	$template_indent{$template_id.$T_Q_list} = $ident_aver;
}

=pob
foreach(sort{$template_indent{$b} <=> $template_indent{$a}} keys %template_indent)
{
	print $_."\t".$template_indent{$_}."\n";
}
=cut

my @keys = sort{$template_indent{$b} <=> $template_indent{$a}} keys %template_indent;
my @model_T= split /\-/, $keys[0];
my $model_T_id = $model_T[0];
my @T_Q_align;
for(my $i=1; $i<=$#model_T; $i++){
	my @tmp = split /\_/, $model_T[$i];
	# print $tmp[0]."--".$tmp[1]."\n";
	push @T_Q_align, [@tmp];
}

# 开始进行序列比对
system "cp $templates/structure/$model_T_id\.pdb $userDirPath";
my @model_T_seqs = glob("$templates/sequence/$model_T_id/*");
my @model_T_chain_ids;
my @query_seqs_name;
my @ali_files;
for(my $i=0; $i<=$#T_Q_align; $i++)
{
	my @tmp = split /\_|\./, $model_T_seqs[ $T_Q_align[$i][0] ];
	my $model_T_chain_id = $tmp[-2];
	push @model_T_chain_ids, $model_T_chain_id;

	my @tmp = split /\/|\./, $query_seqs[ $T_Q_align[$i][1] ];
	my $query_seq_name = $tmp[-2];
	push @query_seqs_name, $query_seq_name;

	system "salign.py -t $userDirPath\/$model_T_id\.pdb -b $model_T_chain_id -e $model_T_chain_id -q $userDirPath\/$query_seq_name\.ali >>$userDirPath\/modeller.log 2>&1";
	push @ali_files, "$userDirPath\/$query_seq_name\-$model_T_id\_$model_T_chain_id\.ali";
}

# 将模板链从PDB文件中提取出来，另存为一个PDB文件
my $model_T_chains_string = join("",@model_T_chain_ids);
system "$modeller_prog/GetPartChains/GetPartChains $userDirPath\/$model_T_id\.pdb $userDirPath $model_T_chains_string";
$model_T_id = $model_T_id."_".$model_T_chains_string;

# 获取模板蛋白中的链顺序
my $model_T_chains = `$modeller_prog/GetChains_id/GetChains_id $userDirPath\/$model_T_id\.pdb`;
my @model_T_chains = split //, $model_T_chains;

my @ali_files_sorted;
foreach(@model_T_chains){
	for(my $i=0; $i<=$#model_T_chain_ids; $i++){
		if($model_T_chain_ids[$i] eq $_){
			push @ali_files_sorted, $ali_files[$i];
			last;
		}
	}
}

@ali_files = @ali_files_sorted;

# 合并比对后的序列到一个文件中
local $/ = ">";
my @seqs_aligned;
my $model_T_start = 1;
my $ali_files_num = 0;
for(my $i=0; $i<=$#ali_files; $i++)
{
	$ali_files_num ++ ;
	open fr_ali_file, "<$ali_files[$i]";
	my $start=1;
	my %gap;
	my $num=1;
	my @seq_aligned;
	while(<fr_ali_file>)
	{
		chomp $_;
		if($_ eq '' || $_ eq "\n"){
			next;
		}
		my @data_split = split /\n/, $_;
		my $seq = join('',@data_split[2..$#data_split]);
		if($num == 1){
			my @label_split = split /\:/, $data_split[1];
			$start = $label_split[2];
			if($ali_files_num == 1){
				$model_T_start = $start;
			}
			while($seq =~ /(\-+)/g){
				my $gap_start = index($seq, $1);
				my $gap_length = length($1);
				$gap{$gap_start} = $gap_length;
				#substr($seq, $gap_start, $gap_length) = '+' x $gap_length;
			}
			#$seq =~ s/(\++)//g;
		}else{
			for my $key (keys %gap){
				#substr($seq, $key, $gap{$key}) = '+' x $gap{$key};
			}
			#$seq =~ s/(\++)//g;
		}
		push @seq_aligned, $seq;
		$num ++;

	}
	chop(@seq_aligned);
	push @seqs_aligned, [@seq_aligned];
}

local $/ = "\n";

open fw_align, ">$userDirPath\/all_seqs_aligned.ali";
print fw_align ">P1; $model_T_id;$template_indent{$keys[0]}\n";
print fw_align "structureX\:$model_T_id\.pdb\:$model_T_start\:$model_T_chain_ids[0]\:\:\:\:\:0.00\:0.00\n";
for(my $i=0; $i<=$#seqs_aligned; $i++)
{
	if($i<$#seqs_aligned){
		print fw_align "$seqs_aligned[$i][0]"."\/";
	}else{
		print fw_align "$seqs_aligned[$i][0]"."\*\n";
	}
}

my $query_label = join("-",@query_seqs_name);
print fw_align ">P1; $query_label\n".
			   "Sequence\:$query_label\:\:\:\:\:\:\:0.00\:0.00\n";
for(my $i=0; $i<=$#seqs_aligned; $i++)
{
	if($i<$#seqs_aligned){
		print fw_align "$seqs_aligned[$i][1]"."\/";
	}else{
		print fw_align "$seqs_aligned[$i][1]"."\*\n";
	}
}

close fw_align;


chdir "$userDirPath";
system "echo all_seqs_aligned.ali $model_num | model-multichain.py >>modeller.log 2>&1";

system("perl -pi -e 's/$year-$mon-$day $hour:$min:$sec\n//' ..\/status.txt");
# system("tar -zcf ..\/$last_dir.tar.gz .");


open fr_logfile, "<modeller.log";
my $flag=0;
while(<fr_logfile>){
	chomp $_;
	if($_ =~ /molpdf/){
		$flag=1;
		next;
	}
	if($flag==1){
		if($_ =~ /------/){
			next;
		}elsif($_ eq ''){
			next;
		}else{
			print fw_config $_."\n";
		}
	}
}
close fw_config;






