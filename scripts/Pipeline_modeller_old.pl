#!/usr/bin/perl
if($#ARGV != 1 )
{
	print "usage: Pipline_modeller.pl [templates set] [query seq set]\n";
	exit;
}

if(-e "./model_build"){
	system "rm -r ./model_build";
}
mkdir "./model_build";

my @query_seq_set = glob("./$ARGV[1]/*");
foreach(@query_seq_set)
{
	# 寻找最相似的模板
	my @templates = glob("./$ARGV[0]/*.fasta");
	my $query = $_;
	my %simis;
	foreach(@templates)
	{
		my @tmp = split /\/|\./, $_;
		my $template_id = $tmp[-2];
		$simis{$template_id} = `./NWalign/src/NWalign $query $_`;
	}

	my @keys = sort {$simis{$b} <=> $simis{$a}} keys %simis;
	my $simi_temp = $keys[0];
	print $simi_temp."\t".$simis{$simi_temp}."\n";
=pob
	for my $key (sort {$simis{$b} <=> $simis{$a}} keys %simis){
		print $key."\t".$simis{$key}."\n";
	}
=cut
	my @simi_temp_info = split /\_/, $simi_temp;
	my $simi_temp_id = $simi_temp_info[0];
	my $simi_temp_alpha = $simi_temp_info[1];
	my $simi_temp_beta = $simi_temp_info[2];

	my @query_info = split /\/|\./, $query;
	mkdir "./model_build/$query_info[-2]";
	my @query_chains = split /\-/, $query_info[-2];
	my $query_alpha = $query_chains[0];
	my $query_beta = $query_chains[1];

	# 分别比对alpha和beta链
	system "cp .\/HLA-DRB\/$simi_temp_id\.pdb ./";
	system "python3 salign.py -t $simi_temp_id\.pdb -b $simi_temp_alpha -e $simi_temp_alpha -q ./query_pir_alpha\/$query_alpha\.ali >out.log";
	system "python3 salign.py -t $simi_temp_id\.pdb -b $simi_temp_beta -e $simi_temp_beta -q ./query_pir_beta\/$query_beta\.ali >out.log";
	system "rm *.pap";

	open fr_alpha_align, "<$query_alpha\-$simi_temp_id\_$simi_temp_alpha\.ali";
	open fr_beta_align, "<$query_beta\-$simi_temp_id\_$simi_temp_beta\.ali";
	open fw_align, ">$query_alpha\_$query_beta\-$simi_temp_id\_$simi_temp_alpha$simi_temp_beta\.ali";

	local $/ = ">";
	my @alpha_align;
	my $alpha_start=1;
	my %alpha_gap;
	my $num=1;
	while(<fr_alpha_align>)
	{
		chomp $_;
		if($_ eq '' || $_ eq "\n"){
			next;
		}
		my @data_split = split /\n/, $_;
		my $alpha_seq = join('',@data_split[2..$#data_split]);
		if($num == 1){
			my @label_split = split /\:/, $data_split[1];
			$alpha_start = $label_split[2];	
			while($alpha_seq =~ /(\-+)/g){
				my $gap_start = index($alpha_seq, $1);
				my $gap_length = length($1);
				$alpha_gap{$gap_start} = $gap_length;
				substr($alpha_seq, $gap_start, $gap_length) = '+' x $gap_length;
			}
			$alpha_seq =~ s/(\++)//g;
		}else{
			for my $key (keys %alpha_gap){
				substr($alpha_seq, $key, $alpha_gap{$key}) = '+' x $alpha_gap{$key};
			}
			print $alpha_seq."\n";
			$alpha_seq =~ s/(\++)//g;
		}

		push @alpha_align, $alpha_seq;
		$num ++;

	}
	chop(@alpha_align);

	my @beta_align;
	my %beta_gap;
	my $num=1;
	while(<fr_beta_align>)
	{
		chomp $_;
		if($_ eq '' || $_ eq "\n"){
			next;
		}
		my @data_split = split /\n/, $_;
		my $beta_seq = join('',@data_split[2..$#data_split]);
		if($num == 1){
			while($beta_seq =~ /(\-+)/g){
				my $gap_start = index($beta_seq, $1);
				my $gap_length = length($1);
				$beta_gap{$gap_start} = $gap_length;
				substr($beta_seq, $gap_start, $gap_length) = '+' x $gap_length;
			}
			$beta_seq =~ s/(\++)//g;
		}else{
			for my $key (keys %beta_gap){
				substr($beta_seq, $key, $beta_gap{$key}) = '+' x $beta_gap{$key};
			}
			print $beta_seq."\n";
			$beta_seq =~ s/(\++)//g;
		}
		push @beta_align, $beta_seq;
		$num++;
	}
	chop(@beta_align);

	print fw_align ">P1;$simi_temp_id\_$simi_temp_alpha$simi_temp_beta;$simis{$simi_temp}\n".
			"structureX\:$simi_temp_id\.pdb\:$alpha_start\:$simi_temp_alpha\:\:$simi_temp_beta\:\:\:0.00\:0.00\n".
			"$alpha_align[0]"."\/"."$beta_align[0]\*\n";

	print fw_align ">P1;$query_alpha\_$query_beta\n".
			"Sequence\:$query_alpha\_$query_beta\:\:\:\:\:\:\:0.00\:0.00\n".
			"$alpha_align[1]"."\/"."$beta_align[1]\*\n";

	system "echo $query_alpha\_$query_beta\-$simi_temp_id\_$simi_temp_alpha$simi_temp_beta\.ali | python3 model-multichain.py >out.log";

	system "mv $query_alpha* $query_beta* $simi_temp_id* ./out.log ./model_build/$query_info[-2] ";
}











