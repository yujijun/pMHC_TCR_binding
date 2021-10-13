open fr, "<$ARGV[0]";

my $sequences = '';
while(<fr>){
	$_ =~ s/\s+//g;
	if($_ eq ""){
		next;
	}
	$sequences .= $_."\n";
}
$sequences =~ s/\n/\//g;
$sequences =~ s/\s+//g;
#print($sequences);

my $flag = 1;
my $seqs_num = 0;
if(substr($sequences,0,1) ne '>'){
	$flag = 0; #没有以“>”开头
	#print("没有以>开头");
}elsif(substr($sequences,1,1) eq '/'){
	$flag = 0; #序列标签为空
	#print("序列标签为空");
}else{
	my @seqs = split />/, $sequences;
	foreach(@seqs){
		if($_ eq ''){
			next;
		}
		my @seq_info = split /\//,$_;
		if($#seq_info == 0){
			$flag = 0; #序列为空
			#print("序列为空");
		}else{
			my $seq_label = $seq_info[0];
			my $seq_content = join '',@seq_info[1..$#seq_info];
			#print($seq_label."--".$seq_content."\n");
			$seqs_num ++;
			if($seq_content =~ /[^A-Za-z]/){
				$flag = 0; #序列内容不规范
				#print("序列内容不规范");
			}
		}
	}
}
print($flag);
