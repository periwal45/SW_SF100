#Program to parse stacker-reader output files in .txt format and generate tab separated files by plate.
#!/usr/bin/perl
use strict;
use warnings;

opendir(DIR, "Data/") || die $!;
my @dir = readdir(DIR);

foreach my $f1(@dir){

	next unless $f1 =~ /(^NT\d+)$/;	#to read all bug folders starting with NT

	my $name1 = $1;

	opendir(NT, "Data/$f1") || die $!;
	
	my @contents = readdir(NT);

	foreach my $f2(@contents){

		next unless $f2 =~ /(^Rep\d).*/;	#to read all Replicate folders starting with Replicate

		my $name2 = $1;

		opendir(REP, "Data/$f1/$f2") || die $!;
		my @rep = readdir(REP);

		foreach my $f3(@rep){

		next unless $f3 =~ /(^plate.*)_NT.*.txt/;

		my $name3 = $1;

		my $file = "Data/$f1/$f2/$f3";
		open my $fh, $file or die $!;
	
		my $outfile = "$name1"."_"."$name2"."_"."$name3".".tab";
		open(OUT, ">Data/$f1/$outfile");

		print "$outfile\n";

		open(FF, ">>file_summary.txt");

		while(my $line = <$fh>){

			# #chomp $line;
			# if($line =~ /^Experiment\sFile\sPath:/){

			# 	print FF "$line\t";
			
			# }elsif($line =~ /^Plate\sNumber/){

			# 	print FF "$line\n";
			# }

			if ($line =~ /^Time\tT/ .. $line =~ /^Results/) {
        		
        		$line =~ s/,/\./g;

        		if($line =~ /^Results/){

        			next;
        		
        		}elsif($line =~ /^0:00:00/){

        			next;
        		
        		}elsif($line =~ /^\s+/){

        			next;

        		}else{

        		print OUT "$line";

        		}
   			 } 
		}	
	}
}
	
}	

#Below code removes the Temperature column from all files and formats time column

opendir(DIR, "Data/") || die $!;
@dir = readdir(DIR);

foreach my $f1(@dir){

	next unless $f1 =~ /^NT\d+$/;	#to read all bug folders starting with NT

	opendir(NT, "Data/$f1") || die $!;
	my @contents = readdir(NT);

	foreach my $file(@contents){

	next unless $file =~ /(^NT.*\.tab)/;

	my $name = $1;

	#print "$name\n";

	my $ff  = `LC_ALL=C cut -f 1,3- Data/$f1/$file`; #default delimiter for cut is tab

	open(OUT, ">Data/$f1/$name") || die $!;
	print OUT $ff;
	close OUT;

	my @f1 = open(FF, "Data/$f1/$file") || die $!;
	@f1=<FF>;

	open(OUT, ">Data/$f1/$name") || die $!;

	for(my $i=0; $i<scalar(@f1); $i++){

		my $line = $f1[$i];
		#print $line;
		chomp $line;
		my $t=$i-1;							#required to set time point to 0 else reading starts from header which is index 0
		$line =~ s/\d+:\d\d:\d\d/$t/g;		#numbers time points starting from 0 onwards
		$line =~ s/,/\t/g;
		print OUT "$line\n";

}
}
}


close DIR;
close NT;
close REP;
close OUT;	
