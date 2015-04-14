#! /usr/bin/perl

use warnings;
use strict;
use Bio::SeqIO;
use File::Temp qw(tempdir);

my $dir = tempdir( CLEANUP => 1);

if (@ARGV != 3){
	die("Syntax: $0 [protein db] [cds db] [file]\n");
}

my $inseq = Bio::SeqIO->new(-file => $ARGV[2], -format => 'fasta');
my $protdb = $ARGV[0];
my $cdsdb = $ARGV[1];

my $name = $ARGV[2];
$name =~ s/\.f.+$//;
my $prot = $name . "_macse_prot.fas";
open(PROT, ">$prot");
my $orfout = $name . "_macse_orf.fas";
open(ORF, ">$orfout");

my $seqtmp = $dir . "/blasttmp.fa";
my $outxtmp = $dir . "/blastxoutmp.txt";
my $outntmp = $dir . "/blastnoutmp.txt";
my $newfile = $dir .  "/seq.fas";
my $aln = $dir . "/aln.fas";
my $newfile2 = $dir . "/seq2.fas";
	

while (my $ori = $inseq->next_seq){
		
	# Print sequence to file
	open (TMP, ">$seqtmp");
	print TMP ">" . $ori->display_id . "\n";
	print TMP $ori->seq . "\n";
	close(TMP);
	
	my $blast = "blastx -db $protdb -query $seqtmp -out $outxtmp -outfmt '6 std qframe qseq sseq' -max_target_seqs 1";
	system($blast);

	
	# Process results
	my @res = process_tab($outxtmp);
		
	my $sta = $res[0];
	my $end = $res[1];
	my $fail = $res[3];
	
	if ($end == 0 || $fail > 0){
		next;
	}
	
	getorthos($ori,$seqtmp,$outntmp,$newfile);
	
	my $seqcount = `grep -c ">" $newfile`;
	if ($seqcount < 3){
		next;
	}
	
	my $ori2 = $ori->trunc($sta,$end);
	open(OUT, ">>$newfile");
	print OUT ">" . $ori2->display_id . "\n";
	print OUT $ori2->seq . "\n";
	close(OUT);
	
	system("java -jar -Xmx2000m /DataStore/Progs/macse/jar_file/macse.jar -prog alignSequences -seq $newfile -out_NT $aln > STDERR");
	
	if (! -e $aln || -z $aln){
		print "MACSE failed to work on " . $ori->display_id . "\n";
		next;
	}

	my $reg;
	my $gogo = 0;
	open(IN, $aln);
	while (<IN>){
		chomp; 
		$_ =~ s/!/N/g;
		my $id = $ori2->display_id;
		if (/^>\Q$id\E/){
			$gogo++;
		}elsif ($gogo > 0){
			$reg = $_;
			$gogo = 0;
		}
	}
	close(IN);
	$reg =~ s/\-//g;
		
	my $orf = $reg;
	$orf =~ s/\s+//g;
	$orf =~ s/-//g;
				
	print ORF ">" . $ori->display_id;
	print ORF " ORF\n";
	print ORF $orf . "\n";
	
	my $seqobj2 = Bio::Seq->new(-seq => $orf);
	$seqobj2 = $seqobj2->translate;
	
	print PROT ">" . $ori->display_id;
	print PROT " PROT\n";
	print PROT $seqobj2->seq . "\n";
}

sub process_tab { # Process BLASTx results in tabular format
	my $outxtmp = $_[0];
	
	my $hid;
	my $sta = 999999999999999;
	my $end = 0;
	my $hsp = 0;
	my $last;
	my $fail = 0;
	my $neg = 0;
	
	open (OUT, "$outxtmp");
	while (my $line = <OUT>){
		chomp $line;
		my @data = split /\t/, $line;
		if ($data[10] > 1e-03){
			last;
		}
		#print $data[10] . "\n";
		$hsp++;
		$hid = $data[1];
		my $qstr = $data[12];
		#print $qstr . "\n";
		if ($hsp == 1){
			$last = $data[12];
		}else{
			#If hit has hsps in different strands, skip this sequence
			if (($last > 0 && $data[12] < 0) || ($last < 0 && $data[12] > 0)){
				$fail++;
			}
		}
		if ($data[12] < 0){
			#print "revcom\n";
			$neg++;
		}
		my $s = $data[6];
		my $e = $data[7];
		if ($e < $s){
			$s = $data[7];
			$e = $data[6];
		}
		if ($s < $sta){
			$sta = $s;
		}
		if ($e > $end){
			$end = $e;
		}
	}
	close(OUT);
	
	return ($sta,$end,$neg,$fail,$hsp,$hid);
}

sub getorthos { # Run BLASTn and process results, pulling back homologous sequences
	my $ori = $_[0];
	my $seqtmp = $_[1];
	my $outntmp = $_[2];
	my $newfile = $_[3];
	
	my $blastn = "blastn -task blastn -db $cdsdb -query $seqtmp -out $outntmp -outfmt '6 std hspnum' -max_target_seqs 3";
	system($blastn);
	
	my @hits;
	open (IN, "$outntmp");
	while (my $line = <IN>){
		chomp $line;
		my @data = split /\t/, $line;
		if ($data[10] < 1e-10){
			$data[1] =~ s/^\w+\|//;
			$data[1] =~ s/\|$//;
			if (!grep/$data[1]/,@hits){
				push (@hits, $data[1]);
				my $getseq;
				if ($data[8] > $data[9]){
					$getseq = "blastdbcmd -entry $data[1] -strand minus -db $cdsdb >> $newfile";
				}else{
					$getseq = "blastdbcmd -entry $data[1] -strand plus -db $cdsdb >> $newfile";
				}
				system($getseq);
			}
		}
	}
	close(IN);
}
