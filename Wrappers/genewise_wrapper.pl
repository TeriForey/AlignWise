#! /usr/bin/perl

use warnings;
use strict;
use Bio::SeqIO;
use File::Temp qw(tempdir);

my $dir = tempdir( CLEANUP => 1);

if (@ARGV != 2){
	die("Syntax: $0 [protein db] [file]\n");
}

my $inseq = Bio::SeqIO->new(-file => $ARGV[1], -format => 'fasta');
my $protdb = $ARGV[0];

my $name = $ARGV[1];
$name =~ s/\.f.+$//;
my $prot = $name . "_gwise_prot.fas";
open(PROT, ">$prot");
my $orfout = $name . "_gwise_orf.fas";
open(ORF, ">$orfout");

my $seqtmp = $dir . "/blasttmp.fa";
my $outxtmp = $dir . "/blastxoutmp.txt";
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

	my $hid = $res[5];
	
	$hid =~ s/^ref\|//;
	$hid =~ s/\|$//;
		
	system("blastdbcmd -entry $hid -db $protdb > $newfile2");
	
	open(OUT, ">$newfile");
	print OUT ">" . $ori->display_id . "\n";
	print OUT $ori->seq . "\n";
	close(OUT);
	
	system("genewise $newfile2 $newfile -alg 333 -cdna -sum -both -silent > $aln");
	
	if (! -e $aln){
		next;
	}
	
	open(IN, $aln);
	my $sep = 0;
	my $regs;
	my $gogo = 0;
	while (<IN>){
		chomp;
		if (/^\/\//){
			$sep++;
			@{$regs->{$sep}} = ();
			$gogo =0;
		}
		if (/^>/){
			$gogo++;
		}elsif ($gogo > 0){
			push(@{$regs->{$sep}},$_);
		}
	}
	close(IN);
	
	$gogo = 0;
	
	foreach my $sep (sort {$a<=>$b} keys %{$regs}){
		my @seqs = @{$regs->{$sep}};
		if (!@seqs){
			next;
		}
		$gogo++;
		
		my $orf = "@seqs";
		$orf =~ s/\s+//g;
		$orf =~ s/-//g;
				
		print ORF ">" . $ori->display_id . "_" . $gogo;
		print ORF " ORF\n";
		print ORF $orf . "\n";
		
		my $seqobj2 = Bio::Seq->new(-seq => $orf);
		$seqobj2 = $seqobj2->translate;
		
		print PROT ">" . $ori->display_id . "_" . $gogo;
		print PROT " PROT\n";
		print PROT $seqobj2->seq . "\n";
	}
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
