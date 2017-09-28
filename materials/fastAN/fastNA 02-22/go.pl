#!/usr/bin/perl -w
use Getopt::Long;  # to parse command line parameters
use strict;
use warnings;

my $line;
my $inFile;
my $outFile;
my $inFileName1;
my $inFileName2;
my $inFileName3;
my $outFileName;

my @files = (  ); 
my $count = 0;
my %sequencesHash = ( );
my $G1NodeID;
my $G2NodeID;
my @values;
my %G1GoHash;
my %G2GoHash;

my %G1EdgeHash;
my %G2EdgeHash;
my %AlignHash;

$inFileName1 = "data/dmela/celeg_gene_to_go.txt";
$inFileName2 = "data/hsapi/dmela_gene_to_go.txt";
#$inFileName3 = "ga/dmela.ga_hsapi.ga.out_name";
$inFileName3 = "celeg.ga_dmela.ga.out_name";
#$inFileName3 = "data/dm-hs/dmela-hsapi-spinalmatch-original7";

#$inFileName1 = "data/celeg/celeg_gene_to_go_experimental.txt";
#$inFileName2 = "data/dmela/dmela_gene_to_go_experimental.txt";
$inFileName1 = "data/celeg/celeg_gene_to_go.txt";
$inFileName2 = "data/dmela/dmela_gene_to_go.txt";
#$inFileName3 = "ga8/celeg.ga_dmela.ga.out_name";
#$inFileName3 = "data/ce-dm/celeg-dmela-spinalmatch-original7";
$inFileName3 = "mi-graal/cd01.aln";


$inFileName1 = "data/celeg/celeg_gene_to_go.txt";
$inFileName2 = "data/hsapi/hsapi_gene_to_go.txt";
#$inFileName3 = "ga/dmela.ga_hsapi.ga.out_name";
$inFileName3 = "celeg.ga_hsapi.ga.out_name";
#$inFileName3 = "data/ce-hs/celeg-hsapi-spinalmatch-original7";

#$inFileName1 = $ARGV[1];
#$inFileName2 = $ARGV[2];
#$inFileName3 = $ARGV[3];

#$outFileName = $inFileName3.".go";

open ( $inFile, '<', $inFileName1 )		   or die "Cannot open the input file";
while ($line = <$inFile>) {
	chomp($line);
	@values = split(" ", $line);
	my $pName = $values[0];
	my $goName;
	if (defined ($values[1])) {
		$goName = $values[1];
		@values = split(/\|/, $pName);
		$pName = $values[0];
		$G1GoHash{$pName} = $goName;
#		print "count: $pName\n";
	}

}
close $inFile;

open ( $inFile, '<', $inFileName2 )		   or die "Cannot open the input file";
while ($line = <$inFile>) {
	chomp($line);
	@values = split(" ", $line);
	my $pName = $values[0];
	my $goName;
	if (defined ($values[1])) {
		$goName = $values[1];
		@values = split(/\|/, $pName);
		$pName = $values[0];
		$G2GoHash{$pName} = $goName;
#		print "count: $pName\n";
	}

}
close $inFile;

open ( $inFile, '<',  $inFileName3 )		   or die "Cannot open the input file";
#open ( $outFile, '>',  $outFileName )		   or die "Cannot open the output file";

$line = <$inFile>;
my @GO1Array;
my @GO2Array;
my @countGOArray;
my $GOC=0;

while ($line = <$inFile>) {
	chomp($line);
	@values = split(" ", $line);
	my $G1pName = $values[0];
	my $G2pName = $values[1];
	my $GOs1 = $G1GoHash{$G1pName};
		if (defined ($GOs1)) {
		@GO1Array = split(/\|/, $GOs1);
#		print "$GOs1\n";
	}
	my $GOs2 = $G2GoHash{$G2pName};
		if (defined ($GOs2)) {
		@GO2Array = split(/\|/, $GOs2);
#		print "$GOs1\n";
	}
	
	my %GO1Map=map{$_ =>1} @GO1Array;
	my %GO2Map=map{$_=>1} @GO2Array;
	my @intersect = grep( $GO1Map{$_}, @GO2Array );
	my $iVal = scalar(@intersect);
	
	my %union = ();
	foreach(@GO1Array,@GO2Array){
    	$union{$_}=1;
	}
	my @union2 = keys %union;
	my $uVal = scalar(@union2);
	if ($uVal == 0) {next;}
	$GOC += $iVal/$uVal;

	

	$countGOArray[$iVal]++;
}
close $inFile;
#close $outFile;

print "\n\tinput: $inFileName3\n\tGOC: $GOC\n\n\t#count\t#GO_terms\n";
for (my $i=1; $i<scalar(@countGOArray); $i++){
	if (defined ($countGOArray[$i])) {
		print "\t$countGOArray[$i]\t$i\n";
	}
}
print "\n";







