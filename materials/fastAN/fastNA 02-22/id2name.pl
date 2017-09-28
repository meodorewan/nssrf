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
my %G1IDNameHash;
my %G2IDNameHash;
my %G1EdgeHash;
my %G2EdgeHash;
my %AlignHash;

$inFileName1 = "data/dmela/dmela.tab.gml";
$inFileName2 = "data/hsapi/hsapi.tab.gml";
$inFileName3 = "ga8/dmela.ga_hsapi.ga.out";

$inFileName1 = "data/scere/scere.tab.gml";
$inFileName2 = "data/dmela/dmela.tab.gml";
$inFileName3 = "ga8/scere.ga_dmela.ga.out";

$inFileName1 = "data/celeg/celeg.tab.gml";
$inFileName2 = "data/hsapi/hsapi.tab.gml";
$inFileName3 = "celeg.ga_hsapi.ga.out";

#$inFileName1 = $ARGV[1];
#$inFileName2 = $ARGV[2];
#$inFileName3 = $ARGV[3];

$outFileName = $inFileName3."_name";

open ( $inFile, '<', $inFileName1 )		   or die "Cannot open the input file";
while ($line = <$inFile>) {
	chomp($line);
	if ($line =~ /node \[/) {
#		print "count: $line\n";
		$line = <$inFile>;
		@values = split(" ", $line);
		my $NodeID = $values[1];
#		print "count: $NodeID\n";
		$line = <$inFile>;
		@values = split("\"", $line);
		my $NodeName = $values[1];
#		print "count: $NodeName\n";
		$G1IDNameHash{$NodeID} = $NodeName;
#		print "$NodeID\t$NodeName\n";
	}
}
close $inFile;

open ( $inFile, '<',  $inFileName2 )		   or die "Cannot open the input file";
while ($line = <$inFile>) {
	chomp($line);
	if ($line =~ /node \[/) {
#		print "count: $line\n";
		$line = <$inFile>;
		@values = split(" ", $line);
		my $NodeID = $values[1];
#		print "count: $NodeID\n";
		$line = <$inFile>;
		@values = split("\"", $line);
		my $NodeName = $values[1];
#		print "count: $NodeName\n";
		$G2IDNameHash{$NodeID} = $NodeName;
#		print "$NodeID\t$NodeName\n";
	}
}
close $inFile;

open ( $inFile, '<',  $inFileName3 )		   or die "Cannot open the input file";
open ( $outFile, '>',  $outFileName )		   or die "Cannot open the output file";

$line = <$inFile>;
print $outFile "!".$line;

while ($line = <$inFile>) {
	chomp($line);
	@values = split(" ", $line);
	my $G1NID = $values[0];
	my $G2NID = $values[1];
	my $G1NName = $G1IDNameHash{$G1NID};
	my $G2NName = $G2IDNameHash{$G2NID};
	print $outFile "$G1NName $G2NName\n";
}

close $inFile;
close $outFile;





