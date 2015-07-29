#!/usr/bin/perl -w
##########MergeResultsAcrossPools.pl########
#Author: Ronak Shah
#Date: 02/24/2014
#Version: 1.0
#Description: MergeSV's across pools
#process.
#############
###Log:
##02/24/2014
#v1.0
use strict;
use Getopt::Long;
use POSIX;
use Tie::IxHash;
my ( $tumortypeFile, $poolList, $outFile );
if (
	@ARGV < 1
	or !GetOptions(
		'tumortypeFile|ttf:s' => \$tumortypeFile,
		'poolList|i:s'        => \$poolList,
		'outFile|o:s'         => \$outFile,

	)
  )
{
	Usage();
}

#Checks for Input
#poolList
if ( ( !$poolList ) or ( !( -f $poolList ) ) ) {
	print
"Please enter the name of the pool list file to be analyzed. See Usage.\n";
	Usage();
	exit(1);
}
else {
	print "Pool List:$poolList\n";
}

#tumortypeFile
if ( ( !$tumortypeFile ) or ( !( -f $tumortypeFile ) ) ) {
	print
	  "Tumortype file not given, tumor type column would not be popultated.\n";
}
else {
	print "TumorType File:$tumortypeFile\n";
}

#outFile
if ( !$outFile ) {
	print "Please Enter the output fle name with location. See Usage.\n";
	Usage();
	exit;
}
else {
	print "Output File:$outFile\n";
}
#####################################
#####################################
#How to use the script.
sub Usage {
	print "Unknow option: @_\n" if (@_);
	print "\nUsage : MergeResultsAcrossPools.pl [options]
        [--poolList|i 							S	 	input Structural Variant[SV] pool list file.]
		[--tumortypeFile|ttf 					S	 	File Containing ID\tTumortype.]	
        [--outFile|o					        S		Name of the Output file.]
      
	\n";
	print "For any question please contact Ronak Shah (shahr2\@mskcc.org)\n";
	print "!~!~! Comments, flames and suggestion are welcome !~!~!\n";
	exit;
}
my $tumortype = "";
my $ttheader  = "";
if ($tumortypeFile) {
	( $ttheader, $tumortype ) = &ReadTumorTypeFile($tumortypeFile);
}
my ($list) = &ReadPoolList($poolList);
&MergeData();
exit;
#####################################
#####################################
#Read TumorType File
sub ReadTumorTypeFile {
	my ($file) = @_;
	my @header;
	tie( my %ttHash, 'Tie::IxHash' );
	open( FH, $file ) || die "ReadTumorTypeFile:Cannot open $file, Error:$!\n";
	while (<FH>) {
		chomp;
		if ( $. == 1 ) {
			@header = split( "\t", $_ );
			shift(@header);
		}
		else {
			my @data     = split( "\t", $_ );
			my $id       = shift(@data);
			my @allIdVal = split( "-", $id );
			my $key      = $allIdVal[0] . "-" . $allIdVal[1];

			#print "$key\n";
			if ( scalar @data <= 6 ) {
				my $totalVal = 6 - ( scalar @data );
				for ( my $i = 0 ; $i < $totalVal ; $i++ ) {
					push( @data, "-" );
				}
			}
			my $value = join( ";", @data );
			$ttHash{$key} = $value;
		}
	}
	close(FH);
	return ( \@header, \%ttHash );
}
#####################################
#####################################
#Read PoolList
sub ReadPoolList {
	my ($file) = @_;
	my @list = ();
	open( FH, $file ) || die "ReadPoolList:Cannot open $file, Error:$!\n";
	while (<FH>) {
		chomp;
		push( @list, $_ );
	}
	close(FH);
	return ( \@list );
}
#####################################
#####################################
#MergeData
sub MergeData {
	my (@poolList)     = @$list;
	my (@ttheaderName) = @$ttheader;
	my %tumortypeInfo  = %$tumortype;
	my $count          = 1;
	open( OFH, ">", $outFile ) || die "Merge:Cannot open $outFile, Error:$!\n";
	foreach my $poolLoc (@poolList) {
		my $poolName   = pop @{ [ split( "/", $poolLoc ) ] };
		my $file       = $poolLoc . "/" . $poolName . "_All_AnnotatedSVs.txt";
		my @allheaders = ();
		my @alldata    = ();
		open( FH, $file ) || die "Merge:Cannot open $file, Error:$!\n";
		print "COUNT=$count\n";
		while (<FH>) {
			chomp;
			if ( $. == 1 ) {
				if ( $count <= 1 ) {
					my @header = split( "\t", $_ );
					@allheaders = ( @header, @ttheaderName, "PoolName" );
					my $header = join( "\t", @allheaders );
					print OFH "$header\n";
				}
				next;
			}
			else {
				my @data     = split( "\t", $_ );
				my @allIdVal = split( "-",  $data[0] );
				my $key      = $allIdVal[0] . "-" . $allIdVal[1];
				my $ttype;
				$ttype = $tumortypeInfo{$key};
				if ( !$ttype ) {
					$ttype = "-;-;-;-;-;-";
				}
				my @ttinfo = split( ";", $ttype );
				@alldata = ( @data, @ttinfo, $poolName );
				my $data = join( "\t", @alldata );
				print OFH "$data\n";

			}
		}
		$count++;
		close(FH);
	}
	close(OFH);
	return;
}
