#!/usr/bin/perl -w
##########FilterSV.pl########
#Author: Ronak Shah
#Date: 01/13/2014
#Version: 1.0
#Description: Merge's Traslocations for a samples & filter's Them
#process.
#############
###Log:
##01/13/2014
#v1.0
use strict;
use Getopt::Long;
use POSIX;

#--This variable holds the current time
my $now = time;
my (
	 $dir,                             $hotspotFile,
	 $tumorId,                         $normalId,
	 $outdir,                          $outputFile,
	 $overallSupportingReadsPE,        $sampleSupportingReadsPE,
	 $overallSupportingReadsSR,        $overallSupportingReadsPE_Hotspot,
	 $sampleSupportingReadsPE_Hotspot, $overallSupportingReadsSR_Hotspot,
	 $overallMAPQ,                     $overallMAPQ_Hotspot,
	 $distanceInTNcall,                $annotationInput
);
if (
	 @ARGV < 1
	 or !GetOptions(
				 'directory|dir:s'                => \$dir,
				 'hotspotFile|hsf:s'              => \$hotspotFile,
				 'tumorId|tid:s'                  => \$tumorId,
				 'normalId|nid:s'                 => \$normalId,
				 'outdir:s'                       => \$outdir,
				 'outputFile|o:s'                 => \$outputFile,
				 'overallSupportingReadsPE|ope:i' => \$overallSupportingReadsPE,
				 'overallSupportingReadsSR|osr:i' => \$overallSupportingReadsSR,
				 'overallSupportingReadsPEforHotspot|opeh:i' =>
				   \$overallSupportingReadsPE_Hotspot,
				 'overallSupportingReadsSRforHotspot|osrh:i' =>
				   \$overallSupportingReadsSR_Hotspot,
				 'overallMAPQ|omq:i'            => \$overallMAPQ,
				 'overallMAPQForHotspot|omqh:i' => \$overallMAPQ_Hotspot,
				 'distanceInTNcall|d:i'         => \$distanceInTNcall
	 )
  )
{
	Usage();
}

#Checks for Input
#dir
if ( ( !$dir ) or ( !( -d $dir ) ) )
{
	print "Please enter the name of the directory to be analyzed. See Usage.\n";
	Usage();
	exit(1);
} else
{
	print "DIR:$dir\n";
}

#hotspotFile
if ( ( !$hotspotFile ) or ( !( -f $hotspotFile ) ) )
{
	print
	  "Please enter the name of the hotspot file to be analyzed. See Usage.\n";
	Usage();
	exit(1);
} else
{
	print "HotspotFile:$hotspotFile\n";
}

#tumorId
if ( !$tumorId )
{
	print "Please enter the id for the tumor bam file. See Usage.\n";
	Usage();
	exit(1);
} else
{
	print "TumorID:$tumorId\n";
}

#normalId
if ( !$normalId )
{
	print "Please enter the id for the tumor bam file. See Usage.\n";
	Usage();
	exit(1);
} else
{
	print "NormalID:$normalId\n";
}

#outdir
if ( !$outdir )
{
	print "Output dir not given default will be used.\n";
	$outdir = getcwd;
	print "OUPUTDIR:$outdir\n";
} else
{
	if ( !( -d $outdir ) )
	{
		print "Output dir given is incorrect or doesnot exists.\n";
		exit(1);
	} else
	{
		print "OUTPUTDIR:$outdir.\n";
	}
}

#outputFile
if ( !$outputFile )
{
	print "Please enter the output file name. See Usage.\n";
	Usage();
	exit(1);
} else
{
	if ( $outputFile =~ /\./ )
	{
		($annotationInput) = $outputFile =~ /(.*)\./;
		$annotationInput = $annotationInput . "_dRangerInput.txt";
	} else
	{
		($annotationInput) = $outputFile;
		$annotationInput = $annotationInput . "_dRangerInput.txt";
	}
}

#overallSupportingReadsPE
if ( !$overallSupportingReadsPE )
{
	print "overallSupportingReadsPE not given default will be used.\n";
	$overallSupportingReadsPE = 10;
	print "overallSupportingReadsPE:$overallSupportingReadsPE\n";
} else
{
	print "overallSupportingReadsPE:$overallSupportingReadsPE.\n";
}

#sampleSupportingReadsPE
if ( !$sampleSupportingReadsPE )
{
	print "sampleSupportingReadsPE not given default will be used.\n";
	$sampleSupportingReadsPE = 2;
	print "sampleSupportingReadsPE:$sampleSupportingReadsPE\n";
} else
{
	print "sampleSupportingReadsPE:$sampleSupportingReadsPE.\n";
}

#sampleSupportingReadsPE_Hotspot
if ( !$sampleSupportingReadsPE_Hotspot )
{
	print "sampleSupportingReadsPE_Hotspot not given default will be used.\n";
	$sampleSupportingReadsPE_Hotspot = 1;
	print "sampleSupportingReadsPE_Hotspot:$sampleSupportingReadsPE_Hotspot\n";
} else
{
	print "sampleSupportingReadsPE:$sampleSupportingReadsPE.\n";
}

#overallSupportingReadsSR
if ( !$overallSupportingReadsSR )
{
	print "overallSupportingReadsSR not given default will be used.\n";
	$overallSupportingReadsSR = 3;
	print "overallSupportingReadsSR:$overallSupportingReadsSR\n";
} else
{
	print "overallSupportingReadsSR:$overallSupportingReadsSR.\n";
}

#overallSupportingReadsPE_Hotspot
if ( !$overallSupportingReadsPE_Hotspot )
{
	print "overallSupportingReadsPE_Hotspot not given default will be used.\n";
	$overallSupportingReadsPE_Hotspot = 5;
	print
	  "overallSupportingReadsPE_Hotspot:$overallSupportingReadsPE_Hotspot\n";
} else
{
	print
	  "overallSupportingReadsPE_Hotspot:$overallSupportingReadsPE_Hotspot.\n";
}

#overallSupportingReadsSR_Hotspot
if ( !$overallSupportingReadsSR_Hotspot )
{
	print "overallSupportingReadsSR_Hotspot not given default will be used.\n";
	$overallSupportingReadsSR_Hotspot = 1;
	print
	  "overallSupportingReadsSR_Hotspot:$overallSupportingReadsSR_Hotspot\n";
} else
{
	print
	  "overallSupportingReadsSR_Hotspot:$overallSupportingReadsSR_Hotspot.\n";
}

#overallMAPQ
if ( !$overallMAPQ )
{
	print "overallMAPQ not given default will be used.\n";
	$overallMAPQ = 10;
	print "overallMAPQ:$overallMAPQ\n";
} else
{
	print "overallMAPQ:$overallMAPQ.\n";
}

#overallMAPQ_Hotspot
if ( !$overallMAPQ_Hotspot )
{
	print "overallMAPQ_Hotspot not given default will be used.\n";
	$overallMAPQ_Hotspot = 5;
	print "overallMAPQ_Hotspot:$overallMAPQ_Hotspot\n";
} else
{
	print "overallMAPQ_Hotspot:$overallMAPQ_Hotspot.\n";
}

#distanceInTNcall
if ( !$distanceInTNcall )
{
	print "distanceInTNcall not given default will be used.\n";
	$distanceInTNcall = 10000;
	print "odistanceInTNcall:$distanceInTNcall\n";
} else
{
	print "distanceInTNcall:$distanceInTNcall.\n";
}
my ($tjmp)         = MergingPEandSRresults($tumorId);
my ($njmp)         = MergingPEandSRresults($normalId);
my ($tFilteredJmp) = FilterResults($tjmp);
my ($nFilteredJmp) = FilterResults($njmp);
CheckOccurenceInNormal( $tFilteredJmp, $nFilteredJmp );

#--Calculate total runtime (current time minus start time)
$now = time - $now;

#--Print runtime
printf( "\n\nTotal running time: %02d:%02d:%02d\n\n",
		int( $now / 3600 ),
		int( ( $now % 3600 ) / 60 ),
		int( $now % 60 ) );
print "\n!!!!Done, Thanks for using the script!!!\n";
exit;
#####################################
#####################################
#How to use the script.
sub Usage
{
	print "Unknow option: @_\n" if (@_);
	print "\nUsage : StructuralVariantFinder.pl [options]
        [--directory|dir 						S	 	Dir containing results for the translocation files.]
        [--hotspotFile|hsf 						S	 	BED4 Format bed file.]
        [--normalId|tid					S		ID for the normal bam file.]
        [--tumorId|tid					S		ID for the tumor bam file.]
        [--outdir						S		Output directory.]
        [--outputFile|o					S		Name of the Output file.]
        [--overallSupportingReadsPE|pe			I		Threshold for number of pair-end reads supporting the SV.]
        [--overallSupportingReadsSR|sr			I		Threshold for number of split reads supporting the SV.]
        [--overallSupportingReadsPE_Hotspot|peh			I		Threshold for number of pair-end reads supporting the SV for hotspot location.]
        [--overallSupportingReadsSR_Hotspot|srh			I		Threshold for number of split reads supporting the SV for hotspot location.]
        [--overallMAPQ|mq	I	Threshold for Mapping Quality.]
        [--overallMAPQ_Hotspot|mqh	I	Threshold for Mapping Quality for hotspot locations.]
        [--distanceInTNcall|d	I	Distance between a tumor call & a normal call.]
	\n";
	print "For any question please contact Ronak Shah (shahr2\@mskcc.org)\n";
	print "!~!~! Comments, flames and suggestion are welcome !~!~!\n";
	exit;
}
#####################################
#####################################
#Read the Hotspot Locations from Bed4 File
sub ReadHotspotFile
{
	my %hotspotLocs = ();
	open( HFH, $hotspotFile )
	  || die "Cannot open file $hotspotFile, Error: $!\n";
	while (<HFH>)
	{
		chomp($_);
		my ( $chr, $start, $end, $gene ) = split( "\t", $_ );

		#print "$chr, $start, $end, $gene\n";
		if ( exists $hotspotLocs{$chr} )
		{
			my $val = $hotspotLocs{$chr};
			$hotspotLocs{$chr} = $val . "#" . $start . ":" . $end . ":" . $gene;
		} else
		{
			$hotspotLocs{$chr} = $start . ":" . $end . ":" . $gene;
		}
	}
	close(HFH);
	return ( \%hotspotLocs );
}
#####################################
#####################################
#Merge All Files to a single record.
sub MergingPEandSRresults
{
	my ($id)          = @_;
	my $jmp           = "$dir/$id\_jmp.txt";
	my $jmpbrkpts     = "$dir/$id\_jmpbrkpts.txt";
	my %jmpData       = &ReadJMP($jmp);
	my %jmpbrkptsData = &ReadJMP($jmpbrkpts);

	#Work on Both PE & SR
	while ( my ( $Aid, $Adata ) = each(%jmpData) )
	{
		my ( $Allchr1, $Allpos1, $Allchr2, $Allpos2, $Allsupport, $Allmapq ) =
		  split( ",", $Adata );
		if ( exists $jmpbrkptsData{$Aid} )
		{
			my $value = $jmpbrkptsData{$Aid};
			my ( $chr1, $pos1, $chr2, $pos2, $support, $mapq ) =
			  split( ",", $value );
			$jmpData{$Aid} =
			  "$chr1,$pos1,$chr2,$pos2,$Allsupport,$Allmapq,$support";
		} else
		{
			$jmpData{$Aid} =
			  "$Allchr1,$Allpos1,$Allchr2,$Allpos2,$Allsupport,$Allmapq,0";
		}
	}
	return ( \%jmpData );
}
#####################################
#####################################
#Check if location is hotspot or not
sub CheckHotspot
{
	my ( $chr1, $pos1, $chr2, $pos2 ) = @_;

	#Read Hotspot Locations
	my $hotspotHash = &ReadHotspotFile();
	my %hotspotLocs = %$hotspotHash;
	my $hotspotFlag;
	if ( exists $hotspotLocs{$chr1} )
	{
		my $hotspotLocations = $hotspotLocs{$chr1};
		if ( $hotspotLocations =~ m/#/ )
		{
			my @multiLocs = split( "#", $hotspotLocations );
			foreach my $hsLocation (@multiLocs)
			{
				my ( $hsStart, $hsEnd, $hsGene ) =
				  split( ":", $hsLocation );
				if ( $pos1 >= $hsStart && $pos1 <= $hsEnd )
				{
					$hotspotFlag = 2;
					last;
				} else
				{
					$hotspotFlag = 1;
				}
			}
		} else
		{
			my ( $hsStart, $hsEnd, $hsGene ) =
			  split( ":", $hotspotLocations );
			if ( $pos1 >= $hsStart && $pos1 <= $hsEnd )
			{
				$hotspotFlag = 2;
			} else
			{
				$hotspotFlag = 1;
			}
		}
	} elsif ( exists $hotspotLocs{$chr2} )
	{
		my $hotspotLocations = $hotspotLocs{$chr2};
		if ( $hotspotLocations =~ m/#/ )
		{
			my @multiLocs = split( "#", $hotspotLocations );
			foreach my $hsLocation (@multiLocs)
			{
				my ( $hsStart, $hsEnd, $hsGene ) =
				  split( ":", $hsLocation );
				if ( $pos2 >= $hsStart && $pos2 <= $hsEnd )
				{
					$hotspotFlag = 2;
					last;
				} else
				{
					$hotspotFlag = 1;
				}
			}
		} else
		{
			my ( $hsStart, $hsEnd, $hsGene ) =
			  split( ":", $hotspotLocations );
			if ( $pos2 >= $hsStart && $pos2 <= $hsEnd )
			{
				$hotspotFlag = 2;
			} else
			{
				$hotspotFlag = 1;
			}
		}
	} else
	{
		$hotspotFlag = 1;
	}
	return ($hotspotFlag);
}
#####################################
#####################################
#Filter Results for each SV
sub FilterResults
{
	my ($jmp)      = @_;
	my %jmpData    = %$jmp;
	my %outjmpData = ();
	while ( my ( $Aid, $Adata ) = each(%jmpData) )
	{
		my ( $chr1, $pos1, $chr2, $pos2, $PEreads, $mapq, $SRreads ) =
		  split( ",", $Adata );

		#Start Checking if its Hotspot or not.
		my $hotspotFlag = CheckHotspot( $chr1, $pos1, $chr2, $pos2 );

		#Check Hotspot and filter
		if ( $hotspotFlag == 2 )
		{
			if (     ( $PEreads >= $overallSupportingReadsPE_Hotspot )
				 and ( $mapq >= $overallMAPQ_Hotspot )
				 and ( $SRreads >= $overallSupportingReadsSR_Hotspot ) )
			{
				$outjmpData{$Aid} = $Adata;

			   #print "$chr1\t$pos1\t$chr2\t$pos2\t$PEreads\t$SRreads\t$mapq\n";
			}
		} else
		{
			if (     ( $PEreads >= $overallSupportingReadsPE )
				 and ( $mapq >= $overallMAPQ )
				 and ( $SRreads >= $overallSupportingReadsSR ) )
			{
				$outjmpData{$Aid} = $Adata;

			   #print "$chr1\t$pos1\t$chr2\t$pos2\t$PEreads\t$SRreads\t$mapq\n";
			}
		}
	}
	return ( \%outjmpData );
}
#####################################
#####################################
#Check if there is evidence in SV
sub CheckOccurenceInNormal
{
	my ( $Tjmp, $Njmp ) = @_;
	my %tjmp                   = %$Tjmp;
	my %njmp                   = %$Njmp;
	my %normalChrCordHash      = ();
	my %normalChrBasedDataHash = ();
	open( OFH, ">", "$outdir/$outputFile" )
	  || die "Cannot open $outdir/$outputFile. Error $!\n";
	print OFH
"Chr1\tTumorPos1\tChr2\tTumorPos2\tDirection\tTumorPE\tTumorSR\tTumorMAPQ\tCanBeSomatic\tNormalPos1\tNormalPos2\tNormalPE\tNormalSR\tNormalMAPQ\n";

	#Open the dranger input File
	open( DFH, ">", "$outdir/$annotationInput" )
	  || die "Cannot open file $outdir/$annotationInput, Error:$!\n";
	print DFH "chr1\tpos1\tstr1\tchr2\tpos2\tstr2\n";

	#Reporcess data for noraml set making two different hashes.
	while ( my ( $Aid, $Adata ) = each(%njmp) )
	{
		my ( $chr1, $pos1, $chr2, $pos2, $PEreads, $mapq, $SRreads ) =
		  split( ",", $Adata );
		if ( exists $normalChrCordHash{"$chr1:$chr2"} )
		{
			my $exisistingPos = $normalChrCordHash{"$chr1:$chr2"};
			$normalChrCordHash{"$chr1:$chr2"} =
			  $exisistingPos . "#" . "$pos1;$pos2";
		} else
		{
			$normalChrCordHash{"$chr1:$chr2"} = "$pos1;$pos2";
		}
		$normalChrBasedDataHash{"$chr1,$pos1,$chr2,$pos2"} =
		  "$PEreads,$mapq,$SRreads";
	}

	#Traverse through the Tumor hash to find Comman SV with Normal.
	while ( my ( $Aid, $Adata ) = each(%tjmp) )
	{
		my ( $chr1, $pos1, $chr2, $pos2, $PEreads, $mapq, $SRreads ) =
		  split( ",", $Adata );
		next if(($chr1 !~ /^(\d{1,2}|X|Y)/) or ($chr2 !~ /^(\d{1,2}|X|Y)/));
		my $convertedChr1 = $chr1;
		my $convertedChr2 = $chr2;
		$convertedChr1 = "23" if($chr1 =~ /^X/);
		$convertedChr1 = "24" if($chr1 =~ /^Y/);
		$convertedChr2 = "23" if($chr2 =~ /^X/);
		$convertedChr2 = "24" if($chr2 =~ /^Y/);
		my @splitId = split( "_", $Aid );
		my $direction = $splitId[1];
		my ( $startCT, $endCT ) = split( "to", $direction );
		my ( $str1, $str2 ) = "";
		$str1 = "0" if ( $startCT == 3 );
		$str2 = "0" if ( $endCT == 3 );
		$str1 = "1" if ( $startCT == 5 );
		$str2 = "1" if ( $endCT == 5 );

		#Start Checking if its Hotspot or not.
		my $hotspotFlag = CheckHotspot( $chr1, $pos1, $chr2, $pos2 );
		if ( exists $normalChrCordHash{"$chr1:$chr2"} )
		{
			my $allNormalPos = $normalChrCordHash{"$chr1:$chr2"};
			if ( $allNormalPos =~ /#/ )
			{
				my @allPos = split( "#", $allNormalPos );
				my $finalCords = "";
				foreach my $pos (@allPos)
				{
					my ( $npos1, $npos2 ) = split( ";", $pos );
					my $npos1_start = $npos1 - $distanceInTNcall;
					my $npos2_start = $npos2 - $distanceInTNcall;
					my $npos1_end   = $npos1 + $distanceInTNcall;
					my $npos2_end   = $npos2 + $distanceInTNcall;
					if (
						(
						   ( $npos1_end >= $pos1 ) and ( $pos1 >= $npos1_start )
						)
						and (     ( $npos2_end >= $pos2 )
							  and ( $pos2 >= $npos2_start ) )
					  )
					{
						$finalCords = $pos;
						last;
					}
				}
				if ( $finalCords eq "" )
				{
					if ( $hotspotFlag == 2 )
					{
						if ( $mapq >= $overallMAPQ_Hotspot )
						{
							print OFH
"$chr1\t$pos1\t$chr2\t$pos2\t$direction\t$PEreads\t$SRreads\t$mapq\tYES\t-\t-\t-\t-\t-\n";
							print DFH
							  "$convertedChr1 \t$pos1\t$str1\t$convertedChr2\t$pos2\t$str2\n";
						}
					} else
					{
						if ( $mapq >= $overallMAPQ )
						{
							print OFH
"$chr1\t$pos1\t$chr2\t$pos2\t$direction\t$PEreads\t$SRreads\t$mapq\tYES\t-\t-\t-\t-\t-\n";
							print DFH
							  "$convertedChr1\t$pos1\t$str1\t$convertedChr2\t$pos2\t$str2\n";
						}
					}
				} else
				{
					my ( $npos1, $npos2 ) = split( ";", $finalCords );
					my $npos1_end = $npos1 + $distanceInTNcall;
					my $npos2_end = $npos2 + $distanceInTNcall;
					my $nValues   =
					  $normalChrBasedDataHash{"$chr1,$npos1,$chr2,$npos2"};
					my ( $nPEreads, $nmapq, $nSRreads ) =
					  split( ",", $nValues );
					if ( $hotspotFlag == 2 )
					{
						if (     #( $nPEreads <= 2 )
							#and 
							( $mapq >= $overallMAPQ_Hotspot ) )
						{
							print OFH
"$chr1\t$pos1\t$chr2\t$pos2\t$direction\t$PEreads\t$SRreads\t$mapq\tNO\t$npos1\t$npos2\t$nPEreads\t$nSRreads\t$nmapq\n";
							print DFH
							  "$convertedChr1\t$pos1\t$str1\t$convertedChr2\t$pos2\t$str2\n";
						}
					} else
					{
						if ( #( $nPEreads <= 1 ) and 
						( $mapq >= $overallMAPQ ) )
						{
							print OFH
"$chr1\t$pos1\t$chr2\t$pos2\t$direction\t$PEreads\t$SRreads\t$mapq\tNO\t$npos1\t$npos2\t$nPEreads\t$nSRreads\t$nmapq\n";
							print DFH
							  "$convertedChr1\t$pos1\t$str1\t$convertedChr2\t$pos2\t$str2\n";
						}
					}
				}
			} else
			{
				my ( $npos1, $npos2 ) = split( ";", $allNormalPos );
				my $npos1_start = $npos1 - $distanceInTNcall;
				my $npos2_start = $npos2 - $distanceInTNcall;
				my $npos1_end   = $npos1 + $distanceInTNcall;
				my $npos2_end   = $npos2 + $distanceInTNcall;
				if ( ( ( $npos1_end >= $pos1 ) and ( $pos1 >= $npos1_start ) )
					 and
					 ( ( $npos2_end >= $pos2 ) and ( $pos2 >= $npos2_start ) ) )
				{

					#print "$npos1\t$npos1_end\t$npos2\t$npos2_end\n";
					my $nValues =
					  $normalChrBasedDataHash{"$chr1,$npos1,$chr2,$npos2"};
					my ( $nPEreads, $nmapq, $nSRreads ) =
					  split( ",", $nValues );
					if ( $hotspotFlag == 2 )
					{
						if (     #( $nPEreads <= 2 )
							 #and 
							 ( $mapq >= $overallMAPQ_Hotspot ) )
						{
							print OFH
"$chr1\t$pos1\t$chr2\t$pos2\t$direction\t$PEreads\t$SRreads\t$mapq\tNO\t$npos1\t$npos2\t$nPEreads\t$nSRreads\t$nmapq\n";
							print DFH
							  "$convertedChr1\t$pos1\t$str1\t$convertedChr2\t$pos2\t$str2\n";
						}
					} else
					{
						if ( #( $nPEreads <= 1 ) and 
						( $mapq >= $overallMAPQ ) )
						{
							print OFH
"$chr1\t$pos1\t$chr2\t$pos2\t$direction\t$PEreads\t$SRreads\t$mapq\tNO\t$npos1\t$npos2\t$nPEreads\t$nSRreads\t$nmapq\n";
							print DFH
							  "$convertedChr1\t$pos1\t$str1\t$convertedChr2\t$pos2\t$str2\n";
						}
					}
				} else
				{
					if ( $mapq >= $overallMAPQ )
					{
						print OFH
"$chr1\t$pos1\t$chr2\t$pos2\t$direction\t$PEreads\t$SRreads\t$mapq\tYES\t-\t-\t-\t-\t-\n";
						print DFH "$convertedChr1\t$pos1\t$str1\t$convertedChr2\t$pos2\t$str2\n";
					}
				}
			}
		} else
		{
			if ( $mapq >= $overallMAPQ )
			{
				print OFH
"$chr1\t$pos1\t$chr2\t$pos2\t$direction\t$PEreads\t$SRreads\t$mapq\tYES\t-\t-\t-\t-\t-\n";
				print DFH "$convertedChr1\t$pos1\t$str1\t$convertedChr2\t$pos2\t$str2\n";
			}
		}
	}
	close(OFH);
	close(DFH);
	return;
}
#####################################
#####################################
#Read the std Traslocation File
sub ReadJMP
{
	my ($file) = @_;
	my %jmpData = ();
	open( FH, $file ) || die "Cannot Open $file. Error:$!\n";
	while (<FH>)
	{
		chomp;
		if ( $_ =~ /Translocation/ )
		{
			my ( $chr1, $pos1, $chr2, $pos2, $support, $mapq, $id ) =
			  split( "\t", $_ );
			$jmpData{$id} = "$chr1,$pos1,$chr2,$pos2,$support,$mapq";
		}
	}
	close(FH);
	return (%jmpData);
}
