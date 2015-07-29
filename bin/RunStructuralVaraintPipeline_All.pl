#!/usr/bin/perl -w
##########StructuralVaraintFinder.pl########
#Author: Ronak Shah
#Date: 07/01/2013
#LastModified: 08/22/2013
#Version: 1.0
#Description: Get In the data and run the mapping and structural variant finding
#process.
#############
###Log:
##07/01/2013
#v1.0
#use lib qw(/ifs/e63data/bergerm1/Resources/Libs/PerlLib/v5.10.1/CPAN/lib64/perl5/auto /home/shahr2/perl/lib/CPAN/lib/perl5/site_perl/5.8.8/ /home/shahr2/perl/lib/CPAN/lib64/perl5/site_perl/5.8.8/x86_64-linux-thread-multi/auto/ /home/shahr2/Software/CREST/current /home/shahr2/perl/lib/CPAN/lib/perl5/5.8.8 /home/shahr2/perl/lib/CPAN/lib/perl5/x86_64-linux-thread-multi /home/shahr2/perl/lib/CPAN/lib/perl5/ /home/shahr2/Software/Vcftools/current/perl /home/shahr2/perl/lib/CPAN/ /home/shahr2/perl/lib/CPAN/lib/ /home/shahr2/perl/lib/CPAN/lib/lib/ /usr/lib64/perl5/site_perl/5.8.8/x86_64-linux-thread-multi /usr/lib/perl5/site_perl/5.8.8 /usr/lib/perl5/site_perl /usr/lib64/perl5/vendor_perl/5.8.8/x86_64-linux-thread-multi /usr/lib/perl5/vendor_perl/5.8.8 /usr/lib/perl5/vendor_perl /usr/lib64/perl5/5.8.8/x86_64-linux-thread-multi /usr/lib/perl5/5.8.8);
use lib qw(/ifs/e63data/bergerm1/Resources/Libs/PerlLib/v5.10.1/CPAN/share/perl5 /home/shahr2/Software/Vcftools/current/perl);

use strict;
use Getopt::Long;
use IO::File;
use List::Util qw(max sum);
use Tie::IxHash;
use Vcf;
use File::Basename;

#--This variable holds the current time
my $now = time;
my ($sampleFile,$stdNormal,$titleFile,$fof,$poolName,$projectName,$process,$datadir,$outdir,$mvFiles,$bwa,$PeSV,$baitIntervalFile,$targetIntervalFile,$targetIntervalBedFile,$refFile,$TMPDIR,$PICARD,$JAVA,$chrSepRef,$samtoolsDir,$dbUser,$dbPass,$dbName,$dbLocStr,$hg,$fastqSource,$barcodeFile,$config_file,$awk,$nprocessors,$ntchr,$GATK,$isTargeted,$refFasta,$bwaDir,$Meerkat,$prog,$blastDir,$RefGeneFile,$RepeatMaskFile,$queue,$PERL,$bwaFlag,$DELLY);

if (@ARGV < 1 or !GetOptions (
            'config|c:s'                        => \$config_file))
	{
		Usage();
	}
#Get Configration File details
my($Location,$Version) = &GetConfigration($config_file,$outdir);
my %location = %$Location;

###Set Variables
$bwaFlag = $location{"BWAFLAG"};
$isTargeted = $location{"isTargeted"};
$TMPDIR = $location{"TMPDIR"};
$PERL = $location{"PERL"};
$JAVA = $location{"JAVA"};
$GATK =  $location{"GATK"};
$refFile =  $location{"Reference"};
$refFasta = $location{"RefFasta"};
$chrSepRef =  $location{"Refseq"};
$PICARD =  $location{"PICARD"};
$baitIntervalFile =  $location{"BaitInterval"};
$targetIntervalFile =  $location{"TargetInterval"};
$targetIntervalBedFile =  $location{"TargetIntervalBed"};
$bwa =  $location{"BWA"};
$bwaDir =  $location{"BWADIR"};
$samtoolsDir = $location{"SAMTOOLSDIR"};
$sampleFile = $location{"SampleFile"};
$titleFile = $location{"TitleFile"};
$PeSV = $location{"PeSVFisher"};
$Meerkat = $location{"Meerkat"};
$DELLY = $location{"Delly"};
$hg = $location{"HumanGenomeVersion"};
$awk =  $location{"awk"};
$barcodeFile =  $location{"barcodeFile"};
$fastqSource = $location{"fastqSource"};
$outdir = $location{"outdir"};
$datadir = $location{"datadir"};
$stdNormal = $location{"stdNormal"};
$fof = $location{"FOF"};
$poolName = $location{"poolName"};
$projectName = $location{"projectName"};
$dbUser = $location{"databaseUserId"};
$dbPass = $location{"databasePassword"};
$dbName = $location{"databaseName"};
$dbLocStr = $location{"databaseLocationString"};
$mvFiles = $location{"moveFiles"};
$process = $location{"Process"};
$ntchr = $location{"NumberOfChromosomes"};
$nprocessors = $location{"NumberOfProcessors"};
$prog = $location{"Program"};
$blastDir = $location{"BlastDir"};
$RepeatMaskFile = $location{"RepeatMaskFile"};
$RefGeneFile = $location{"RefGeneFile"};
$queue = $location{"SGE_QUEUE"};


my $PrintConfigFile = "RunConfigration.txt";
open (VERSION,">","$outdir/$PrintConfigFile") || die "Cannot open $outdir/$PrintConfigFile, Reason:$!\n";
#Prin Version of tools and Files used
print VERSION "Tools|Files\tVersion\n";
while(my($tools_files,$version) = each(%$Version))
{
    print VERSION "$tools_files\t$version\n";
}
close(VERSION);

if((! $sampleFile) or (! $titleFile))
{
    print "Please provide the sample information as well as the title file name.See Usage\n";
    Usage();
    exit;
}
if(! $projectName)
{
    print"Please enter the project name. See Usage\n";
    Usage();
    exit;
}
if(! $poolName)
{
    print"Please enter the pool name. See Usage\n";
    Usage();
    exit;
}
if(! $mvFiles)
{
    print"Folders will be created and Files will be moved.\n";
    $mvFiles = 1;
}
else
{
    if($mvFiles == 1)
    {
	print"Folders will be created and Files will be moved.\n";
    }
    if($mvFiles == 2)
    {
	print"Folders will not be created and Files will not be moved.\n";
    }
}
if(! $bwaFlag)
{
    print"Default:BWA ALN will be used for alignment.\n";
    $bwaFlag = 1;
}
else
{
    if($bwaFlag == 1)
    {
	print"BWA ALN will be used for alignment.\n";
    }
    if($bwaFlag == 2)
    {
	print"BWA MEM will be used for alignment\n";
    }
}
if(! $stdNormal)
{
    $stdNormal = "NULL";
}
else
{
    print "Starndard Normal is given: $stdNormal\n";
}
if(! $datadir)
{
    print "Please enter the directory that contains the data to be processed.See Usage\n";
    Usage();
    exit;
}
if (! $outdir)
{
   print "Please enter the directory where to write the data while/after being processed.See Usage\n";
    Usage();
    exit;
}
else
{
    print "Results will be written in $outdir\n";
}
if(! $fastqSource)
{
    print "Assume fastq files came from GCL.\n";
    $fastqSource = "GCL";
}
else
{
    if($fastqSource ne "GCL" && $fastqSource ne "Path")
    {
	 print "Please indicate fastqSource. See Usage\n";
	 Usage();
	 exit;
     }
}
if($barcodeFile)
{
    print"The barcode file in use is $barcodeFile.\n";
}
else
{
    $barcodeFile = "/home/shahr2/Scripts/All/barcodekey.txt";
    print"The barcode file in use is $barcodeFile.\n";
}
if($config_file)
{
    print"The configration file in use is $config_file.\n";
}
else
{
    $config_file = "/home/shahr2/Scripts/All/StrVarConfig_v1.0_v4.txt";
    print"The configration file in use is $config_file.\n";
}
if(! $TMPDIR)
{
    print "Path to temporary directory is not give and default will be used.\n";
    $TMPDIR = "/ifs/e63data/bergerm1/Resources/TMPDIR";
}
else
{
    print "TMPDIR=$TMPDIR\n";
}
if(! $PICARD)
{
    print "Path to samtools executables is not give and default will be used.\n";
    $samtoolsDir = "/home/shahr2/Software/Samtools/samtools-0.1.17";
}
else
{
    print "SAMTOOLSDIR=$samtoolsDir\n";
}
if(! $PICARD)
{
    print "Path to Picard's executables is not give and default will be used.\n";
    $PICARD = "/ifs/e63data/bergerm1/Resources/Tools/Picard/current/picard-tools-1.79";
}
else
{
    print "PICARD=$PICARD\n";
}
if(! $bwa)
{
    print "Path to BWA MEM executables is not give and default will be used.\n";
    $bwa = "/ifs/e63data/bergerm1/Resources/Tools/common_bin/bwa";
}
else
{
    print "BWA=$bwa\n";
}
if(! $PeSV)
{
    print "Path to PeSV source folder is not give and default will be used.\n";
    $PeSV = "/ifs/e63data/bergerm1/Resources/Tools/PeSVFisher/PeSVFisher-0.93";
}
else
{
    print "PeSV=$PeSV\n";
}
if(! $baitIntervalFile)
{
    print "Bait Interval file is not give and default will be used.\n";
    $baitIntervalFile = "/home/shahr2/References/v5_hg19_picard_baits.interval_list";
}
else
{
    print "BAIT_INTERVAL=$baitIntervalFile\n";
}
if(! $targetIntervalFile)
{
    print "Target Interval file is not give and default will be used.\n";
    $targetIntervalFile = "/home/shahr2/References/v5_hg19_picard_targets.interval_list";
}
else
{
    print "Target_INTERVAL=$targetIntervalFile\n";
}
if(! $refFile)
{
    print "Reference Sequence file is not give and default will be used.\n";
    $refFile = "/home/shahr2/References/Homo_sapiens_assembly19.fasta";
}
else
{
    print "Reference_Sequence=$refFile\n";
}
if(! $chrSepRef)
{
    print "Reference Sequence file is not give and default will be used.\n";
    $chrSepRef = "/ifs/e63data/bergerm1/Resources/PubData/hg19_chrs_seprated";
}
else
{
    print "Chr_Sep_Reference_Folder=$chrSepRef\n";
}
if(! $nprocessors)
{
    print "Number of processors to use not given default value will be used\n";
    $nprocessors = 1;
}
{
    print "Number of processor = $nprocessors\n";
}
if(! $ntchr)
{
    print "Number of chromosomes to use not given default value will be used\n";
    $ntchr = 25;
}
{
    print "Number of chromosomes = $ntchr\n";
}
if(! $JAVA)
{
    print "Path to java executables is not give and default will be used.\n";
    $JAVA = "/usr/java/latest/bin/java";
}
else
{
    print "JAVA=$JAVA\n";
}
if(! $awk)
{
    print "Path to awk executables is not give and default will be used.\n";
    $JAVA = "/usr/bin/awk";
}
else
{
    print "AWK=$awk\n";
}
if(! $dbUser)
{
   print "MySQL db user name not found defult user id will be used\n";
   $dbUser="dmp_prod";
}
else
{
    print "DB User=$dbUser\n";
}
if(! $dbPass)
{
   print "MySQL db user password not found defult user pasword will be used\n";
   $dbPass="TRAIN-734dev";
}
if(! $dbName)
{
   print "MySQL db name not found defult database will be used\n";
   $dbName="rnd_dmp";
}
else
{
    print "DB Name = $dbName\n";
}
if(! $dbLocStr)
{
   print "MySQL db location not found defult database location will be used\n";
   $dbLocStr="unagi.cbio.private";
}
else
{
    print "DB Location = $dbLocStr\n";
}

#Read Sample File
my($fcId,$lane,$sampleId,$sampleRef,$index,$description,$control,$recipe,$operator,$sampleProject) = &ReadSampleFile($sampleFile,$projectName,$outdir);
#Read Title File
my($barcode,$pool,$titleSampleId,$collabId,$patientId,$class,$sampleType,$inputNg,$libraryYeild,$poolInput,$baitVersion) = &ReadTitleFile($titleFile,$outdir);
tie (my %classPerBarcode, 'Tie::IxHash');
for(my $i = 0; $i < scalar(@$barcode); $i++)
{
    #print "$$barcode[$i] => $$class[$i]\n";
    $classPerBarcode{$$barcode[$i] . "_" . $$titleSampleId[$i]} = $$class[$i];
}
# Make a dummy wait file to run after process
&MakeCSH($outdir);
#Split porces to know how many to run
my @allProcess = split (",",$process);
my $numberOfProcess = scalar (@allProcess);

#Call Subroutine according to number of process.
if($numberOfProcess == 1)
{
    &Select1(\@allProcess,$barcodeFile,$titleFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId,$collabId,$patientId,$class,$sampleType,$inputNg,$datadir,$outdir,$fof,$stdNormal);
}
elsif($numberOfProcess == 2)
{
    &Select2(\@allProcess,$barcodeFile,$titleFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId,$collabId,$patientId,$class,$sampleType,$inputNg,$datadir,$outdir,$fof,$stdNormal);
}
elsif($numberOfProcess == 3)
{
    &Select3(\@allProcess,$barcodeFile,$titleFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId,$collabId,$patientId,$class,$sampleType,$inputNg,$datadir,$outdir,$fof,$stdNormal);
}
elsif($numberOfProcess == 4)
{
    &Select4(\@allProcess,$barcodeFile,$titleFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId,$collabId,$patientId,$class,$sampleType,$inputNg,$datadir,$outdir,$fof,$stdNormal);
}
elsif($numberOfProcess == 5)
{
    &Select5(\@allProcess,$barcodeFile,$titleFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId,$collabId,$patientId,$class,$sampleType,$inputNg,$datadir,$outdir,$fof,$stdNormal);
}
elsif($numberOfProcess == 6)
{
    &Select6(\@allProcess,$barcodeFile,$titleFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId,$collabId,$patientId,$class,$sampleType,$inputNg,$datadir,$outdir,$fof,$stdNormal);
}
else
{
    print"Please specify the process number correctly. See Usage\n";
}
#--Calculate total runtime (current time minus start time) 
$now = time - $now;

#--Print runtime 
printf("\n\nTotal running time: %02d:%02d:%02d\n\n", int($now / 3600), int(($now % 3600) / 60), int($now % 60));

print "\n!!!!Done, Thanks for using the script!!!\n";
exit;

#####################################
#####################################
#How to use the script.

sub Usage
{
	print "Unknow option: @_\n" if (@_);

	print "\nUsage : StructuralVariantFinder.pl [options]
        [--config|c                        S Path to configration file(required;Template:/home/shahr2/Scripts/All/StrVarConfig_(script_version)_(bait_version.txt)]
        Inside Config File:
        >Location:
        SampleFile                         S csv file describing details about the sample (required and submit with full path)
        TitleFile                          S tab-delimited title file for the samples (required and submit with full path)
        stdNormal                          S file to be used as standard normal #full path of bam file, *.bai file to be located in same folder)
        list                               S Name of the files with there path (fof:Paired files;one per line;one after another) (optional)
        projectName                        S Name of the project(required,e.g:Colons).
        poolName                           S Name of the pool(required,e.g:Colon5P1).]
        Process                            S 1 => MergingFastq 2 => Run Mapping. 3 => Run HsMetrics. 4 => Call Structural Variants. 5 => Filter Structural Varaints 6 => Annotate Structural Variants\"1,2,3,4,5,6\" for all or in valid combination.]
        datadir                            S Path where all the files to be processed are located. (required)
        outdir                             S Path where all the output files will be written. (required)
        barcodeFile                        S tab-delimited barcode file describing name of the barcode(optional)
        BWA                                S Path to bwa mem program. (optional;default:/ifs/e63data/bergerm1/Resources/Tools/common_bin/bwa)
        PICARD                             S Path to picard tools (optional;default:/ifs/e63data/bergerm1/Resources/Tools/Picard/current/picard-tools-1.79)
        PeSVFisher                         S Path to PeSVfisher source folder.(optional;default:/ifs/e63data/bergerm1/Resources/Tools/PeSVFisher/PeSVFisher-0.93)
        SAMTOOLSDIR                        S Path to samtools source folder.(optional;default:/home/shahr2/Software/Samtools/samtools-0.1.17)
        BaitInterval                       S Path to baits interval file. (optional;default:/home/shahr2/References/v5_hg19_picard_baits.interval_list)
        TargetInterval                     S Path to targets interval file. (optional;default:/home/shahr2/References/v5_hg19_picard_targets.interval_list)
        Reference                          S Path to genome reference file. (optional;default:/home/shahr2/References/Homo_sapiens_assembly19.fasta)
        RefChrSep                          S Path to genome reference chromosome seprated folder (optional;default:/ifs/e63data/bergerm1/Resources/PubData/hg19_chrs_seprated)
        TMPDIR                             S Path to temporary directory. (optional;default:/ifs/e63data/bergerm1/Resources/TMPDIR)
        JAVA                               S Path to java executables. (optional;default:/usr/java/latest/bin/java)
        moveFiles                          I 2 => Skip Making Folders & Moving Files. 1 => Make Folders and Move Files. (default:1,optional)
        NumberOfChromosomes                I Number of chromosomes on which analysis needs to be done.eg: 25 => 1-22,X,Y,MT (default:25,optional)
        NumberOfProcessors                 I Number of processors to use for analysis (default:1,optional)
        FOF                                S Name of the files with there path. (fof:Paired files;one per line;one after another) (optional)
        databaseUserId                     S MySQL DB user id. (optional;default:dmp_prod)
        databasePassword                   S MySQL DB password. (optional;default:password of dmp_prod)
        databaseName                       S MySQL DB name. (optional;default:rnd_dmp)
        databaseLocationString             S MySQL DB location. (optional;default:unagi.cbio.private)
        awk                                S Location to awk executables. (optional;default:/user/bin/awk).
	\n";

	print "For any question please contact Ronak Shah (shahr2\@mskcc.org)\n";
	print "!~!~! Comments, flames and suggestion are welcome !~!~!\n";
	exit;
}

###################################################
###################################################
#--GET CONFIGRATION DETAIL
sub GetConfigration
 {
     my($config_file,$outdir) = @_;
     my @data = ();
     tie (my %location,'Tie::IxHash');
     tie (my %version,'Tie::IxHash');
     open(CONFIG,"$config_file") || die "Cannot open CONFIGFile: $config_file,$!\n";
     while(my $text = <CONFIG>)
     {
	 next if($text =~ /^#/);
	 chomp($text);
	 #Location
	 if($text =~ />Location/)
	 {
	     while(my $locLine = <CONFIG>)
	     {
		#Version
		if($locLine =~ />Version/)
		{
		    while(my $verLine = <CONFIG>)
		    {
			@data = split("=",$verLine);
			$data[0] =~ s/\s//g;
			$data[1] =~ s/\s//g;
			$version{$data[0]} = $data[1];
		    }
		}
		else
		{
		    @data = split(":",$locLine);
		    $data[0] =~ s/\s//g;
		    $data[1] =~ s/\s//g;
		    $location{$data[0]} = $data[1];
		}
	    }
	 }
     }
     close(CONFIG);
     return(\%location,\%version);
}

#####################################
#####################################
#Run only single process
sub Select1
{
    my($allProcess,$barcodeFile,$titleFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId,$collabId,$patientId,$class,$sampleType,$inputNg,$datadir,$outdir,$fof,$stdNormal) = @_;
    #print "F:$fof\n";
    my @process = @$allProcess;
    my $parseFilenames = "";
    if($process[0] == 1)
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
    }
    elsif($process[0] == 2)
    {
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof);
    }
    elsif($process[0] == 3)
    {
        $parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif($process[0] == 4)
    {
        $parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof);
    } 
    elsif($process[0] == 5)
    {
        $parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif($process[0] == 6)
    {
        $parseFilenames = &AnnotateStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    } 
    else
    {
        print"The process number entered does not exists, See Usage.\n";
        exit;
    }
    return;

}

#####################################
#####################################
#Run two process
sub Select2
{
    my ($allProcess,$barcodeFile,$titleFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId,$collabId,$patientId,$class,$sampleType,$inputNg,$datadir,$outdir,$fof,$stdNormal) = @_;
    my @process = @$allProcess;
    my $parseFilenames = "";
    if(($process[0] == 1) and ($process[1] == 2))
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof);
    }
    elsif(($process[0] == 1) and ($process[1] == 3))
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 1) and ($process[1] == 4))
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 1) and ($process[1] == 5))
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 1) and ($process[1] == 6))
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &AnnotateStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 2) and ($process[1] == 3))
    {
	$parseFilenames = &DoMapping($parseFilenames,$outdir,$fof);
        $parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 2) and ($process[1] == 4))
    {
	$parseFilenames = &DoMapping($parseFilenames,$outdir,$fof);
        $parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 2) and ($process[1] == 5))
    {
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 2) and ($process[1] == 6))
    {
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &AnnotateStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 3) and ($process[1] == 4))
    {
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 3) and ($process[1] == 5))
    {
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 3) and ($process[1] == 6))
    {
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &AnnotateStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 4) and ($process[1] == 5))
    {
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 4) and ($process[1] == 6))
    {
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof);
	$parseFilenames = &AnnotateStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 5) and ($process[1] == 6))
    {
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &AnnotateStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    else
    {
	print "Please See Usage and Input proper Options\n";
	Usage();
	exit;
    }
    return;
}
#####################################
#####################################
#Run three process
sub Select3
{
    my ($allProcess,$barcodeFile,$titleFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId,$collabId,$patientId,$class,$sampleType,$inputNg,$datadir,$outdir,$fof,$stdNormal) = @_;

    my @process = @$allProcess;
    my $parseFilenames = "";
    if(($process[0] == 1) and ($process[1] == 2) and ($process[2] == 3))
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof);
        $parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 1) and ($process[1] == 2) and ($process[2] == 4))
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 1) and ($process[1] == 2) and ($process[2] == 5))
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 1) and ($process[1] == 2) and ($process[2] == 6))
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &AnnotateStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 1) and ($process[1] == 3) and ($process[2] == 4))
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof);
    } 
    elsif(($process[0] == 1) and ($process[1] == 3) and ($process[2] == 5))
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 1) and ($process[1] == 3) and ($process[2] == 6))
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &AnnotateStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    } 
    elsif(($process[0] == 1) and ($process[1] == 4) and ($process[2] == 5))
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 1) and ($process[1] == 4) and ($process[2] == 6))
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &AnnotateStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 1) and ($process[1] == 5) and ($process[2] == 6))
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &AnnotateStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 2) and ($process[1] == 3) and ($process[2] == 4))
    {
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof);
    } 
    elsif(($process[0] == 2) and ($process[1] == 3) and ($process[2] == 5))
    {
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 2) and ($process[1] == 3) and ($process[2] == 6))
    {
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &AnnotateStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    } 
    elsif(($process[0] == 2) and ($process[1] == 4) and ($process[2] == 5))
    {
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    } 
    elsif(($process[0] == 2) and ($process[1] == 4) and ($process[2] == 6))
    {
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    } 
    elsif(($process[0] == 3) and ($process[1] == 4) and ($process[2] == 5))
    {
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &AnnotateStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    } 
    elsif(($process[0] == 3) and ($process[1] == 4) and ($process[2] == 6))
    {
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &AnnotateStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    } 
    elsif(($process[0] == 4) and ($process[1] == 5) and ($process[2] == 6))
    {
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &AnnotateStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    else
    {
	print "Please See Usage and Input proper Options\n";
	Usage();
	exit;
    }
    return;
}
#####################################
#####################################
#Run four process
sub Select4
{
    my ($allProcess,$barcodeFile,$titleFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId,$collabId,$patientId,$class,$sampleType,$inputNg,$datadir,$outdir,$fof,$stdNormal) = @_;

    my @process = @$allProcess;
    my $parseFilenames = "";

    if(($process[0] == 1) and ($process[1] == 2) and ($process[2] == 3) and ($process[3] == 4))
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
    }
    elsif(($process[0] == 1) and ($process[1] == 2) and ($process[2] == 3) and ($process[3] == 5))
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 1) and ($process[1] == 2) and ($process[2] == 4) and ($process[3] == 5))
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 1) and ($process[1] == 2) and ($process[2] == 4) and ($process[3] == 6))
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &AnnotateStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 1) and ($process[1] == 3) and ($process[2] == 4) and ($process[3] == 5))
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 1) and ($process[1] == 3) and ($process[2] == 4) and ($process[3] == 6))
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &AnnotateStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 1) and ($process[1] == 3) and ($process[2] == 5) and ($process[3] == 6))
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &AnnotateStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 1) and ($process[1] == 4) and ($process[2] == 5) and ($process[3] == 6))
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &AnnotateStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 2) and ($process[1] == 3) and ($process[2] == 4) and ($process[3] == 5))
    {
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 2) and ($process[1] == 3) and ($process[2] == 5) and ($process[3] == 6))
    {
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &AnnotateStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 2) and ($process[1] == 4) and ($process[2] == 5) and ($process[3] == 6))
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &AnnotateStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 3) and ($process[1] == 4) and ($process[2] == 5) and ($process[3] == 6))
    {
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &AnnotateStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    else
    {
	print "Please See Usage and Input proper Options\n";
	Usage();
	exit;
    }
    return;
}
#####################################
#####################################
#Run five process
sub Select5
{
    my ($allProcess,$barcodeFile,$titleFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId,$collabId,$patientId,$class,$sampleType,$inputNg,$datadir,$outdir,$fof,$stdNormal) = @_;

    my @process = @$allProcess;
    my $parseFilenames = "";
    if(($process[0] == 1) and ($process[1] == 2) and ($process[2] == 3) and ($process[3] == 4) and ($process[4] == 5))
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof); 
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 1) and ($process[1] == 2) and ($process[2] == 3) and ($process[3] == 4) and ($process[4] == 6))
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof);
	$parseFilenames = &AnnotateStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    elsif(($process[0] == 2) and ($process[1] == 3) and ($process[2] == 4) and ($process[3] == 5) and ($process[4] == 6))
    {
	$parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof);
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &AnnotateStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    else
    {
	print "Please See Usage and Input proper Options\n";
	Usage();
	exit;
    }
    return;
}
#####################################
#####################################
#Run six process
sub Select6
{
    my ($allProcess,$barcodeFile,$titleFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId,$collabId,$patientId,$class,$sampleType,$inputNg,$datadir,$outdir,$fof,$stdNormal) = @_;

    my @process = @$allProcess;
    my $parseFilenames = "";
    if(($process[0] == 1) and ($process[1] == 2) and ($process[2] == 3) and ($process[3] == 4) and ($process[4] == 5) and ($process[5] == 6))
    {
	$parseFilenames = &MergeDataFromDirectory($datadir,$outdir,$barcodeFile,$lane,$sampleId,$index,$barcode,$pool,$titleSampleId);
        $parseFilenames = &DoMapping($parseFilenames,$outdir,$fof); 
	$parseFilenames = &CalcHsMetrics($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &CallStructuralVariants($parseFilenames,$stdNormal,$patientId,$datadir,$outdir,$fof);
	$parseFilenames = &FilterStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
	$parseFilenames = &AnnotateStructuralVariants($parseFilenames,$datadir,$outdir,$fof);
    }
    else
    {
	print "Please See Usage and Input proper Options\n";
	Usage();
	exit;
    }
    return;
}
###################################################
###################################################
#--Make Notification file
sub MakeCSH
 {
        my($outdir) = @_;
        my $filename = $outdir . "/Notify.csh";
	if(! -e $filename)
	{
	    my $ntmp = new IO::File(">$filename");
	    print $ntmp "#!/bin/csh\n";
	    print $ntmp "#Notification File\n";
	    print $ntmp "echo"," This is Done","\n";
	    $ntmp->close();
	    `chmod +x $filename`;
	}
	else
	{
	    print "$filename exists and wont be created.\n";
	}
        return;
}
###################################################
###################################################
#--Check Error Files

sub CheckErrorFiles
{
	my($outdir,@errorfilenames) = @_;
	print "\nChecking Error Files.\n";
	foreach my $efile(@errorfilenames)
	{
	    next if($efile eq "NULL");
	    if(-e "$outdir/$efile")
	    {
		if((-s "$outdir/$efile") != 0 )
		{
		    print "Please Check $efile; Something went wrong here.\n";
		}
		else
		{
		   `rm $outdir/$efile`;
		    next;
		}
	    }
	    else
	    {
		print "Please check, no error file was created named $efile; Something went wrong here. \n";
	    }
	}
	return;
}
###################################################
###################################################
#--Waiting for the process to finish

sub WaitToFinish
{
	my($outdir,@waitfilenames) = @_;
	print "Waiting for the Process to finish...\n";
	foreach my $wfile(@waitfilenames)
	{
	    next if($wfile eq "NULL");
	    wait while(! -e "$outdir/$wfile");
	    #print "$outdir/$wfile\n";
	    while(-e "$outdir/$wfile")
	    {
		#print "$outdir/$wfile\n";
		open(FH,"<","$outdir/$wfile");
		while(<FH>)
		{
		    if($_ =~ /This is Done/ig)
		    {
			#print "\nFinished: $wfile\n";
			last;
		    }
		    else
		    {
			wait;
		    }
		}
		last;
		}
	    close(FH);
	}
	foreach my $wfile(@waitfilenames)
	{
	   next if($wfile eq "NULL");
	   `rm $outdir/$wfile`;
	}
	return;
}
###################################################
###################################################
#--Make array of file of files list from the outdir
 sub GetNames
 {
 	my($fof,$outdir) = @_;
	my (@filenames) = ();
	open(FOF,"$outdir/$fof") || die "Cannot open ListFile: \"$outdir/$fof\"\n";
	while(<FOF>)
	{
		$_ =~ s/\s//g;
		my $filename =  pop @{[split("/",$_)]};
		push(@filenames,$filename);	
	}
	return(@filenames);
}

###################################################
###################################################
#--Make Pairs of the files.

sub MAKEPAIRS
{
	my($filenames,$outdir) = @_;
	my @names = @$filenames;
	my $count = scalar @names;
	my (@newnames) = ();
	if($count%2 != 0)
	{
		print STDERR "\nOdd number of files given, please check Input file.\n";
		exit;
	}
	else
	{
		for(my $i =0; $i < scalar (@names); $i+=2)
		{
			chomp($names[$i]);
			chomp($names[$i+1]);
			push(@newnames,"$names[$i],$names[$i+1]");
		}
	}
	return(@newnames);
}
#####################################
#####################################
#Read data related to samples as well as barcodes.

sub ReadSampleFile
{
    my($sampleFile,$projectName,$outdir) = @_;
    my (@fcId,@lane,@sampleId,@sampleRef,@index,@description,@control,@recipe,@operator,@sampleProject) = ();
    my $sampleFileName = "";
    if($sampleFile =~ /\//)
    {
	$sampleFileName = pop @{[split("/",$sampleFile)]};
    }
    else
    {
	$sampleFileName = $sampleFile;
    }
    open(SAMPLEFILE, $sampleFile) || die "Cannot open SAMPLEFILE:$sampleFile,$!\n";
    while(<SAMPLEFILE>)
    {
	next if $. == 1;
	my @dataCols = split(",",$_);
	if($dataCols[0]){push (@fcId,$dataCols[0]);}
	if($dataCols[1]){push (@lane,$dataCols[1]);}
	if($dataCols[2]){push (@sampleId,$dataCols[2]);}
	if($dataCols[3]){push (@sampleRef,$dataCols[3]);}
	if($dataCols[4]){push (@index,$dataCols[4]);}
        if($dataCols[5]){push (@description,$dataCols[5]);}
        if($dataCols[6]){push (@control,$dataCols[6]);}
        if($dataCols[7]){push (@recipe,$dataCols[7]);}
	if($dataCols[8]){push (@operator,$dataCols[8])}
	if($dataCols[9]){push (@sampleProject,$dataCols[8]);}
    }
    close(SAMPLEFILE);
    if(! -e "$outdir/$sampleFileName")
    {
	`cp $sampleFile $outdir/$sampleFileName`;
    }
    return(\@fcId,\@lane,\@sampleId,\@sampleRef,\@index,\@description,\@control,\@recipe,\@operator,\@sampleProject);
}

#####################################
#####################################
#Read data related to samples as well as barcodes from title file.

sub ReadTitleFile
{
    my($titleFile,$outdir) = @_;
    my @barcode = ();
    my @pool = ();
    my @sampleId = ();
    my @collabId = ();
    my @patientId = ();
    my @class = ();
    my @sampleType = ();
    my @inputNg = ();
    my @libraryYeild = ();
    my @poolInput = ();
    my @baitVersion = ();
    my @fof = ();
    my @newfof = ();

    open(TFH,$titleFile)||die"Cannot open file TitleFile:$titleFile, $!\n";
    while(<TFH>)
    {
	next if($. == 1);
	my @dataCols = split("\t",$_);
	my @newDatacols = grep(s/\s*$//g, @dataCols);#remove whitespace if any
	push(@barcode,$newDatacols[0]);
	push(@pool,$newDatacols[1]);
	push(@sampleId,$newDatacols[2]);
	push(@collabId,$newDatacols[3]);
	push(@patientId,$newDatacols[4]);
	push(@class,$newDatacols[5]);
	push(@sampleType,$newDatacols[6]);
	push(@inputNg,$newDatacols[7]);
	push(@libraryYeild,$newDatacols[8]);
	push(@poolInput,$newDatacols[9]);
	push(@baitVersion,$newDatacols[10]);
    }
    close(TFH);
    my $poolName = $pool[0];
    my $newtitleFileName = $poolName . "_title.txt";
    if(! -e "$outdir/$newtitleFileName")
    {
	`cp $titleFile $outdir/$newtitleFileName`;
    }
    return(\@barcode,\@pool,\@sampleId,\@collabId,\@patientId,\@class,\@sampleType,\@inputNg,\@libraryYeild,\@poolInput,\@baitVersion);

}
#####################################
#####################################
#sort by barcode name:

sub lowestNumber
{
    my $files = shift;
    my @filenames = split(",",$files);
    my ($number) = $filenames[0] =~ m/.*_bc(\d{1,2})_.*/g;
    return $number;
}
#####################################
#####################################
#Merge data from reading data from the directory

sub MergeDataFromDirectory
{
    my($datadir,$outdir,$barcodeFile,$lane,$sampleId,$indexSeq,$titleBarcode,$titlePool,$titleSampleId) = @_;

    my @lane = @$lane;
    my @sampleId = @$sampleId;
    my @index = @$indexSeq;
    my @titleBarcode = @$titleBarcode;
    my @titlePool = @$titlePool;
    my @titleSampleId = @$titleSampleId;
    my %barcodes = ();
    my %indexHash = ();
    my %titleInfo = ();
    my @notifyNames = ();
    my @checkErrorFiles = ();
    my $newIndex;
    my $name;
    my $Null = "NULL";
    my @parseFilenames = ();
    my $now = time;
    open(BARCODEFILE, $barcodeFile) || die "Cannot open BARCODEFILE:$barcodeFile,$!\n";
    while(<BARCODEFILE>)
    {
	next if ($. == 1);
	my @dataCols = split ("\t",$_);
	$dataCols[0] =~ s/\s//g; 
	$dataCols[1] =~ s/\s//g;
	$barcodes{$dataCols[0]} = $dataCols[1];
	$indexHash{$dataCols[1]} = $dataCols[0];
    }
    close(BARCODEFILE);

    print "Running merge jobs on SGE at ". localtime() ."\n";

    if($fastqSource eq "Path")
    {
	foreach my $i (0 .. $#titleBarcode)
	{
	    my $newIndex = $titleBarcode[$i];
	    my $name = $titleSampleId[$i] . "_" . $titleBarcode[$i] ."_". $titlePool[$i];
	    my $read1ListName = "";
	    my $read2ListName = "";
	    foreach my $j (0 .. $#sampleId)
	    {
		if($index[$j] eq $indexHash{$titleBarcode[$i]}){
		    $read1ListName .= $datadir . "/" . $sampleId[$j] . "_" . $index[$j] . "_L00" . $lane[$j] . "_R1_001.fastq.gz ";
		    $read2ListName .= $datadir . "/" . $sampleId[$j] . "_" . $index[$j] . "_L00" . $lane[$j] . "_R2_001.fastq.gz ";
		}
	    }
	    my $read1Name = $outdir . "/" . $name . "_L000_R1_mrg.fastq.gz";
	    my $read2Name = $outdir . "/" .  $name . "_L000_R2_mrg.fastq.gz";
	    if((-e $read1Name) and ((-s $read1Name) != 0) and (-e $read2Name) and ((-s $read2Name) !=0))
	    {
		print "Files:\n$read1Name\n$read2Name\n they exists and process will not run to merge files.\n";
		push(@notifyNames,$Null);
		push(@checkErrorFiles,$Null);
		push(@notifyNames, $Null);
		push(@checkErrorFiles,$Null);
		push(@parseFilenames,"$read1Name,$read2Name");
	    }
	    else
	    {
		#Read1
		`qsub -q all.q -V -wd $outdir -N MergeRead1.$newIndex.$i.$$ -l h_vmem=8G,virtual_free=8G -pe smp 1 -e MergeRead1.$newIndex.$i.$$.err -o /dev/null -b y "/bin/zcat $read1ListName | gzip > $read1Name"`;
		`qsub -q all.q -V -wd $outdir -hold_jid MergeRead1.$newIndex.$i.$$ -N NotifyMR.Read1.$i.$$ -l h_vmem=2G,virtual_free=2G -pe smp 1 -e /dev/null -o NotifyMR.Read1.$i.$$.stat -b y "$outdir/Notify.csh"`;
		#Read2
		`qsub -q all.q -V -wd $outdir -N MergeRead2.$newIndex.$i.$$ -l h_vmem=8G,virtual_free=8G -pe smp 1 -e MergeRead2.$newIndex.$i.$$.err -o /dev/null -b y "/bin/zcat $read2ListName | gzip > $read2Name"`;
		`qsub -q all.q -V -wd $outdir -hold_jid MergeRead2.$newIndex.$i.$$ -N NotifyMR.Read2.$i.$$ -l h_vmem=2G,virtual_free=2G -pe smp 1 -e /dev/null -o NotifyMR.Read2.$i.$$.stat -b y "$outdir/Notify.csh"`;
		push(@notifyNames, "NotifyMR.Read1.$i.$$.stat");
		push(@checkErrorFiles,"MergeRead1.$newIndex.$i.$$.err");
		push(@notifyNames, "NotifyMR.Read2.$i.$$.stat");
		push(@checkErrorFiles,"MergeRead2.$newIndex.$i.$$.err");
		push(@parseFilenames,"$read1Name,$read2Name");
	    }
	}
	&WaitToFinish($outdir,@notifyNames);
	&CheckErrorFiles($outdir,@checkErrorFiles);
	$now = time - $now;
	print "Finished running merge jobs on SGE at ". localtime() ."\n";
	printf("Total running time: %02d:%02d:%02d\n\n", int($now / 3600), int(($now % 3600) / 60), int($now % 60));
	my(@sortedparseFilenames) = sort {lowestNumber($a) <=>  lowestNumber($b)} @parseFilenames;
	return(\@sortedparseFilenames);

    }
    else
    {

	for(my $i=0; $i < scalar(@titleBarcode); $i++)
	{
	    $titleInfo{$titleBarcode[$i]} = $titleSampleId[$i] . "_" . $titleBarcode[$i] ."_". $titlePool[$i];
	}

	for(my $sampleNum = 0; $sampleNum < scalar(@sampleId); $sampleNum++)
	{
	    my $read1ListName = $datadir . "/" . $sampleId[$sampleNum] . "_" . $index[$sampleNum] . "_L00" . $lane[$sampleNum] . "_R1_*.fastq.gz";
	    my $read2ListName = $datadir . "/" . $sampleId[$sampleNum] . "_" . $index[$sampleNum] . "_L00" . $lane[$sampleNum] . "_R2_*.fastq.gz";

	    if(exists $barcodes{$index[$sampleNum]})
	    {
		$newIndex = $barcodes{$index[$sampleNum]};
		if(exists $titleInfo{$newIndex})
		{
		    $name = $titleInfo{$newIndex};
		}
		else
		{
		    print "The barcode $newIndex doesnot exists in the title file. Cannot move ahead. Please check and rerun.\n";
		    exit;
		}
	    }
	    else
	    {
		print "The barcode sequence $barcodes{$index[$sampleNum]} does not exists in barcode file. Cannot move ahead. Please check and rerun.\n";
		exit;
	    }
	    my $read1Name = $outdir . "/" . $name . "_L00" . $lane[$sampleNum] . "_R1_mrg.fastq.gz";
	    my $read2Name = $outdir . "/" .  $name . "_L00" . $lane[$sampleNum] . "_R2_mrg.fastq.gz";
            #Run the qsub command to merge the files.
	    #Read1
	    if((-e $read1Name) and ((-s $read1Name) != 0) and (-e $read2Name) and ((-s $read2Name) !=0))
	    {
		print "Files:\n$read1Name\n$read2Name\n they exists and process will not run to merge files.\n";
		push(@notifyNames,$Null);
		push(@checkErrorFiles,$Null);
		push(@notifyNames, $Null);
		push(@checkErrorFiles,$Null);
		push(@parseFilenames,"$read1Name,$read2Name");
		next;
	    }
	    else
	    {
		#Read1
		`qsub -q all.q -V -wd $outdir -N MergeRead1.$newIndex.$sampleNum.$$ -l h_vmem=8G,virtual_free=8G -pe smp 1 -e MergeRead1.$newIndex.$sampleNum.$$.err -o /dev/null -b y "/bin/zcat $read1ListName | gzip > $read1Name"`;
		`qsub -q all.q -V -wd $outdir -hold_jid MergeRead1.$newIndex.$sampleNum.$$ -N NotifyMR.Read1.$sampleNum.$$ -l h_vmem=2G,virtual_free=2G -pe smp 1 -e /dev/null -o NotifyMR.Read1.$sampleNum.$$.stat -b y "$outdir/Notify.csh"`;
		#Read2
		`qsub -q all.q -V -wd $outdir -N MergeRead2.$newIndex.$sampleNum.$$ -l h_vmem=8G,virtual_free=8G -pe smp 1 -e MergeRead2.$newIndex.$sampleNum.$$.err -o /dev/null -b y "/bin/zcat $read2ListName | gzip > $read2Name"`;
		`qsub -q all.q -V -wd $outdir -hold_jid MergeRead2.$newIndex.$sampleNum.$$ -N NotifyMR.Read2.$sampleNum.$$ -l h_vmem=2G,virtual_free=2G -pe smp 1 -e /dev/null -o NotifyMR.Read2.$sampleNum.$$.stat -b y "$outdir/Notify.csh"`;
		push(@notifyNames, "NotifyMR.Read1.$sampleNum.$$.stat");
		push(@checkErrorFiles,"MergeRead1.$newIndex.$sampleNum.$$.err");
		push(@notifyNames, "NotifyMR.Read2.$sampleNum.$$.stat");
		push(@checkErrorFiles,"MergeRead2.$newIndex.$sampleNum.$$.err");
		push(@parseFilenames,"$read1Name,$read2Name");
	    }
	}
	&WaitToFinish($outdir,@notifyNames);
	&CheckErrorFiles($outdir,@checkErrorFiles);
	$now = time - $now;
	print "Finished running merge jobs on SGE at ". localtime() ."\n";
	printf("Total running time: %02d:%02d:%02d\n\n", int($now / 3600), int(($now % 3600) / 60), int($now % 60));
	my(@sortedparseFilenames) = sort {lowestNumber($a) <=>  lowestNumber($b)} @parseFilenames;
	return(\@sortedparseFilenames);
    }
}

#####################################
#####################################
#Do Mapping which includes:
#Run BWA mem
#Sort SAM

sub DoMapping
{
    my ($filenames,$outdir,$fof) = @_;
    my (@names,@notifyNames,@SAFilenames,@SamFilenames,@sortedBamFilenames,@MarkDuplicatesBamFilenames) = ();
    if($filenames){(@names) = @$filenames;} 
    if((scalar(@names) == 0)and($fof))
    {
	my @fnames = &GetNames($fof,$outdir);
	@names = &MAKEPAIRS(\@fnames,$outdir);
    }
    if ($bwaFlag == 1) {
    
        my(@sortedparseFilenames) = sort {lowestNumber($a) <=>  lowestNumber($b)} @names;
        @names = @sortedparseFilenames;
        for(my $i = 0; $i < scalar(@names); $i++){
            my($file1,$file2) = split(",",$names[$i]);
            #print "$file1\n$file2\n";
            my($SA1Filename,$notifyname1) = &RunBwaAln($file1,$outdir,$i,"read1");
            my($SA2Filename,$notifyname2) = &RunBwaAln($file2,$outdir,$i,"read2");
            push(@notifyNames,$notifyname1);
            push(@notifyNames,$notifyname2);
            push(@SAFilenames,"$SA1Filename,$SA2Filename");
        }
        #waiting for bwa aln to finish
        &WaitToFinish($outdir,@notifyNames);
        print "Finished running bwa aln jobs on SGE\n";
        $now = time - $now;
        printf("Total BWA ALN run time: %02d:%02d:%02d\n\n", int($now / 3600), int(($now % 3600) / 60), int($now % 60));
   
        #Running BWA sampe
        $now = time;
        print "Started runing bwa sampe jobs on SGE\n";
        @notifyNames = ();
        for(my $i = 0; $i < scalar(@SAFilenames); $i++)
            {
                my($clippedfile1,$clippedfile2) = split(",",$names[$i]);
                my($SAfile1,$SAfile2) = split(",",$SAFilenames[$i]);
                #print "$SAfile1\n$SAfile2\n";
                my ($samFile,$notifyname) = &RunBwaSampe($clippedfile1,$clippedfile2,$SAfile1,$SAfile2,$outdir,$i);
                push(@notifyNames,$notifyname);
                $samFile = basename($samFile);
                push(@SamFilenames,$samFile);
            }
        #waiting for bwa sampe to finish
        &WaitToFinish($outdir,@notifyNames);
        print "Finished running bwa sampe jobs on SGE\n";
        $now = time - $now;
        printf("Total BWA SAMPE run time: %02d:%02d:%02d\n\n", int($now / 3600), int(($now % 3600) / 60), int($now % 60));
    }
    
    if ( $bwaFlag == 2) {
        
        #Running BwaMem
        $now = time;
        print "Started runing bwa mem jobs on SGE at ". localtime() ."\n";
        @notifyNames = ();
        for(my $i = 0; $i < scalar(@names); $i++)
            {
                my($file1,$file2) = split(",",$names[$i]);
                #print "$file1\n$file2\n";
                my($samFilename,$notifyname) = &RunBwaMem($file1,$file2,$outdir,$i);
                push(@notifyNames,$notifyname);
                push(@SamFilenames,"$samFilename");
            }
        #waiting for bwa aln to finish
        &WaitToFinish($outdir,@notifyNames);
        print "Finished running bwa mem jobs on SGE at ". localtime() ."\n";
        $now = time - $now;
        printf("Total running time: %02d:%02d:%02d\n\n", int($now / 3600), int(($now % 3600) / 60), int($now % 60));
    }
    
    
    #Run Sort Sam
    $now = time;
    print "Started running Sort Sam jobs on SGE at ". localtime() ."\n";
    @notifyNames = ();
    for(my $i = 0; $i < scalar(@SamFilenames); $i++)
    {
	my($sortedBamFile,$notifyname) = &RunSortSam($SamFilenames[$i],$outdir,$i);
	push(@notifyNames,$notifyname);
	push(@sortedBamFilenames,$sortedBamFile);
    }
    #waiting for sort sam to finish
    &WaitToFinish($outdir,@notifyNames);
    $now = time - $now;
    print "Finished running Sort Sam jobs on SGE at ". localtime() ."\n"; 
    #Run MarkDuplicates
    $now = time;
    print "Started running Mark Duplicates jobs on SGE at ". localtime() ."\n";
    @notifyNames = ();
    for(my $i = 0; $i < scalar(@sortedBamFilenames); $i++)
    {
	my($MarkDuplicatesBamFile,$notifyname) = &RunMarkDuplicates($sortedBamFilenames[$i],$outdir,$i);
	push(@notifyNames,$notifyname);
	push(@MarkDuplicatesBamFilenames,$MarkDuplicatesBamFile);
    }
    #waiting for Mark Duplicates to finish
    &WaitToFinish($outdir,@notifyNames);
    $now = time - $now;
    print "Finished running Mark Duplicates jobs on SGE at ". localtime() ."\n"; 
    printf("Total running time: %02d:%02d:%02d\n\n", int($now / 3600), int(($now % 3600) / 60), int($now % 60));
    return (\@MarkDuplicatesBamFilenames);
}
#####################################
#####################################
#BWA Find Suffix Array(SA) Co-ordinates

sub RunBwaAln
{
   my($file,$outdir,$id,$readType) = @_;
   my($basename) = $file =~ /(.*)\.fastq.gz/;
   my($outFilename) = "$basename" . ".sai";
   if((-e "$outFilename") and ((-s "$outFilename") != 0))
   {
       print"Files:\n$outFilename\n they exists and process will not run to make \".sai\" file.\n";
       return("$outFilename",'NULL');
   }
   else
   {
       eval
       {
	   `qsub -q all.q -V -wd $outdir -N bwaAln.$readType.$id.$$ -o bwaAln.$readType.$id.$$.stdout -e bwaAln.$readType.$id.$$.stderr -l h_vmem=3G,virtual_free=3G -pe smp 5 -b y "$bwa aln -t $nprocessors -l 40 -k 2 -f $outFilename $refFile $file"`;
	   `qsub -q all.q -V -wd $outdir -hold_jid bwaAln.$readType.$id.$$ -N NotifyBwaAln.$id.$$ -e NotifyBwaAln.$readType.$id.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1 -o NotifyBwaAln.$readType.$id.$$.stat -b y "$outdir/Notify.csh"`;
       };
       if($@)
       {
	  print "BWA ALN:Job Submission Failed, Error:$@\n";
	  exit(1);
       }
   }
    return("$outFilename","NotifyBwaAln.$readType.$id.$$.stat");
}
#####################################
#####################################
#BWA SAMPE align SA.

sub RunBwaSampe
{
   my($clippedfile1,$clippedfile2,$SAfile1,$SAfile2,$outdir,$id) = @_;
   my($basename) = $clippedfile1 =~ /(.*)_R1.*\.fastq.gz/;
   my $outFilename = "$basename" . "_mrg_cl_aln.sam";
   if($basename =~ /\//)
   {
       $basename = basename($basename);
   }
   my @sampleDetails = split("_bc",$basename);
   my $sampleId = $sampleDetails[0];
   my ($barcode) = $basename =~ /.*_(bc\d{1,2})_.*/;
   my ($pool) = $basename =~ /.*_bc\d{1,2}_(.*)_L\d{1,3}_.*/;
    if((-e "$outFilename") and ((-s "$outFilename") != 0))
   {
       print "Files:\n$outFilename\n they exists and process will not run to make \"_aln.sam\" file.\n";
       return("$outFilename",'NULL');
   }
   else
   {
       eval
       {
	   #`qsub -q $queue -V -wd $outdir -N bwaSampe.$id.$$ -o bwaSampe.$id.$$.stdout -e bwaSampe.$id.$$.stderr -b y "$BWA sampe -r \'\@RG\tID:$basename\tLB:$id\tSM:$sampleId\tPL:Illumina\tPU:$barcode\tCN:BergerLab_MSKCC\' -f $outFilename $Reference $SAfile1 $SAfile2 $clippedfile1 $clippedfile2"`;
	   `qsub -q all.q -V -wd $outdir -N bwaSampe.$id.$$ -o bwaSampe.$id.$$.stdout -e bwaSampe.$id.$$.stderr -l h_vmem=6G,virtual_free=6G -pe smp 1 -b y "$bwa sampe -N 100 -f $outFilename $refFile $SAfile1 $SAfile2 $clippedfile1 $clippedfile2"`;
	   `qsub -q all.q -V -wd $outdir -hold_jid bwaSampe.$id.$$ -N NotifyBwaSampe.$id.$$ -e NotifyBwaSampe.$id.$$.stderr -l h_vmem=2G,virtual_free=2G -pe smp 1 -o NotifyBwaSampe.$id.$$.stat -b y "$outdir/Notify.csh"`;
       };
       if($@)
       {
	   print "BWA SAMPE:Job Submission Failed, Error:$@\n";
	   exit(1);
       }
   }
   return("$outFilename","NotifyBwaSampe.$id.$$.stat");
}

#####################################
#####################################
#BWA MEM to align fastq.

sub RunBwaMem{
    my($fastq1,$fastq2,$outdir,$id) = @_;
    my($basename) = $fastq1 =~ /(.*)_R1.*\.fastq.gz/;
    my $outFilename = "$basename" . "_mrg_cl_aln.sam";
    if ($basename =~ /\//) {
        $basename = basename($basename);
    }
    my @sampleDetails = split("_bc",$basename);
    my $sampleId = $sampleDetails[0];
    my ($barcode) = $basename =~ /.*_(bc\d{1,2})_.*/;
    my ($pool) = $basename =~ /.*bc\d{1,2}_(.*)_L\d{1,3}_.*/;
    if ((-e "$outFilename") and ((-s "$outFilename") != 0)) {
        print "Files:\n$outFilename\n they exists and process will not run to make \"_aln.sam\" file.\n";
        return("$outFilename",'NULL');
    } else {
        `qsub -q all.q -wd $outdir -N bwaMem.$id.$$ -l h_vmem=6G,virtual_free=6G -pe smp $nprocessors -o $outFilename -e /dev/null -b y "$bwa mem -t 4 -PM -R \'\@RG\tID:$basename\tLB:$id\tSM:$sampleId\tPL:Illumina\tPU:$barcode\tCN:BergerLab_MSKCC\' $refFile $fastq1 $fastq2"`;
        `qsub -q all.q -V -wd $outdir -hold_jid bwaMem.$id.$$ -N NotifyBwaMem.$id.$$ -l h_vmem=2G,virtual_free=2G -pe smp 1 -e /dev/null -o NotifyBwaMem.$id.$$.stat -b y "$outdir/Notify.csh"`;
    }
    return("$outFilename","NotifyBwaMem.$id.$$.stat");
}

#####################################
#####################################
#Sort Sam file

sub RunSortSam
{
   my($samFile,$outdir,$id) = @_; 
   my($basename) = $samFile =~ /(.*)_mrg_cl_aln.sam/;
   my $outFilename = $samFile;
   $outFilename =~ s/\.sam/_srt\.bam/g;
   my @sampleDetails = split("_bc",$basename);
   my $sampleId = $sampleDetails[0];
   my ($barcode) = $basename =~ /.*_(bc\d{1,2})_.*/;
   my ($pool) = $basename =~ /.*_bc\d{1,2}_(.*)_L\d{1,3}_.*/;
   my $platform = "Illumina";
   if((-e "$outFilename") and ((-s "$outFilename") != 0))
   {
       print "Files:\n$outFilename\n they exists and process will not run to make \"_srt.bam\" file.\n";
       return("$outFilename",'NULL');
   }
   else
   {
       eval
       {
	   `qsub -q all.q -V -wd $outdir -N SortSam.$id.$$ -o SortSam.$id.$$.stdout -e SortSam.$id.$$.stderr -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y "$JAVA -Xmx4g -jar $PICARD/AddOrReplaceReadGroups.jar I=$samFile O=$outFilename SO=coordinate RGID=$basename RGLB=$id RGPL=$platform RGPU=$barcode RGSM=$sampleId RGCN=MSKCC TMP_DIR=$TMPDIR COMPRESSION_LEVEL=0 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT"`;
	   `qsub -q all.q -V -wd $outdir -hold_jid SortSam.$id.$$ -N NotifySortSam.$id.$$ -e NotifySortSam.$id.$$.stderr -o NotifySortSam.$id.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
       };
       if($@)
       {
	   print "SortSam:Job Submission Failed, Error:$@\n";
	   exit(1);
       }
   }
   return("$outFilename","NotifySortSam.$id.$$.stat");
}

=begin
    #####################################
    #####################################
    #Sort Sam file

    sub RunSortSam
    {
        my($samFile,$outdir,$id) = @_;
        my $outFilename = $samFile;
        $outFilename =~ s/\.sam/_srt\.bam/;
        if ((-e "$outFilename") and ((-s "$outFilename") != 0)) {
            print "Files:\n$outFilename\n they exists and process will not run to make \"_srt.bam\" file.\n";
            return("$outFilename",'NULL');
        } else {
            `qsub -q all.q -wd $outdir -N SortSam.$id.$$ -l h_vmem=8G,virtual_free=8G -pe smp 1 -o /dev/null -e /dev/null -b y "$JAVA -Xmx4g -jar $PICARD/SortSam.jar I=$samFile O=$outFilename SO=coordinate TMP_DIR=$TMPDIR COMPRESSION_LEVEL=0 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT"`;
            `qsub -q all.q -V -wd $outdir -hold_jid SortSam.$id.$$ -N NotifySortSam.$id.$$ -l h_vmem=2G,virtual_free=2G -pe smp 1 -e /dev/null -o NotifySortSam.$id.$$.stat -b y "$outdir/Notify.csh"`;
        }
        return("$outFilename","NotifySortSam.$id.$$.stat");
    }
=cut

#####################################
#####################################
#Mark Duplicates in Bam

sub RunMarkDuplicates
{
   my($bamFile,$outdir,$id) = @_;
   my $outFilename = $bamFile;
   my $metricsFilename = $bamFile;
   $outFilename =~ s/\.bam/_MD\.bam/g;
   $metricsFilename =~ s/\.bam/_MD\.metrics/g;
   if((-e "$outFilename") and ((-s "$outFilename") != 0))
   {
       print "Files:\n$outFilename\n they exists and process will not run to make \"_MD.bam\" file.\n";
       return("$outFilename",'NULL');
   }
   else
   {
       eval
       {
	   `qsub -q all.q -V -wd $outdir -N MD.$id.$$  -o /dev/null -e /dev/null -l h_vmem=8G,virtual_free=8G -pe smp 1 -b y "$JAVA -Xmx4g -jar $PICARD/MarkDuplicates.jar I=$bamFile O=$outFilename ASSUME_SORTED=true METRICS_FILE=$metricsFilename TMP_DIR=$TMPDIR COMPRESSION_LEVEL=0 CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT"`;
	   `qsub -q all.q -V -wd $outdir -hold_jid MD.$id.$$ -N NotifyMD.$id.$$ -e NotifyMD.$id.$$.stderr -o NotifyMD.$id.$$.stat -l h_vmem=2G,virtual_free=2G -pe smp 1 -b y "$outdir/Notify.csh"`;
       };
       if($@)
       {
	   print "MarkDuplicates:Job Submission Failed, Error:$@\n";
	   exit(1);
       }
   }
   return("$outFilename","NotifyMD.$id.$$.stat");
}
#####################################
#####################################
#This will calculate and compile metrics for BAM files:
sub CalcHsMetrics
{
    my($filenames,$datadir,$outdir,$fof) = @_;
    my @names = ();
    if($filenames){(@names) = @$filenames;}
    if((scalar(@names) == 0)and($fof)){@names = &GetNames($fof,$outdir);}
    my @notifyNames = ();
    my(@sortedparseFilenames) = sort {lowestNumber($a) <=>  lowestNumber($b)} @names;
    @names = @sortedparseFilenames;
    ##################
    #Calculate HsMetrics
    $now = time;
    print "Started running metrics calculation jobs on SGE at ". localtime() ."\n";
    for(my $i = 0; $i < scalar(@names); $i++)
    {
	my $waitFileNames = &RunHsMetrics($names[$i],$outdir,$i);
	foreach my $waitName (@$waitFileNames)
	{
	    push(@notifyNames,$waitName);
	}
    }
    #waiting for metrics calculations to finish
    &WaitToFinish($outdir,@notifyNames);
    $now = time - $now;
    print "Finished running metrics calculation jobs on SGE at ". localtime() ."\n";
    printf("Total running time: %02d:%02d:%02d\n\n", int($now / 3600), int(($now % 3600) / 60), int($now % 60));

  return(\@names);

}
#####################################
#####################################
#Run Picard HsMetrics
sub RunHsMetrics
{
    my($bamFile,$outdir,$id) = @_;
    my($basename) = $bamFile =~ /(.*)\.bam/;
    my $HSmetricsFilename = $basename . ".HSmetrics.txt";
    my @notifynames = ();
     #Calculate Hybrid Selection specific metrics
    if((-e "$HSmetricsFilename") and ((-s "$HSmetricsFilename") != 0))
    {
	print "Files:\n$HSmetricsFilename\n they exists and process will not run to make \".HSmetrics.txt\" file.\n";
	push(@notifynames,"NULL");
    }
    else
    {
	`qsub -q all.q -wd $outdir -N HSmetrics.$id.$$ -l h_vmem=8G,virtual_free=8G -pe smp 1 -o /dev/null -e /dev/null -b y "$JAVA -Xmx4g -jar $PICARD/CalculateHsMetrics.jar I=$bamFile O=$HSmetricsFilename BI=$baitIntervalFile TI=$targetIntervalFile REFERENCE_SEQUENCE=$refFile TMP_DIR=$TMPDIR VALIDATION_STRINGENCY=LENIENT"`;
       `qsub -q all.q -V -wd $outdir -hold_jid HSmetrics.$id.$$ -N NotifyHSmetrics.$id.$$ -l h_vmem=2G,virtual_free=2G -pe smp 1 -e /dev/null -o NotifyHSmetrics.$id.$$.stat -b y "$outdir/Notify.csh"`;
	push(@notifynames,"NotifyHSmetrics.$id.$$.stat");
    }
    return(\@notifynames);
}
#####################################
#####################################
#This will help to call:
#Somatic SVs: PeSV
sub CallStructuralVariants
{
    my($filenames,$gStdNormal,$patientId,$datadir,$outdir,$fof) = @_;
    my @names = ();
    if($filenames){(@names) = @$filenames;}
    print "F:$fof\n";
    if((scalar(@names) == 0)and($fof)){@names = &GetNames($fof,$outdir);}
    my @notifyNames = ();
    tie (my %groupedFilenames, 'Tie::IxHash');
    my @somaticSVfiles = ();
    my @CoveragePerSample = ();
    my @meanCoverageValues = ();
    tie (my %coverageForNormal, 'Tie::IxHash'); 
    tie (my %NormalPerFile, 'Tie::IxHash');
    my $standardNormal;
    my $now = time;
    my(@sortedparseFilenames) = sort {lowestNumber($a) <=>  lowestNumber($b)} @names;
    @names = @sortedparseFilenames; 
    my ($poolName) = $names[0] =~ /.*_bc\d{1,2}_(.*)_L\d{1,3}.*/;
    my $NormalUsed = $poolName . "_NormalUsedInSVcalling.txt";
    my @fileNames = ();
    
    if ($names[0]=~/\//) {
        foreach ( @names) {
            my $filename =  pop @{[split("/",$_)]};
            push(@fileNames,$filename);
        }
        @names = @fileNames;
    }
   
    #Call Somatic SVs
    print "Started running Somatics Variant jobs on SGE at ". localtime() ."\n";
    my $fCount = 0;
    #Group the files
    foreach my $file (@names)
    {
        print "$file\n";
        
	if(exists $groupedFilenames{@$patientId[$fCount]})
	{
            #print "$fCount:@$patientId[$fCount]\n";
	    my $files =  $groupedFilenames{@$patientId[$fCount]};
	    $files = "$files" . ",$file";
	    $groupedFilenames{@$patientId[$fCount]} = "$files";
	}
	else
	{
           # print "$fCount:@$patientId[$fCount]\n";
	    $groupedFilenames{@$patientId[$fCount]} = "$file";
	}
	$fCount++;
    }
    
    #Get Mean Coverage from HSmetrics file.
    foreach my $file (@names)
    {
	#print "$file\n";
	my($fileBarcode) = $file =~ /.*_(bc\d{1,2})_.*/;
	my($fileSampleId) = $file =~ /(.*)_bc\d{1,2}_/;
	#print $fileBarcode . "_" . $fileSampleId, "\n";
	my $fileClass = $classPerBarcode{$fileBarcode . "_" . $fileSampleId};
        #print "$fileClass\n";
        
	if($fileClass =~ m/Normal/i)
	{
	    my $HSmetricsFile = $file;
	    $HSmetricsFile =~ s/\.bam/\.HSmetrics\.txt/g;
	    print "HS:$HSmetricsFile\n";
	    open(FH,"$outdir/$HSmetricsFile") or die "Cannot Open HSmetricsFile:$outdir/$HSmetricsFile, $!\n";
	    my $meanCov;
	    my $CovForFile;
	    while(<FH>)
	    {
		next until($_ =~ /^BAIT_SET/);
		while(<FH>)
		{
		    next if(($_ =~ /^BAIT_SET/) or ($_ =~ /^\s$/));
		    my(@values) = split("\t",$_);
		    #print "MeanCOv:$values[21]\n";
		    $CovForFile = $values[21];
		    $meanCov = $values[21];
		}
	    }
	    close(FH);
	    $coverageForNormal{$file} = $CovForFile;
	    push (@CoveragePerSample, $meanCov);
	}
	else
	{
	    next;
	}
    }
    
    #Get file that will be used as standard normal
    my $maxCoverage = max @CoveragePerSample;
    #print "MAX:$maxCoverage\n";
    while((my $key, my $value) = each (%coverageForNormal))
    {
	if($value == $maxCoverage)
	{
	    $standardNormal = $key;
	}
	else
	{
	    next;
	}
    }
    if(! $standardNormal)
    {
	$standardNormal = $gStdNormal;
    }
    #print "SN:$standardNormal\n";
    #Running Mutect and Somatic Indel Caller
    my $count = 0;
    while((my $key, my $value) = each (%groupedFilenames))
    {
	my @files = split(",",$value);
	# Section of Normal
	my @normalSamples = ();
	tie (my %coverageForSampleNormals, 'Tie::IxHash'); 
	my @CoverageForMultipleNormal = ();
	my $normal;
	foreach my $file (@files)
	{
	    my ($fileBarcode) = $file =~ /.*_(bc\d{1,2})_.*/;
	    my($fileSampleId) = $file =~ /(.*)_bc\d{1,2}_/;
	    my $fileClass = $classPerBarcode{$fileBarcode . "_" . $fileSampleId};
	    if ($fileClass =~ m/Normal/i)
	    {
		push (@normalSamples, $file)
	    }
	}
	foreach my $file (@normalSamples)
	{
	    my $HSmetricsFile = $file;
	    $HSmetricsFile =~ s/\.bam/\.HSmetrics\.txt/g;
	    #print "HS:$HSmetricsFile\n";
	    open(FH,"$outdir/$HSmetricsFile") or die "Cannot Open HSmetricsFile:$outdir/$HSmetricsFile, $!\n";
	    my $meanCov;
	    my $CovForFile;
	    while(<FH>)
	    {
		next until($_ =~ /^BAIT_SET/);
		while(<FH>)
		{
		    next if(($_ =~ /^BAIT_SET/) or ($_ =~ /^\s$/));
		    my(@values) = split("\t",$_);
		    #print "MeanCOv:$values[21]\n";
		    $CovForFile = $values[21];
		    $meanCov = $values[21];
		}
		$coverageForSampleNormals{$file} = $CovForFile;
		push (@CoverageForMultipleNormal, $meanCov);
	    }
	    close(FH);
	}
	#Get file that will be used as normal
	my $maxCoverage = max @CoverageForMultipleNormal;
	#print "MAX:$maxCoverage\n";
	if (scalar @normalSamples > 1)
	{
	    while((my $key, my $value) = each (%coverageForSampleNormals))
	    {
		if(($value == $maxCoverage) and ($value >= 50))
		{
		    $normal = $key;
		}
		else
		{
		    $normal = $standardNormal;
		}
	    }
	}
	else
	{
	    if(scalar @normalSamples == 1)
	    {
		my $coverage = $coverageForSampleNormals{$normalSamples[0]};
		if($coverage >= 50)
		{

		    $normal = $normalSamples[0];
		}
		else
		{
		    $normal = $standardNormal;
		}
	    }
	    else
	    {
		$normal = $standardNormal;
	    }
	}
	print "Normal : $normal\n";
	foreach my $file (@files)
	{
	    #print "Final1:$file\n";
	    #next if($file =~ m/.*(N\d|N|N-\d{1,3}ng|N\d-\d{1,3}ng|N-\d)_bc\d{1,2}.*/); 
	    my ($fileBarcode) = $file =~ /.*_(bc\d{1,2})_.*/;
	    my($fileSampleId) = $file =~ /(.*)_bc\d{1,2}_/;
	    my $fileClass = $classPerBarcode{$fileBarcode . "_" . $fileSampleId};
	    next if ($fileClass =~ m/Normal/i);
	    print "Final2:T->$file\nN->$normal\n\n";
	    my ($tFileId) = $file =~ /(.*)_bc\d{1,2}_/;
	    my ($nFileId) =  $normal =~ /(.*)_bc\d{1,2}_/;
	    $NormalPerFile{$tFileId} = $nFileId;
            my($waitFileNames,$dirctory);
            
            if ($prog =~ /PeSV/i) {
                ($waitFileNames) = &RunPeSVfisher($normal,$file,$nFileId,$tFileId,$PeSV,$outdir,$count);
            }
            if ($prog =~ /Meerkat/i) {
                ($waitFileNames,$dirctory) = &RunMeerkat($normal,$file,$nFileId,$tFileId,$outdir,$count);
	    }
            if ($prog =~ /Delly/i) {
                ($waitFileNames,$dirctory) = &RunDelly($normal,$file,$nFileId,$tFileId,$outdir,$count);
	    }
            foreach my $waitName (@$waitFileNames)
	    {
		push(@notifyNames,$waitName);
	    }
	    $count++;
	}
    }
    &WaitToFinish($outdir,@notifyNames);

=begin
           my $somaticSVfile = MergeSVs($outdir,$distanceSV,$oriSV,$orderSV,$chrOriSV,$chrOrderSV);
	    push(@somaticSVfiles,$somaticSVfile);
    open (NFH,">","$outdir/$NormalUsed") || die "Cannot open NormalUsedinSVFile:$outdir/$NormalUsed;$!\n";
    while(my($key,$value) = each (%NormalPerFile))
    {
	print NFH "$key\t$value\n";
    }
    close(NFH);
=cut

    $now = time - $now;
    print "Finished running Germline and Somatic Variant jobs on SGE at ". localtime() ."\n";
    printf("Total running time: %02d:%02d:%02d\n\n", int($now / 3600), int(($now % 3600) / 60), int($now % 60));

    return (\@names);
}

#######################################
#######################################
#Run PeSVfisher
sub RunPeSVfisher
{   
    my ($normal,$tumor,$nId,$tId,$PeSV,$outdir,$count) = @_;
    my @waitFilenames = ();
    #Changing the filenames so as to have smaller file names for PeSV.
    my $nNormalId = $nId;
    $nNormalId =~ s/-/_/g;
    my $tTumorId = $tId;
    $tTumorId =~ s/-/_/g;
    #get the index files
    chomp($normal);chomp($tumor);
    my $normalBai = $normal;chop($normalBai);$normalBai = $normalBai . "i";
    my $tumorBai = $tumor;chop($tumorBai);$tumorBai = $tumorBai . "i";
    print "$nNormalId:$tTumorId\n";
    #PeSV folder for Tumor only
    my $pesvTFolderName = "PeSV_" . $tId . "_$count" . "_T"; 
    #PeSV folder for Tumor Normal Analysis
    my $pesvTNFolderName = "PeSV_" . $tId . "_$count" . "_TN";
    my $PeSVtfolder = "$outdir/$pesvTFolderName";
    my $PeSVtnfolder = "$outdir/$pesvTNFolderName";
    print "$PeSVtfolder;$PeSVtnfolder\n";
    #Link executables from source dir in tumor dir
    if(! -d $PeSVtfolder)
    {
	`mkdir $PeSVtfolder`;
	`ln -s $PeSV/* $PeSVtfolder/.`;
    }
    else
    {
	#`rm -rf $PeSVtfolder`;
	#`mkdir $PeSVtfolder`;
	print "$PeSVtfolder folder exists and thus anlysis for this folder will not be done.\n";
	return();
    }
    #Link executables from source dir in tumor normal
    if(! -d $PeSVtnfolder)
    {
	`mkdir $PeSVtnfolder`;
	`ln -s $PeSV/* $PeSVtnfolder/.`;
    }
    else
    {
	#`rm -rf $PeSVtnfolder`;
	#`mkdir $PeSVtnfolder`;
	print "$PeSVtnfolder folder exists and thus anlysis for this folder will not be done.\n";
	return();
    }
    #Link tumor bam file in tumor dir
    if(-e "$PeSVtfolder/$tTumorId.bam")
    {
	print "$PeSVtfolder/$tTumorId.bam exits and wont be softlinked\n";
    }
    else
    {
	`ln -s $outdir/$tumor $PeSVtfolder/$tTumorId.bam`;
	`ln -s $outdir/$tumorBai $PeSVtfolder/$tTumorId.bai`;
    }
    #Link tumor and normal bam file in tumor normal dir
    if(-e "$PeSVtnfolder/$nNormalId.bam")
    {
	print "$PeSVtnfolder/$nNormalId.bam exits and wont be softlinked\n";
    }
    else
    {
	`ln -s $outdir/$normal $PeSVtnfolder/$nNormalId.bam`;
	`ln -s $outdir/$normalBai $PeSVtnfolder/$nNormalId.bai`;
    }
    if(-e "$PeSVtnfolder/$tTumorId.bam")
    {
	print "$PeSVtnfolder/$tTumorId.bam exits and wont be softlinked\n";
    }
    else
    {
	`ln -s $outdir/$tumor $PeSVtnfolder/$tTumorId.bam`;
	`ln -s $outdir/$tumorBai $PeSVtnfolder/$tTumorId.bai`;
    }
    #getting the configration file ready for both versions.
    my $pesvSomConfigTemplate = "$PeSVtnfolder/config/Template.som.cfg";
    my $pesvNsomConfigTemplate = "$PeSVtfolder/config/Template.Nonsom.cfg";
    my $pesvSomConfig = "$PeSVtnfolder/$pesvTNFolderName.cfg";
    my $pesvNsomConfig = "$PeSVtfolder/$pesvTFolderName.cfg";
    tie (my %somHash, 'Tie::IxHash');
    tie (my %NsomHash, 'Tie::IxHash');
    %somHash = &ReadPeSVTemplate($pesvSomConfigTemplate);
    %NsomHash = &ReadPeSVTemplate($pesvNsomConfigTemplate);
    #populate somatic & non somatic hash
    while (my($key,$value) = each (%somHash))
    {
	if($key eq "ntchr")
	{
	    $somHash{$key} = $ntchr;
	    $NsomHash{$key} = $ntchr;
	}
	if($key eq "nprocess")
	{
	    $somHash{$key} = $nprocessors;
	    $NsomHash{$key} = $nprocessors;
	}
	if($key eq "ndofcprocess")
	{
	    $somHash{$key} = $nprocessors;
	    $NsomHash{$key} = $nprocessors;
	}
	if($key eq "basedir")
	{
	    $somHash{$key} = $PeSVtnfolder . "/";
	    $NsomHash{$key} = $PeSVtfolder . "/";
	}
	if($key eq "samtoolsdir")
	{
	    $somHash{$key} = "samtools";
	    $NsomHash{$key} = "samtools";
	}
	if($key eq "awk")
	{
	    $somHash{$key} = $awk;
	    $NsomHash{$key} = $awk;
	}
	if($key eq "refaprefix")
	{
	    $somHash{$key} = $chrSepRef . "/"; 
	    $NsomHash{$key} = $chrSepRef . "/";
	}
	if($key eq "hg")
	{
	    $somHash{$key} = $hg;
	    $NsomHash{$key} = $hg;
	}
	if($key eq "istargeted")
	{
            if ($isTargeted == 1) {
                $somHash{$key} = 1;
                $NsomHash{$key} = 1;
            }
            else {
                $somHash{$key} = 0;
                $NsomHash{$key} = 0;
            }
                
	}
	if($key eq "tfile")
	{
	    $somHash{$key} = $targetIntervalBedFile;
	    $NsomHash{$key} = $targetIntervalBedFile;
	}
	if($key eq "datadirbam")
	{
	    $somHash{$key} = $PeSVtnfolder . "/"; 
	    $NsomHash{$key} = $PeSVtfolder . "/";
	}
	if($key eq "prefix")
	{
	     $NsomHash{$key} = $tTumorId;
	}
	if($key eq "control")
	{
	    $somHash{$key} = $nNormalId;
	}
	if($key eq "disease")
	{
	    $somHash{$key} = $tTumorId;
	}
	if($key eq "login")
	{
	     $somHash{$key} = $dbUser;
	     $NsomHash{$key} = $dbUser;
	}
	if($key eq "passwd")
	{
	     $somHash{$key} = $dbPass;
	     $NsomHash{$key} = $dbPass;
	}
	if($key eq "database")
	{
	     $somHash{$key} = $dbName;
	     $NsomHash{$key} = $dbName;
	}
	if($key eq "conn")
	{
	     $somHash{$key} = $dbLocStr;
	     $NsomHash{$key} = $dbLocStr;
	}
    }
    #print som hash
    open (SFH,">",$pesvSomConfig) or die "Cannot open PeSVsomaticConfigFile:$pesvSomConfig,$!\n";
    while(my($key,$value) = each (%somHash))
    {
	print SFH "$key=$value\n";
    }
    close(SFH);
    open (NSFH,">",$pesvNsomConfig) or die "Cannot open PeSVnon-somaticConfigFile:$pesvNsomConfig,$!\n";
    while(my($key,$value) = each (%NsomHash))
    {
	print NSFH "$key=$value\n";
    }
    close(NSFH);
    my $path = $ENV{PATH};
    $ENV{PATH} = $path . ":" . $samtoolsDir;
    my $PYTHON3="/ifs/e63data/bergerm1/Resources/SupportTools/Python/installs/Python-3.3.0/bin/python3";
    #Run Somatic Mode
    `qsub -q all.q -v PSVFDIR=$PeSVtnfolder,PSVFPYTHON=$PYTHON3 -V -wd $PeSVtnfolder -N PeSVsom.$count.$$ -o PeSVsom.$count.$$.log -l h_vmem=5G,virtual_free=5G -pe smp $nprocessors -b y "$PeSVtnfolder/bin/PeSVFisher $pesvSomConfig"`;
    `qsub -q all.q -wd $PeSVtnfolder -hold_jid PeSVsom.$count.$$ -N NotifyPeSVsom.$count.$$ -l h_vmem=2G,virtual_free=2G -pe smp 1 -e /dev/null -o NotifyPeSVsom.$count.$$.stat -b y "$outdir/Notify.csh"`;
    `qsub -q all.q -v PSVFDIR=$PeSVtfolder,PSVFPYTHON=$PYTHON3 -V -wd $PeSVtfolder -N PeSVnsom.$count.$$ -o PeSVnsom.$count.$$.log -l h_vmem=5G,virtual_free=5G -pe smp $nprocessors -b y "$PeSVtfolder/bin/PeSVFisher $pesvNsomConfig"`;
    `qsub -q all.q -wd $PeSVtfolder -hold_jid PeSVnsom.$count.$$ -N NotifyPeSVnsom.$count.$$ -l h_vmem=2G,virtual_free=2G -pe smp 1 -e /dev/null -o NotifyPeSVnsom.$count.$$.stat -b y "$outdir/Notify.csh"`;
    push(@waitFilenames,"NotifyPeSVsom.$count.$$.stat");
    push(@waitFilenames,"NotifyPeSVnsom.$count.$$.stat");
    return(\@waitFilenames);
}
#######################################
#######################################
#Run Meerkat
sub RunMeerkat
{   
    my ($normal,$tumor,$nId,$tId,$outdir,$count) = @_;
    my @waitFilenames = ();
    my $meerkatOutdir = $outdir . "/MeerkatOutput";
    my $sampleNormalOutput = $meerkatOutdir ."/". $nId . "-" . $count;
    my $sampleTumorOutput = $meerkatOutdir . "/" . $tId . "-" . $count;
    my ($tFlag,$nFlag);
    
    #Make Ouput dir
    if (! (-d "$meerkatOutdir")) {
        `mkdir $meerkatOutdir`;
    }
    else {
        if ($count == 0) {
            warn "$meerkatOutdir exists !!\n";
        }
        
    }
    #for making link to files
    my $lnNormal = $normal;
    chop($lnNormal);
    my $lnTumor = $tumor;
    chop($lnTumor);
    #Make Normal Sample Output Dir
    if (!(-d "$sampleNormalOutput")) {
        `mkdir $sampleNormalOutput`;
        `ln -s $outdir/$lnNormal* $sampleNormalOutput/`;
        `ln -s $outdir/$lnTumor* $sampleNormalOutput/`;
    }
    else {
        warn "$sampleNormalOutput exists !!\n";
        $tFlag = 1;
    }
    #Make Tumor Sample Output Dir
    if (!(-d "$sampleTumorOutput")) {
        `mkdir $sampleTumorOutput`;
        `ln -s $outdir/$lnTumor* $sampleTumorOutput/`;
        `ln -s $outdir/$lnNormal* $sampleTumorOutput/`;
    }
    else {
        warn "$sampleTumorOutput exists !!\n";
        $nFlag = 1;
    }
    my $mem;
    #Assign Memory
    if ($nprocessors >= 2 ) {
        $mem = "4G";
        
    }
    else {
        $mem = "10G";
        
    }
    #Assign Queue
    my  $runQueue = $queue;
    #Notify CMD
    my $notify_cmd = "$outdir/Notify.csh";
    #Normal CMD values for Preprocess
    my $N1_cmd = "$PERL $Meerkat/scripts/pre_process.pl -b $normal -I $refFile -A $refFile.fai -S $samtoolsDir -W $bwaDir -k 1500 -s 20 -q 15 -t $nprocessors";
    my $N1_jname = "preprocess_$nId.$$.$count";
    my $N1_stdout = $N1_jname . ".stdout";
    my $N1_stderr = $N1_jname . ".stderr";

    #Normal CMD values for Meerkat
    my $N2_cmd;
    my $N2_jname;
    my $N2_stdout;
    my $N2_stderr;
    
    if ($bwaFlag == 1) {
        $N2_cmd = "$PERL $Meerkat/scripts/meerkat.pl -b $normal -F $refFasta -S $samtoolsDir -W $bwaDir -B $blastDir -s 20 -p 3 -o 1 -t $nprocessors";
        $N2_jname = "meerkat_$nId.$$.$count";
        $N2_stdout = $N2_jname . ".stdout";
        $N2_stderr = $N2_jname . ".stderr";
    }
    if ($bwaFlag == 2) {
        
        $N2_cmd = "$PERL $Meerkat/scripts/meerkat.pl -b $normal -F $refFasta -S $samtoolsDir -W $bwaDir -B $blastDir -u 1 -s 20 -p 3 -o 1 -t $nprocessors";
        $N2_jname = "meerkat_$nId.$$.$count";
        $N2_stdout = $N2_jname . ".stdout";
        $N2_stderr = $N2_jname . ".stderr";
    }
    #Tumor CMD values for Preprocess
    my $T1_cmd = "$PERL $Meerkat/scripts/pre_process.pl -b $tumor -I $refFile -A $refFile.fai -S $samtoolsDir -W $bwaDir -k 1500 -s 20 -q 15 -t $nprocessors";
    my $T1_jname = "preprocess_$tId.$$.$count";
    my $T1_stdout = $T1_jname . ".stdout";
    my $T1_stderr = $T1_jname . ".stderr";

    #Tumor CMD values for Meerkat
    my $T2_cmd;
    my $T2_jname;
    my $T2_stdout;
    my $T2_stderr;
    if ($bwaFlag == 1) {
        $T2_cmd = "$PERL $Meerkat/scripts/meerkat.pl -b $tumor -F $refFasta -S $samtoolsDir -W $bwaDir -B $blastDir -s 20 -p 3 -o 1 -t $nprocessors";
        $T2_jname = "meerkat_$tId.$$.$count";
        $T2_stdout = $T2_jname . ".stdout";
        $T2_stderr = $T2_jname . ".stderr";
    }
    if ($bwaFlag == 2) {
        $T2_cmd = "$PERL $Meerkat/scripts/meerkat.pl -b $tumor -F $refFasta -S $samtoolsDir -W $bwaDir -B $blastDir -u 1 -s 20 -p 3 -o 1 -t $nprocessors";
        $T2_jname = "meerkat_$tId.$$.$count";
        $T2_stdout = $T2_jname . ".stdout";
        $T2_stderr = $T2_jname . ".stderr";
    }
    
    #Tumor CMD values for Mechanism
    my $T3_cmd = "$PERL $Meerkat/scripts/mechanism.pl -b $tumor -R $RepeatMaskFile";
    my $T3_jname = "mechanism_$tId.$$.$count";
    my $T3_stdout = $T3_jname . ".stdout";
    my $T3_stderr = $T3_jname . ".stderr";
    #Define file names
    
    my($variantName) = $tumor =~ /(.*)\.bam/;
    my $variantFile = $variantName . ".variants";
    my $variantFileA = $variantName . "_A.variants";
    my $variantFileB = $variantName . "_B.variants";
    my $variantFileC = $variantName . "_C.variants";
    my $variantFileD = $variantName . "_D.variants";
    my $variantFileE = $variantName . "_E.variants";
    my $variantFileF = $variantName . "_F.variants";
    my $variantFileG = $variantName . "_G.variants";
    #Tumor CMD values for sv_somtic filter
    #name of folder contains all *.discord files from normal genomes to filter
    #germline events, recommended to filter against all normal genomes from one tumor type
    my $Ta_cmd = "$PERL $Meerkat/scripts/somatic_sv.pl -R $RepeatMaskFile -i $variantFile -o $variantFileA -F $sampleNormalOutput/";
    my $Ta_jname = "somaticSVa_$tId.$$.$count";
    my $Ta_stdout = $Ta_jname . ".stdout";
    my $Ta_stderr = $Ta_jname . ".stderr";
    #Tumor CMD values for sv_somtic filter by total number of discordant read pairs in matched normal genome
    my $Tb_cmd = "$PERL $Meerkat/scripts/somatic_sv.pl -S $samtoolsDir -R $RepeatMaskFile -i $variantFileA -o $variantFileB -n 1 -b $normal";
    my $Tb_jname = "somaticSVb_$tId.$$.$count";
    my $Tb_stdout = $Tb_jname . ".stdout";
    my $Tb_stderr = $Tb_jname . ".stderr";
    #Tumor CMD values for sv_somtic filter by non-uniq mapped reads in matched normal genome, only works for BWA aligned bam with XT tag
    my $Tc_cmd = "$PERL $Meerkat/scripts/somatic_sv.pl -S $samtoolsDir -R $RepeatMaskFile -i $variantFileB -o $variantFileC -u 1 -b $normal";
    my $Tc_jname = "somaticSVc_$tId.$$.$count";
    my $Tc_stdout = $Tc_jname . ".stdout";
    my $Tc_stderr = $Tc_jname . ".stderr";
    #Tumor CMD values for sv_somtic filter by soft-clipped reads in matched normal genome
    my $Td_cmd = "$PERL $Meerkat/scripts/somatic_sv.pl -S $samtoolsDir -R $RepeatMaskFile -i $variantFileC -o $variantFileD -f 1 -b $normal";
    my $Td_jname = "somaticSVd_$tId.$$.$count";
    my $Td_stdout = $Td_jname . ".stdout";
    my $Td_stderr = $Td_jname . ".stderr";
    #Tumor CMD values for sv_somtic filter by total number of discordant read pairs in tumor genome, if certain breakpoint has
    #too many discordant read pairs supporting different events, it probably is an artifact
    my $Te_cmd = "$PERL $Meerkat/scripts/somatic_sv.pl -S $samtoolsDir -R $RepeatMaskFile -i $variantFileD -o $variantFileE -e 1 -B $tumor";
    my $Te_jname = "somaticSVe_$tId.$$.$count";
    my $Te_stdout = $Te_jname . ".stdout";
    my $Te_stderr = $Te_jname . ".stderr";
    #Tumor CMD values for sv_somtic
    #filter by number of supporting discordant read pairs, default 3
    #filter by number of supporting split reads, default 1
    #filter by sum of supporting discordant read pairs and supporting split reads, default 6
    my $Tf_cmd = "$PERL $Meerkat/scripts/somatic_sv.pl -S $samtoolsDir -R $RepeatMaskFile -i $variantFileE -o $variantFileF -z 1";
    my $Tf_jname = "somaticSVf_$tId.$$.$count";
    my $Tf_stdout = $Tf_jname . ".stdout";
    my $Tf_stderr = $Tf_jname . ".stderr";
    #Tumor CMD values for sv_somtic
    #max homology allowed for deletion and intra-chr events
    #max homology allowed for inter chromosomal translocation events
    my $Tg_cmd = "$PERL $Meerkat/scripts/somatic_sv.pl -S $samtoolsDir -R $RepeatMaskFile -i $variantFileF -o $variantFileG -d 40-t 20";
    my $Tg_jname = "somaticSVg_$tId.$$.$count";
    my $Tg_stdout = $Tg_jname . ".stdout";
    my $Tg_stderr = $Tg_jname . ".stderr";
    #Tumor CMD for annotation
    my $Tfusions_cmd = "$PERL $Meerkat/scripts/fusions.pl -i $variantFileG -G $RefGeneFile";
    my $Tfusions_jname = "somaticSVg_$tId.$$.$count";
    my $Tfusions_stdout = $Tfusions_jname . ".stdout";
    my $Tfusions_stderr = $Tfusions_jname . ".stderr";
    #Notify Tumor CMD values
    my $Tnotify_hjname = $Tfusions_jname;
    my $Tnotify_jname = "NotifyMeerkat.$tId.$$.$count";
    my $Tnotify_stdout = $Tnotify_jname . ".stat";
    my $Tnotify_stderr = $Tnotify_jname . ".stderr";

    #Launch only if folder does not exists
    if (($tFlag) and ($nFlag)) {
   
        if (($tFlag == 1) and ($nFlag == 1)) {
            print "Resuts for tumor & normal sample exists. Thus Meerkat would not be ran\n";
            push(@waitFilenames,"NULL");
            return(\@waitFilenames,"$sampleTumorOutput,$sampleNormalOutput");
        }
        
    }

    #Run Normal Samples
    if (! $nFlag) {
        
        &launchQsub ($N1_cmd,$sampleNormalOutput,$mem,$N1_stdout,$N1_stderr,$nprocessors,$runQueue,$N1_jname,"Null");
        &launchQsub ($N2_cmd,$sampleNormalOutput,$mem,$N2_stdout,$N2_stderr,$nprocessors,$runQueue,$N2_jname,$N1_jname);
    }
    else {
        print "Resuts for normal sample exists. Thus Meerkat would not be ran\n";
        
    }
    #Run Tumor Samples
    if (! $tFlag) {
    
        
        &launchQsub ($T1_cmd,$sampleTumorOutput,$mem,$T1_stdout,$T1_stderr,$nprocessors,$runQueue,$T1_jname,"Null");
        if (! $nFlag) {
        
            &launchQsub ($T2_cmd,$sampleTumorOutput,$mem,$T2_stdout,$T2_stderr,$nprocessors,$runQueue,$T2_jname,"$T1_jname,$N2_jname");
            
        }
        else {
            
            &launchQsub ($T2_cmd,$sampleTumorOutput,$mem,$T2_stdout,$T2_stderr,$nprocessors,$runQueue,$T2_jname,$T1_jname);
            
        }
       
        
        &launchQsub ($T3_cmd,$sampleTumorOutput,$mem,$T3_stdout,$T3_stderr,$nprocessors,$runQueue,$T3_jname,$T2_jname);
        
        &launchQsub ($Ta_cmd,$sampleTumorOutput,"5G",$Ta_stdout,$Ta_stderr,"1",$runQueue,$Ta_jname,$T3_jname);
        
        &launchQsub ($Tb_cmd,$sampleTumorOutput,"5G",$Tb_stdout,$Tb_stderr,"1",$runQueue,$Tb_jname,$Ta_jname);

        &launchQsub ($Tc_cmd,$sampleTumorOutput,"5G",$Tc_stdout,$Tc_stderr,"1",$runQueue,$Tc_jname,$Tb_jname);
                
        &launchQsub ($Td_cmd,$sampleTumorOutput,"5G",$Td_stdout,$Td_stderr,"1",$runQueue,$Td_jname,$Tc_jname);
        
        &launchQsub ($Te_cmd,$sampleTumorOutput,"5G",$Te_stdout,$Te_stderr,"1",$runQueue,$Te_jname,$Td_jname);
        
        &launchQsub ($Tf_cmd,$sampleTumorOutput,"5G",$Tf_stdout,$Tf_stderr,"1",$runQueue,$Tf_jname,$Te_jname);
        
        &launchQsub ($Tg_cmd,$sampleTumorOutput,"5G",$Tg_stdout,$Tg_stderr,"1",$runQueue,$Tg_jname,$Tf_jname);
        
        &launchQsub ($Tfusions_cmd,$sampleTumorOutput,"5G",$Tfusions_stdout,$Tfusions_stderr,"1",$runQueue,$Tfusions_jname,$Tg_jname);
        
        &launchQsub ($notify_cmd,$outdir,$mem,$Tnotify_stdout,$Tnotify_stderr,"1",$runQueue,$Tnotify_jname,$Tnotify_hjname);
        push(@waitFilenames,$Tnotify_stdout);
        
    }
    else {

        print "Resuts for tumor sample exists. Thus Meerkat would not be ran\n";
        
    }
    
    return(\@waitFilenames,"$sampleTumorOutput,$sampleNormalOutput");

}

#######################################
#######################################
#Run Delly
sub RunDelly {
    
    my ($normal,$tumor,$nId,$tId,$outdir,$count) = @_;
    my @waitFilenames = ();
    my $dellyOutdir = $outdir . "/DellyOutput";
    my $sampleNormalOutput = $dellyOutdir ."/". $nId . "-" . $count;
    my $sampleTumorOutput = $dellyOutdir . "/" . $tId . "-" . $count;
    my ($tFlag,$nFlag);
    
    #Make Ouput dir
    if (! (-d "$dellyOutdir")) {
        `mkdir $dellyOutdir`;
    }
    else {
        if ($count == 0) {
            warn "$dellyOutdir exists !!\n";
        }
    }
    #for making link to files
    my $lnNormal = $normal;
    chop($lnNormal);
    my $lnTumor = $tumor;
    chop($lnTumor);
    #Make Normal Sample Output Dir
    if (!(-d "$sampleNormalOutput")) {
        `mkdir $sampleNormalOutput`;
        `ln -s $outdir/$lnNormal* $sampleNormalOutput/`;
        `ln -s $outdir/$lnTumor* $sampleNormalOutput/`;
    }
    else {
        warn "$sampleNormalOutput exists !!\n";
        $tFlag = 1;
    }
    #Make Tumor Sample Output Dir
    if (!(-d "$sampleTumorOutput")) {
        `mkdir $sampleTumorOutput`;
        `ln -s $outdir/$lnTumor* $sampleTumorOutput/`;
        `ln -s $outdir/$lnNormal* $sampleTumorOutput/`;
    }
    else {
        warn "$sampleTumorOutput exists !!\n";
        $nFlag = 1;
    }
    my $mem;
    #Assign Memory
    if ($nprocessors >= 2 ) {
        $mem = "4G";
        
    }
    else {
        $mem = "10G";
        
    }
    #Assign Queue
    my  $runQueue = $queue;
    #Notify CMD
    my $notify_cmd = "$outdir/Notify.csh";
    #Delly Normal CMD
    my $dellyN_cmd = "$DELLY/delly/delly -p -g $refFile -i $nId -o $nId\_del.txt -b $nId\_delbrkpts.txt $normal";
    my $dellyN_jname = "delly_$nId.$$.$count";
    my $dellyN_stdout =  $dellyN_jname . ".stdout";
    my $dellyN_stderr =  $dellyN_jname . ".stderr";
    #Delly Tumor CMD
    my $dellyT_cmd = "$DELLY/delly/delly -p -g $refFile -i $tId -o $tId\_del.txt -b $tId\_delbrkpts.txt $tumor";
    my $dellyT_jname = "delly_$tId.$$.$count";
    my $dellyT_stdout =  $dellyT_jname . ".stdout";
    my $dellyT_stderr =  $dellyT_jname . ".stderr";
    #Duppy Normal CMD
    my $duppyN_cmd = "$DELLY/duppy/duppy -p -g $refFile -i $nId -o $nId\_dup.txt -b $nId\_dupbrkpts.txt $normal";
    my $duppyN_jname = "duppy_$nId.$$.$count";
    my $duppyN_stdout =  $duppyN_jname . ".stdout";
    my $duppyN_stderr =  $duppyN_jname . ".stderr";
    #Duppy Tumor CMD
    my $duppyT_cmd = "$DELLY/duppy/duppy -p -g $refFile -i $tId -o $tId\_dup.txt -b $tId\_dupbrkpts.txt $tumor";
    my $duppyT_jname = "duppy_$tId.$$.$count";
    my $duppyT_stdout =  $duppyT_jname . ".stdout";
    my $duppyT_stderr =  $duppyT_jname . ".stderr";
    #Invy Normal CMD
    my $invyN_cmd = "$DELLY/invy/invy -p -g $refFile -i $nId -o $nId\_inv.txt -r $nId\_invmerged.txt -b $nId\_invbrkpts.txt -k $nId\_invbrkptsmerged.txt $normal";
    my $invyN_jname = "invy_$nId.$$.$count";
    my $invyN_stdout =  $invyN_jname . ".stdout";
    my $invyN_stderr =  $invyN_jname . ".stderr";
    #Invy Tumor CMD
    my $invyT_cmd = "$DELLY/invy/invy -p -g $refFile -i $tId -o $tId\_inv.txt -r $tId\_invmerged.txt -b $tId\_invbrkpts.txt -k $tId\_invbrkptsmerged.txt $tumor";
    my $invyT_jname = "invy_$tId.$$.$count";
    my $invyT_stdout =  $invyT_jname . ".stdout";
    my $invyT_stderr =  $invyT_jname . ".stderr";
    #Jumpy Normal CMD
    my $jumpyN_cmd = "$DELLY/jumpy/jumpy -p -g $refFile -i $nId -o $nId\_jmp.txt -r $nId\_jmpmerged.txt -b $nId\_jmpbrkpts.txt -k $nId\_jmpbrkptsmerged.txt $normal";
    my $jumpyN_jname = "jumpy_$nId.$$.$count";
    my $jumpyN_stdout =  $jumpyN_jname . ".stdout";
    my $jumpyN_stderr =  $jumpyN_jname . ".stderr";
    #Invy Tumor CMD
    my $jumpyT_cmd = "$DELLY/jumpy/jumpy -p -g $refFile -i $tId -o $tId\_jmp.txt -r $tId\_jmpmerged.txt -b $tId\_jmpbrkpts.txt -k $tId\_jmpbrkptsmerged.txt $tumor";
    my $jumpyT_jname = "jumpy_$tId.$$.$count";
    my $jumpyT_stdout =  $jumpyT_jname . ".stdout";
    my $jumpyT_stderr =  $jumpyT_jname . ".stderr";
    #Notify Normal CMD values
    my $notifyN_hjname = "$dellyN_jname,$duppyN_jname,$invyN_jname,$jumpyN_jname";
    my $notifyN_jname = "NotifyDelly.$nId.$$.$count";
    my $notifyN_stdout = $notifyN_jname . ".stat";
    my $notifyN_stderr = $notifyN_jname . ".stderr";
    #Notify Tumor CMD values
    my $notifyT_hjname = "$dellyT_jname,$duppyT_jname,$invyT_jname,$jumpyT_jname";
    my $notifyT_jname = "NotifyDelly.$tId.$$.$count";
    my $notifyT_stdout = $notifyT_jname . ".stat";
    my $notifyT_stderr = $notifyT_jname . ".stderr";
    #Launch only if folder does not exists
    if (($tFlag) and ($nFlag)) {
   
        if (($tFlag == 1) and ($nFlag == 1)) {
            print "Resuts for tumor & normal sample exists. Thus Meerkat would not be ran\n";
            push(@waitFilenames,"NULL");
            return(\@waitFilenames,"$sampleTumorOutput,$sampleNormalOutput");
        }
        
    }
    if (! $nFlag) {
        &launchQsub($dellyN_cmd,$sampleNormalOutput,"10G",$dellyN_stdout,$dellyN_stderr,"1",$runQueue,$dellyN_jname,"Null");
        &launchQsub($duppyN_cmd,$sampleNormalOutput,"10G",$duppyN_stdout,$duppyN_stderr,"1",$runQueue,$duppyN_jname,"Null");
        &launchQsub($invyN_cmd,$sampleNormalOutput,"10G",$invyN_stdout,$invyN_stderr,"1",$runQueue,$invyN_jname,"Null");
        &launchQsub($jumpyN_cmd,$sampleNormalOutput,"10G",$jumpyN_stdout,$jumpyN_stderr,"1",$runQueue,$jumpyN_jname,"Null");
        &launchQsub($notify_cmd,$outdir,"2G",$notifyN_stdout,$notifyN_stderr,"1",$runQueue,$notifyN_jname,$notifyN_hjname);
        push(@waitFilenames,$notifyN_stdout);
    }
    else {
        
         print "Resuts for normal:$nId sample exists. Thus Delly would not be ran\n";
         push(@waitFilenames,"NULL");
            
    }
    if (! $tFlag) {
        &launchQsub($dellyT_cmd,$sampleTumorOutput,"10G",$dellyT_stdout,$dellyT_stderr,"1",$runQueue,$dellyT_jname,"Null");
        &launchQsub($duppyT_cmd,$sampleTumorOutput,"10G",$duppyT_stdout,$duppyT_stderr,"1",$runQueue,$duppyT_jname,"Null");
        &launchQsub($invyT_cmd,$sampleTumorOutput,"10G",$invyT_stdout,$invyT_stderr,"1",$runQueue,$invyT_jname,"Null");
        &launchQsub($jumpyT_cmd,$sampleTumorOutput,"10G",$jumpyT_stdout,$jumpyT_stderr,"1",$runQueue,$jumpyT_jname,"Null");
        &launchQsub($notify_cmd,$outdir,"2G",$notifyT_stdout,$notifyT_stderr,"1",$runQueue,$notifyT_jname,$notifyT_hjname);
        push(@waitFilenames,$notifyT_stdout);
    }
    else {
        
         print "Resuts for Tumor:$tId sample exists. Thus Delly would not be ran\n";
         push(@waitFilenames,"NULL");
            
    }
    
    return(\@waitFilenames,"$sampleTumorOutput,$sampleNormalOutput");
    
    
}

#######################################
#######################################
#Run the the cmd as qsub
sub launchQsub{
    my ($cmd,$outdir,$mem,$stdout,$stderr,$processors,$queue,$jobname,$holdjobname) = @_;
    
    #Run Job with hold job id
    if ($holdjobname ne "Null") {

        `qsub -q $queue -V -wd $outdir -hold_jid $holdjobname -N $jobname -o $stdout -e $stderr -l h_vmem=$mem,virtual_free=$mem -pe smp $processors -b y $cmd`;
    }

    #Run Jobs without hold job Id
    if ($holdjobname eq "Null") {
        
        `qsub -q $queue -V -wd $outdir -N $jobname -o $stdout -e $stderr -l h_vmem=$mem,virtual_free=$mem -pe smp $processors -b y $cmd`;
    }
    
    
    return;
    
}
 
#####################################
#####################################
#Read Template PeSV file

sub ReadPeSVTemplate
{
    my($pesvTemplate) = @_;
    tie (my %pesvConfigHash, 'Tie::IxHash');
    open(FH,$pesvTemplate) or die "Cannot open $pesvTemplate,$!\n";
    while(<FH>)
    {
	chomp($_);
	my @data = split("=",$_);
	$pesvConfigHash{$data[0]} = $data[1];
    }
    close(FH);
    return(%pesvConfigHash);
}

#####################################
#####################################
#This will help to Filter:
#Somatic SVs
sub FilterStructuralVariants
{
    return;
}
#####################################
#####################################
#This will help to Annotate:
#Somatic SVs
sub AnnotateStructuralVariants
{
    return;
}
