#!/usr/bin/perl -w
use Parallel::ForkManager;
use Data::Dumper;
use POSIX; 
use Getopt::Long;

my %param =(
'FASTQDIR'=>'fastq',
'GENOME' => undef,
'FASTALIST' => undef,
'PROJECTNAME' => './',
'THREADS' => 1,
'JOBS' => 1,
'METRICS'=>0,
'ALIGN'=>0,
'LIBCOMPLEX'=>1,
'FILTER'=>0
);

my $help;

usage() if ( @ARGV < 1 or
          ! GetOptions('fastalist=s' => \$param{'FASTALIST'}, 'fastqdir:s' => \$param{'FASTQDIR'}, 'project=s' => \$param{'PROJECTNAME'}, 'jobs:i' => \$param{'JOBS'},'threads:i' => \$param{'THREADS'}, 'genome:s' => \$param{'GENOME'},'align:i' =>  \$param{'ALIGN'},'metrics:i' =>  \$param{'METRICS'},'libcomplex:i' =>  \$param{'LIBCOMPLEX'}, 'filter:i' =>  \$param{'FILTER'})
		);
 
sub usage
{
  print "Unknown option: @_\n" if ( @_ );
  print "usage: program [--fastalist FASTAFILE list ] [--fastqdir fastq ] [--project PROJECTNAME ] [--jobs 1] [--threads 1] [--genome] [--filter 1]  [--metrics 0]  [--libcomplex 1]\n";
  exit;
}


$param{'ALNDIR'} ="$param{'PROJECTNAME'}/STAR";
$param{'RESULTSDIR'} ="$param{'PROJECTNAME'}/results";
$param{'DATADIR'} ="$param{'PROJECTNAME'}/data";



# make directory with projectname and subdirectories gsnap, cufflinks and genecounts
#system("mkdir $param{'PROJECTNAME'}") unless (-d $param{'PROJECTNAME'});

system("mkdir $param{'PROJECTNAME'}") unless (-d $param{'PROJECTNAME'});
system("mkdir $param{'ALNDIR'}") unless (-d $param{'ALNDIR'});

system("mkdir $param{'RESULTSDIR'}") unless (-d $param{'RESULTSDIR'});
system("mkdir $param{'DATADIR'}") unless (-d $param{'DATADIR'});

#move the list of fasta files to the project dir to be used later
system("cp $param{'FASTALIST'} $param{'DATADIR'}/samples.csv");

#open($LOG, '>', "$param{'PROJECTNAME'}/log");


#open file with filename which is value of the key FASTALIST
open FILE, $param{'FASTALIST'} or die 
<FILE>;
my @keys = map{chomp;$_} split(',',<FILE>);


#create an array samples
my %samples;


#splitting each line based on ',' and storing in an array @r
#pushing the reference of this array in another array @samples

while(<FILE>){
	chomp;
	#my @r = split(',');
	my %h;
	@h{@keys}=split(',');
	my $sample_name = $h{'sample'};
		
	if($samples{$sample_name}{'fastq'}){
		$samples{$sample_name}{'fastq'}= "$samples{$sample_name}{'fastq'},fastq/$h{'fastq'}"; 
	}else{
		$samples{$sample_name}{'fastq'}= "fastq/$h{'fastq'}"; 
	}	

}


#run parallel jobs
my $pm=new Parallel::ForkManager($param{'JOBS'});


##########################  first passs #####################################


foreach (keys %samples)
{
	$pm->start and next;
	
	if($param{'ALIGN'}){
		my $STARcmd = "STAR --runThreadN $param{'THREADS'}  --limitBAMsortRAM  10000000000 --runMode alignReads --genomeDir ~/NGSshare/".$param{'GENOME'}."_STAR --outSAMtype BAM SortedByCoordinate --clip3pAdapterSeq AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --alignIntronMax 1 --alignEndsType EndToEnd --outReadsUnmapped Fastx --outSAMattributes NH HI AS NM MD --readFilesCommand zcat --readFilesIn $samples{$_}->{'fastq'}  --outFileNamePrefix $param{'ALNDIR'}/$_";
		print $STARcmd,"\n";		
		system($STARcmd);
		system("samtools index $param{'ALNDIR'}/".$_."Aligned.sortedByCoord.out.bam");
		}

	if($param{'METRICS'}){
		my $metricscmd = "java -jar /opt/picard/picard.jar CollectInsertSizeMetrics I= $param{'ALNDIR'}/".$_."Aligned.sortedByCoord.out.bam O=$param{'ALNDIR'}/".$_.".metricsfull H=$param{'ALNDIR'}/".$_.".metrics.pdf";
		print $metricscmd,"\n";
		system($metricscmd);
		#system("grep PF_BASES -A1 $param{'ALNDIR'}/".$_.".metricsfull > $param{'ALNDIR'}/".$_.".metrics");
		}

	if($param{'LIBCOMPLEX'}){
		my $libcomplexcmd = "java -jar /opt/picard/picard.jar EstimateLibraryComplexity I= $param{'ALNDIR'}/".$_."_filter.bam O=$param{'ALNDIR'}/".$_.".libcomplexfull";
		print $libcomplexcmd,"\n";
		system($libcomplexcmd);
		system("grep LIBRARY -A1 $param{'ALNDIR'}/".$_.".libcomplexfull > $param{'ALNDIR'}/".$_.".libcomplex");
		}

	if($param{'FILTER'}){
		my $FilterSingle = "samtools view -b -q 255 $param{'ALNDIR'}/$_\Aligned.sortedByCoord.out.bam > $param{'ALNDIR'}/$_\_single.bam";
		print $FilterSingle,"\n";
		system($FilterSingle);
	

		my $dupCmd = "java -jar /opt/picard/picard.jar MarkDuplicates I=$param{'ALNDIR'}/$_\_single.bam O=$param{'ALNDIR'}/$_\_dup.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$param{'ALNDIR'}/$_\_dup.metrics";
		print $dupCmd, "\n";
		system($dupCmd);

		my $filterCmd = "bamtools filter -isDuplicate false -in $param{'ALNDIR'}/$_\_dup.bam -out $param{'ALNDIR'}/$_\_filter.bam";
		print $filterCmd, "\n";
		system($filterCmd);


		my $samstatsCmd = "samtools flagstat $param{'ALNDIR'}/$_\_filter.bam > $param{'ALNDIR'}/$_\_filter.flagstat";
		print $samstatsCmd, "\n";
		system($samstatsCmd)		
	}


	print "$_ completed\n"; 
	$pm->finish;
	#exit;
}

$pm->wait_all_children;

	$Summ_cmd = "~/ngs/bin/parseSTARLog.pl --dir $param{'ALNDIR'} > $param{'RESULTSDIR'}/STAR_summmary.csv";
	print $Summ_cmd,"\n";
	system($Summ_cmd);



print "Run complete\n";

