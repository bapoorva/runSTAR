#!/usr/bin/perl -w
use Parallel::ForkManager;
use Data::Dumper;
use POSIX; 
use Getopt::Long;
use feature "switch";

my %param =(
'FASTQDIR'=>'fastq',
'GENOME' => undef,
'FASTALIST' => undef,
'PROJECTNAME' => './',
'THREADS' => 1,
'ADAPTER' =>'illumina',
'JOBS' => 1,
'MERGE' =>1,
'DEDUP' =>0,
'METRICS'=>0,
'ALIGN'=>1,
'LIBCOMPLEX'=>0,
'TIN'=>0
);

my $help;

usage() if ( @ARGV < 1 or
          ! GetOptions('fastalist=s' => \$param{'FASTALIST'}, 'fastqdir:s' => \$param{'FASTQDIR'}, 'project=s' => \$param{'PROJECTNAME'}, 'adapter=s' => \$param{'ADAPTER'}, 'jobs:i' => \$param{'JOBS'},'threads:i' => \$param{'THREADS'}, 'genome:s' => \$param{'GENOME'},'align:i' =>  \$param{'ALIGN'}, 'mergefiles:i' =>  \$param{'MERGE'},'dedup:i' =>  \$param{'DEDUP'},'metrics:i' =>  \$param{'METRICS'},'libcomplex:i' =>  \$param{'LIBCOMPLEX'},'tin:i' =>  \$param{'TIN'} )
		);
 
sub usage
{
  print "Unknown option: @_\n" if ( @_ );
  print "usage: program [--fastalist FASTAFILE list ] [--fastqdir fastq ] [--project PROJECTNAME ] [--jobs 1] [--threads 1] [--genome hg19|mm9|mm10 ] [--adapter illumina|nextera] [--mergefiles 1] [--dedup 0] [--metrics 0] [--libcomplex 0] [--align 1] [--tin 0]\n";
  exit;
}


### adapter Seq
my $adapterseq='';
if ($param{'ADAPTER'} eq 'illumina'){
$adapterseq='AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC';
}elsif($param{'ADAPTER'} eq 'nextera'){
$adapterseq='CTGTCTCTTATACACATCT';
}


$param{'ALNDIR'} ="$param{'PROJECTNAME'}/STAR";
$param{'RESULTSDIR'} ="$param{'PROJECTNAME'}/results";
$param{'DATADIR'} ="$param{'PROJECTNAME'}/data";

########## use proper reFlat file for metrics #####################
my $REFFLAT; 
    for ($param{'GENOME'}) {
        when (/hg19/) { $REFFLAT = '/project/labprojects/share/hg19_data/Homo_sapiens.GRCh37.75.refflat' }
        when (/mm9/) { $REFFLAT = '/project/labprojects/share/mm9_data/refFlat.txt' }
        default       { $REFFLAT='' }
    }




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

	my $sample_name;
	$sample_name = $h{'SAMP_name'};
	$samples{$sample_name}{'mate1'}= "$param{'FASTQDIR'}/".$h{'mate1'}; 
	$samples{$sample_name}{'mate2'}= "$param{'FASTQDIR'}/".$h{'mate2'}; 

}


#run parallel jobs
my $pm=new Parallel::ForkManager($param{'JOBS'});


##########################  first passs #####################################


foreach (keys %samples)
{
	$pm->start and next;
		if($param{'ALIGN'}){
	#my $STARcmd = "STAR --runThreadN $param{'THREADS'}  --limitBAMsortRAM  10000000000 --runMode alignReads --genomeDir /project/labprojects/share/".$param{'GENOME'}."_STAR --outSAMtype BAM SortedByCoordinate --clip3pAdapterSeq AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --outReadsUnmapped Fastx --quantMode GeneCounts --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04   --alignIntronMin 20   --alignIntronMax 1000000  --alignSJoverhangMin 8   --alignSJDBoverhangMin 1 --sjdbScore 1 --alignMatesGapMax 1000000 --readFilesCommand zcat --readFilesIn $samples{$_}->{'mate1'} $samples{$_}->{'mate2'} --outFileNamePrefix $param{'ALNDIR'}/$_ --outWigType wiggle --outWigStrand Stranded --outWigNorm RPM";
	my $STARcmd = "STAR --runThreadN $param{'THREADS'}  --limitBAMsortRAM  10000000000 --runMode alignReads --genomeDir /project/labprojects/share/".$param{'GENOME'}."_STAR --outSAMtype BAM SortedByCoordinate --clip3pAdapterSeq $adapterseq --outReadsUnmapped Fastx --quantMode GeneCounts --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04   --alignIntronMin 20   --alignIntronMax 1000000  --alignSJoverhangMin 8   --alignSJDBoverhangMin 1 --sjdbScore 1 --alignMatesGapMax 1000000 --readFilesCommand zcat --readFilesIn $samples{$_}->{'mate1'} $samples{$_}->{'mate2'} --outFileNamePrefix $param{'ALNDIR'}/$_";

		print $STARcmd,"\n";		
		system($STARcmd);
		system("samtools index $param{'ALNDIR'}/".$_."Aligned.sortedByCoord.out.bam");
}

### dedup the samples
	if($param{'DEDUP'}){
		my $mrkDupcmd = "java -jar /opt/picard/picard.jar MarkDuplicates I= $param{'ALNDIR'}/".$_."Aligned.sortedByCoord.out.bam O= $param{'ALNDIR'}/".$_."_dedupped.bam  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$param{'ALNDIR'}/".$_.".mrkdupfull";
		print $mrkDupcmd,"\n";
		system($mrkDupcmd);
		system("grep LIBRARY -A1 $param{'ALNDIR'}/".$_.".mrkdupfull > $param{'ALNDIR'}/".$_.".mrkdup");

	}

	if($param{'METRICS'}){
		my $metricscmd = "java -jar /opt/picard/picard.jar CollectRnaSeqMetrics I= $param{'ALNDIR'}/".$_."Aligned.sortedByCoord.out.bam O=$param{'ALNDIR'}/".$_.".metricsfull REF_FLAT=$REFFLAT STRAND=SECOND_READ_TRANSCRIPTION_STRAND";
		print $metricscmd,"\n";
		system($metricscmd);
		system("grep PF_BASES -A1 $param{'ALNDIR'}/".$_.".metricsfull > $param{'ALNDIR'}/".$_.".metrics");
		}

	if($param{'LIBCOMPLEX'}){
		my $libcomplexcmd = "java -jar /opt/picard/picard.jar EstimateLibraryComplexity I= $param{'ALNDIR'}/".$_."Aligned.sortedByCoord.out.bam O=$param{'ALNDIR'}/".$_.".libcomplexfull";
		print $libcomplexcmd,"\n";
		system($libcomplexcmd);
		system("grep LIBRARY -A1 $param{'ALNDIR'}/".$_.".libcomplexfull > $param{'ALNDIR'}/".$_.".libcomplex");
		}

	if($param{'TIN'}){
		my $tincmd = "tin.py -i $param{'ALNDIR'}/".$_."Aligned.sortedByCoord.out.bam -r /project/labprojects/share/hg19_data/hg19_gencode17.bed";
		print $tincmd,"\n";
		system($tincmd);
		}


	print "$_ completed\n"; 
	$pm->finish;
	#exit;
}

$pm->wait_all_children;

	$Summ_cmd = "~/NGS/bin/parseSTARLog.pl --dir $param{'ALNDIR'} > $param{'RESULTSDIR'}/STAR_summmary.csv";
	print $Summ_cmd,"\n";
	system($Summ_cmd);

	$R_cmd = "Rscript --vanilla ~/NGS/bin/createsummary.R $param{'PROJECTNAME'}";
	print $R_cmd,"\n";
	system($R_cmd);
	


print "Run complete\n";

