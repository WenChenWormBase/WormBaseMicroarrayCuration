#!/usr/bin/perl -w

use strict;
use Ace;


#-----------------Find all SMD probes from citace-----------
print "Connecting to database...";
my $db = Ace->connect(-path => '/home/citace/citace',  -program => '/usr/local/bin/tace') || die print "Connection failure: ", Ace->error;
print "done\n";

my $query='QUERY FIND Microarray_results Microarray = Affymetrix_C.elegans_Genome_Array';
my @tmp=$db->find($query);

print scalar @tmp, " Affymetrix (GPL200) Microarray_results objects found\n";

my %os=();
foreach (@tmp) {
    $os{$_}=1;
}
#-------------All probes found ---------------------------------


#-------------print out ace and txt files. --------------------
my $ref = "WBPaper000XXXXX";
open (OUT1, ">WBPaper000XXXXXCMinus.ace") || die "cannot open $!\n"; #data for CitaceMinus
open (OUT2, ">WBPaper000XXXXXMrData.txt") || die "cannot open $!\n"; #data for mr_data
print OUT1 "\n";

my @exp;
my $gsm;
my $i = 0; #experiment id
my @stuff;
my $line;
my %PrintCondition;
my ($condition, $tmp_length, $prep, $geno); 
my $exp_rem="XXXXX. Frequency values reported here are back-transformed to the linear scale, i.e. values are 2^VALUE. Presence: PA = P(present), PS = M(intermediate), NP = A(absent)."; 


open (IN1, "SampleTitle.txt") || die "cannot open $!\n";
while ($line = <IN1>) {
    chomp($line);
    next unless ($line ne "");
    ($stuff[0], $stuff[1]) = split /\s+/, $line;
    $exp[$i] = $stuff[1];
    $gsm = $stuff[0];

    @tmp = split /__/, $stuff[1];
    $condition = $tmp[0];

    if ($condition =~ //) {
	$geno = "";
    } elsif  ($condition =~ //) {
	$geno = "";
    } elsif  ($condition =~ /nth-1/) {
	$geno = "";
    } elsif  ($condition =~ /xpa-1/) {
	$geno = "";
    }

    $prep = "";

    #print experiment and condition info. 
    print OUT1 "\nMicroarray_experiment : \"$ref:$exp[$i]\"\n";
    print OUT1 "Species\t\"Caenorhabditis elegans\"\n";
    print OUT1 "Microarray\t\"Affymetrix_C.elegans_Genome_Array\"\n";
    print OUT1 "Microarray_sample\t\"$ref:$condition\"\n";
    print OUT1 "Remark\t\"GEO record: $gsm\"\n";
    print OUT1 "Remark\t\"$exp_rem\"\n";        
    print OUT1 "Reference\t\"$ref\"\n\n";
    #print OUT1 "\n";
    if ($PrintCondition{$condition}) {
	#Already printed, skip
    } else {
	print OUT1 "Condition : \"$ref:$condition\"\n";

	print OUT1 "Sex\t\"Hermaphrodite\"\n";
	print OUT1 "Food\t\"OP50\"\n";
	print OUT1 "Species\t\"Caenorhabditis elegans\"\n";
	print OUT1 "Temperature\t\"20\"\n";
	print OUT1 "Tissue\t\"\"\n";
	print OUT1 "Life_stage\t\"\" \/\/\n";

	print OUT1 "Genotype\t\"$geno\"\n";
	print OUT1 "Preparation\t\"$prep\"\n";
	#print OUT1 "Remark\t\"$remark\"\n";
	print OUT1 "Reference\t\"$ref\"\n";
	print OUT1 "\n";
	$PrintCondition{$condition} = 1;
    }
    #Done printing. 
    $i++;
}
close (IN1);


#Done printing experiments and conditions.

my @inputFile;
$inputFile[0] = "GSEXXXXX_family.soft";
my @NumExp;
$NumExp[0] = 1; #total number of replicates
my @totalColumn;
$totalColumn[0] = ?;
my @dataHash;

my $c = 0;
$i = -1;
#my $tmp_length;
my ($mr_result, $valid_probe);
my %NotFound;

while ($c < 1) {
    open (IN1, "$inputFile[$c]") || die "cannot open $!\n";
    while ($i == -1) {#Look for the beginning of sample table		
	$line = <IN1>;
	if ($line =~ /sample_table_begin/){
	    $i++;
	}
    }

    #parsing sample tables. 
    while ($line = <IN1>) {
	chomp($line);
	if ($line =~ /sample_table_begin/) {
	    $i++;
	}   
	@tmp=split(/\t/, $line);
	$tmp_length = @tmp;
	next unless ($tmp_length == $totalColumn[$c]);

	#Check if probe is valid. 
	$mr_result = $tmp[0];	
	if ($os{$mr_result}) {
		$valid_probe = 1;
	} else {	    
		    unless (($tmp[0] =~ /^ID_REF/) || ($tmp[0] =~ /^AFFX/) || ($NotFound{$tmp[0]})) {
			print "Cannot find probe for $tmp[0].\n";
			$NotFound{$tmp[0]} = 1;
		    }
		    $valid_probe = 0;
	}
	#Done checking

	#print file for MrData. 
	if ($valid_probe == 1) {
	    #print table for mr_data;
	    $dataHash[0]= "\\N"; #A_vs_B_log_ratio
	    $dataHash[1]= "\\N"; #A_vs_B_SD
	    if (($tmp[1] =~ /NA/i) || ($tmp[1] =~ /Null/i)) {
		$dataHash[2]= "\\N";
	    } else {
		$dataHash[2]= 2**$tmp[1]; #frequency
	    }
	    $dataHash[3]= $NumExp[$c]; #Number_of_experiments
	    $dataHash[4]= "\\N"; #Confidence_level
	    $dataHash[5]= "\\N"; #P_value
	    #$dataHash[6]= "\\N"; #presence

	    if ($tmp[2] eq "P") {
	       $dataHash[6]= "PA"; #presence
	    } elsif ($tmp[2] eq "A") {
	       $dataHash[6]= "NP"; #presence
	    } elsif ($tmp[2] eq "M") {
	       $dataHash[6]= "PS"; #presence
	    } else {
	       $dataHash[6]= "\\N"; #presence
	    }

	    print OUT2 "$mr_result\t$ref:$exp[$i]\t$ref";
	    foreach (@dataHash) {
	    	print OUT2 "\t$_";
	    }
	    print OUT2 "\n";

	}
	#Done printing.
    }
    $c++;
    close (IN1);
}

#Print Expression_cluster based on Affy ID, data are not entered to Microarray_experiment objects. 
#If data need to be entered to Microarray_experiment objects, see the parsing script for WBPaper00034757.  

$c = 0;
@inputFile = ("TableS.csv"); 
my ($Cluster, $ClusDes);
while ($c < 1) {
    open (IN2, "$inputFile[$c]") || die "cannot open $!\n";
    while ($line = <IN2>) {
	chomp($line);
	if ($line =~ /Cluster/) {
	    ($stuff[3], $Cluster, $ClusDes) = split /:/, $line;	
	    print OUT1 "\nExpression_cluster : \"$ref:$Cluster\"\n";
	    print OUT1 "Reference\t\"$ref\"\n";
	    print OUT1 "Description\t\"$ClusDes\"\n";
	    print OUT1 "Algorithm\t\"\"\n";  
	    print OUT1 "Species\t\"Caenorhabditis elegans\"\n";
	    print OUT1 "Life_stage\t\"\" \/\/\n";
	    print OUT1 "GO_term\t\"\" \/\/\n";
	    print OUT1 "Anatomy_term\t\"\" \/\/\n";
	    print OUT1 "WBProcess\t\"\" \/\/\n";	    
	    print OUT1 "Regulated_by_gene\t\"\" \/\/\n";
	    #print OUT1 "Regulated_by_treatment\t\"Bacteria: \"\n";
	    print OUT1 "Regulated_by_molecule\t\"\" \/\/\n";
	    #print OUT1 "Based_on_WB_Release\n";
	    foreach (@exp) {
		print OUT1 "Microarray_experiment\t\"$ref:$_\"\n";
	    }	    
	} elsif ($line ne "") {
	    if ($os{$line}) {
		print OUT1 "Microarray_results\t\"$line\"\n";
	    } else {
		print "Cannot find oligo set for $line.\n";
	    }   
	}
    }
    $c++;
    close (IN2);
}
#Done printing. 

close (OUT1); 
close (OUT2); 

