#Tell Perl which modules to use
use strict;
use warnings;
use Getopt::Long;
use File::Find::Rule;
use Bio::Index::Fasta;

#Declare help and fasta as variables
my $help;
my $fasta;

#Write a usage message
my $usage = "Name:\n get_closest_relatives.pl\n\nContact:\n Ryan Cook <stxrc24\@nottingham.ac.uk>\n\nUsage:\n  perl get_closest_relatives.pl --query <mysequence.fasta> --inphared <directory> [optional arguments]\n\nMandatory arguments:\n  --query       -q <s>   Query fasta sequence to be compared against database.\n  --inphared    -i <s>   Location of inphared output directory. This is required to find the MASH sketch database and associated metadata.\n\nOptional arguments:\n  --help        -h <s>   This help.\n  --max         -m <n>   Maximum number of results to be outputted. Default is no maximum.\n  --threads     -t <n>   Number of MASH threads for parallelism. Default is 10.\n  --distance    -d <n>   Maximum MASH distance to be included. Default is 0.2. Distances over 0.1 should be interpreted with caution.\n  --pvalue      -p <n>   Maximum P-value cutoff to be included. Default is 0.001.\n  --fasta       -f       Use this flag for fasta sequences of matches to be written to output files (one file per query). By default, fasta sequences aren't written to output.\n\nIf you find this script useful, please cite our pre-print:\n\nCook et al. 2021. INfrastructure for a PHAge REference Database: Identification of large-scale biases in the current collection of phage genomes.\nbioRxiv doi: 10.1101/2021.05.01.442102\n";

#Get commandline options
GetOptions(

    #the number of cpus to use in Prokka as a commandline option
    'query=s'    => \(my $query = "jaboogapalooga"),
    
    #The name of the inphared directory, this is to find the MASH sketch database and associated metadata
    'inphared=s'  => \(my $inphared = "jaboogapalooga"),
    
    #help menu
    'help' => \$help,
    
    #Maximum number of results to be outputted, default is (basically) infinite
    'max=i'    => \(my $max = 10000000000000),
    
    #Threads for parallelism, default is 10
    'threads=i'    => \(my $threads = 10),
    
    #Maximum MASH distance, default is 0.2
    'distance=s'    => \(my $distance = 0.2),
    
    #Maximum MASH p-value, default is 0.001
    'pvalue=s'    => \(my $pvalue = 0.001),
    
    #The name of the fasta output file for sequences to be written to, if one is required
    'fasta'  => \$fasta,

);

#Just in case...
chomp $query;
chomp $inphared;
chomp $max;
chomp $threads;
chomp $distance;
chomp $pvalue;

#If any of the arguments are empty, die and print usage statement!
if (($query =~ m/^jaboogapalooga$/) || ($inphared =~ m/^jaboogapalooga$/)) {
    
    die $usage;
};

#If user specified help, let's die and print usage statement!
if ($help) {
    
    die $usage;
};

#Use the inphared directory path to get the MASH sketch database and associated metadata
#First, let's remove a slash from the end if there is one
$inphared =~ s/\/$//gi;

#Let's get the name of the database
my @database_files = File::Find::Rule->file()
                            ->name('*_genomes.fa.msh')
                            ->in("$inphared");

my $database = shift @database_files;
chomp $database;

#Let's get the name of the metadata file
my @metadata_files = File::Find::Rule->file()
                            ->name('*_data.tsv')
                            ->in("$inphared");

my $metadata = shift @metadata_files;
chomp $metadata;

#Let's get the name of the sequence database, in case user wants seqs
my @sequencedb_files = File::Find::Rule->file()
                            ->name('*_genomes.db')
                            ->in("$inphared");

my $sequencedb = shift @sequencedb_files;
chomp $sequencedb;

#Open the metadata file
open FILE1,"<$metadata";

#Produce hash to store metadata
my %data;

#Read metadata file line by line and store data in the hash we just made
while (my $line = <FILE1>) {

    chomp $line;
    my @cells = split /\t/, $line;
    my $key = shift @cells;
    my $value = join("\t", @cells);
    $data{$key} = $value;
}

#Convert input file to SeqIO object
my $inseq = Bio::SeqIO->new(-file   => "$query",
                            -format => 'Fasta', 
                           );

#Print headers to screen
print "Query\tHit\tMash_Distance\tIdentity(%)\tMatching_Kmers\tDescription\tGenus\tSubfamily\tFamily\tmolGC(%)\tLength\tHost\n";

#Loop for grabbing sequences one at a time and putting them through the process
while (my $sequence = $inseq->next_seq) {

    #Declare an array for storing accessions if the user wants seqs
    my @accession_array;
    
    #Grab the sequence ID while you're at it...
    my $id = $sequence->primary_id;
    
    #Name the temporary file after sequence ID...
    my $temp = join('', $id, "temp_seqfile.fa");

    #Creates a file to temporarily store sequences to be used in mash
    my $seq_out = Bio::SeqIO->new( -file   => ">$temp",
                                    -format => 'Fasta',
                                    );

    #Store sequence object temporarily in file, to be used in mash
    $seq_out->write_seq($sequence);

    #First make a sketch of the query for rapid lookup
    my $sketch_command = "mash.2 sketch -s 1000 -k 21 $temp".' >/dev/null 2>/dev/null';
    system("$sketch_command");

    #Grab the name of the new sketch file
    my $mash = $temp; 
    chomp $mash;
    $mash = $mash.".msh";

    #Run MASH to find nearest relatives to query, sort and write to temporary output file
    my $mash_command = "mash.2 dist -p $threads -v $pvalue -d $distance $database $mash | sort -k3 -g | head -n $max >xtemp.txt";
    system("$mash_command");
    
    #Delete the sketch of the query
    unlink $mash;

    #Open temporary MASH results file (see subroutine at bottom of script)
    open RES,"<xtemp.txt" or die;

    #Read it line by line
    while (my $line = <RES>) {

        #Always important for potential trailing strings
        chomp $line;

        #Split on tab characters
        my @data = split(/\t/,$line);

        #Assign useful MASH output fields to variables
        my $accession = $data[0]; 
        my $query = $data[1];
        my $prob = $data[3];
        my $mash_dist = $data[2];
        my $kmers = $data[4];
    
        #Let's exclude RefSeq duplicates for now
        if ($accession !~/NC_*/) { 
        
            #Lookup the relevant metadata using the accession and split on tabs into an array!
            my @entry =split(/\t/,$data{"$accession"});
            my $description = $entry[0]; 
            my $length = $entry[2];
            my $molgc = $entry[4];
            my $host = $entry[13];
            my $genus = $entry[15];
            my $subfamily = $entry[16];
            my $family = $entry[17];
    
            #Always important for potential trailing strings
            chomp $description;
            chomp $length;
            chomp $molgc;
            chomp $host;
            chomp $genus;
            chomp $subfamily;
            chomp $family;
    
            #Calculate percentage ID from MASH distance
            my $per_id = (1-$mash_dist)*100;
    
            #Print useful data to screen
            print "$id\t$accession\t$mash_dist\t$per_id\t$kmers\t$description\t$genus\t$subfamily\t$family\t$molgc\t$length\t$host\n";
        
            #Has user specified a fasta output file?
            if ($fasta) {

                #Push the accessions to an array
                push(@accession_array, $accession);
            }
        }
    }

    #Remove the temporary MASH output file
    unlink ("xtemp.txt");
    unlink ("$temp");

    #Has user specified a fasta output file?
    if ($fasta) {
        
        #Define a file for matches
        my $matches = join('', $id, "_matches.fa");
        
        #Make a BioSeq object with their specified name
        my $out = Bio::SeqIO->new(-format => 'Fasta',
                             -file => ">$matches");
                             
        #For each accession in the array
        foreach my $retrieval (@accession_array) {
    
            #Let's start grabbing the seqs!
            #Tell PERL where the index is (this is the DB of sequences you made with another script)
            my $inx = Bio::Index::Fasta->new(-filename => $sequencedb);
            #Grab the sequence based upon the ID
            my $seq = $inx->fetch($retrieval); # Returns Bio::Seq object
            #Write it to the output file
            $out->write_seq($seq);
        }
    }
}