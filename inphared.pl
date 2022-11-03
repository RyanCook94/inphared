#!/usr/bin/perl

#Give Perl the modules to use. If one of these is not installed, the script will fail. However, Perl should print to the screen which one is missing.
use Getopt::Long;
use POSIX qw(strftime);
use Time::localtime;
use File::Find;
use Bio::SeqIO;
use Bio::Index::Fasta;
use List::MoreUtils qw(uniq);
use v5.23;
use strict;
use warnings;

#Date string to time stamp on databases created 
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date = strftime"%e%b%Y", gmtime;

#Tidy up the date
chomp $date;
$date =~ s/ //gi;

#Declare help as a variable
my $help;

#Get commandline options
GetOptions(

    #the number of cpus to use in Prokka as a commandline option
    'cpus=i'    => \(my $cpu = 8),
    
    #The name of an output directory, if not, one will be produced automatically
    'outdir=s'  => \(my $outdir = "inphared_"."$date"),
    
    #the use of an exlcusion file, if desired
    'exclusion=s'   => \(my $exclusion = 0),
    
    #the use of phrogs for annotation, if desired
    'PHROG=s'   => \(my $phrog = 0),
    
    #help menu
    'help' => \$help,  

);

#Remove a slash from the end of the output directory if the user has specified one
$outdir =~ s/(\/)$//;

#If user specified help, let's die and print something!
if ($help) {
    
    die "Name:\n INPHARED: INfrastructre for a PHAge REference Database v1.8\n\nContact:\n Ryan Cook <stxrc24\@nottingham.ac.uk>\n\nUsage:\n  perl inphared.pl [options]\n\nOptions:\n  --help        -h   This help.\n  --exclusion   -e   Pipe-delimited text file of accessions that should be excluded. This is recommended.\n  --outdir      -o   Output directory for files to be written to. Default is inphared_"."$date"." (changes daily).\n  --cpus        -c   Number of cpus to be used in Prokka gene calling. Default is 8.\n  --PHROG       -P   Absolute path to phrog HMM database for annotations (can be found on millardlab website). If not specified, Prokka will be run with the --noanno flag (no annotations, just gene calling).\n";
}

#Say hello to the user
say "INPHARED: INfrastructre for a PHAge REference Database v1.8\n";

#Get full paths to mash, efetch, esearch, efilter and prokka. If one of these is not installed and available in PATH, the script will fail. It will tell you which it cannot find
say "Searching for dependencies required for this script to run.\n";

#First MASH
say "Searching for MASH...";

#Declare this as a variable
my $mash_path;

#Search for mash
my $tool_mash = "mash";
$mash_path = `which $tool_mash`;
chomp $mash_path;

#Check to see if it found it
if ($mash_path eq "") {

    #If we didn't find it, let's try mash.2
    my $tool_mash2 = "mash.2";
    $mash_path = `which $tool_mash2`;
    chomp $mash_path;
    
    }

die "$tool_mash was not found. Need to install mash in PATH.\n" unless ($mash_path);
say "MASH has been found: $mash_path\n";

#Then efetch
say "Searching for efetch...";
my $tool_efetch = "efetch";
my $efetch_path = `which $tool_efetch`;
chomp $efetch_path;
die "$tool_efetch was not found. Need to install efetch in PATH.\n" unless ($efetch_path);
say "efetch has been found: $efetch_path\n";

#Then esearch
say "Searching for esearch...";
my $tool_esearch = "esearch";
my $esearch_path = `which $tool_esearch`;
chomp $esearch_path;
die "$tool_esearch was not found. Need to install esearch in PATH.\n" unless ($esearch_path);
say "esearch has been found: $esearch_path\n";

#Then efilter
say "Searching for efilter...";
my $tool_efilter = "efilter";
my $efilter_path = `which $tool_efilter`;
chomp $efilter_path;
die "$tool_efilter was not found. Need to install efilter in PATH.\n" unless ($efilter_path);
say "efilter has been found: $efilter_path\n";

#And finally prokka
say "Searching for prokka...";
my $tool_prokka = "prokka";
my $prokka_path = `which $tool_prokka`;
chomp $prokka_path;
die "$tool_prokka was not found. Need to install prokka in PATH.\n" unless ($prokka_path);
say "prokka has been found: $prokka_path\n";

#If we've made it this far, then all dependencies should be installed
say "All dependencies have been found, let's roll!\n";

#Check to see if a GenomesDB directory already exists
say "Searching for GenomesDB directory..."; 

#Perform check
if (! -d "GenomesDB") {
    
    #Couldn't find it...
    say "GenomesDB not found, so we're going to make it!\n";
    system("mkdir GenomesDB");
    
    } else { 
    
    #Found it!
    say "Found GenomesDB, all good.\n";
};

#Give list of accessions for "megaphages" that may not ordinarily be downloaded
my $meta = "BK014188|BK013603|BK014009|BK013977|BK013744|BK013754|KT462700|KT462701|BK013715|BK013911|BK013882|BK013443|BK014194|BK013654|BK013891|BK013480|BK013490|BK013910|BK014035|BK013987|BK013827|BK014032|BK014031|BK013870|BK013499|BK013511|BK013792|BK013396|BK014165|BK013591|BK013689|BK013379|BK013690|BK014142|BK013413|BK014001|BK013454|BK013452|BK014018|BK013787|BK014000|BK013951|BK013534|BK013927|BK013848|BK013455|BK013831|BK013906|BK014119|BK013639|BK014120|BK013737|BK014122|BK014189|BK013485|BK013473|BK014083|BK014028|BK013533|BK013980|BK014023|BK013747|BK013772|BK013806|BK013748|BK013821|BK013678|BK013438|BK013925|BK013945|BK013666|BK013923|BK013439|BK014118|BK013919|BK014176|BK013820|BK014027|BK013686|BK013523|BK013786|BK013570|BK013931|BK013989|BK013947|BK013811|BK013561|BK013407|BK013889|BK013477|BK013546|BK013553|BK013969|BK013926|BK013466|BK013917|BK013790|BK013540|BK013447|BK013384|BK013825|BK013733|BK013595|BK014137|BK013668|BK014090|BK013988|BK013419|BK014121|BK014041|BK013572|BK013920|BK014014|BK013844|BK013858|BK013970|BK013691|BK014139|BK014133|BK013890|BK013484|BK013512|BK013516|BK013809|BK014043|BK013954|BK013851|BK014007|BK014024|BK014193|BK013712|BK013599|BK014144|BK013486|BK013857|BK014030|BK013563|BK013694|BK014064|BK013446|BK013538|BK013864|BK014114|BK013641|BK014152|BK013968|BK013872|BK013389|BK013878|BK013943|BK013993|BK013914|BK013671|BK013596|BK013734|BK014080|BK014196|BK013465|BK013834|BK013887|BK013476|BK013830|BK013932|BK013481|BK013937|BK013567|BK013423|BK014216|BK014161|BK014199|BK013658|BK013601|BK013798|BK013897|BK013659|BK013963|BK013426|BK013950|BK013953|BK013424|BK013871|BK013769|BK013557|BK013504|BK013488|BK013467|BK013556|BK013836|BK013441|BK013510|BK013552|BK013385|BK013527|BK014026|BK013955|BK014148|BK013985|BK013964|BK013451|BK014138|BK013801|BK013894|BK014184|BK014178|BK013588|BK014179|BK014048|BK013578|BK013997|BK013393|BK013942|BK013470|BK013500|BK013837|BK013394|BK013541|BK013939|BK013402|BK013568|BK013459|BK013881|BK013559|BK013999|BK013457|BK013550|BK013430|BK013507|BK013829|BK013410|BK013506|BK013566|BK013879|BK013400|BK013479|BK014005|BK013867|BK013404|BK013381|BK013501|BK013828|BK013794|BK013974|BK013814|BK014004|BK013817|BK013391|BK013388|BK013383|BK013560|BK013775|BK013450|BK013884|BK013813|BK013930|BK013796|BK013432|BK013991|BK013940|BK013962|BK013965|BK013408|BK013468|BK013779|BK013558|BK013449|BK013956|BK013976|BK013909|BK013458|BK014029|BK014042|BK013896|BK013803|BK013982|BK013482|BK013461|BK013928|BK013397|BK013839|BK013866|BK013411|BK013973|BK013808|BK013386|BK013895|BK013845|BK014022|BK014015|BK014033|BK013494|BK013832|BK013417|BK013409|BK013960|BK013495|BK013493|BK013705|BK013924|BK013435|BK014045|BK013535|BK013565|BK014012|BK013903|BK013816|BK013463|BK013804|BK013843|BK014010|BK013632|BK013998|BK013826|BK013600|BK013399|BK014084|BK013704|BK014071|BK014135|BK014168|BK014091|BK013944|BK013756|BK013873|BK014213|BK013571|BK013789|BK013412|BK013979|BK013642|BK013853|BK013763|BK013938|BK013886|BK013908|BK013774|BK013474|BK014019|BK013371|BK013542|BK013425|BK013475|BK013401|BK013576|BK013418|BK013776|BK014044|BK013865|BK014102|BK013992|BK013652|BK014210|BK014131|BK013647|BK014108|BK013723|BK013771|BK013517|BK014008|BK013431|BK013941|BK013536|BK013835|BK013984|BK014020|BK013532|BK013815|BK013520|BK013496|BK013805|BK013971|BK013905|BK013731|BK013589|BK013648|BK013841|BK013529|BK013922|BK013530|BK013509|BK013788|BK013370|BK013525|BK013508|BK013810|BK013868|BK013544|BK013797|BK014187|BK013978|BK013850|BK013874|BK013387|BK013860|BK013912|BK014034|BK013994|BK014013|BK013842|BK013898|BK013907|BK014126|BK013515|BK013869|BK013782|BK013398|BK013904|BK014200|BK013670|BK014117|BK014107|BK013861|BK013555|BK013537|BK014025|BK014065|BK013952|BK014038|BK013818|BK013543|BK014124|BK013587|BK014163|BK013586|BK013380|BK013740|BK013590|BK013584|BK013518|BK013681|BK014181|BK013634|BK013460|BK013428|BK013633|BK013876|BK013359|BK013703|BK014162|BK013730|BK013675|BK013657|BK013742|BK014059|BK014053|BK013741|BK013360|BK013575|BK013362|KT462698|BK013638|BK014052|BK013717|BK013422|BK014097|BK013364|BK013372|BK013693|BK013680|BK014076|BK014057|BK013478|BK013966|BK013579|BK014109|BK014074|BK013663|BK014055|BK014159|BK014154|BK013604|BK014078|BK013795|BK013444|BK013793|BK013462|BK014169|BK013621|BK013673|BK013710|BK014141|BK014167|BK014130|BK013709|BK014166|BK013716|BK013644|BK014003|BK013746|BK013738|BK014125|BK013750|BK013762|BK014204|BK013627|BK014145|BK014150|BK013367|BK013645|BK013777|BK014021|BK013472|BK014006|BK013735|BK013582|BK013687|BK013656|BK014218|BK014201|BK013852|BK013736|BK013757|BK014087|BK013683|BK013649|BK013958|BK013390|BK013892|BK013610|BK014103|BK013729|BK013382|BK013718|BK013913|BK013791|BK013933|BK013732|BK013724|BK013698|BK013573|BK014073|BK014068|BK014058|BK014183|BK013613|BK014143|BK014086|BK014202|BK013660|BK013725|BK013577|BK013807|BK013719|BK014147|BK014132|BK014082|BK013415|BK013497|BK013406|BK013833|BK013959|BK014075|BK013800|BK014217|BK013551|BK013505|BK013519|BK013946|BK013513|BK013883|BK013355|BK013672|BK013916|BK014051|BK013822|BK013990|BK013503|BK013759|BK013608|BK013662|BK014115|BK013592|BK014113|BK013366|BK013487|BK013464|BK013436|BK013961|BK014077|LT821717|BK013679|BK013651|BK013628|BK014128|BK013623|BK013996|BK013569|BK013489|BK013983|BK013824|BK013491|BK014170|BK013356|BK014172|BK014197|BK013614|BK013626|BK014191|BK013605|BK013395|BK013863|BK014089|BK014085|BK014190|BK013665|BK013812|BK013799|BK013612|BK014160|BK013593|BK013655|BK014070|BK014066|BK013521|BK013899|BK013445|BK013838|BK013893|BK013606|BK014208|BK014164|BK014099|BK013745|BK014177|BK013531|BK013986|BK013434|BK014180|BK013706|BK013674|BK013854|BK013900|BK013676|BK014049|BK013957|BK013374|BK013637|BK013699|BK013784|BK014105|BK013753|BK013682|BK013846|BK013427|BK013819|BK013783|BK014158|BK014101|BK014063|BK013646|BK014146|BK013376|BK014116|BK014061|BK013617|BK013598|BK013581|BK013619|BK014195|BK014054|BK014046|BK014060|BK014095|BK013700|BK013378|BK014198|BK013669|BK014215|BK013758|BK014207|BK013615|BK013721|BK013767|BK014205|BK013585|BK014056|BK013429|BK013616|BK013739|BK014096|BK013405|BK014079|BK013877|BK013847|BK013677|BK014088|BK013702|BK013692|BK013640|BK014050|BK013363|BK013752|BK013361|BK013369|BK013802|BK013622|BK013849|BK013921|BK014155|BK013554|BK013743|BK013635|BK013643|BK013625|BK013618|BK013358|BK014067|BK013365|BK013901|BK013433|BK013528|BK013823|BK014036|BK013859|BK013778|BK013972|BK013856|BK014016|BK013650|BK014221|BK013722|BK013594|BK014174|BK014106|BK013620|BK013929|BK013442|BK013547|BK013357|BK013580|BK013708|BK013707|BK013636|BK013995|BK013437|BK014040|BK014111|BK013936|BK013697|BK013453|BK013440|BK013377|BK013714|BK014069|BK014062|BK013597|BK013667|BK014134|BK013624|BK013602|BK013701|BK013696|BK013583|BK014153|BK014186|BK014182|BK014220|BK014047|BK013720|BK014175|BK014157|BK013548|BK014100|BK013766|BK013688|BK014127|BK014037|BK013483|BK013549|BK013918|BK013773|BK013749|BK013862|BK013471|BK013368|BK014173|BK014092|BK013685|BK013611|BK013751|BK014219|BK013949|BK013492|BK013403|BK014104|BK014149|BK014209|BK014094|BK013761|BK013416|BK013967|BK013607|BK014017|BK013781|BK013456|BK013526|BK013514|BK013564|BK014002|BK014185|BK014156|BK013661|BK014093|BK014214|BK013764|BK014203|BK014151|BK013629|BK013373|BK013770|BK014011|BK013935|BK013948|BK014039|BK013902|BK013975|BK014212|BK013765|BK013726|BK014072|BK013880|BK014171|BK013711|BK013840|BK013855|BK013915|BK013414|BK013785|BK013885|BK014206|BK013695|BK013664|BK013630|BK014123|BK014129|BK013755|BK014211|BK013609|BK014081|BK013420|BK013684|BK013768|BK013981|BK013888|BK013502|BK013875|BK013934|BK013392|BK013524|BK014112|BK013728|BK013574|BK013448|BK013760|BK013653|BK013469|BK013780|BK013522|BK014140|BK014136|BK013545|BK013727|BK013539|BK013562|BK014110|BK013375|BK014098|BK013631|BK013713|BK014192|BK013498|BK013421|MF361639|MG592455|MG592402|MG592394|MG592427|MG592615|JQ995537|MH616963|MH552500|AP013538|AP013544|KT997847|AP013545|AP013543|GU943031|KX397280|AP013549|HQ633062|AP013547|AP013546|GU943068|AP013539|GU943054|AP013542|AP013541|AP013548|AP013540|AC171169|AY576273|MF417929|KK213166|AE006468|AY319521|CP038625|LR797314|CP023680|CP023686|CP006891|CP014526|CP019275|CP052639|HM208303|MW175890|GU339467|EF489910|AM040673|MN335248|HM066936|AF059242|DQ836063|JX045650|JX625144|JX045649|AF334111|LR756501|CACVHF010000000|LR745206|LR745208|LR756500|LR756501|LR756502|LR756503|LR756504|LR756508|LR756511|MK250015|MK250016|MK250017|MK250018|MK250019|MK250020|MK250021|MK250022|MK250023|MK250024|MK250025|MK250026|MK250027|MK250028|MK250029";
chomp $meta;

#Produce the output directory if it doesn't already exist
if (! -d "$outdir") {

    #make it!
    system("mkdir $outdir");
    }

#Let the user know the details of the output directory
say "Output will be written to $outdir. This can be changed using the -o flag in the commandline. If not specified, the default is inphared_"."$date".".\n";

#Set name of phage database output. Note, these are the raw pre-filtered Genbank entries for all downloaded genomes
my $phagedb = $outdir."/".$date."_phages_downloaded_from_genbank.gb"; 

#Is an exclusion list being used?
my $exclusion_list = "peanuts";
if($exclusion) {
    
    #Open and read it
    open EXC, "<$exclusion" or die;
    $exclusion_list = <EXC>;
    chomp $exclusion_list;
    say "Exclusion list specified as $exclusion\nExclusion lists can be specified using the -e flag. This should be a pipe-delimited file of unwanted accessions.\n";

} else {

    #Give a warning
    say "No exclusion list specified.\nWARNING some eroneous genomes will not be excluded.\nExclusion lists can be specified using the -e flag. This should be a pipe-delimited file of unwanted accessions.\n";
}

#Set name of log file and open it for writing. This will contain the excluded genome IDs
my $log = $outdir."/".$date.".log";
open LOG, ">$log";

#Create vConTACT2 input files for writing
#Declare name for mapping file
my $mapping = $outdir."/".$date."_vConTACT2_gene_to_genome.csv";

#Check to see if a vConTACT2 mapping file already exists
say "Searching for vConTACT2 mapping file...";

#Perform check
if (! -f "$mapping") {
    
    #Couldn't find it...
    say "$mapping not found, so we're going to make it!\n";
    
    #Open mapping file for writing
    open(MAP,">$mapping") or die;
    
    #Give it some headers
    print MAP "protein_id,contig_id,keywords\n";

    } else { 
    
    #Found it!
    say "Found $mapping, all good.\n";
};

##############################################################

#Create an empty hash to store hex codes for the annotation files
my %hex_hash;

#Store Unclassified and Unspecified as pale grey. Change this hex code if you'd rather them a different colour
$hex_hash{"Unclassified"} = "#D3D3D3";
$hex_hash{"Unspecified"} = "#D3D3D3";

#Create empty hashes to store accession/genus data, accession/sub-family data, accession/family data, accession/host data, and accession/description data
my %genus_hash;
my %sub_family_hash;
my %family_hash;
my %lowest_taxa_hash;
my %host_hash;
my %description_hash;
my %length_hash;

#Create empty arrays to store accessions, genera, sub-family, family and host data
my @accession_array;
my @genus_array;
my @sub_family_array;
my @family_array;
my @lowest_taxa_array;
my @host_array;

#Create vConTACT2 annotation files for writing

#Declare name for GENUS annotation file
my $vcontact_genus = $outdir."/".$date."_vConTACT2_genus_annotations.tsv";

#Check to see if file already exists
say "Searching for $vcontact_genus file...";

#Perform check
if (! -f "$vcontact_genus") {
    
    #Couldn't find it...
    say "$vcontact_genus not found, so we're going to make it!\n";
    
    #Open GENUS file for writing
    open(VGENUS,">$vcontact_genus") or die;
    
    #Give it some headers
    print VGENUS "Node_ID\tDescription\tGenus\tColour (Genus)\n";

    } else { 
    
    #Found it!
    say "Found $vcontact_genus, all good.\n";
};

#Declare name for SUB-FAMILY annotation file
my $vcontact_sub = $outdir."/".$date."_vConTACT2_subfamily_annotations.tsv";

#Check to see if file already exists
say "Searching for $vcontact_sub file...";

#Perform check
if (! -f "$vcontact_sub") {
    
    #Couldn't find it...
    say "$vcontact_sub not found, so we're going to make it!\n";
    
    #Open SUB-FAMILY file for writing
    open(VSUB,">$vcontact_sub") or die;
    
    #Give it some headers
    print VSUB "Node_ID\tDescription\tSub-family\tColour (Sub-family)\n";

    } else { 
    
    #Found it!
    say "Found $vcontact_sub, all good.\n";
};

#Declare name for FAMILY annotation file
my $vcontact_family = $outdir."/".$date."_vConTACT2_family_annotations.tsv";

#Check to see if file already exists
say "Searching for $vcontact_family file...";

#Perform check
if (! -f "$vcontact_family") {
    
    #Couldn't find it...
    say "$vcontact_family not found, so we're going to make it!\n";
    
    #Open FAMILY file for writing
    open(VFAMILY,">$vcontact_family") or die;
    
    #Give it some headers
    print VFAMILY "Node_ID\tDescription\tFamily\tColour (Family)\n";

    } else { 
    
    #Found it!
    say "Found $vcontact_family, all good.\n";
};

#Declare name for LOWEST TAXA annotation file
my $vcontact_lowest_taxa = $outdir."/".$date."_vConTACT2_lowest_taxa_annotations.tsv";

#Check to see if file already exists
say "Searching for $vcontact_lowest_taxa file...";

#Perform check
if (! -f "$vcontact_lowest_taxa") {
    
    #Couldn't find it...
    say "$vcontact_lowest_taxa not found, so we're going to make it!\n";
    
    #Open LOWEST TAXA file for writing
    open(VLOWEST,">$vcontact_lowest_taxa") or die;
    
    #Give it some headers
    print VLOWEST "Node_ID\tDescription\tLowest Taxa\tColour (Lowest taxa)\n";

    } else { 
    
    #Found it!
    say "Found $vcontact_lowest_taxa, all good.\n";
};

#Declare name for HOST annotation file
my $vcontact_host = $outdir."/".$date."_vConTACT2_host_annotations.tsv";

#Check to see if file already exists
say "Searching for $vcontact_host file...";

#Perform check
if (! -f "$vcontact_host") {
    
    #Couldn't find it...
    say "$vcontact_host not found, so we're going to make it!\n";
    
    #Open HOST file for writing
    open(VHOST,">$vcontact_host") or die;
    
    #Give it some headers
    print VHOST "Node_ID\tDescription\tHost\tColour (Host)\n";

    } else { 
    
    #Found it!
    say "Found $vcontact_host, all good.\n";
};

##############################################################

#Show user the date
say "Today is: $date";

#Check if the database already exists. This is the file which contains the raw genbank files
say "Checking for database named: $phagedb";

#Perform check
if (-e "$phagedb") { 
    
    #It already exists, so move on. No point re-downloading the same sequences. If you think your database is incomplete, please delete the files and try a fresh run with the script
    say "$phagedb already exists. If you think it's incomplete, please delete output files and try a clean run.\n"; 
    
    } else { 

    #If it doesn't exist, then download the genomes! Current search parameters download genomes with the PHG identifier and minimum/maximum length cutoffs. These are also filtered later on in the script
    say "$phagedb doesn't yet exist, so let's make it!\n";
    say "This step may take some. Please free free to put the kettle on...\n";
    system ("$esearch_path -db nucleotide -query \"$meta\" | $efilter_path -query \"1417:800000 [SLEN] \" | $efetch_path -format gb >$phagedb");
    system ("$esearch_path -db nucleotide -query \"gbdiv_PHG\"[prop] | $efilter_path -query \"1417:800000 [SLEN] \" | $efetch_path -format gb >>$phagedb"); 
};

#Create output file with relevant extension for tab separated data relating to phage genomes 
my $output = $outdir."/".$date."_data.tsv";

#And make one that excludes RefSeq sequences...
my $output_without_refseq = $outdir."/".$date."_data_excluding_refseq.tsv";

#Check to see if this .tsv file already exists
say "Searching for $output";

#Perform check
if (-f $output) { 
    
    #It already exists
    say "Output file $output already exists";
    
    } else { 
    
    #It doesn't exist yet, so let's make it!
    say "$output doesn't yet exist. Will now process...\n";
    
    #Let the user know how many cpus will be used in Prokka
    say "Will use $cpu cpus in Prokka pipeline. This can be changed using the -c flag in the commandline. If not specified, the default is 8.\n";

    my @input;
    push @input, $output;
    push @input, $phagedb;
    push @input, $date;

    #run subroutine 
    &filter_genomes(@input);
  
};

#Create fastadb for convenient lookup 
&create_fastadb($date);

#########################################################

##############CREATE ITOL ANNOTATION FILES###############

#Declare name for GENUS file
my $itol_genus = $outdir."/".$date."_itol_genus_annotations.txt";

#Check to see if this file already exists
say "Searching for $itol_genus file...";

#Perform check
if (! -f "$itol_genus") {
    
    #Couldn't find it...
    say "$itol_genus not found, so we're going to make it!\n";
    
    #Open annotation file for writing
    open(ITOL_GENUS,">$itol_genus") or die;
    
    } else { 
    
    #Found it!
    say "Found $itol_genus, all good.\n";
};

#Sort array alphabetically and grab unique values only
my @sorted_genus = sort @genus_array;
@sorted_genus = uniq(@sorted_genus);

#Create an array to store the genus colour codes
my @genus_code_array;

#Add information to genus annotation file
print ITOL_GENUS "#This annotation file produces a coloured ring around ITOL trees which shows viral genera of nodes.\n\nDATASET_COLORSTRIP\n\nSEPARATOR TAB\n\nDATASET_LABEL\tViral Genera\n\nCOLOR\t#ff0000\n\nCOLOR_BRANCHES\t0\n\nLEGEND_TITLE\tViral Genera\n\n";

#Add the legend labels
print ITOL_GENUS "LEGEND_LABELS\t";

foreach my $uniq_genus (@sorted_genus) {

    #Trailing strings
    chomp $uniq_genus;
    
    #Lookup the hex colour code that matches the genus
    my $genus_code = $hex_hash{"$uniq_genus"};
    
    #Push this into an array
    push @genus_code_array, $genus_code;
    
    #Print them all
    print ITOL_GENUS "$uniq_genus\t"; 
}

#Add the legend colours
print ITOL_GENUS "\n\nLEGEND_COLORS\t";

foreach my $genus_codes (@genus_code_array) {

    #Trailing strings
    chomp $genus_codes;
    
    #Print them all
    print ITOL_GENUS "$genus_codes\t";
}

#Add the legend shapes
print ITOL_GENUS "\n\nLEGEND_SHAPES\t";

#This makes sure the number of shapes is right
foreach my $genus_codes (@genus_code_array) {

    #Trailing strings
    chomp $genus_codes;
    
    #Print them all
    print ITOL_GENUS "2\t";
}

#Add the optional settings
print ITOL_GENUS "\n\nSTRIP_WIDTH\t15\n\nDATA\n\n";

#Now to get the individual data!
#Go through accessions one at a time, lookup genus, lookup hexcode. Print in tsv. V simple!
foreach my $accession (@accession_array) {

    #Trailing strings
    chomp $accession;
    
    #Lookup the genus that matches the accession
    my $genus_lookup = $genus_hash{"$accession"};
    
    #Lookup the hex colour code that matches the genus
    my $hex_code = $hex_hash{"$genus_lookup"};
    
    print ITOL_GENUS "$accession\t$hex_code\t$genus_lookup\n";
}

#####################################################

#Declare name for SUB-FAMILY file
my $itol_sub = $outdir."/".$date."_itol_subfamily_annotations.txt";

#Check to see if this file already exists
say "Searching for $itol_sub file...";

#Perform check
if (! -f "$itol_sub") {
    
    #Couldn't find it...
    say "$itol_sub not found, so we're going to make it!\n";
    
    #Open annotation file for writing
    open(ITOL_SUB,">$itol_sub") or die;
    
    } else { 
    
    #Found it!
    say "Found $itol_sub, all good.\n";
};

#Sort array alphabetically and grab unique values only
my @sorted_sub = sort @sub_family_array;
@sorted_sub = uniq(@sorted_sub);

#Create an array to store the sub-family colour codes
my @sub_code_array;

#Add information to sub-family annotation file
print ITOL_SUB "#This annotation file produces a coloured ring around ITOL trees which shows viral sub-family of nodes.\n\nDATASET_COLORSTRIP\n\nSEPARATOR TAB\n\nDATASET_LABEL\tViral Sub-family\n\nCOLOR\t#ff0000\n\nCOLOR_BRANCHES\t0\n\nLEGEND_TITLE\tViral Sub-family\n\n";

#Add the legend labels
print ITOL_SUB "LEGEND_LABELS\t";

foreach my $uniq_sub (@sorted_sub) {

    #Trailing strings
    chomp $uniq_sub;
    
    #Lookup the hex colour code that matches the sub-family
    my $sub_code = $hex_hash{"$uniq_sub"};
    
    #Push this into an array
    push @sub_code_array, $sub_code;
    
    #Print them all
    print ITOL_SUB "$uniq_sub\t"; 
}

#Add the legend colours
print ITOL_SUB "\n\nLEGEND_COLORS\t";

foreach my $sub_codes (@sub_code_array) {

    #Trailing strings
    chomp $sub_codes;
    
    #Print them all
    print ITOL_SUB "$sub_codes\t";
}

#Add the legend shapes
print ITOL_SUB "\n\nLEGEND_SHAPES\t";

#This makes sure the number of shapes is right
foreach my $sub_codes (@sub_code_array) {

    #Trailing strings
    chomp $sub_codes;
    
    #Print them all
    print ITOL_SUB "2\t";
}

#Add the optional settings
print ITOL_SUB "\n\nSTRIP_WIDTH\t15\n\nDATA\n\n";

#Now to get the individual data!
#Go through accessions one at a time, lookup sub-family, lookup hexcode. Print in tsv. V simple!
foreach my $accession (@accession_array) {

    #Trailing strings
    chomp $accession;
    
    #Lookup the sub-family that matches the accession
    my $sub_lookup = $sub_family_hash{"$accession"};
    
    #Lookup the hex colour code that matches the sub-family
    my $hex_code = $hex_hash{"$sub_lookup"};
    
    print ITOL_SUB "$accession\t$hex_code\t$sub_lookup\n";
}

#########################################################

#Declare name for FAMILY file
my $itol_family = $outdir."/".$date."_itol_family_annotations.txt";

#Check to see if this file already exists
say "Searching for $itol_family file...";

#Perform check
if (! -f "$itol_family") {
    
    #Couldn't find it...
    say "$itol_family not found, so we're going to make it!\n";
    
    #Open annotation file for writing
    open(ITOL_FAM,">$itol_family") or die;
    
    } else { 
    
    #Found it!
    say "Found $itol_family, all good.\n";
};

#Sort array alphabetically and grab unique values only
my @sorted_family = sort @family_array;
@sorted_family = uniq(@sorted_family);

#Create an array to store the family colour codes
my @family_code_array;

#Add information to family annotation file
print ITOL_FAM "#This annotation file produces a coloured ring around ITOL trees which shows viral family of nodes.\n\nDATASET_COLORSTRIP\n\nSEPARATOR TAB\n\nDATASET_LABEL\tViral Family\n\nCOLOR\t#ff0000\n\nCOLOR_BRANCHES\t0\n\nLEGEND_TITLE\tViral Family\n\n";

#Add the legend labels
print ITOL_FAM "LEGEND_LABELS\t";

foreach my $uniq_family (@sorted_family) {

    #Trailing strings
    chomp $uniq_family;
    
    #Lookup the hex colour code that matches the family
    my $family_code = $hex_hash{"$uniq_family"};
    
    #Push this into an array
    push @family_code_array, $family_code;
    
    #Print them all
    print ITOL_FAM "$uniq_family\t"; 
}

#Add the legend colours
print ITOL_FAM "\n\nLEGEND_COLORS\t";

foreach my $family_codes (@family_code_array) {

    #Trailing strings
    chomp $family_codes;
    
    #Print them all
    print ITOL_FAM "$family_codes\t";
}

#Add the legend shapes
print ITOL_FAM "\n\nLEGEND_SHAPES\t";

#This makes sure the number of shapes is right
foreach my $family_codes (@family_code_array) {

    #Trailing strings
    chomp $family_codes;
    
    #Print them all
    print ITOL_FAM "2\t";
}

#Add the optional settings
print ITOL_FAM "\n\nSTRIP_WIDTH\t15\n\nDATA\n\n";

#Now to get the individual data!
#Go through accessions one at a time, lookup family, lookup hexcode. Print in tsv. V simple!
foreach my $accession (@accession_array) {

    #Trailing strings
    chomp $accession;
    
    #Lookup the sub-family that matches the accession
    my $family_lookup = $family_hash{"$accession"};
    
    #Lookup the hex colour code that matches the family
    my $hex_code = $hex_hash{"$family_lookup"};
    
    print ITOL_FAM "$accession\t$hex_code\t$family_lookup\n";
}

#########################################################

#Declare name for LOWEST TAXA file
my $itol_lowest_taxa = $outdir."/".$date."_itol_lowest_taxa_annotations.txt";

#Check to see if this file already exists
say "Searching for $itol_lowest_taxa file...";

#Perform check
if (! -f "$itol_lowest_taxa") {
    
    #Couldn't find it...
    say "$itol_lowest_taxa not found, so we're going to make it!\n";
    
    #Open annotation file for writing
    open(ITOL_LOW,">$itol_lowest_taxa") or die;
    
    } else { 
    
    #Found it!
    say "Found $itol_lowest_taxa, all good.\n";
};

#Sort array alphabetically and grab unique values only
my @sorted_lowest_taxa = sort @lowest_taxa_array;
@sorted_lowest_taxa = uniq(@sorted_lowest_taxa);

#Create an array to store the lowest taxa colour codes
my @lowest_taxa_code_array;

#Add information to lowest taxa annotation file
print ITOL_LOW "#This annotation file produces a coloured ring around ITOL trees which shows the lowest viral taxa of nodes.\n\nDATASET_COLORSTRIP\n\nSEPARATOR TAB\n\nDATASET_LABEL\tLowest Viral Taxa\n\nCOLOR\t#ff0000\n\nCOLOR_BRANCHES\t0\n\nLEGEND_TITLE\tLowest Viral Taxa\n\n";

#Add the legend labels
print ITOL_LOW "LEGEND_LABELS\t";

foreach my $uniq_lowest_taxa (@sorted_lowest_taxa) {

    #Trailing strings
    chomp $uniq_lowest_taxa;
    
    #Lookup the hex colour code that matches the lowest taxa
    my $lowest_taxa_code = $hex_hash{"$uniq_lowest_taxa"};
    
    #Push this into an array
    push @lowest_taxa_code_array, $lowest_taxa_code;
    
    #Print them all
    print ITOL_LOW "$uniq_lowest_taxa\t"; 
}

#Add the legend colours
print ITOL_LOW "\n\nLEGEND_COLORS\t";

foreach my $lowest_taxa_codes (@lowest_taxa_code_array) {

    #Trailing strings
    chomp $lowest_taxa_codes;
    
    #Print them all
    print ITOL_LOW "$lowest_taxa_codes\t";
}

#Add the legend shapes
print ITOL_LOW "\n\nLEGEND_SHAPES\t";

#This makes sure the number of shapes is right
foreach my $lowest_taxa_codes (@lowest_taxa_code_array) {

    #Trailing strings
    chomp $lowest_taxa_codes;
    
    #Print them all
    print ITOL_LOW "2\t";
}

#Add the optional settings
print ITOL_LOW "\n\nSTRIP_WIDTH\t15\n\nDATA\n\n";

#Now to get the individual data!
#Go through accessions one at a time, lookup lowest taxa, lookup hexcode. Print in tsv. V simple!
foreach my $accession (@accession_array) {

    #Trailing strings
    chomp $accession;
    
    #Lookup the sub-family that matches the accession
    my $lowest_taxa_lookup = $lowest_taxa_hash{"$accession"};
    
    #Lookup the hex colour code that matches the lowest taxa
    my $hex_code = $hex_hash{"$lowest_taxa_lookup"};
    
    print ITOL_LOW "$accession\t$hex_code\t$lowest_taxa_lookup\n";
}

#########################################################

#Declare name for HOST file
my $itol_host = $outdir."/".$date."_itol_host_annotations.txt";

#Check to see if this file already exists
say "Searching for $itol_family file...";

#Perform check
if (! -f "$itol_host") {
    
    #Couldn't find it...
    say "$itol_host not found, so we're going to make it!\n";
    
    #Open annotation file for writing
    open(ITOL_HOST,">$itol_host") or die;
    
    } else { 
    
    #Found it!
    say "Found $itol_host, all good.\n";
};

#Sort array alphabetically and grab unique values only
my @sorted_host = sort @host_array;
@sorted_host = uniq(@sorted_host);

#Create an array to store the host colour codes
my @host_code_array;

#Add information to family annotation file
print ITOL_HOST "#This annotation file produces a coloured ring around ITOL trees which shows bacterial host of nodes.\n\nDATASET_COLORSTRIP\n\nSEPARATOR TAB\n\nDATASET_LABEL\tBacterial Host\n\nCOLOR\t#ff0000\n\nCOLOR_BRANCHES\t0\n\nLEGEND_TITLE\tBacterial Host\n\n";

#Add the legend labels
print ITOL_HOST "LEGEND_LABELS\t";

foreach my $uniq_host (@sorted_host) {

    #Trailing strings
    chomp $uniq_host;
    
    #Lookup the hex colour code that matches the host
    my $host_code = $hex_hash{"$uniq_host"};
    
    #Push this into an array
    push @host_code_array, $host_code;
    
    #Print them all
    print ITOL_HOST "$uniq_host\t"; 
}

#Add the legend colours
print ITOL_HOST "\n\nLEGEND_COLORS\t";

foreach my $host_codes (@host_code_array) {

    #Trailing strings
    chomp $host_codes;
    
    #Print them all
    print ITOL_HOST "$host_codes\t";
}

#Add the legend shapes
print ITOL_HOST "\n\nLEGEND_SHAPES\t";

#This makes sure the number of shapes is right
foreach my $host_codes (@host_code_array) {

    #Trailing strings
    chomp $host_codes;
    
    #Print them all
    print ITOL_HOST "2\t";
}

#Add the optional settings
print ITOL_HOST "\n\nSTRIP_WIDTH\t15\n\nDATA\n\n";

#Now to get the individual data!
#Go through accessions one at a time, lookup host, lookup hexcode. Print in tsv. V simple!
foreach my $accession (@accession_array) {

    #Trailing strings
    chomp $accession;
    
    #Lookup the sub-family that matches the accession
    my $host_lookup = $host_hash{"$accession"};
    
    #Lookup the hex colour code that matches the host
    my $hex_code = $hex_hash{"$host_lookup"};
    
    print ITOL_HOST "$accession\t$hex_code\t$host_lookup\n";
}

#########################################################

#Declare name for NODE LABELS file
my $itol_node = $outdir."/".$date."_itol_node_label_annotations.txt";

#Check to see if this file already exists
say "Searching for $itol_node file...";

#Perform check
if (! -f "$itol_node") {
    
    #Couldn't find it...
    say "$itol_node not found, so we're going to make it!\n";
    
    #Open annotation file for writing
    open(ITOL_NODE,">$itol_node") or die;
    
    } else { 
    
    #Found it!
    say "Found $itol_node, all good.\n";
};

#Print useful data to the file
print ITOL_NODE "#This annotation file adds more informative labels for nodes, replacing the accession number with the virus' name.\n\nLABELS\n\nSEPARATOR TAB\n\nDATA\n\n";

#Now to get the individual data!
#Go through accessions one at a time, lookup description. Print in tsv. V simple!
foreach my $accession (@accession_array) {

    #Trailing strings
    chomp $accession;
    
    #Lookup the description that matches the accession
    my $description_lookup = $description_hash{"$accession"};
    
    print ITOL_NODE "$accession\t$description_lookup ($accession)\n";
}

#########################################################

#Declare name for LENGTH BAR file
my $itol_length = $outdir."/".$date."_itol_length_annotations.txt";

#Check to see if this file already exists
say "Searching for $itol_length file...";

#Perform check
if (! -f "$itol_length") {
    
    #Couldn't find it...
    say "$itol_length not found, so we're going to make it!\n";
    
    #Open annotation file for writing
    open(ITOL_LEN,">$itol_length") or die;
    
    } else { 
    
    #Found it!
    say "Found $itol_length, all good.\n";
};

#Print useful data to the file
print ITOL_LEN "#This annotation file adds bars which indicate the length of the viral genome.\n\nDATASET_SIMPLEBAR\n\nSEPARATOR TAB\n\nDATASET_LABEL\tGenome Length (KB)\n\nCOLOR\t#3186f5\n\nWIDTH\t100\n\nDATA\n\n";

#Now to get the individual data!
#Go through accessions one at a time, lookup length. Print in tsv. V simple!
foreach my $accession (@accession_array) {

    #Trailing strings
    chomp $accession;
    
    #Lookup the description that matches the accession
    my $length_lookup = $length_hash{"$accession"};
    
    #Change it to KB
    $length_lookup = $length_lookup / 1000;
    
    print ITOL_LEN "$accession\t$length_lookup\n";
}

#########################################################

############CREATE_MASH_DB###############################

#Give name for MASH database
my $mash = $outdir."/".$date."_genomes.fa.msh"; 

#Check if it already exists first
if (-f $mash) { 

    #If it exists, then don't make it again
    say "MASH database $mash already exists. Re-run to make a new one.";
    
    } else {
    
    #If it doesn't, then make it!
    say "Will create MASH database $mash";
    
    #Pass process to Perl subroutine
    &create_mash($date);
    
    #Tell user what's going on
    say "Finished sketching MASH database $mash\n";
};

#########################################################

########Extract RefSeq Sequences Only####################

#Name output file for RefSeq sequences only
my $refseq = $outdir."/".$date."_refseq_genomes.fa";

#Tell user what you're doing
say "Extracting REFSEQ genomes...";

#Check if the RefSeq output file already exists
if (-f $refseq) { 
    
    #If it does, don't repeat it
    say "REFSEQ genomes have already been extracted. Re-run to create new set of RefSeq sequences.\n";
    
    } else {

    #If not, then let's grab them
    my $refseq_out = Bio::SeqIO->new( -file   => ">$refseq",
                               -format => 'fasta',
                             );

    #Open the tab delimited file for reading
    open REF,"<$output" or die "couldn't open tab delimited file for reading\n";

    #Read it line by line
    while (my $line = <REF>) {
 
        #Split on tab characters
        my @line =split(/\t/,$line) ;
        
        #If the accession begins with NC_, you know it's a RefSeq genome
        if ($line[0] =~ m/NC_*/ ) {
            
            #Store these in an array
            my @temp;
            push @temp,$refseq_out; 
            push @temp,$line[0];
            push @temp,$date;
            
            #Run the extract_fasta subroutine to grab the relevant sequences and write to output
            &extract_fasta(@temp);
        }
    }
    
    #Tell User what's happening
    say "Finished extracting RefSeq sequences and written to $refseq\n";
}

#This is currently the endpoint of the script
say "Finished performing analysis. Thank you for using INPHARED.pl.";

#########################################################

#############EXTRACT_FASTA SUBROUTINE####################

#Run Perl subroutine extract_fasta
sub extract_fasta { 
    
    #Declare variables
    my $date = $_[2]; 
    my $out = $_[0];
    my $id = $_[1];
    my $fastadb = $outdir."/".$date."_genomes.db";
    my $inx = Bio::Index::Fasta->new(-filename => "$fastadb");
    
    #Returns Bio::Seq object
    my $seq = $inx->fetch($id);
    
    #Write the sequence to output
    $out->write_seq($seq);
}

#########################################################

############CREATE_MASH SUBROUTINE#######################

#Create MASH database in Perl subroutine
sub create_mash {
    
    #Declare the variables
    my $date = $_[0];
    my $fasta = $outdir."/".$date."_genomes.fa";
    my $run = "$mash_path sketch -i -s 1000 $fasta"; 
    
    #Pass the command to system
    system($run);
}

#########################################################

############CREATE_FASTADB SUBROUTINE####################

#Create fastadb in perl subroutine
sub create_fastadb { 

    #Declare variables
    my $date = $_[0];
    my $fastadb = $outdir."/".$date."_genomes.db";
    my $in = $outdir."/".$date."_genomes.fa";
    
    #Tell user what you're doing
    say "Creating fasta index $fastadb for easy lookup of sequences...";
    
    #Check if fastadb already exists
    if (-f $fastadb) { 
        
        #If it does exist, let user know
        say "$fastadb already exists, delete it to run again.\n";
        
        } else { 
        
        #If it doesn't exist, then make it!
        my $Index_File_Name = $fastadb or die $!;
        my $inx = Bio::Index::Fasta->new(-filename => $Index_File_Name,
                                                    -write_flag => 1);
        $inx->make_index($in);
        
        #Print it for the user
        say "$fastadb has been produced.\n";
    }
}

#########################################################

#############FILTER_GENOMES SUBROUTINE###################

#Run Perl subroutine to filter genomes and grab useful information about them
sub filter_genomes { 

    #Assign variables
    my $output = $_[0];
    my $phagedb = $_[1];
    my $date = $_[2];
    
    #Produce table for upload
    my $html = "$outdir"."/"."$date"."_millardlab_website_table.txt";
    
    #open table for upload 
    open(HTML,">$html") or die;
    
    #Print useful headers into table for upload
    print HTML "Accession\tDescription\tGenome Length (KB)\tmolGC (%)\tGenus\tSub-family\tFamily\tHost\n";
    
    #Open TSV file for writing useful phage information
    open(TAXA,">$output") or die;
    
    #Open TSV file for writing useful phage information EXCLUDING REFSEQ SEQUENCES
    open(TAXA_without_refseq,">$output_without_refseq") or die;

    #Print useful headers into TSV file
    print TAXA "Accession\tDescription\tClassification\tGenome Length (bp)\tJumbophage\tmolGC (%)\tMolecule\tModification Date\tNumber CDS\tPositive Strand (%)\tNegative Strand (%)\tCoding Capacity (%)\tLow Coding Capacity Warning\ttRNAs\tHost\tLowest Taxa\tGenus\tSub-family\tFamily\tOrder\tClass\tPhylum\tKingdom\tRealm\tBaltimore Group\tGenbank Division\tIsolation Host (beware inconsistent and nonsense values)\n";
    
    #And the one excluding refseq...
    print TAXA_without_refseq "Accession\tDescription\tClassification\tGenome Length (bp)\tJumbophage\tmolGC (%)\tMolecule\tModification Date\tNumber CDS\tPositive Strand (%)\tNegative Strand (%)\tCoding Capacity(%)\tLow Coding Capacity Warning\ttRNAs\tHost\tLowest Taxa\tGenus\tSub-family\tFamily\tOrder\tClass\tPhylum\tKingdom\tRealm\tBaltimore Group\tGenbank Division\tIsolation Host (beware inconsistent and nonsense values)\n";

    #Read genomes in from the phage DB
    say "Reading $phagedb for filtering...\n";

    #Read it in as a Bio::SeqIO object...
    my $seq_in = Bio::SeqIO->new( -file   => "$phagedb",
                                -format => "genbank",
                                );

    #Produce new fasta file to write filtered genomes into
    my $seq_out = Bio::SeqIO->new( -format => 'fasta',
                                -file   => ">$outdir"."/"."$date"."_genomes.fa",
                                );
                                
    #Produce new fasta file to write filtered genomes into
    my $exc_refseq_fasta = $outdir."/".$date."_genomes_excluding_refseq.fa";
    my $seq_out_excluding_refseq = Bio::SeqIO->new( -format => 'fasta',
                                -file   => ">$exc_refseq_fasta",
                                );
                                
    #Define a name for the vConTACT2 proteins file
    my $proteins_file = $outdir."/".$date."_vConTACT2_proteins.faa";
    
    #Produce new fasta file to write proteins into (excluding RefSeq)
    my $protein = Bio::SeqIO->new( -format => 'fasta',
                                -file   => ">$proteins_file",
                                );
    
    #Go through the phageDB Genbank file one entry at a time
    while (my $inseq_obj = $seq_in->next_seq) {

        #Grab and chomp useful information
        my $desc = $inseq_obj->description;
        my $acc =  $inseq_obj->accession_number;
        chomp $desc;
        chomp $acc;
        
        #Does the accession already exist?
        unless ($acc ~~ @accession_array) {

        #Now time for the actual filtering!
        #If the accession is in our exclusion list, print that accession in the log file so user can access it
        if ($acc =~ m/$exclusion_list/) { 
            print LOG "$acc \texclusion\n";
        };

        #If the description contains "Complete" and "Genome", is a valid sequence object and isn't in the exclusion list...
        if ((($desc =~ m/complete/i) && ($desc =~ /Genome/i) && ($inseq_obj->length > 100) && ($inseq_obj->seq) && ($acc !~ m/$exclusion_list/)) ||
        #OR... The sequence is over 10kb, a valid sequence object and isn't in the exclusion list...
        (($inseq_obj->length > 10000) && ($inseq_obj->seq) && ($acc !~ m/$exclusion_list/)) ||
        #OR... The description contains "Complete" and "Sequence", sequence is over 10kb, is a valid sequence object and isn't in the exclusion list...
        (($inseq_obj->length > 10000) && ($desc =~ m/complete/i) && ($desc =~ /sequence/i) && ($inseq_obj->seq) && ($acc !~ m/$exclusion_list/)) ||
        #OR... The accession is one of these manually curated complete phage genomes that may otherwise be missed (i.e. they're small)
        ($acc =~ m/NC_002180|NC_008355|NC_002180|NC_019920|NC_019922|NC_027984|NC_028993|NC_028994|NC_028651|NC_025471|NC_026013|NC_025444|NC_029032|NC_028998|NC_004173|NC_004174|NC_004175|NC_004171|NC_004170|NC_003301|NC_003300|NC_003299|NC_003714|NC_003715|NC_003716|NC_023586|NC_021866|NC_020869|NC_026612|NC_026613|NC_024711|NC_026582|JN377905/) && ($inseq_obj->seq) && ($acc !~ m/$exclusion_list/)) {

            #If it passed filtering, we can write sequence to output fasta file!
            $seq_out->write_seq($inseq_obj);
            
            #Show user the sequence which passed filtering
            say "Filtering: $acc \t $desc\n";

            #set variables to print out 
            my $primary_id = $inseq_obj->display_id;
            my $division = $inseq_obj->division;
            my $length = $inseq_obj->length;
            my $jumbo;
            my $species_string = $inseq_obj->species->node_name;
            my @classification = $inseq_obj->species->classification;
            my @dates = $inseq_obj->get_dates;
            my $molecule = $inseq_obj->molecule;
            my $seq = $inseq_obj->seq;
            my $gc = $seq =~ tr/gcGC//;
            my $pergc = $gc / $length * 100;
            $pergc = sprintf("%.3f", $pergc);
            
            #Get classification as a scalar
            my $scalar_classification = join(" ", @classification);
            
            ########################################################
            
            #Declare a variable for isolation host
            my $isolation_host = "Unspecified";
            
            #Get the features from the original GB file
            for my $gb_feature ($inseq_obj->get_SeqFeatures) {
            
                #Find the source tag
                for my $tag ($gb_feature->get_all_tags) {
                
                    #Always important
                    chomp $tag;
                    
                    #Is the tag supposedly the host?
                    if ($tag =~ m/^host$|^lab_host&/) {
                    
                        #Get the value of it
                        my @values = $gb_feature->get_tag_values("$tag");
                        
                        my $val = $values[0];
                        chomp $val;
                        
                        $isolation_host = $val;
                    }
                    
                }
            
            }
            
            ########################################################
            
            #Determine whether phage would be classified as a "Jumbophage"
            if ($length >= 200000) {
            
                #It is jumbo
                $jumbo = "TRUE";
            
            } else {
            
                #It is not jumbo
                $jumbo = "FALSE";
            }

            #Define name for individual sequence file
            my $seqfile = "GenomesDB/".$acc.".fasta";
            
            #Check if sequence exists in GenomesDB directory
            if (-f $seqfile) { 

                #If it exists, then don't make it again
                say "$seqfile already exists, no need to repeat it!";
    
                } else {
    
                #If it doesn't, then make it!
                say "$seqfile does not yet exist. So let's write it!";
                
                #Define a Bio::SeqIO object for writing
                my $ind_seq_out = Bio::SeqIO->new( -format => 'fasta',
                                        -file   => ">$seqfile",
                                            );
                
                #Write sequence to the file!
                $ind_seq_out->write_seq($inseq_obj);
            };
            
            #Now to call genes on the individual sequence files to produce useful output in a standardised format
            
            #Give the name of the Prokka output directory
            my $prokka_directory = "GenomesDB/".$acc;
            
            #Check if the Prokka output directory for each sequence exists
            if (-d $prokka_directory) { 

                #If it exists, then don't make it again
                say "$prokka_directory has already been through Prokka, no need to repeat it!\n";
    
                } else {
    
                #If it doesn't, then make it!
                say "Now processing $prokka_directory with Prokka.\n";
                
                #Unless it's one of the accessions that matches the code15 phages, annotate normally
                unless($acc =~ m/MK250015|MK250016|MK250017|MK250018|MK250019|MK250020|MK250021|MK250022|MK250023|MK250024|MK250025|MK250026|MK250027|MK250028|MK250029/) {
                
                    #Let's work out if the PHROG flag has been used
                    if($phrog) {
                    
                        #Define and run the Prokka command (Note that the number of cpus is set as a commandline argument, with a default value of 8)
                        my $prokka = "$prokka_path --quiet --hmms $phrog --outdir $prokka_directory --prefix $acc --locustag $acc $seqfile --cpus $cpu";
                        system("$prokka");
                        
                        } else {
                        
                        #Define and run the Prokka command (Note that the number of cpus is set as a commandline argument, with a default value of 8)
                        my $prokka = "$prokka_path --noanno --outdir $prokka_directory --prefix $acc --locustag $acc $seqfile --cpus $cpu";
                        system("$prokka");
                        
                        }
                    
                    } else {
                    
                    #Let's work out if the PHROG flag has been used
                    if($phrog) {
                    
                        #And if it is one, we'll annotate using code15
                        my $prokka = "$prokka_path --quiet --hmms $phrog --outdir $prokka_directory --prefix $acc --locustag $acc $seqfile --cpus $cpu --gcode 15";
                        system("$prokka");
                        
                        } else {
                        
                        #And if it is one, we'll annotate using code15
                        my $prokka = "$prokka_path --quiet --noanno --outdir $prokka_directory --prefix $acc --locustag $acc $seqfile --cpus $cpu --gcode 15";
                        system("$prokka");
                        
                        }
                    
                    }
            };
            
            #Use the Prokka .gbf files to grab useful data
            
            #Define the annotation file I'm using to grab this information from
            my $gbf = "GenomesDB/".$acc."/".$acc.".gbf";
            
            if (-f $gbf) {
            
                #Do nothing
                
            } else {
            
                #Then gbf becomes something a bit different
                $gbf = "GenomesDB/".$acc."/".$acc.".gbk";
            }
            
            #Read it as a Bio::SeqIO object
            my $gbf_object = Bio::SeqIO->new(-file => "$gbf");         
            my $gbf_seq_object = $gbf_object->next_seq;
            
            #Set up some useful counters starting at zero         
            my $plus = 0; #Number of CDS on plus strand
            my $minus = 0; #Number of CDS on minus strand
            my $coding_count = 0; #Sum length of coding features
            my $trna_count = 0; #Number of tRNAs
            my @seen; #This array stores unique barcodes for each feature, to ensure they're not counted twice when calculating coding capacity
            
            #Take the features and read them as a feature object one by one
            for my $feat_object ($gbf_seq_object->get_SeqFeatures) {
                
                #If it's a CDS (this is for strand bias)...
                if ($feat_object->primary_tag eq 'CDS') {
                    
                    #Get the strand and add to the relevant counter
                    my $strand = $feat_object->strand();
                    if ($strand > 0) {$plus++};
                    if ($strand < 0) {$minus++};              
                    }
                    
                #If it's a tRNA (this is for tRNA count)...
                if ($feat_object->primary_tag eq 'tRNA') {
                
                    #Add to the relevant counter
                    $trna_count++
                    }
                
                #If it's anything other than source (this is for coding capacity)...
                if ($feat_object->primary_tag eq 'source') {
                    
                    #Do nothing, it's the source
                    } else {
                    
                    #This isn't source. So we want it. Start by defining variables for it's length and coordinates
                    my $feat_length;
                    my $start;
                    my $end;
                    
                    #Some features are split over the ends of the genome, which makes them slightly trickier to deal with. This identifies them
                    if ($feat_object->location->isa('Bio::Location::SplitLocationI')) {
                        
                        #Assign a variable for it's length (taken from two sub-sequences added together)
                        my $splitty = 0;
                
                        #This takes both parts of the split feature
                        for my $location ($feat_object->location->sub_Location ) {
                            
                            #Assign variables for coordinates
                            my $split_start = $location->start();
                            my $split_end = $location->end();
                            chomp $split_start;
                            chomp $split_end;
                            my $split_length = $split_end - $split_start;
                            $splitty = $splitty + $split_length;
                            $feat_length = $splitty;
                            $start = "NA";
                            $end = "NA";
                            }
                    
                    #Else statement for features which are NOT split
                    } else {
            
                        #Grab it's coordinates and length
                        $feat_length = $feat_object->length();
                        $start = $feat_object->start();
                        $end = $feat_object->end();
                        chomp $feat_length;
                        chomp $start;
                        chomp $end;
                        }
                    
                    #This barcode system acts as a fail-safe way to make sure that features aren't duplicated when counting coding capacity (i.e. same feature listed as gene and CDS if user changes the Prokka command)
                    my $barcode = "$start"."_"."$end"."_"."$feat_length";
                    chomp $barcode;
                    
                    #Turn the seen array into a hash
                    my %seen;
                    $seen{$_}++ for (@seen);
                    
                    #Has the barcode already been seen? (i.e. searching for duplicated features)
                    unless(exists($seen{$barcode})) {
                    
                        #If not, add it to array and add to counter!
                        push @seen, $barcode;
                        $coding_count = $coding_count + $feat_length;
                        }         
                    }
                }
            
            #Add both strand counters together for total number of CDS
            my $total_CDS = $plus + $minus;
            
            #Calculate strand usage as a percentage
            my $per_plus = $plus / $total_CDS * 100;
            my $per_minus = $minus / $total_CDS * 100;
                        
            #Now let's express as percentage
            my $coding_capacity = ($coding_count/$length)*100;
            
            #Let's add a warning flag if this is less than 50
            my $coding_warning;
            
            if($coding_capacity < 50) {
            
                $coding_warning = "WARNING";
            
            } else {
                
                $coding_warning = "NA";
            }
            
            #The script will grab host data from the phage description where possible (although some nonsense hosts will inevitably come through)
            my $host = "Unspecified";
            my $pos;
            
            #This finds the word that precludes "prophage" and takes it as the host name. Obviously not perfect, but best way to automate this!
            if ($species_string =~ m/ prophage/ ) {
                $host = substr($species_string,0,index($species_string,"prophage"));
                $host =~ m/([A-Z]*)/i;
                $host = $1;
                $host =~ s/ |\t|\n//gi;
                }
                
            #This finds the word that precludes "bacteriophage" and takes it as the host name. Obviously not perfect, but best way to automate this!
            if ($species_string =~ m/ bacteriophage/ ) {
                $host = substr($species_string,0,index($species_string,"phage"));
                $host =~ m/([A-Z]*)/i;
                $host = $1;
                $host =~ s/ |\t|\n//gi;
                } 
            
            #This finds the word that precludes "phage" and takes it as the host name. Obviously not perfect, but best way to automate this!
            if ($species_string =~ m/ phage/ ) {
                $host = substr($species_string,0,index($species_string,"phage"));
                $host =~ m/([A-Z]*)/i;
                $host = $1;
                $host =~ s/ |\t|\n//gi;
                }
                
            #A small number are capitalised
            if ($species_string =~ m/ Phage/ ) {
                $host = substr($species_string,0,index($species_string,"Phage"));
                $host =~ m/([A-Z]*)/i;
                $host = $1;
                $host =~ s/ |\t|\n//gi;
                }

            #Same as above, but for the word "virus"
            if ($species_string =~ m/ virus/ ) {
                $host = substr($species_string,0,index($species_string,"virus"));
                $pos = index($host," ");
                $host = substr($host,0,$pos+1);
                $host =~ s/ |\t|\n//gi;
                }
                
            #A few manual ones...
            if ($species_string =~ m/Cyanophage S-|Cyanophage KBS-|Cyanophage Syn/) {
            
                $host = "Synechococcus";
            }
            
            if ($species_string =~ m/Cyanophage P-|Cyanophage 9515-10a|Cyanophage MED4-117|Cyanophage NATL|Cyanophage PSS2|Cyanophage SS120-1/) {
            
                $host = "Prochlorococcus";
            }
            
            if ($species_string =~ m/Enterobacteriaphage/) {
            
                $host = "Enterobacteria";
            }
            
            if ($species_string =~ m/Mycobacteriophage/) {
            
                $host = "Mycobacterium";
            }
            
            if ($species_string =~ m/Sulfolobus/) {
            
                $host = "Sulfolobus";
            }
            
            if ($species_string =~ m/Vibriophage/) {
            
                $host = "Vibrio";
            }
            
            if (($species_string =~ m/Enterobacteria/) && ($scalar_classification =~ m/Escherichia/)) {
            
                $host = "Escherichia";
            }
            
            #Manual list of nonsense hosts which are changed to unspecified
            if ($host =~ m/^Methylophilaceae$|^Idiomarinaceae$|^Chlorobiaceae$|^Alteromonadaceae$|^abalone$|^ANMV-1$|^Archaeal$|^Deep$|^Halocynthia$|^His$|^His2$|^IAS$|^Los$|^Marine$|^Salicola$|^Stx$|^Hardygib1$|^CrAss-like$|^Proteobacteria$|^Aquificae$|^Erysiphe$|^Methylophilales$|^Sulfolobales$|^TM$|^uncultured$/i) {
            
                #Change it to Unspecified
                $host = "Unspecified";
            }
            
            #And some other edits for hosts (mostly typos)
            $host =~ s/Enterobacterial/Enterobacteria/gi;
            $host =~ s/Enterococus/Enterococcus/gi;
            $host =~ s/Entercoccus/Enterococcus/gi;
            $host =~ s/Eschericha/Escherichia/gi;
            $host =~ s/Panteoa/Pantoea/gi;
            $host =~ s/Pseudomonad/Pseudomonas/gi;
            $host =~ s/Baciilus/Bacillus/gi;
            $host =~ s/Kelbsiella/Klebsiella/gi;
            $host =~ s/Pectinobacterium/Pectobacterium/gi;
            $host =~ s/Stretomyces/Streptomyces/gi;
            $host =~ s/Synechoccus/Synechococcus/gi;
            
            #This script will grab the lowest taxa rank available from classification where possible (although some nonsense taxa will inevitably come through)
            my $lowest_taxa = "Unclassified";
            my $genus = "Unclassified";
            my $sub_family = "Unclassified";
            my $family = "Unclassified";
            my $order = "Unclassified";
            my $class = "Unclassified";
            my $phylum = "Unclassified";
            my $kingdom = "Unclassified";
            my $realm = "Unclassified";
            my $baltimore = "Unclassified";
            
            #Go through each word in the classification, one at a time (left to right)
            foreach my $preword (@classification) {
                
                #Trailing characters
                chomp $preword;
                
                #Split the word on whitespace. Some classifications go into the array strangely. This helps with that
                my @word = split ' ', $preword;
                
                #Go through each word
                foreach my $word (@word) {
                
                    #Potential trailing strings
                    chomp $word;
                
                    #Grab the genus where possible (these are words that end in virus but aren't only the word virus)
                    if($word =~ m/^[A-Z][a-z]*virus$/ && $word !~ m/^virus$/i) {
                
                        #This is my genus
                        $genus = $word;
                    }
                
                    #Grab the subfamily where possible (these are words that end in virinae)
                    if($word =~ m/^[A-Z][a-z]*virinae$/) {
                
                        #This is my sub-family
                        $sub_family = $word;
                    }
                
                    #Grab the family where possible (these are words that end in viridae)
                    if($word =~ m/^[A-Z][a-z]*viridae$/) {
                
                        #This is my family
                        $family = $word;
                    }
                    
                    #Grab the order where possible (these are words that end in virales)
                    if($word =~ m/^[A-Z][a-z]*virales$/) {
                
                        #This is my phylum
                        $order = $word;
                    }
                    
                    #Grab the class where possible (these are words that end in viricetes)
                    if($word =~ m/^[A-Z][a-z]*viricetes$/) {
                
                        #This is my phylum
                        $class = $word;
                    }
                    
                    #Grab the phylum where possible (these are words that end in viricota, is important for baltimore group)
                    if($word =~ m/^[A-Z][a-z]*viricota$/) {
                
                        #This is my phylum
                        $phylum = $word;
                    }
                    
                    #Grab the kingdom where possible (these are words that end in virae)
                    if($word =~ m/^[A-Z][a-z]*virae$/) {
                
                        #This is my phylum
                        $kingdom = $word;
                    }
                    
                    #Grab the realm where possible (these are words that end in viria)
                    if($word =~ m/^[A-Z][a-z]*viria$/) {
                
                        #This is my realm
                        $realm = $word;
                    }
                }
            }
            
            #Iteratively grab the lowest level of taxonomy
            #First with genus
            if($genus !~ m/^Unclassified$/) {
            
                #Then assign it...
                $lowest_taxa = $genus;
            }
            
            #Then sub-family
            if($genus =~ m/^Unclassified$/ && $sub_family !~ m/^Unclassified$/) {
            
                #Then assign it...
                $lowest_taxa = $sub_family;
            }
            
            #And finally family
            if($genus =~ m/^Unclassified$/ && $sub_family =~ m/^Unclassified$/ && $family !~ m/^Unclassified$/) {
            
                #Then assign it...
                $lowest_taxa = $family;
            }
            
            #Let's use the realm to try and grab the Baltimore group
            if($phylum =~ m/^Dividoviricota$|^Preplasmiviricota$|^Taleaviricota$|^Uroviricota$/) {
            
                #It must be Group I
                $baltimore = "Group I";
            
            }
            
            if($phylum =~ m/^Hofneiviricota$|^Phixviricota$|^Saleviricota$/) {
            
                #It must be Group I
                $baltimore = "Group II";
            
            }
            
            if($phylum =~ m/^Duplornaviricota$/) {
            
                #It must be Group I
                $baltimore = "Group III";
            
            }
            
            if($phylum =~ m/^Lenarviricota$/) {
            
                #It must be Group I
                $baltimore = "Group IV";
            
            }
            
            #Print useful information to TSV file
            print TAXA "$primary_id\t$species_string\t@classification\t$length\t$jumbo\t$pergc\t$molecule\t@dates\t$total_CDS\t$per_plus\t$per_minus\t$coding_capacity\t$coding_warning\t$trna_count\t$host\t$lowest_taxa\t$genus\t$sub_family\t$family\t$order\t$class\t$phylum\t$kingdom\t$realm\t$baltimore\t$division\t$isolation_host\n";
                        
            #Sort formating for html links
            my $postfix =  "\">$primary_id</a>";
            my $prefix = "<a href=\"https://www.ncbi.nlm.nih.gov/nuccore/";
            my $url = $prefix.$primary_id.$postfix;
            
            #Get the values for html table
            my $length_kb = $length / 1000;
            my $html_length = int($length_kb + 0.5);
            my $html_pergc = int($pergc + 0.5);
            
            print HTML "$url\t$species_string\t$html_length\t$html_pergc\t$genus\t$sub_family\t$family\t$host\n";
            
            #Declare a variable for the genus hex code
            my $genus_hex;
            my $sub_family_hex;
            my $family_hex;
            my $host_hex;
            my @genus_hex_chars = ('0'..'9', 'A'..'F');
            my $genus_hex_len = 6;
            my @sub_hex_chars = ('0'..'9', 'A'..'F');
            my $sub_hex_len = 6;
            my @fam_hex_chars = ('0'..'9', 'A'..'F');
            my $fam_hex_len = 6;
            my @host_hex_chars = ('0'..'9', 'A'..'F');
            my $host_hex_len = 6;
            
            #Check if the viral genus has been seen before
            if(exists($hex_hash{$genus})) {
            
                #Then look up the relevant hexcode as a hash-key lookup
                $genus_hex = $hex_hash{"$genus"};
            
            } else {
            
                #And if not, let's make a new hex code for it!
                while($genus_hex_len--){$genus_hex .= $genus_hex_chars[rand @genus_hex_chars]};
                
                #And add a hash at the start...
                $genus_hex = "#".$genus_hex;
                
                #And add the genus and hex as a hash-key pair
                $hex_hash{$genus} = $genus_hex;
            }
            
            #Check if the viral sub-family has been seen before
            if(exists($hex_hash{$sub_family})) {
            
                #Then look up the relevant hexcode as a hash-key lookup
                $sub_family_hex = $hex_hash{"$sub_family"};
            
            } else {
            
                #And if not, let's make a new hex code for it!
                while($sub_hex_len--){$sub_family_hex .= $sub_hex_chars[rand @sub_hex_chars]};
                
                #And add a hash at the start...
                $sub_family_hex = "#".$sub_family_hex;
                
                #And add the sub-family and hex as a hash-key pair
                $hex_hash{$sub_family} = $sub_family_hex;
            }
            
            #Check if the viral family has been seen before
            if(exists($hex_hash{$family})) {
            
                #Then look up the relevant hexcode as a hash-key lookup
                $family_hex = $hex_hash{"$family"};
            
            } else {
            
                #And if not, let's make a new hex code for it!
                while($fam_hex_len--){$family_hex .= $fam_hex_chars[rand @fam_hex_chars]};
                
                #And add a hash at the start...
                $family_hex = "#".$family_hex;
                
                #And add the family and hex as a hash-key pair
                $hex_hash{$family} = $family_hex;
            }
            
            #Check if the host has been seen before
            if(exists($hex_hash{$host})) {
            
                #Then look up the relevant hexcode as a hash-key lookup
                $host_hex = $hex_hash{"$host"};
            
            } else {
            
                #And if not, let's make a new hex code for it!
                while($host_hex_len--){$host_hex .= $host_hex_chars[rand @host_hex_chars]};
                
                #And add a hash at the start...
                $host_hex = "#".$host_hex;
                
                #And add the host and hex as a hash-key pair
                $hex_hash{$host} = $host_hex;
            }
            
            #Store the genus, sub-family, family, lowest taxa, host data, length, and description in the hashes we made earlier (these are useful for the annotation files)
            $genus_hash{$primary_id} = $genus;
            $sub_family_hash{$primary_id} = $sub_family;
            $family_hash{$primary_id} = $family;
            $lowest_taxa_hash{$primary_id} = $lowest_taxa;
            $host_hash{$primary_id} = $host;
            $description_hash{$primary_id} = $species_string;
            $length_hash{$primary_id} = $length;
            
            #Add the accession, genus, sub-family, family, lowest taxa, and hosts to their respective arrays
            push @accession_array, $primary_id;
            push @genus_array, $genus;
            push @sub_family_array, $sub_family;
            push @family_array, $family;
            push @lowest_taxa_array, $lowest_taxa;
            push @host_array, $host;
            
            #Take those that are not RefSeq and write them to their own file, their own .tsv file, and produce useful vConTACT2 files for them (including annotations). The non-RefSeq sequences are separated here, as to avoid duplicates in the useful files
            if ($primary_id !~ m/NC_*/ ) {
            
                #Write sequences to fasta file that excludes refseq sequences
                $seq_out_excluding_refseq->write_seq($inseq_obj);
                
                #Write these to the tsv file that does not include RefSeq sequences
                print TAXA_without_refseq "$primary_id\t$species_string\t@classification\t$length\t$jumbo\t$pergc\t$molecule\t@dates\t$total_CDS\t$per_plus\t$per_minus\t$coding_capacity\t$coding_warning\t$trna_count\t$host\t$lowest_taxa\t$genus\t$sub_family\t$family\t$order\t$class\t$phylum\t$kingdom\t$realm\t$baltimore\t$division\t$isolation_host\n";
                
                #Grab whatever hex code was used for the lowest taxa
                my $lowest_taxa_hex = $hex_hash{"$lowest_taxa"};
                
                #Print relevant information to the vConTACT2 annotation files
                print VGENUS "$primary_id\t$species_string\t$genus\t$genus_hex\n";
                print VSUB "$primary_id\t$species_string\t$sub_family\t$sub_family_hex\n";
                print VFAMILY "$primary_id\t$species_string\t$family\t$family_hex\n";
                print VLOWEST "$primary_id\t$species_string\t$lowest_taxa\t$lowest_taxa_hex\n";
                print VHOST "$primary_id\t$species_string\t$host\t$host_hex\n";
                
                #Give the path to the amino acid fasta file. This is to write vConTACT2 input files
                my $faa = "GenomesDB/".$acc."/".$acc.".faa";
            
                #Read it as a Bio::SeqIO object
                my $faa_object = Bio::SeqIO->new(-file => "$faa");         
                
                #Read them one at a time
                while (my $faa_seq_object = $faa_object->next_seq) {
                
                    #Get the primary ID
                    my $protein_id = $faa_seq_object->display_id;
                    
                    #Print this to the mapping file
                    print MAP "$protein_id,$primary_id,none\n";
                    
                    #Write the amino acid sequence to faa file
                    $protein->write_seq($faa_seq_object);                   
                }
            }
        }
    }
    }
}
