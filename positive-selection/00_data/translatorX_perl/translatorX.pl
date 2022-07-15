#!/usr/bin/perl -w
#Usage: perl translatorx_vLocal.pl -i borrame.ntseqs.fasta -o trx.out -p M -t F -w 1 -c 5

#Citation: Abascal F, Zardoya R, Telford MJ (2010)
#TranslatorX server: multiple alignment of nucleotide sequences guided by amino acid information 
#Nucleic Acids Research

use strict;
use Cwd;

my $old_input_rec_sep = $/;

sub print_in_H_and_exit {
	my $fh  = shift;
	my $msg = shift;
	print $fh "<font color=\"red\">$msg</font>\n";
	print STDERR $msg;
	exit;
}

sub tryToSolveCodonUncertainty {
	my $codon   = shift;
	my $p_code  = shift;
	my $fh      = shift;
	my %correspondences;
	$correspondences{"A"} = ["A"];
	$correspondences{"T"} = ["T"];
	$correspondences{"C"} = ["C"];
	$correspondences{"G"} = ["G"];
	$correspondences{"U"} = ["T"];
	$correspondences{"M"} = ["A","C"];
	$correspondences{"R"} = ["A","G"];
	$correspondences{"W"} = ["A","T"];
	$correspondences{"S"} = ["C","G"];
	$correspondences{"Y"} = ["C","T"];
	$correspondences{"K"} = ["G","T"];
	$correspondences{"B"} = ["C","G","T"];
	$correspondences{"V"} = ["A","C","G"];
	$correspondences{"H"} = ["A","C","T"];
	$correspondences{"D"} = ["A","G","T"];
	$correspondences{"N"} = ["A","C","G","T"];
	$correspondences{"X"} = ["A","C","G","T"];
	my($f,$s,$t) = (split(//,$codon))[0,1,2];
	foreach my $l ( $f, $s, $t) {
		if(!exists($correspondences{$l})) {
			print $fh    "<font color=\"red\">Error: codon $codon is not in IUPAC code ($l?). Substituting $l by'N'</font><br>\n";
			print STDERR "Error: codon $codon is not in IUPAC code ($l?). Substituting $l by'N'\n";
		}
	}
	if($f !~ /[ATCGUMRWSYKBVHDNX]/)	{	$f = "N";	}
	if($s !~ /[ATCGUMRWSYKBVHDNX]/)	{	$s = "N";	}
	if($t !~ /[ATCGUMRWSYKBVHDNX]/)	{	$t = "N";	}
	my(@codons,@codons2,@codons3);
	foreach my $c ( @{$correspondences{$f}} ) {
		push(@codons,$c);
		foreach my $alt_codons ( @codons ) {
			foreach my $c2 ( @{$correspondences{$s}} ) {
				push(@codons2,$alt_codons.$c2);
				foreach my $alt_codons2 ( @codons2 ) {
					foreach my $c3 ( @{$correspondences{$t}} ) {
						push(@codons3,$alt_codons.$c2.$c3);
					}
				}
			}
		}
	}
	my %aas;
	foreach my $c ( @codons3 ) {
		$aas{$p_code->{$c}}++;
	}
	if(scalar(keys(%aas)) > 1) {
		print $fh    "<font color=\"red\">Warning: Unable to disambiguate codon $codon. It can be translated as ",scalar(keys(%aas))," different amino acids [",keys(%aas),"]</font><br>\n";
		print STDERR "Warning: Unable to disambiguate codon $codon. It can be translated as ",scalar(keys(%aas))," different amino acids [",keys(%aas),"]\n";
		return "X";
	} else {
		return $p_code->{$codons3[0]};
	}
}

sub usage {
	my $error = shift;
	print STDERR " TranslatorX v1.1 (march 2010)\n\n";
	print STDERR " Error: $error\n\n";
	print STDERR "             Abascal F, Zardoya R, Telford MJ (2010)\n";
	print STDERR "   TranslatorX server: multiple alignment of nucleotide sequences\n";
	print STDERR "                guided by amino acid information.\n";
	print STDERR "               Nucleic Acids Research  (in press).\n\n\n";
	print STDERR " Basic usage:    perl translatorx.pl -i nt_file\n";
	print STDERR " Advanced usage: perl translatorx.pl -i nt_file -a aa_file -j gc_file\n";
	print STDERR "  -i: the file containing the nucleotide sequences in FASTA format (Required) \n";
	print STDERR "  -o: output file (Optional). Default: \"translatorx_res\".\n";
	print STDERR "  -a: file containing the amino acid sequence alignment (Optional)\n";
	print STDERR "  -p: program to build the multiple alignment (Optional). Available options are:\n";
	print STDERR "      M/C/F/P, standing for Muscle, Clustalw, Prank, and maFft\n";
	print STDERR "      Default: Muscle\n";
	print STDERR "  -c: genetic code to translate the sequences (Optional). Available options are:\n";
	print STDERR "      1   Standard\n";
	print STDERR "      2   Vertebrate Mitochondrial\n";
	print STDERR "      3   Yeast Mitochondrial\n";
	print STDERR "      4   Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma\n";
	print STDERR "      5   Invertebrate Mitochondrial\n";
	print STDERR "      6   Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear\n";
	print STDERR "      9   Echinoderm Mitochondrial; Flatworm Mitochondrial\n";
	print STDERR "      10  Euplotid Nuclear\n";
	print STDERR "      11  Bacterial and Plant Plastid\n";
	print STDERR "      12  Alternative Yeast Nuclear\n";
	print STDERR "      13  Ascidian Mitochondrial\n";
	print STDERR "      14  Alternative Flatworm Mitochondrial\n";
	print STDERR "      15  Blepharisma Macronuclear\n";
	print STDERR "      16  Chlorophycean Mitochondrial\n";
	print STDERR "      21  Trematode Mitochondrial\n";
	print STDERR "      22  Scenedesmus obliquus Mitochondrial\n";
	print STDERR "      23  Thraustochytrium Mitochondrial\n";
	print STDERR "      100 Ancestral Arthropod Mitochondrial Code (AGG=K)\n";
	print STDERR "      101 Hemichordate Mitochondrial\n";
	print STDERR "      Default: 1 (Standard code)\n";
	print STDERR "  -j: file containing alternative genetic codes for each taxon. (Optional)\n";
	print STDERR "      File format: Taxon\tgenetic_code[number]\n";
	print STDERR "  -g: parameters for GBlocks (Optional)\n";
	print STDERR "      Example: -g \"-b2 x -b3 x -b4 x...\"\n"; #POR HACER
	print STDERR "  -t: F/T. Guess reading frame (Optional)[default F]\n";
	print STDERR "\n";
	print STDERR "\n";
	exit;
}


my @aas = ("D", "E", "C", "M", "K", "R", "S", "T", "F", "Y", "N", "Q", "G", "L", "V", "I", "A", "W", "H", "P");
my %aa2colors;
$aa2colors{'D'} = "#A60A0A";		$aa2colors{'E'} = "#F61A2A";		$aa2colors{'C'} = "#E6E600";
$aa2colors{'M'} = "#E6E162";		$aa2colors{'K'} = "#444AFF";		$aa2colors{'R'} = "#296AFF";
$aa2colors{'S'} = "#FA9600";		$aa2colors{'T'} = "#BA8600";		$aa2colors{'F'} = "#8282AA";
$aa2colors{'Y'} = "#5252AA";		$aa2colors{'N'} = "#00AAAA";		$aa2colors{'Q'} = "#00DFDF";
$aa2colors{'G'} = "#EBEBEB";		$aa2colors{'L'} = "#1F921F";		$aa2colors{'V'} = "#7FA27F";
$aa2colors{'I'} = "#5F825F";		$aa2colors{'A'} = "#C8C8C8";		$aa2colors{'W'} = "#B45AB4";
$aa2colors{'H'} = "#8282D2";		$aa2colors{'P'} = "#DC9682";		$aa2colors{'X'} = "#FFFFFF";
$aa2colors{'*'} = "#FFFFFF";		$aa2colors{'-'} = "#FFFFFF";		$aa2colors{'.'} = "#FFFFFF";
$aa2colors{'?'} = "#FFFFFF";


my %params = @ARGV;
my %valid_params;
$valid_params{"-i"} = 1;	$valid_params{"-o"} = 1;	$valid_params{"-a"} = 1;	
$valid_params{"-p"} = 1;	$valid_params{"-c"} = 1;	$valid_params{"-j"} = 1;
$valid_params{"-w"} = 1;	$valid_params{"-g"} = 1;	$valid_params{"-t"} = 1;

foreach my $param ( keys %params ) {
	if(!exists($valid_params{$param})) {
		&usage("Option -$param not recognized");
	}
	if(!defined($params{$param})) {
		&usage("Option -$param requires an accompanying value");
	}
}


my $from_web_server = 0;
if(exists($params{"-w"})) {
	$from_web_server = 1;
} 


if(!exists($params{"-i"})) {	&usage("Option -i not defined") ;	}
my $nt_file = $params{"-i"};
if(!-e $nt_file) 
{	&usage("File $nt_file does not exist");		}


my $program = "muscle";
my $aa_file = "";
if(!exists($params{"-a"})) {	
	if(!exists($params{"-p"})) 			{	$program = "muscle";		}
	elsif ($params{"-p"} =~ /^[Cc]/)	{	$program = "clustalw";}
	elsif ($params{"-p"} =~ /^[Tt]/)	{	$program = "t_coffee";}
	elsif ($params{"-p"} =~ /^[Pp]/)	{	$program = "prank";	}
	elsif ($params{"-p"} =~ /^[Ff]/)	{	$program = "mafft";	}
	else 								{	$program = "muscle";		}
} else 									{	$program = "";			
										$aa_file = $params{"-a"};	
										if(!-e $aa_file) 
										{	&usage("File $aa_file does not exist");		}
										}

my $gblocks_params = "";
if(exists($params{"-g"})) {
	$gblocks_params = $params{"-g"};
}


my $geneticcode = 1;
my $geneticcode_file = "";
if(!exists($params{"-j"})) {
	if(!exists($params{"-c"})) 	{	$geneticcode = 1;		}
	else						{	$geneticcode = $params{"-c"};		}
} else							{	$geneticcode = 0;						#0 means "variable"	
									$geneticcode_file = $params{"-j"};	}

my $guessRF = 0; #POR HACER: ponerlo a 0 por defecto.
if(exists($params{"-t"})) {
	if($params{"-t"} =~ /[Tt]/) {
		$guessRF = 1;
	}
}

my $out_file = "";
if(exists($params{"-o"}))		{	$out_file = $params{"-o"};	}
else							{	$out_file = "translatorx_res";		}


my $html_file = $out_file.".html";
open(H, ">$html_file") || print_in_H_and_exit(\*H, "Unable to write to the html file: $html_file\n");
print H <<HTML_HEADER;
<html>
<head>
	<style type='text/css'>


div#rounded-box2 {margin:10 auto;padding:0;border-radius:7px;-moz-border-radius: 7px; -webkit-border-radius : 7px;}  
div#rounded-box2 p {margin:0;padding:10px;}  	  
	  
	.elementoVisible {display:block;}
	.elementoOculto {display:none;}
	.linkContraido {
		cursor: pointer;
		background: url(EC_contraido.png) no-repeat;
		color: #444444;
		font-size: 10px;
		font-weight: bold;
		font-family: Verdana,Georgia,Helvetica;
	}
	.linkExpandido {
		cursor: pointer;
		background: url(EC_expandido.png) no-repeat;
		color: #444444;
		font-size: 10px;
		font-weight: bold;
		font-family: Verdana,Georgia,Helvetica;
	}
	TD, TR, BODY
	{
		font-size: 12px;
		font-family: Verdana, Arial, Helvetica, Helv, sans-serif;
	}
	textarea
	{
  	border:1px solid #999999;
	  width:98%;
  	margin:5px 0;
	  padding:1%;
	}
	.rjt
	{
 		color : #3333BB;
 		font-size : 12px; 
		font-weight: normal;
		text-decoration : underline;
		font-family: Verdana, Arial, Helvetica, Helv, sans-serif;
	}	
	

	
	</style>
	
	<script type='text/JavaScript'>
	function desplegarContraer(cual,desde) {
		var elElemento=document.getElementById(cual);
		if(elElemento.className == 'elementoVisible') {
		elElemento.className = 'elementoOculto';
		desde.className = 'linkContraido';
		} else {
		elElemento.className = 'elementoVisible';
		desde.className = 'linkExpandido';
		}
	}
	function expandir(cual,desde) {
		var elElemento=document.getElementById(cual);
		elElemento.className = 'elementoVisible';
	}
	function contraer(cual,desde) {
		var elElemento=document.getElementById(cual);
		elElemento.className = 'elementoOculto';
	}
	function expandirlink(cual,desde) {
		var elElemento=document.getElementById(cual);
		elElemento.className = 'linkExpandido';
	}
	function contraerlink(cual,desde) {
		var elElemento=document.getElementById(cual);
		elElemento.className = 'linkContraido';
	}
	function alterlinks(cual,elotro,desde) {
		var elElemento=document.getElementById(cual);
		var elOtroElemento =document.getElementById(elotro);
		if(elElemento.className == 'linkContraido') {
		elElemento.className = 'linkExpandido';
		elOtroElemento.className = 'linkContraido';
		} else {
		elElemento.className = 'linkContraido';
		elOtroElemento.className = 'linkExpandido';
		}
	}
	function send2file (nombre_elemento) {
		var ta = document.getElementById(nombre_elemento);
		var taval = ta.value;
		var doc = document.open("text/save");
		doc.writeln(taval);
		doc.close(); 
	}
	function send2file_webserver (file) {
		window.location = "/cgi-bin/sequences.fasta?file="+file;
	}
	function goBackToTranslatorx () {
		window.location.replace("/index_v5.html");
	}
	</script>
	<title>TranslatorX results</title>
</head>

<body>

<table cellspacing=0 style="border : thin solid; border-color:#CCCCCC" width="100%" bgcolor="#CCCCCC"><tr><td>
<table width="100%" border=0 cellspacing=5 bgcolor="#EEEEEE"><tr><td>




HTML_HEADER


#POR HACER: poner T R A N S L A T O R - X con letras de colores aleatorios. Seguro que queda chuli

my $range = 19;
my %yausados;
my @fname = ("T","R","A","N","S","L","A","T","O","R","-","X");
my $fcounter = 0;
#print H "<input type=\"button\" value=\"Back\" onclick=\"history.go(-2)\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<input type=\"button\" value=\"Help\" onclick=\"window.open('/translatorx_help.html','TranslatorX help','width=640,height=600,scrollbars=1')\"><br>\n";
print H "<input type=\"button\" value=\"Back\" onclick=\"goBackToTranslatorx();\">&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<input type=\"button\" value=\"Help\" onclick=\"window.open('/translatorx_help.html','TranslatorX help','width=640,height=600,scrollbars=1')\"><br>\n";
print H "<table border=0 cellspacing=10><tr valign=\"center\"><td><table cellpadding=5><tr>\n";
foreach my $letter ( @fname ) { 
	my $random_number = int(rand($range)+1);
	while(exists($yausados{$random_number})) {	$random_number = int(rand($range)+1);	}
	$yausados{$random_number} = 1;
	print H "<td bgcolor=\"$aa2colors{$aas[$random_number]}\"><b>$letter</b></td>\n";
	$fcounter++;
	if($fcounter == 4) {
		print H "</tr><tr>\n"; 
		$fcounter = 0;
	}
}
print H "</tr></table><br></td><td>\n";
print H "<b style=\"font-size: 20px;\">TranslatorX results</b><br>\n";
print H "<font color=\"#804444\" style=\"font-size:10px;\">Powered by ReadSeq, GBlocks, Jalview, Muscle, Clustalw, MAFFT, T-Coffee and PRank</font>\n";
print H "</td></tr></table>\n";
if($from_web_server) {
	my $tmp_file = $nt_file;
	$tmp_file = (split(/\//,$tmp_file))[-1];
	print H "<a href=\"$tmp_file\">Check if Readseq successfully converted your original alignment to fasta format</a>.<br>\n";
	if($program && $program ne "") {
		$tmp_file = $out_file;
		$tmp_file = (split(/\//,$tmp_file))[-1];
		print H "<a href=\"$tmp_file.aaseqs.fasta\">Amino acid translation of the nucleotide sequences</a>.<br>\n";
		print H "<a href=\"$tmp_file.$program.log\">Check <b style=\"color:green;\">$program</b>"."'s log</a>.<br>\n";
	}
}


my(	%gccontent_sp_1st, 	%gccontent_sp_2nd, 	%gccontent_sp_3rd,
	%gccontent_1st, 	%gccontent_2nd, 	%gccontent_3rd,
	%gccontent_sp, 		%gccontent					);
my(	%clean_gccontent_sp_1st, 	%clean_gccontent_sp_2nd, 	%clean_gccontent_sp_3rd,
	%clean_gccontent_1st, 		%clean_gccontent_2nd, 		%clean_gccontent_3rd,
	%clean_gccontent_sp, 		%clean_gccontent			);


# Variables.
my @nt_seq;
my $nt_seq;
my $nt_name;
my %nt_nameseq;
my @nt_spnames;
my $code;
my $codehash;
my %startingCodons;
my %namecode;
my $namecode_name;
my $namecode_code;
my @tempseqarray;
my $temptriplet;
my @tempaminoacidseq;
my $pauser;
my $clustaloutfile;
my $aa_name;
my @aa_seq;
my $aa_seq;
my %aa_nameseq;
my @aa_spnames;
my $shifted_nts;
my $fastaoutfile;
my $file1; 
my $clustalinfile;
my $pirfirstline;
my $pir_seq;
my @pir_seq;
my $pir_name;
my %pir_nameseq;
my @pir_spnames;
my $muscleinfile;
my $muscleoutfile;
my $samecodeYN; 
my $spcode; 
my $namecodefile; 
my $xcounter; 
my $discardline;
my $aa_posn; 
my @aa_array; 
my @nt_array; 
my $aa_length;
my @finalnt_array; 
my $finalaoutfile;
my $firstline;
my $secondline;
my $t_coffeeoutname;
my $filetype;
my $function;
my $tempnt =();
my $tempaa =();
my $tempntlength = ();
my $tempaalength = ();

print STDERR "Path: $0\n";
my $path = $0;
$path =~ s/\/.[^\/]*?$/\//;
$path = $path. "/";
print STDERR "Path: $path\n";
$path ="";


my(%code_names, %starting_codes, %codes);

#Code: 1 Standard
$code_names{1} = "Standard";
$starting_codes{1}->{"ATG"} = 1;
$starting_codes{1}->{"CTG"} = 1;
$starting_codes{1}->{"TTG"} = 1;

$codes{1}->{"AAA"} = "K";	$codes{1}->{"AAC"} = "N";	$codes{1}->{"AAG"} = "K";	$codes{1}->{"AAT"} = "N";	
$codes{1}->{"ACA"} = "T";	$codes{1}->{"ACC"} = "T";	$codes{1}->{"ACG"} = "T";	$codes{1}->{"ACT"} = "T";	
$codes{1}->{"AGA"} = "R";	$codes{1}->{"AGC"} = "S";	$codes{1}->{"AGG"} = "R";	$codes{1}->{"AGT"} = "S";	
$codes{1}->{"ATA"} = "I";	$codes{1}->{"ATC"} = "I";	$codes{1}->{"ATG"} = "M";	$codes{1}->{"ATT"} = "I";	
$codes{1}->{"CAA"} = "Q";	$codes{1}->{"CAC"} = "H";	$codes{1}->{"CAG"} = "Q";	$codes{1}->{"CAT"} = "H";	
$codes{1}->{"CCA"} = "P";	$codes{1}->{"CCC"} = "P";	$codes{1}->{"CCG"} = "P";	$codes{1}->{"CCT"} = "P";	
$codes{1}->{"CGA"} = "R";	$codes{1}->{"CGC"} = "R";	$codes{1}->{"CGG"} = "R";	$codes{1}->{"CGT"} = "R";	
$codes{1}->{"CTA"} = "L";	$codes{1}->{"CTC"} = "L";	$codes{1}->{"CTG"} = "L";	$codes{1}->{"CTT"} = "L";	
$codes{1}->{"GAA"} = "E";	$codes{1}->{"GAC"} = "D";	$codes{1}->{"GAG"} = "E";	$codes{1}->{"GAT"} = "D";	
$codes{1}->{"GCA"} = "A";	$codes{1}->{"GCC"} = "A";	$codes{1}->{"GCG"} = "A";	$codes{1}->{"GCT"} = "A";	
$codes{1}->{"GGA"} = "G";	$codes{1}->{"GGC"} = "G";	$codes{1}->{"GGG"} = "G";	$codes{1}->{"GGT"} = "G";	
$codes{1}->{"GTA"} = "V";	$codes{1}->{"GTC"} = "V";	$codes{1}->{"GTG"} = "V";	$codes{1}->{"GTT"} = "V";	
$codes{1}->{"TAA"} = "X";	$codes{1}->{"TAC"} = "Y";	$codes{1}->{"TAG"} = "X";	$codes{1}->{"TAT"} = "Y";	
$codes{1}->{"TCA"} = "S";	$codes{1}->{"TCC"} = "S";	$codes{1}->{"TCG"} = "S";	$codes{1}->{"TCT"} = "S";	
$codes{1}->{"TGA"} = "X";	$codes{1}->{"TGC"} = "C";	$codes{1}->{"TGG"} = "W";	$codes{1}->{"TGT"} = "C";	
$codes{1}->{"TTA"} = "L";	$codes{1}->{"TTC"} = "F";	$codes{1}->{"TTG"} = "L";	$codes{1}->{"TTT"} = "F";	

#Code: 2 Vertebrate Mitochondrial
$code_names{2} = "Vertebrate Mitochondrial";
$starting_codes{2}->{"ATA"} = 1;
$starting_codes{2}->{"ATC"} = 1;
$starting_codes{2}->{"ATG"} = 1;
$starting_codes{2}->{"ATT"} = 1;
$starting_codes{2}->{"GTG"} = 1;

$codes{2}->{"AAA"} = "K";	$codes{2}->{"AAC"} = "N";	$codes{2}->{"AAG"} = "K";	$codes{2}->{"AAT"} = "N";	
$codes{2}->{"ACA"} = "T";	$codes{2}->{"ACC"} = "T";	$codes{2}->{"ACG"} = "T";	$codes{2}->{"ACT"} = "T";	
$codes{2}->{"AGA"} = "X";	$codes{2}->{"AGC"} = "S";	$codes{2}->{"AGG"} = "X";	$codes{2}->{"AGT"} = "S";	
$codes{2}->{"ATA"} = "M";	$codes{2}->{"ATC"} = "I";	$codes{2}->{"ATG"} = "M";	$codes{2}->{"ATT"} = "I";	
$codes{2}->{"CAA"} = "Q";	$codes{2}->{"CAC"} = "H";	$codes{2}->{"CAG"} = "Q";	$codes{2}->{"CAT"} = "H";	
$codes{2}->{"CCA"} = "P";	$codes{2}->{"CCC"} = "P";	$codes{2}->{"CCG"} = "P";	$codes{2}->{"CCT"} = "P";	
$codes{2}->{"CGA"} = "R";	$codes{2}->{"CGC"} = "R";	$codes{2}->{"CGG"} = "R";	$codes{2}->{"CGT"} = "R";	
$codes{2}->{"CTA"} = "L";	$codes{2}->{"CTC"} = "L";	$codes{2}->{"CTG"} = "L";	$codes{2}->{"CTT"} = "L";	
$codes{2}->{"GAA"} = "E";	$codes{2}->{"GAC"} = "D";	$codes{2}->{"GAG"} = "E";	$codes{2}->{"GAT"} = "D";	
$codes{2}->{"GCA"} = "A";	$codes{2}->{"GCC"} = "A";	$codes{2}->{"GCG"} = "A";	$codes{2}->{"GCT"} = "A";	
$codes{2}->{"GGA"} = "G";	$codes{2}->{"GGC"} = "G";	$codes{2}->{"GGG"} = "G";	$codes{2}->{"GGT"} = "G";	
$codes{2}->{"GTA"} = "V";	$codes{2}->{"GTC"} = "V";	$codes{2}->{"GTG"} = "V";	$codes{2}->{"GTT"} = "V";	
$codes{2}->{"TAA"} = "X";	$codes{2}->{"TAC"} = "Y";	$codes{2}->{"TAG"} = "X";	$codes{2}->{"TAT"} = "Y";	
$codes{2}->{"TCA"} = "S";	$codes{2}->{"TCC"} = "S";	$codes{2}->{"TCG"} = "S";	$codes{2}->{"TCT"} = "S";	
$codes{2}->{"TGA"} = "W";	$codes{2}->{"TGC"} = "C";	$codes{2}->{"TGG"} = "W";	$codes{2}->{"TGT"} = "C";	
$codes{2}->{"TTA"} = "L";	$codes{2}->{"TTC"} = "F";	$codes{2}->{"TTG"} = "L";	$codes{2}->{"TTT"} = "F";	

#Code: 3 Yeast Mitochondrial
$code_names{3} = "Yeast Mitochondrial";
$starting_codes{3}->{"ATA"} = 1;
$starting_codes{3}->{"ATG"} = 1;

$codes{3}->{"AAA"} = "K";	$codes{3}->{"AAC"} = "N";	$codes{3}->{"AAG"} = "K";	$codes{3}->{"AAT"} = "N";	
$codes{3}->{"ACA"} = "T";	$codes{3}->{"ACC"} = "T";	$codes{3}->{"ACG"} = "T";	$codes{3}->{"ACT"} = "T";	
$codes{3}->{"AGA"} = "R";	$codes{3}->{"AGC"} = "S";	$codes{3}->{"AGG"} = "R";	$codes{3}->{"AGT"} = "S";	
$codes{3}->{"ATA"} = "M";	$codes{3}->{"ATC"} = "I";	$codes{3}->{"ATG"} = "M";	$codes{3}->{"ATT"} = "I";	
$codes{3}->{"CAA"} = "Q";	$codes{3}->{"CAC"} = "H";	$codes{3}->{"CAG"} = "Q";	$codes{3}->{"CAT"} = "H";	
$codes{3}->{"CCA"} = "P";	$codes{3}->{"CCC"} = "P";	$codes{3}->{"CCG"} = "P";	$codes{3}->{"CCT"} = "P";	
$codes{3}->{"CGA"} = "R";	$codes{3}->{"CGC"} = "R";	$codes{3}->{"CGG"} = "R";	$codes{3}->{"CGT"} = "R";	
$codes{3}->{"CTA"} = "T";	$codes{3}->{"CTC"} = "T";	$codes{3}->{"CTG"} = "T";	$codes{3}->{"CTT"} = "T";	
$codes{3}->{"GAA"} = "E";	$codes{3}->{"GAC"} = "D";	$codes{3}->{"GAG"} = "E";	$codes{3}->{"GAT"} = "D";	
$codes{3}->{"GCA"} = "A";	$codes{3}->{"GCC"} = "A";	$codes{3}->{"GCG"} = "A";	$codes{3}->{"GCT"} = "A";	
$codes{3}->{"GGA"} = "G";	$codes{3}->{"GGC"} = "G";	$codes{3}->{"GGG"} = "G";	$codes{3}->{"GGT"} = "G";	
$codes{3}->{"GTA"} = "V";	$codes{3}->{"GTC"} = "V";	$codes{3}->{"GTG"} = "V";	$codes{3}->{"GTT"} = "V";	
$codes{3}->{"TAA"} = "X";	$codes{3}->{"TAC"} = "Y";	$codes{3}->{"TAG"} = "X";	$codes{3}->{"TAT"} = "Y";	
$codes{3}->{"TCA"} = "S";	$codes{3}->{"TCC"} = "S";	$codes{3}->{"TCG"} = "S";	$codes{3}->{"TCT"} = "S";	
$codes{3}->{"TGA"} = "W";	$codes{3}->{"TGC"} = "C";	$codes{3}->{"TGG"} = "W";	$codes{3}->{"TGT"} = "C";	
$codes{3}->{"TTA"} = "L";	$codes{3}->{"TTC"} = "F";	$codes{3}->{"TTG"} = "L";	$codes{3}->{"TTT"} = "F";	

#Code: 4 Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma
$code_names{4} = "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma";
$starting_codes{4}->{"ATA"} = 1;
$starting_codes{4}->{"ATC"} = 1;
$starting_codes{4}->{"ATG"} = 1;
$starting_codes{4}->{"ATT"} = 1;
$starting_codes{4}->{"CTG"} = 1;
$starting_codes{4}->{"GTG"} = 1;
$starting_codes{4}->{"TTA"} = 1;
$starting_codes{4}->{"TTG"} = 1;

$codes{4}->{"AAA"} = "K";	$codes{4}->{"AAC"} = "N";	$codes{4}->{"AAG"} = "K";	$codes{4}->{"AAT"} = "N";	
$codes{4}->{"ACA"} = "T";	$codes{4}->{"ACC"} = "T";	$codes{4}->{"ACG"} = "T";	$codes{4}->{"ACT"} = "T";	
$codes{4}->{"AGA"} = "R";	$codes{4}->{"AGC"} = "S";	$codes{4}->{"AGG"} = "R";	$codes{4}->{"AGT"} = "S";	
$codes{4}->{"ATA"} = "I";	$codes{4}->{"ATC"} = "I";	$codes{4}->{"ATG"} = "M";	$codes{4}->{"ATT"} = "I";	
$codes{4}->{"CAA"} = "Q";	$codes{4}->{"CAC"} = "H";	$codes{4}->{"CAG"} = "Q";	$codes{4}->{"CAT"} = "H";	
$codes{4}->{"CCA"} = "P";	$codes{4}->{"CCC"} = "P";	$codes{4}->{"CCG"} = "P";	$codes{4}->{"CCT"} = "P";	
$codes{4}->{"CGA"} = "R";	$codes{4}->{"CGC"} = "R";	$codes{4}->{"CGG"} = "R";	$codes{4}->{"CGT"} = "R";	
$codes{4}->{"CTA"} = "L";	$codes{4}->{"CTC"} = "L";	$codes{4}->{"CTG"} = "L";	$codes{4}->{"CTT"} = "L";	
$codes{4}->{"GAA"} = "E";	$codes{4}->{"GAC"} = "D";	$codes{4}->{"GAG"} = "E";	$codes{4}->{"GAT"} = "D";	
$codes{4}->{"GCA"} = "A";	$codes{4}->{"GCC"} = "A";	$codes{4}->{"GCG"} = "A";	$codes{4}->{"GCT"} = "A";	
$codes{4}->{"GGA"} = "G";	$codes{4}->{"GGC"} = "G";	$codes{4}->{"GGG"} = "G";	$codes{4}->{"GGT"} = "G";	
$codes{4}->{"GTA"} = "V";	$codes{4}->{"GTC"} = "V";	$codes{4}->{"GTG"} = "V";	$codes{4}->{"GTT"} = "V";	
$codes{4}->{"TAA"} = "X";	$codes{4}->{"TAC"} = "Y";	$codes{4}->{"TAG"} = "X";	$codes{4}->{"TAT"} = "Y";	
$codes{4}->{"TCA"} = "S";	$codes{4}->{"TCC"} = "S";	$codes{4}->{"TCG"} = "S";	$codes{4}->{"TCT"} = "S";	
$codes{4}->{"TGA"} = "W";	$codes{4}->{"TGC"} = "C";	$codes{4}->{"TGG"} = "W";	$codes{4}->{"TGT"} = "C";	
$codes{4}->{"TTA"} = "L";	$codes{4}->{"TTC"} = "F";	$codes{4}->{"TTG"} = "L";	$codes{4}->{"TTT"} = "F";	

#Code: 5 Invertebrate Mitochondrial
$code_names{5} = "Invertebrate Mitochondrial";
$starting_codes{5}->{"ATA"} = 1;
$starting_codes{5}->{"ATC"} = 1;
$starting_codes{5}->{"ATG"} = 1;
$starting_codes{5}->{"ATT"} = 1;
$starting_codes{5}->{"GTG"} = 1;
$starting_codes{5}->{"TTG"} = 1;

$codes{5}->{"AAA"} = "K";	$codes{5}->{"AAC"} = "N";	$codes{5}->{"AAG"} = "K";	$codes{5}->{"AAT"} = "N";	
$codes{5}->{"ACA"} = "T";	$codes{5}->{"ACC"} = "T";	$codes{5}->{"ACG"} = "T";	$codes{5}->{"ACT"} = "T";	
$codes{5}->{"AGA"} = "S";	$codes{5}->{"AGC"} = "S";	$codes{5}->{"AGG"} = "S";	$codes{5}->{"AGT"} = "S";	
$codes{5}->{"ATA"} = "M";	$codes{5}->{"ATC"} = "I";	$codes{5}->{"ATG"} = "M";	$codes{5}->{"ATT"} = "I";	
$codes{5}->{"CAA"} = "Q";	$codes{5}->{"CAC"} = "H";	$codes{5}->{"CAG"} = "Q";	$codes{5}->{"CAT"} = "H";	
$codes{5}->{"CCA"} = "P";	$codes{5}->{"CCC"} = "P";	$codes{5}->{"CCG"} = "P";	$codes{5}->{"CCT"} = "P";	
$codes{5}->{"CGA"} = "R";	$codes{5}->{"CGC"} = "R";	$codes{5}->{"CGG"} = "R";	$codes{5}->{"CGT"} = "R";	
$codes{5}->{"CTA"} = "L";	$codes{5}->{"CTC"} = "L";	$codes{5}->{"CTG"} = "L";	$codes{5}->{"CTT"} = "L";	
$codes{5}->{"GAA"} = "E";	$codes{5}->{"GAC"} = "D";	$codes{5}->{"GAG"} = "E";	$codes{5}->{"GAT"} = "D";	
$codes{5}->{"GCA"} = "A";	$codes{5}->{"GCC"} = "A";	$codes{5}->{"GCG"} = "A";	$codes{5}->{"GCT"} = "A";	
$codes{5}->{"GGA"} = "G";	$codes{5}->{"GGC"} = "G";	$codes{5}->{"GGG"} = "G";	$codes{5}->{"GGT"} = "G";	
$codes{5}->{"GTA"} = "V";	$codes{5}->{"GTC"} = "V";	$codes{5}->{"GTG"} = "V";	$codes{5}->{"GTT"} = "V";	
$codes{5}->{"TAA"} = "X";	$codes{5}->{"TAC"} = "Y";	$codes{5}->{"TAG"} = "X";	$codes{5}->{"TAT"} = "Y";	
$codes{5}->{"TCA"} = "S";	$codes{5}->{"TCC"} = "S";	$codes{5}->{"TCG"} = "S";	$codes{5}->{"TCT"} = "S";	
$codes{5}->{"TGA"} = "W";	$codes{5}->{"TGC"} = "C";	$codes{5}->{"TGG"} = "W";	$codes{5}->{"TGT"} = "C";	
$codes{5}->{"TTA"} = "L";	$codes{5}->{"TTC"} = "F";	$codes{5}->{"TTG"} = "L";	$codes{5}->{"TTT"} = "F";	

#Code: 6 Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear
$code_names{6} = "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear";
$starting_codes{6}->{"ATG"} = 1;

$codes{6}->{"AAA"} = "K";	$codes{6}->{"AAC"} = "N";	$codes{6}->{"AAG"} = "K";	$codes{6}->{"AAT"} = "N";	
$codes{6}->{"ACA"} = "T";	$codes{6}->{"ACC"} = "T";	$codes{6}->{"ACG"} = "T";	$codes{6}->{"ACT"} = "T";	
$codes{6}->{"AGA"} = "R";	$codes{6}->{"AGC"} = "S";	$codes{6}->{"AGG"} = "R";	$codes{6}->{"AGT"} = "S";	
$codes{6}->{"ATA"} = "I";	$codes{6}->{"ATC"} = "I";	$codes{6}->{"ATG"} = "M";	$codes{6}->{"ATT"} = "I";	
$codes{6}->{"CAA"} = "Q";	$codes{6}->{"CAC"} = "H";	$codes{6}->{"CAG"} = "Q";	$codes{6}->{"CAT"} = "H";	
$codes{6}->{"CCA"} = "P";	$codes{6}->{"CCC"} = "P";	$codes{6}->{"CCG"} = "P";	$codes{6}->{"CCT"} = "P";	
$codes{6}->{"CGA"} = "R";	$codes{6}->{"CGC"} = "R";	$codes{6}->{"CGG"} = "R";	$codes{6}->{"CGT"} = "R";	
$codes{6}->{"CTA"} = "L";	$codes{6}->{"CTC"} = "L";	$codes{6}->{"CTG"} = "L";	$codes{6}->{"CTT"} = "L";	
$codes{6}->{"GAA"} = "E";	$codes{6}->{"GAC"} = "D";	$codes{6}->{"GAG"} = "E";	$codes{6}->{"GAT"} = "D";	
$codes{6}->{"GCA"} = "A";	$codes{6}->{"GCC"} = "A";	$codes{6}->{"GCG"} = "A";	$codes{6}->{"GCT"} = "A";	
$codes{6}->{"GGA"} = "G";	$codes{6}->{"GGC"} = "G";	$codes{6}->{"GGG"} = "G";	$codes{6}->{"GGT"} = "G";	
$codes{6}->{"GTA"} = "V";	$codes{6}->{"GTC"} = "V";	$codes{6}->{"GTG"} = "V";	$codes{6}->{"GTT"} = "V";	
$codes{6}->{"TAA"} = "Q";	$codes{6}->{"TAC"} = "Y";	$codes{6}->{"TAG"} = "Q";	$codes{6}->{"TAT"} = "Y";	
$codes{6}->{"TCA"} = "S";	$codes{6}->{"TCC"} = "S";	$codes{6}->{"TCG"} = "S";	$codes{6}->{"TCT"} = "S";	
$codes{6}->{"TGA"} = "X";	$codes{6}->{"TGC"} = "C";	$codes{6}->{"TGG"} = "W";	$codes{6}->{"TGT"} = "C";	
$codes{6}->{"TTA"} = "L";	$codes{6}->{"TTC"} = "F";	$codes{6}->{"TTG"} = "L";	$codes{6}->{"TTT"} = "F";	

#Code: 9 Echinoderm Mitochondrial; Flatworm Mitochondrial
$code_names{9} = "Echinoderm Mitochondrial; Flatworm Mitochondrial";
$starting_codes{9}->{"ATG"} = 1;
$starting_codes{9}->{"GTG"} = 1;

$codes{9}->{"AAA"} = "N";	$codes{9}->{"AAC"} = "N";	$codes{9}->{"AAG"} = "K";	$codes{9}->{"AAT"} = "N";	
$codes{9}->{"ACA"} = "T";	$codes{9}->{"ACC"} = "T";	$codes{9}->{"ACG"} = "T";	$codes{9}->{"ACT"} = "T";	
$codes{9}->{"AGA"} = "S";	$codes{9}->{"AGC"} = "S";	$codes{9}->{"AGG"} = "S";	$codes{9}->{"AGT"} = "S";	
$codes{9}->{"ATA"} = "I";	$codes{9}->{"ATC"} = "I";	$codes{9}->{"ATG"} = "M";	$codes{9}->{"ATT"} = "I";	
$codes{9}->{"CAA"} = "Q";	$codes{9}->{"CAC"} = "H";	$codes{9}->{"CAG"} = "Q";	$codes{9}->{"CAT"} = "H";	
$codes{9}->{"CCA"} = "P";	$codes{9}->{"CCC"} = "P";	$codes{9}->{"CCG"} = "P";	$codes{9}->{"CCT"} = "P";	
$codes{9}->{"CGA"} = "R";	$codes{9}->{"CGC"} = "R";	$codes{9}->{"CGG"} = "R";	$codes{9}->{"CGT"} = "R";	
$codes{9}->{"CTA"} = "L";	$codes{9}->{"CTC"} = "L";	$codes{9}->{"CTG"} = "L";	$codes{9}->{"CTT"} = "L";	
$codes{9}->{"GAA"} = "E";	$codes{9}->{"GAC"} = "D";	$codes{9}->{"GAG"} = "E";	$codes{9}->{"GAT"} = "D";	
$codes{9}->{"GCA"} = "A";	$codes{9}->{"GCC"} = "A";	$codes{9}->{"GCG"} = "A";	$codes{9}->{"GCT"} = "A";	
$codes{9}->{"GGA"} = "G";	$codes{9}->{"GGC"} = "G";	$codes{9}->{"GGG"} = "G";	$codes{9}->{"GGT"} = "G";	
$codes{9}->{"GTA"} = "V";	$codes{9}->{"GTC"} = "V";	$codes{9}->{"GTG"} = "V";	$codes{9}->{"GTT"} = "V";	
$codes{9}->{"TAA"} = "X";	$codes{9}->{"TAC"} = "Y";	$codes{9}->{"TAG"} = "X";	$codes{9}->{"TAT"} = "Y";	
$codes{9}->{"TCA"} = "S";	$codes{9}->{"TCC"} = "S";	$codes{9}->{"TCG"} = "S";	$codes{9}->{"TCT"} = "S";	
$codes{9}->{"TGA"} = "W";	$codes{9}->{"TGC"} = "C";	$codes{9}->{"TGG"} = "W";	$codes{9}->{"TGT"} = "C";	
$codes{9}->{"TTA"} = "L";	$codes{9}->{"TTC"} = "F";	$codes{9}->{"TTG"} = "L";	$codes{9}->{"TTT"} = "F";	

#Code: 10 Euplotid Nuclear
$code_names{10} = "Euplotid Nuclear";
$starting_codes{10}->{"ATG"} = 1;

$codes{10}->{"AAA"} = "K";	$codes{10}->{"AAC"} = "N";	$codes{10}->{"AAG"} = "K";	$codes{10}->{"AAT"} = "N";	
$codes{10}->{"ACA"} = "T";	$codes{10}->{"ACC"} = "T";	$codes{10}->{"ACG"} = "T";	$codes{10}->{"ACT"} = "T";	
$codes{10}->{"AGA"} = "R";	$codes{10}->{"AGC"} = "S";	$codes{10}->{"AGG"} = "R";	$codes{10}->{"AGT"} = "S";	
$codes{10}->{"ATA"} = "I";	$codes{10}->{"ATC"} = "I";	$codes{10}->{"ATG"} = "M";	$codes{10}->{"ATT"} = "I";	
$codes{10}->{"CAA"} = "Q";	$codes{10}->{"CAC"} = "H";	$codes{10}->{"CAG"} = "Q";	$codes{10}->{"CAT"} = "H";	
$codes{10}->{"CCA"} = "P";	$codes{10}->{"CCC"} = "P";	$codes{10}->{"CCG"} = "P";	$codes{10}->{"CCT"} = "P";	
$codes{10}->{"CGA"} = "R";	$codes{10}->{"CGC"} = "R";	$codes{10}->{"CGG"} = "R";	$codes{10}->{"CGT"} = "R";	
$codes{10}->{"CTA"} = "L";	$codes{10}->{"CTC"} = "L";	$codes{10}->{"CTG"} = "L";	$codes{10}->{"CTT"} = "L";	
$codes{10}->{"GAA"} = "E";	$codes{10}->{"GAC"} = "D";	$codes{10}->{"GAG"} = "E";	$codes{10}->{"GAT"} = "D";	
$codes{10}->{"GCA"} = "A";	$codes{10}->{"GCC"} = "A";	$codes{10}->{"GCG"} = "A";	$codes{10}->{"GCT"} = "A";	
$codes{10}->{"GGA"} = "G";	$codes{10}->{"GGC"} = "G";	$codes{10}->{"GGG"} = "G";	$codes{10}->{"GGT"} = "G";	
$codes{10}->{"GTA"} = "V";	$codes{10}->{"GTC"} = "V";	$codes{10}->{"GTG"} = "V";	$codes{10}->{"GTT"} = "V";	
$codes{10}->{"TAA"} = "X";	$codes{10}->{"TAC"} = "Y";	$codes{10}->{"TAG"} = "X";	$codes{10}->{"TAT"} = "Y";	
$codes{10}->{"TCA"} = "S";	$codes{10}->{"TCC"} = "S";	$codes{10}->{"TCG"} = "S";	$codes{10}->{"TCT"} = "S";	
$codes{10}->{"TGA"} = "C";	$codes{10}->{"TGC"} = "C";	$codes{10}->{"TGG"} = "W";	$codes{10}->{"TGT"} = "C";	
$codes{10}->{"TTA"} = "L";	$codes{10}->{"TTC"} = "F";	$codes{10}->{"TTG"} = "L";	$codes{10}->{"TTT"} = "F";	

#Code: 11 Bacterial and Plant Plastid
$code_names{11} = "Bacterial and Plant Plastid";
$starting_codes{11}->{"ATA"} = 1;
$starting_codes{11}->{"ATC"} = 1;
$starting_codes{11}->{"ATG"} = 1;
$starting_codes{11}->{"ATT"} = 1;
$starting_codes{11}->{"CTG"} = 1;
$starting_codes{11}->{"GTG"} = 1;
$starting_codes{11}->{"TTG"} = 1;

$codes{11}->{"AAA"} = "K";	$codes{11}->{"AAC"} = "N";	$codes{11}->{"AAG"} = "K";	$codes{11}->{"AAT"} = "N";	
$codes{11}->{"ACA"} = "T";	$codes{11}->{"ACC"} = "T";	$codes{11}->{"ACG"} = "T";	$codes{11}->{"ACT"} = "T";	
$codes{11}->{"AGA"} = "R";	$codes{11}->{"AGC"} = "S";	$codes{11}->{"AGG"} = "R";	$codes{11}->{"AGT"} = "S";	
$codes{11}->{"ATA"} = "I";	$codes{11}->{"ATC"} = "I";	$codes{11}->{"ATG"} = "M";	$codes{11}->{"ATT"} = "I";	
$codes{11}->{"CAA"} = "Q";	$codes{11}->{"CAC"} = "H";	$codes{11}->{"CAG"} = "Q";	$codes{11}->{"CAT"} = "H";	
$codes{11}->{"CCA"} = "P";	$codes{11}->{"CCC"} = "P";	$codes{11}->{"CCG"} = "P";	$codes{11}->{"CCT"} = "P";	
$codes{11}->{"CGA"} = "R";	$codes{11}->{"CGC"} = "R";	$codes{11}->{"CGG"} = "R";	$codes{11}->{"CGT"} = "R";	
$codes{11}->{"CTA"} = "L";	$codes{11}->{"CTC"} = "L";	$codes{11}->{"CTG"} = "L";	$codes{11}->{"CTT"} = "L";	
$codes{11}->{"GAA"} = "E";	$codes{11}->{"GAC"} = "D";	$codes{11}->{"GAG"} = "E";	$codes{11}->{"GAT"} = "D";	
$codes{11}->{"GCA"} = "A";	$codes{11}->{"GCC"} = "A";	$codes{11}->{"GCG"} = "A";	$codes{11}->{"GCT"} = "A";	
$codes{11}->{"GGA"} = "G";	$codes{11}->{"GGC"} = "G";	$codes{11}->{"GGG"} = "G";	$codes{11}->{"GGT"} = "G";	
$codes{11}->{"GTA"} = "V";	$codes{11}->{"GTC"} = "V";	$codes{11}->{"GTG"} = "V";	$codes{11}->{"GTT"} = "V";	
$codes{11}->{"TAA"} = "X";	$codes{11}->{"TAC"} = "Y";	$codes{11}->{"TAG"} = "X";	$codes{11}->{"TAT"} = "Y";	
$codes{11}->{"TCA"} = "S";	$codes{11}->{"TCC"} = "S";	$codes{11}->{"TCG"} = "S";	$codes{11}->{"TCT"} = "S";	
$codes{11}->{"TGA"} = "X";	$codes{11}->{"TGC"} = "C";	$codes{11}->{"TGG"} = "W";	$codes{11}->{"TGT"} = "C";	
$codes{11}->{"TTA"} = "L";	$codes{11}->{"TTC"} = "F";	$codes{11}->{"TTG"} = "L";	$codes{11}->{"TTT"} = "F";	

#Code: 12 Alternative Yeast Nuclear
$code_names{12} = "Alternative Yeast Nuclear";
$starting_codes{12}->{"ATG"} = 1;
$starting_codes{12}->{"CTG"} = 1;

$codes{12}->{"AAA"} = "K";	$codes{12}->{"AAC"} = "N";	$codes{12}->{"AAG"} = "K";	$codes{12}->{"AAT"} = "N";	
$codes{12}->{"ACA"} = "T";	$codes{12}->{"ACC"} = "T";	$codes{12}->{"ACG"} = "T";	$codes{12}->{"ACT"} = "T";	
$codes{12}->{"AGA"} = "R";	$codes{12}->{"AGC"} = "S";	$codes{12}->{"AGG"} = "R";	$codes{12}->{"AGT"} = "S";	
$codes{12}->{"ATA"} = "I";	$codes{12}->{"ATC"} = "I";	$codes{12}->{"ATG"} = "M";	$codes{12}->{"ATT"} = "I";	
$codes{12}->{"CAA"} = "Q";	$codes{12}->{"CAC"} = "H";	$codes{12}->{"CAG"} = "Q";	$codes{12}->{"CAT"} = "H";	
$codes{12}->{"CCA"} = "P";	$codes{12}->{"CCC"} = "P";	$codes{12}->{"CCG"} = "P";	$codes{12}->{"CCT"} = "P";	
$codes{12}->{"CGA"} = "R";	$codes{12}->{"CGC"} = "R";	$codes{12}->{"CGG"} = "R";	$codes{12}->{"CGT"} = "R";	
$codes{12}->{"CTA"} = "L";	$codes{12}->{"CTC"} = "L";	$codes{12}->{"CTG"} = "S";	$codes{12}->{"CTT"} = "L";	
$codes{12}->{"GAA"} = "E";	$codes{12}->{"GAC"} = "D";	$codes{12}->{"GAG"} = "E";	$codes{12}->{"GAT"} = "D";	
$codes{12}->{"GCA"} = "A";	$codes{12}->{"GCC"} = "A";	$codes{12}->{"GCG"} = "A";	$codes{12}->{"GCT"} = "A";	
$codes{12}->{"GGA"} = "G";	$codes{12}->{"GGC"} = "G";	$codes{12}->{"GGG"} = "G";	$codes{12}->{"GGT"} = "G";	
$codes{12}->{"GTA"} = "V";	$codes{12}->{"GTC"} = "V";	$codes{12}->{"GTG"} = "V";	$codes{12}->{"GTT"} = "V";	
$codes{12}->{"TAA"} = "X";	$codes{12}->{"TAC"} = "Y";	$codes{12}->{"TAG"} = "X";	$codes{12}->{"TAT"} = "Y";	
$codes{12}->{"TCA"} = "S";	$codes{12}->{"TCC"} = "S";	$codes{12}->{"TCG"} = "S";	$codes{12}->{"TCT"} = "S";	
$codes{12}->{"TGA"} = "X";	$codes{12}->{"TGC"} = "C";	$codes{12}->{"TGG"} = "W";	$codes{12}->{"TGT"} = "C";	
$codes{12}->{"TTA"} = "L";	$codes{12}->{"TTC"} = "F";	$codes{12}->{"TTG"} = "L";	$codes{12}->{"TTT"} = "F";	

#Code: 13 Ascidian Mitochondrial
$code_names{13} = "Ascidian Mitochondrial";
$starting_codes{13}->{"ATA"} = 1;
$starting_codes{13}->{"ATG"} = 1;
$starting_codes{13}->{"GTG"} = 1;
$starting_codes{13}->{"TTG"} = 1;

$codes{13}->{"AAA"} = "K";	$codes{13}->{"AAC"} = "N";	$codes{13}->{"AAG"} = "K";	$codes{13}->{"AAT"} = "N";	
$codes{13}->{"ACA"} = "T";	$codes{13}->{"ACC"} = "T";	$codes{13}->{"ACG"} = "T";	$codes{13}->{"ACT"} = "T";	
$codes{13}->{"AGA"} = "G";	$codes{13}->{"AGC"} = "S";	$codes{13}->{"AGG"} = "G";	$codes{13}->{"AGT"} = "S";	
$codes{13}->{"ATA"} = "M";	$codes{13}->{"ATC"} = "I";	$codes{13}->{"ATG"} = "M";	$codes{13}->{"ATT"} = "I";	
$codes{13}->{"CAA"} = "Q";	$codes{13}->{"CAC"} = "H";	$codes{13}->{"CAG"} = "Q";	$codes{13}->{"CAT"} = "H";	
$codes{13}->{"CCA"} = "P";	$codes{13}->{"CCC"} = "P";	$codes{13}->{"CCG"} = "P";	$codes{13}->{"CCT"} = "P";	
$codes{13}->{"CGA"} = "R";	$codes{13}->{"CGC"} = "R";	$codes{13}->{"CGG"} = "R";	$codes{13}->{"CGT"} = "R";	
$codes{13}->{"CTA"} = "L";	$codes{13}->{"CTC"} = "L";	$codes{13}->{"CTG"} = "L";	$codes{13}->{"CTT"} = "L";	
$codes{13}->{"GAA"} = "E";	$codes{13}->{"GAC"} = "D";	$codes{13}->{"GAG"} = "E";	$codes{13}->{"GAT"} = "D";	
$codes{13}->{"GCA"} = "A";	$codes{13}->{"GCC"} = "A";	$codes{13}->{"GCG"} = "A";	$codes{13}->{"GCT"} = "A";	
$codes{13}->{"GGA"} = "G";	$codes{13}->{"GGC"} = "G";	$codes{13}->{"GGG"} = "G";	$codes{13}->{"GGT"} = "G";	
$codes{13}->{"GTA"} = "V";	$codes{13}->{"GTC"} = "V";	$codes{13}->{"GTG"} = "V";	$codes{13}->{"GTT"} = "V";	
$codes{13}->{"TAA"} = "X";	$codes{13}->{"TAC"} = "Y";	$codes{13}->{"TAG"} = "X";	$codes{13}->{"TAT"} = "Y";	
$codes{13}->{"TCA"} = "S";	$codes{13}->{"TCC"} = "S";	$codes{13}->{"TCG"} = "S";	$codes{13}->{"TCT"} = "S";	
$codes{13}->{"TGA"} = "W";	$codes{13}->{"TGC"} = "C";	$codes{13}->{"TGG"} = "W";	$codes{13}->{"TGT"} = "C";	
$codes{13}->{"TTA"} = "L";	$codes{13}->{"TTC"} = "F";	$codes{13}->{"TTG"} = "L";	$codes{13}->{"TTT"} = "F";	

#Code: 14 Alternative Flatworm Mitochondrial
$code_names{14} = "Alternative Flatworm Mitochondrial";
$starting_codes{14}->{"ATG"} = 1;

$codes{14}->{"AAA"} = "N";	$codes{14}->{"AAC"} = "N";	$codes{14}->{"AAG"} = "K";	$codes{14}->{"AAT"} = "N";	
$codes{14}->{"ACA"} = "T";	$codes{14}->{"ACC"} = "T";	$codes{14}->{"ACG"} = "T";	$codes{14}->{"ACT"} = "T";	
$codes{14}->{"AGA"} = "S";	$codes{14}->{"AGC"} = "S";	$codes{14}->{"AGG"} = "S";	$codes{14}->{"AGT"} = "S";	
$codes{14}->{"ATA"} = "I";	$codes{14}->{"ATC"} = "I";	$codes{14}->{"ATG"} = "M";	$codes{14}->{"ATT"} = "I";	
$codes{14}->{"CAA"} = "Q";	$codes{14}->{"CAC"} = "H";	$codes{14}->{"CAG"} = "Q";	$codes{14}->{"CAT"} = "H";	
$codes{14}->{"CCA"} = "P";	$codes{14}->{"CCC"} = "P";	$codes{14}->{"CCG"} = "P";	$codes{14}->{"CCT"} = "P";	
$codes{14}->{"CGA"} = "R";	$codes{14}->{"CGC"} = "R";	$codes{14}->{"CGG"} = "R";	$codes{14}->{"CGT"} = "R";	
$codes{14}->{"CTA"} = "L";	$codes{14}->{"CTC"} = "L";	$codes{14}->{"CTG"} = "L";	$codes{14}->{"CTT"} = "L";	
$codes{14}->{"GAA"} = "E";	$codes{14}->{"GAC"} = "D";	$codes{14}->{"GAG"} = "E";	$codes{14}->{"GAT"} = "D";	
$codes{14}->{"GCA"} = "A";	$codes{14}->{"GCC"} = "A";	$codes{14}->{"GCG"} = "A";	$codes{14}->{"GCT"} = "A";	
$codes{14}->{"GGA"} = "G";	$codes{14}->{"GGC"} = "G";	$codes{14}->{"GGG"} = "G";	$codes{14}->{"GGT"} = "G";	
$codes{14}->{"GTA"} = "V";	$codes{14}->{"GTC"} = "V";	$codes{14}->{"GTG"} = "V";	$codes{14}->{"GTT"} = "V";	
$codes{14}->{"TAA"} = "Y";	$codes{14}->{"TAC"} = "Y";	$codes{14}->{"TAG"} = "X";	$codes{14}->{"TAT"} = "Y";	
$codes{14}->{"TCA"} = "S";	$codes{14}->{"TCC"} = "S";	$codes{14}->{"TCG"} = "S";	$codes{14}->{"TCT"} = "S";	
$codes{14}->{"TGA"} = "W";	$codes{14}->{"TGC"} = "C";	$codes{14}->{"TGG"} = "W";	$codes{14}->{"TGT"} = "C";	
$codes{14}->{"TTA"} = "L";	$codes{14}->{"TTC"} = "F";	$codes{14}->{"TTG"} = "L";	$codes{14}->{"TTT"} = "F";	

#Code: 15 Blepharisma Macronuclear
$code_names{15} = "Blepharisma Macronuclear";
$starting_codes{15}->{"ATG"} = 1;

$codes{15}->{"AAA"} = "K";	$codes{15}->{"AAC"} = "N";	$codes{15}->{"AAG"} = "K";	$codes{15}->{"AAT"} = "N";	
$codes{15}->{"ACA"} = "T";	$codes{15}->{"ACC"} = "T";	$codes{15}->{"ACG"} = "T";	$codes{15}->{"ACT"} = "T";	
$codes{15}->{"AGA"} = "R";	$codes{15}->{"AGC"} = "S";	$codes{15}->{"AGG"} = "R";	$codes{15}->{"AGT"} = "S";	
$codes{15}->{"ATA"} = "I";	$codes{15}->{"ATC"} = "I";	$codes{15}->{"ATG"} = "M";	$codes{15}->{"ATT"} = "I";	
$codes{15}->{"CAA"} = "Q";	$codes{15}->{"CAC"} = "H";	$codes{15}->{"CAG"} = "Q";	$codes{15}->{"CAT"} = "H";	
$codes{15}->{"CCA"} = "P";	$codes{15}->{"CCC"} = "P";	$codes{15}->{"CCG"} = "P";	$codes{15}->{"CCT"} = "P";	
$codes{15}->{"CGA"} = "R";	$codes{15}->{"CGC"} = "R";	$codes{15}->{"CGG"} = "R";	$codes{15}->{"CGT"} = "R";	
$codes{15}->{"CTA"} = "L";	$codes{15}->{"CTC"} = "L";	$codes{15}->{"CTG"} = "L";	$codes{15}->{"CTT"} = "L";	
$codes{15}->{"GAA"} = "E";	$codes{15}->{"GAC"} = "D";	$codes{15}->{"GAG"} = "E";	$codes{15}->{"GAT"} = "D";	
$codes{15}->{"GCA"} = "A";	$codes{15}->{"GCC"} = "A";	$codes{15}->{"GCG"} = "A";	$codes{15}->{"GCT"} = "A";	
$codes{15}->{"GGA"} = "G";	$codes{15}->{"GGC"} = "G";	$codes{15}->{"GGG"} = "G";	$codes{15}->{"GGT"} = "G";	
$codes{15}->{"GTA"} = "V";	$codes{15}->{"GTC"} = "V";	$codes{15}->{"GTG"} = "V";	$codes{15}->{"GTT"} = "V";	
$codes{15}->{"TAA"} = "X";	$codes{15}->{"TAC"} = "Y";	$codes{15}->{"TAG"} = "Q";	$codes{15}->{"TAT"} = "Y";	
$codes{15}->{"TCA"} = "S";	$codes{15}->{"TCC"} = "S";	$codes{15}->{"TCG"} = "S";	$codes{15}->{"TCT"} = "S";	
$codes{15}->{"TGA"} = "X";	$codes{15}->{"TGC"} = "C";	$codes{15}->{"TGG"} = "W";	$codes{15}->{"TGT"} = "C";	
$codes{15}->{"TTA"} = "L";	$codes{15}->{"TTC"} = "F";	$codes{15}->{"TTG"} = "L";	$codes{15}->{"TTT"} = "F";	

#Code: 16 Chlorophycean Mitochondrial
$code_names{16} = "Chlorophycean Mitochondrial";
$starting_codes{16}->{"ATG"} = 1;

$codes{16}->{"AAA"} = "K";	$codes{16}->{"AAC"} = "N";	$codes{16}->{"AAG"} = "K";	$codes{16}->{"AAT"} = "N";	
$codes{16}->{"ACA"} = "T";	$codes{16}->{"ACC"} = "T";	$codes{16}->{"ACG"} = "T";	$codes{16}->{"ACT"} = "T";	
$codes{16}->{"AGA"} = "R";	$codes{16}->{"AGC"} = "S";	$codes{16}->{"AGG"} = "R";	$codes{16}->{"AGT"} = "S";	
$codes{16}->{"ATA"} = "I";	$codes{16}->{"ATC"} = "I";	$codes{16}->{"ATG"} = "M";	$codes{16}->{"ATT"} = "I";	
$codes{16}->{"CAA"} = "Q";	$codes{16}->{"CAC"} = "H";	$codes{16}->{"CAG"} = "Q";	$codes{16}->{"CAT"} = "H";	
$codes{16}->{"CCA"} = "P";	$codes{16}->{"CCC"} = "P";	$codes{16}->{"CCG"} = "P";	$codes{16}->{"CCT"} = "P";	
$codes{16}->{"CGA"} = "R";	$codes{16}->{"CGC"} = "R";	$codes{16}->{"CGG"} = "R";	$codes{16}->{"CGT"} = "R";	
$codes{16}->{"CTA"} = "L";	$codes{16}->{"CTC"} = "L";	$codes{16}->{"CTG"} = "L";	$codes{16}->{"CTT"} = "L";	
$codes{16}->{"GAA"} = "E";	$codes{16}->{"GAC"} = "D";	$codes{16}->{"GAG"} = "E";	$codes{16}->{"GAT"} = "D";	
$codes{16}->{"GCA"} = "A";	$codes{16}->{"GCC"} = "A";	$codes{16}->{"GCG"} = "A";	$codes{16}->{"GCT"} = "A";	
$codes{16}->{"GGA"} = "G";	$codes{16}->{"GGC"} = "G";	$codes{16}->{"GGG"} = "G";	$codes{16}->{"GGT"} = "G";	
$codes{16}->{"GTA"} = "V";	$codes{16}->{"GTC"} = "V";	$codes{16}->{"GTG"} = "V";	$codes{16}->{"GTT"} = "V";	
$codes{16}->{"TAA"} = "X";	$codes{16}->{"TAC"} = "Y";	$codes{16}->{"TAG"} = "L";	$codes{16}->{"TAT"} = "Y";	
$codes{16}->{"TCA"} = "S";	$codes{16}->{"TCC"} = "S";	$codes{16}->{"TCG"} = "S";	$codes{16}->{"TCT"} = "S";	
$codes{16}->{"TGA"} = "X";	$codes{16}->{"TGC"} = "C";	$codes{16}->{"TGG"} = "W";	$codes{16}->{"TGT"} = "C";	
$codes{16}->{"TTA"} = "L";	$codes{16}->{"TTC"} = "F";	$codes{16}->{"TTG"} = "L";	$codes{16}->{"TTT"} = "F";	

#Code: 21 Trematode Mitochondrial
$code_names{21} = "Trematode Mitochondrial";
$starting_codes{21}->{"ATG"} = 1;
$starting_codes{21}->{"GTG"} = 1;

$codes{21}->{"AAA"} = "N";	$codes{21}->{"AAC"} = "N";	$codes{21}->{"AAG"} = "K";	$codes{21}->{"AAT"} = "N";	
$codes{21}->{"ACA"} = "T";	$codes{21}->{"ACC"} = "T";	$codes{21}->{"ACG"} = "T";	$codes{21}->{"ACT"} = "T";	
$codes{21}->{"AGA"} = "S";	$codes{21}->{"AGC"} = "S";	$codes{21}->{"AGG"} = "S";	$codes{21}->{"AGT"} = "S";	
$codes{21}->{"ATA"} = "M";	$codes{21}->{"ATC"} = "I";	$codes{21}->{"ATG"} = "M";	$codes{21}->{"ATT"} = "I";	
$codes{21}->{"CAA"} = "Q";	$codes{21}->{"CAC"} = "H";	$codes{21}->{"CAG"} = "Q";	$codes{21}->{"CAT"} = "H";	
$codes{21}->{"CCA"} = "P";	$codes{21}->{"CCC"} = "P";	$codes{21}->{"CCG"} = "P";	$codes{21}->{"CCT"} = "P";	
$codes{21}->{"CGA"} = "R";	$codes{21}->{"CGC"} = "R";	$codes{21}->{"CGG"} = "R";	$codes{21}->{"CGT"} = "R";	
$codes{21}->{"CTA"} = "L";	$codes{21}->{"CTC"} = "L";	$codes{21}->{"CTG"} = "L";	$codes{21}->{"CTT"} = "L";	
$codes{21}->{"GAA"} = "E";	$codes{21}->{"GAC"} = "D";	$codes{21}->{"GAG"} = "E";	$codes{21}->{"GAT"} = "D";	
$codes{21}->{"GCA"} = "A";	$codes{21}->{"GCC"} = "A";	$codes{21}->{"GCG"} = "A";	$codes{21}->{"GCT"} = "A";	
$codes{21}->{"GGA"} = "G";	$codes{21}->{"GGC"} = "G";	$codes{21}->{"GGG"} = "G";	$codes{21}->{"GGT"} = "G";	
$codes{21}->{"GTA"} = "V";	$codes{21}->{"GTC"} = "V";	$codes{21}->{"GTG"} = "V";	$codes{21}->{"GTT"} = "V";	
$codes{21}->{"TAA"} = "X";	$codes{21}->{"TAC"} = "Y";	$codes{21}->{"TAG"} = "X";	$codes{21}->{"TAT"} = "Y";	
$codes{21}->{"TCA"} = "S";	$codes{21}->{"TCC"} = "S";	$codes{21}->{"TCG"} = "S";	$codes{21}->{"TCT"} = "S";	
$codes{21}->{"TGA"} = "W";	$codes{21}->{"TGC"} = "C";	$codes{21}->{"TGG"} = "W";	$codes{21}->{"TGT"} = "C";	
$codes{21}->{"TTA"} = "L";	$codes{21}->{"TTC"} = "F";	$codes{21}->{"TTG"} = "L";	$codes{21}->{"TTT"} = "F";	

#Code: 22 Scenedesmus obliquus Mitochondrial
$code_names{22} = "Scenedesmus obliquus Mitochondrial";
$starting_codes{22}->{"ATG"} = 1;

$codes{22}->{"AAA"} = "K";	$codes{22}->{"AAC"} = "N";	$codes{22}->{"AAG"} = "K";	$codes{22}->{"AAT"} = "N";	
$codes{22}->{"ACA"} = "T";	$codes{22}->{"ACC"} = "T";	$codes{22}->{"ACG"} = "T";	$codes{22}->{"ACT"} = "T";	
$codes{22}->{"AGA"} = "R";	$codes{22}->{"AGC"} = "S";	$codes{22}->{"AGG"} = "R";	$codes{22}->{"AGT"} = "S";	
$codes{22}->{"ATA"} = "I";	$codes{22}->{"ATC"} = "I";	$codes{22}->{"ATG"} = "M";	$codes{22}->{"ATT"} = "I";	
$codes{22}->{"CAA"} = "Q";	$codes{22}->{"CAC"} = "H";	$codes{22}->{"CAG"} = "Q";	$codes{22}->{"CAT"} = "H";	
$codes{22}->{"CCA"} = "P";	$codes{22}->{"CCC"} = "P";	$codes{22}->{"CCG"} = "P";	$codes{22}->{"CCT"} = "P";	
$codes{22}->{"CGA"} = "R";	$codes{22}->{"CGC"} = "R";	$codes{22}->{"CGG"} = "R";	$codes{22}->{"CGT"} = "R";	
$codes{22}->{"CTA"} = "L";	$codes{22}->{"CTC"} = "L";	$codes{22}->{"CTG"} = "L";	$codes{22}->{"CTT"} = "L";	
$codes{22}->{"GAA"} = "E";	$codes{22}->{"GAC"} = "D";	$codes{22}->{"GAG"} = "E";	$codes{22}->{"GAT"} = "D";	
$codes{22}->{"GCA"} = "A";	$codes{22}->{"GCC"} = "A";	$codes{22}->{"GCG"} = "A";	$codes{22}->{"GCT"} = "A";	
$codes{22}->{"GGA"} = "G";	$codes{22}->{"GGC"} = "G";	$codes{22}->{"GGG"} = "G";	$codes{22}->{"GGT"} = "G";	
$codes{22}->{"GTA"} = "V";	$codes{22}->{"GTC"} = "V";	$codes{22}->{"GTG"} = "V";	$codes{22}->{"GTT"} = "V";	
$codes{22}->{"TAA"} = "X";	$codes{22}->{"TAC"} = "Y";	$codes{22}->{"TAG"} = "L";	$codes{22}->{"TAT"} = "Y";	
$codes{22}->{"TCA"} = "X";	$codes{22}->{"TCC"} = "S";	$codes{22}->{"TCG"} = "S";	$codes{22}->{"TCT"} = "S";	
$codes{22}->{"TGA"} = "X";	$codes{22}->{"TGC"} = "C";	$codes{22}->{"TGG"} = "W";	$codes{22}->{"TGT"} = "C";	
$codes{22}->{"TTA"} = "L";	$codes{22}->{"TTC"} = "F";	$codes{22}->{"TTG"} = "L";	$codes{22}->{"TTT"} = "F";	

#Code: 23 Thraustochytrium Mitochondrial
$code_names{23} = "Thraustochytrium Mitochondrial";
$starting_codes{23}->{"ATG"} = 1;
$starting_codes{23}->{"ATT"} = 1;
$starting_codes{23}->{"GTG"} = 1;

$codes{23}->{"AAA"} = "K";	$codes{23}->{"AAC"} = "N";	$codes{23}->{"AAG"} = "K";	$codes{23}->{"AAT"} = "N";	
$codes{23}->{"ACA"} = "T";	$codes{23}->{"ACC"} = "T";	$codes{23}->{"ACG"} = "T";	$codes{23}->{"ACT"} = "T";	
$codes{23}->{"AGA"} = "R";	$codes{23}->{"AGC"} = "S";	$codes{23}->{"AGG"} = "R";	$codes{23}->{"AGT"} = "S";	
$codes{23}->{"ATA"} = "I";	$codes{23}->{"ATC"} = "I";	$codes{23}->{"ATG"} = "M";	$codes{23}->{"ATT"} = "I";	
$codes{23}->{"CAA"} = "Q";	$codes{23}->{"CAC"} = "H";	$codes{23}->{"CAG"} = "Q";	$codes{23}->{"CAT"} = "H";	
$codes{23}->{"CCA"} = "P";	$codes{23}->{"CCC"} = "P";	$codes{23}->{"CCG"} = "P";	$codes{23}->{"CCT"} = "P";	
$codes{23}->{"CGA"} = "R";	$codes{23}->{"CGC"} = "R";	$codes{23}->{"CGG"} = "R";	$codes{23}->{"CGT"} = "R";	
$codes{23}->{"CTA"} = "L";	$codes{23}->{"CTC"} = "L";	$codes{23}->{"CTG"} = "L";	$codes{23}->{"CTT"} = "L";	
$codes{23}->{"GAA"} = "E";	$codes{23}->{"GAC"} = "D";	$codes{23}->{"GAG"} = "E";	$codes{23}->{"GAT"} = "D";	
$codes{23}->{"GCA"} = "A";	$codes{23}->{"GCC"} = "A";	$codes{23}->{"GCG"} = "A";	$codes{23}->{"GCT"} = "A";	
$codes{23}->{"GGA"} = "G";	$codes{23}->{"GGC"} = "G";	$codes{23}->{"GGG"} = "G";	$codes{23}->{"GGT"} = "G";	
$codes{23}->{"GTA"} = "V";	$codes{23}->{"GTC"} = "V";	$codes{23}->{"GTG"} = "V";	$codes{23}->{"GTT"} = "V";	
$codes{23}->{"TAA"} = "X";	$codes{23}->{"TAC"} = "Y";	$codes{23}->{"TAG"} = "X";	$codes{23}->{"TAT"} = "Y";	
$codes{23}->{"TCA"} = "S";	$codes{23}->{"TCC"} = "S";	$codes{23}->{"TCG"} = "S";	$codes{23}->{"TCT"} = "S";	
$codes{23}->{"TGA"} = "X";	$codes{23}->{"TGC"} = "C";	$codes{23}->{"TGG"} = "W";	$codes{23}->{"TGT"} = "C";	
$codes{23}->{"TTA"} = "X";	$codes{23}->{"TTC"} = "F";	$codes{23}->{"TTG"} = "L";	$codes{23}->{"TTT"} = "F";	


#Our own codes:
#Code: 100 Ancestral Arthropod Mitochondrial Code (AGG=K)
$code_names{100} = "Ancestral Arthropod Mitochondrial Code (AGG=K)";
$starting_codes{100}->{"ATA"} = 1;
$starting_codes{100}->{"ATC"} = 1;
$starting_codes{100}->{"ATG"} = 1;
$starting_codes{100}->{"ATT"} = 1;
$starting_codes{100}->{"GTG"} = 1;
$starting_codes{100}->{"TTG"} = 1;

$codes{100}->{"AAA"} = "K";	$codes{100}->{"AAC"} = "N";	$codes{100}->{"AAG"} = "K";	$codes{100}->{"AAT"} = "N";	
$codes{100}->{"ACA"} = "T";	$codes{100}->{"ACC"} = "T";	$codes{100}->{"ACG"} = "T";	$codes{100}->{"ACT"} = "T";	
$codes{100}->{"AGA"} = "S";	$codes{100}->{"AGC"} = "S";	$codes{100}->{"AGG"} = "K";	$codes{100}->{"AGT"} = "S";	
$codes{100}->{"ATA"} = "M";	$codes{100}->{"ATC"} = "I";	$codes{100}->{"ATG"} = "M";	$codes{100}->{"ATT"} = "I";	
$codes{100}->{"CAA"} = "Q";	$codes{100}->{"CAC"} = "H";	$codes{100}->{"CAG"} = "Q";	$codes{100}->{"CAT"} = "H";	
$codes{100}->{"CCA"} = "P";	$codes{100}->{"CCC"} = "P";	$codes{100}->{"CCG"} = "P";	$codes{100}->{"CCT"} = "P";	
$codes{100}->{"CGA"} = "R";	$codes{100}->{"CGC"} = "R";	$codes{100}->{"CGG"} = "R";	$codes{100}->{"CGT"} = "R";	
$codes{100}->{"CTA"} = "L";	$codes{100}->{"CTC"} = "L";	$codes{100}->{"CTG"} = "L";	$codes{100}->{"CTT"} = "L";	
$codes{100}->{"GAA"} = "E";	$codes{100}->{"GAC"} = "D";	$codes{100}->{"GAG"} = "E";	$codes{100}->{"GAT"} = "D";	
$codes{100}->{"GCA"} = "A";	$codes{100}->{"GCC"} = "A";	$codes{100}->{"GCG"} = "A";	$codes{100}->{"GCT"} = "A";	
$codes{100}->{"GGA"} = "G";	$codes{100}->{"GGC"} = "G";	$codes{100}->{"GGG"} = "G";	$codes{100}->{"GGT"} = "G";	
$codes{100}->{"GTA"} = "V";	$codes{100}->{"GTC"} = "V";	$codes{100}->{"GTG"} = "V";	$codes{100}->{"GTT"} = "V";	
$codes{100}->{"TAA"} = "X";	$codes{100}->{"TAC"} = "Y";	$codes{100}->{"TAG"} = "X";	$codes{100}->{"TAT"} = "Y";	
$codes{100}->{"TCA"} = "S";	$codes{100}->{"TCC"} = "S";	$codes{100}->{"TCG"} = "S";	$codes{100}->{"TCT"} = "S";	
$codes{100}->{"TGA"} = "W";	$codes{100}->{"TGC"} = "C";	$codes{100}->{"TGG"} = "W";	$codes{100}->{"TGT"} = "C";	
$codes{100}->{"TTA"} = "L";	$codes{100}->{"TTC"} = "F";	$codes{100}->{"TTG"} = "L";	$codes{100}->{"TTT"} = "F";	


#Code: 101 Hemichordate Mitochondrial
$code_names{101} = "Hemichordate Mitochondrial";

$starting_codes{101}->{"ATA"} = 1;
$starting_codes{101}->{"ATC"} = 1;
$starting_codes{101}->{"ATG"} = 1;
$starting_codes{101}->{"ATT"} = 1;
$starting_codes{101}->{"GTG"} = 1;
$starting_codes{101}->{"TTG"} = 1;


my %tmp  = (	"TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" => "L",
				"TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S",
				"TAT" => "Y", "TAC" => "Y", "TAA" => "X", "TAG" => "X",
				"TGT" => "C", "TGC" => "C", "TGA" => "W", "TGG" => "W",
				"CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L",
				"CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
				"CAT" => "H", "CAC" => "H", "CAA" => "Q", "CAG" => "Q",
				"CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R",
				"ATT" => "I", "ATC" => "I", "ATA" => "I", "ATG" => "M",
				"ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
				"AAT" => "N", "AAC" => "N", "AAA" => "K", "AAG" => "K",
				"AGT" => "S", "AGC" => "S", "AGA" => "S", "AGG" => "S",
				"GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
				"GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
				"GAT" => "D", "GAC" => "D", "GAA" => "E", "GAG" => "E",
				"GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G", "---" => "-");
$codes{101} = \%tmp;

foreach my $code ( keys %codes ) {		$codes{$code}->{"---"} = "-";	}





	
# #########################################################################
# 
# Load up unaligned nucleotide file
# 
# 
# #########################################################################


# Read in nucleotide file and put it into a hash keyed by the species name	

open (FILE1, $nt_file)  || print_in_H_and_exit(\*H, "Unable to write to nt_file: $nt_file\n"); # put dcse file into hash

$/ = '>'; # alter record separator for fasta/nbrf file

while (<FILE1>) {
	chomp;
	if ($_ =~ m/\w/i) {
		if ($_ =~ m/\w\w;/ and $_ =~ /\*/i) { # This means it is an nbrf/pir file
			($firstline, $nt_name, @nt_seq) = split "\n"; # name into $name, several seq lines into @Sequence
		}
		else { # its a FASTA file
			($nt_name, @nt_seq) = split "\n"; # name into $name, several seq lines into @Sequence
		}
		$nt_seq = join ('',@nt_seq); #now join up elements in @n_seq to create $n_seq.
		$nt_name =~ s/[ \t\r]+$//g;
		$nt_name =~ s/\s/_/g;
		$nt_name =~ s/[\(\)]//g;# getting rid of brackets
		$nt_seq =~ s/\s//g;
		$nt_name =~ s/[\[\]]//g;# not getting rid of square brackets
		$nt_nameseq{$nt_name} = $nt_seq;
		push @nt_spnames, $nt_name;
	}
	else {
		next;
	}
}
	
	

# #########################################################################
# 
# 				Decide how to use the software
# 
# 
# #########################################################################


if($program eq "")	{	$function = 2;	}
else 				{	$function = 1;	}


# #########################################################################
# 
# 			Translate nucleotides, align the amino acids and 
# 			then go to the end bit where they are recombined
# 
# #########################################################################

if ($function == 1) {

	#print "\n\nAligned nucleotides will have suffix .fasta.\nClustal output files will have suffixes .aln .dnd and .pir (aligned aa seqs)\nT_coffee output files will have suffixes .dnd and .pir (aligned aa seqs)\n\nOutput file prefix....\n";
	$clustalinfile = $out_file. ".aaseqs"; #BY FEDE

	
	
	# Find out what genetic code is to be used
	
	#print "Do your taxa share the same genetic code? (Y/N)\n";
	if($geneticcode != 0)	{	$samecodeYN = "Y";	}
	else 					{	$samecodeYN = "N";	}
	
	if ($samecodeYN eq "Y") {
		$codehash = \%{$codes{$geneticcode}};
		if($geneticcode > 23 && $geneticcode != 100 && $geneticcode != 101)
		{		print STDERR "Genetic code $geneticcode is not a valid choice\n";
				print H      "<font color=\"red\">Genetic code $geneticcode is not a valid choice</font>\n";
				exit 1;						}
		if($geneticcode == 7 || $geneticcode == 8 || $geneticcode == 17 || $geneticcode == 18 || $geneticcode == 19 || $geneticcode == 20)
		{		print STDERR "Genetic code $geneticcode is not a valid choice\n";
				print H      "<font color=\"red\"><b>Genetic code $geneticcode is not a valid choice. Exiting...</b></font>\n";
				exit 1;						}
		%startingCodons = %{$starting_codes{$geneticcode}};
		if($geneticcode > 23 && $geneticcode != 100 && $geneticcode != 101)
		{		print STDERR "Genetic code $geneticcode is not a valid choice\n";
				print H      "<font color=\"red\"><b>Genetic code $geneticcode is not a valid choice. Exiting...</b></font>\n";
				exit 1;						}
	} else {
	# If codes vary between taxa they can be downloaded from a file or input one by one.
		$spcode = 0;
		$namecodefile = $geneticcode_file;
		$/ = $old_input_rec_sep;
		open (NAMECODEFILE, $namecodefile) || print_in_H_and_exit(\*H,  "Unable to open namecodefile: $namecodefile\n"); 
		while (<NAMECODEFILE>) {
			chomp;
			($namecode_name, $namecode_code) = (split(/\t/,$_))[0,1];
			$namecode_name =~ s/[ \t\r\n]+$//g;
			$namecode_name =~ s/\s/_/g;
			$namecode_name =~ s/[\(\)]//g;
			$namecode_code =~ s/\s//g;
			$namecode{$namecode_name} = $namecode_code;
		}
		close NAMECODEFILE;
	}
	
	
# Open up the file to be sent to clustal.  
# This has to have each nucleotide translated into amino acids according to its appropriate code.
	$muscleinfile=$out_file . '.aaseqs.fasta';
	$xcounter =0;
	open (CLUSTALOUT, ">$clustalinfile");
	open (MUSCLEOUT, ">$muscleinfile");
	foreach (@nt_spnames) {
			if ($samecodeYN eq "N") {
				if (! $namecode{$_}){
					print STDERR "No genetic code defined for $_. Assuming the Standard Genetic Code\n"; 
					print H      "<font color=\"orange\">No genetic code was defined for $_. Assuming the Standard Genetic Code.</font><br>\n"; 
					$namecode{$_} = 1;
					#exit 1;
				}
				$code = $namecode{$_};
				if($code > 23 && $code != 100 && $code != 101)
				{		print STDERR "Genetic code $code is not a valid choice\n";
						print H      "<font color=\"red\"><b>Genetic code $code is not a valid choice. Exiting...</b></font>\n";
						exit 1;				}
				if($code == 7 || $code == 8 || $code == 17 || $code == 18 || $code == 19 || $code == 20)
				{		print STDERR "Genetic code $code is not a valid choice\n";
						print H      "<font color=\"red\"><b>Genetic code $code is not a valid choice. Exiting...</b></font>\n";
						exit 1;				}
				$codehash = \%{$codes{$code}};
				%startingCodons = %{$starting_codes{$code}};
			}

		#NEW!!!! (added by Fede)
		if($guessRF == 1) {
			my $min_stop_codons = 9999999999;
			my $bestrf = 999;
			my $best_seq = "";
			my $stops_in_first = 9999;
			my $tmpseq = $nt_nameseq{$_};
			foreach my $rf (1,2,3) {
				my $stop_counter = 0;
				my @tmpnt = $tmpseq =~ /(.{3})/g;
				foreach my $triplet ( @tmpnt ) {
					$triplet = uc($triplet);
					if($$codehash{$triplet} && $$codehash{$triplet} =~ /x/i) {
						$stop_counter++;
					}
				}
				if($rf == 1) {	$stops_in_first = $stop_counter;	}
				print STDERR "$_, $rf: $stop_counter\n";
				if($stop_counter < $min_stop_codons) {
					$bestrf = $rf;
					$min_stop_codons = $stop_counter;
					$best_seq = $tmpseq;
				}
				$tmpseq =~ s/^.//;#we remove first base.
				print STDERR "First seven chars: ".substr($tmpseq,0,7)."\n";
			}
			$nt_nameseq{$_} = $best_seq;
			if($bestrf != 1) {
				print H "<font color='brown'>$_: reading frame +$bestrf better than default rf +1 ($min_stop_codons vs $stops_in_first stop codons)</font><br>\n";
			}
		}
		#end of NEW.


		#POR HACER: aqu sera un buen sitio para calcular el contenido GC original!!
		my $is_first_codon = 1;
		my @tempseqarray;
		$temptriplet = '';
		my @tempaminoacidseq;
		@tempseqarray = $nt_nameseq{$_} =~ /(.{3})/g;
		foreach  $temptriplet (@tempseqarray) {
			$temptriplet = (uc $temptriplet);
			if ($$codehash{$temptriplet}) {
				if ($$codehash{$temptriplet} =~/x/i) { #stop codon
					$xcounter++;
				}
				if($is_first_codon == 1 && exists($startingCodons{$temptriplet})) {
					push  (@tempaminoacidseq, "M");
				} else {
					push  (@tempaminoacidseq, $$codehash{$temptriplet});
				}
			} else {
				my $aa_ = &tryToSolveCodonUncertainty($temptriplet,$codehash,\*H);
				if($aa_ ne "X") {
					push(@tempaminoacidseq, $aa_);
					print H "<font color=\"orange\">Codon $temptriplet successfully disambiguated as '$aa_'</font><br>\n";
				} else {
					push (@tempaminoacidseq, "X");
					print H "<font color=\"red\">Warning: codon $temptriplet translated as 'X'</font><br>\n";
				}
			}
			$is_first_codon = 0;
		}
		if ($xcounter >1) {
			print H "<font color=\"red\"><b>Error: Species $_ has $xcounter termination codons. Are you sure you selected the proper genetic code? Alternatively, is the sequence in the correct reading frame?</b></font><br>\n";
			print STDERR "in file $nt_file Species $_ has $xcounter termination codons.  Suggest check correct code used and that in frame\n"
		}
		print CLUSTALOUT">P1;$_\n$_\n",@tempaminoacidseq,"*\n";
		$xcounter =0;
		print MUSCLEOUT">$_\n",@tempaminoacidseq, "\n";
	}			
	
	close CLUSTALOUT;
	close MUSCLEOUT;
	close FILE1;
	
	
	
	$clustaloutfile = $out_file . '.aa_ali.pir';
	$muscleoutfile  = $out_file . '.aa_ali.fasta';
	
	
	
	$pauser = $program;
	chomp $pauser;
	
# Make a file that will send instructions to Clustal with appropriate name of input file
	
	print STDERR ">>>>>>pauser is: $pauser\n";
	
	if ($pauser eq "clustalw") {
	
		open (CLUSTALINPUT, ">$out_file".".clustal.input");
		print CLUSTALINPUT "1\n$clustalinfile\n2\n9\n1\n2\n\n1\n\n\nx\n\nx\n\nx\n\n";
		close CLUSTALINPUT;
		
		# Start up CLustalW and run it according to the instruction file just made.
		
		print STDERR "Executing clustalw...";
		open(LOG, ">$out_file.$pauser.log") || die "Cannot create file outfile.pauser.log: $out_file.$pauser.log\n";
		my $pid = open(PH, "$path"."clustalw < $out_file".".clustal.input 2>&1 |");              # or with an open pipe
		while(<PH>) {
			print STDERR $_;
			print LOG $_;
		}
		close(LOG);
		#system "$path/clustalw < $out_file".".clustal.input";
		print STDERR "...done!\n";
		print STDERR "$path"."clustalw < $out_file".".clustal.input\n";
		print STDERR "Renaming : ", $out_file.".pir", " to ", $clustaloutfile, "\n";
		rename($out_file.".pir", $clustaloutfile);
		rename($out_file.".dnd", "$out_file.aa_ali.dnd");
	}
	
	elsif ($pauser eq "t_coffee") { #POR HACER: cambiar esto. Era t-coffee
		if(!exists($ENV{'HOME'})) {
			$ENV{'HOME'}  = "$path"."tmp";
			$ENV{'PATH'}  = "$path"."tmp:/usr/local/bin:/usr/bin:$ENV{'PATH'}";
			$ENV{'HOME_4_TCOFFEE'} = "$path"."tmp";
			$ENV{'VAR'}   = "$path"."tmp";
		}
		$t_coffeeoutname=$out_file . '.aa_ali.pir';
		my $dir = getcwd;
		chdir("$path/tmp");
		open(LOG, ">$out_file.$pauser.log") || die "Cannot create file outfile.pauser.log: $out_file.$pauser.log\n";
		my $pid = open(PH, "/usr/local/bin/t_coffee -infile=$muscleinfile -output=pir_aln -outfile=$clustaloutfile 2>&1 |");              # or with an open pipe
		while(<PH>) {
			print STDERR $_;
			print LOG $_;
		}
		close(LOG);
		#system "/usr/local/bin/t_coffee -infile=$muscleinfile -output=pir_aln -outfile=$clustaloutfile";#-outfile=$t_coffeeoutname";
		print STDERR "/usr/local/bin/t_coffee -infile=$muscleinfile -output=pir_aln -outfile=$clustaloutfile\n";#-outfile=$t_coffeeoutname";
		chdir($dir);
	}
	
	elsif ($pauser eq "muscle") {
		$muscleoutfile=$out_file . '.aa_ali.fasta';	
		open(LOG, ">$out_file.$pauser.log") || die "Cannot create file outfile.pauser.log: $out_file.$pauser.log\n";
		my $pid = open(PH, "$path"."muscle -in $muscleinfile -out $muscleoutfile 2>&1 |");              # or with an open pipe
		while(<PH>) {
			print STDERR $_;
			print LOG $_;
		}
		close(LOG);
		#system "$path/muscle -in $muscleinfile -out $muscleoutfile";
		print STDERR "$path"."muscle -in $muscleinfile -out $muscleoutfile\n";
	}
	
	elsif ($pauser eq "mafft") {
		$muscleoutfile=$out_file . '.aa_ali.fasta';	
		#print STDERR "/usr/local/bin/mafft --auto $muscleinfile > $muscleoutfile\n"; # commented SAC 210713
		print STDERR "mafft --auto $muscleinfile > $muscleoutfile\n";
		open(LOG, ">$out_file.$pauser.log") || die "Cannot create file outfile.pauser.log: $out_file.$pauser.log\n";
		#my $pid = open(PH, "/usr/local/bin/mafft --auto $muscleinfile 2>&1 1>$muscleoutfile |");              # or with an open pipe # commented SAC 210713
		my $pid = open(PH, "mafft --auto $muscleinfile 2>&1 1>$muscleoutfile |");              # or with an open pipe
		while(<PH>) {
			print STDERR $_;
			print LOG $_;
		}
		close(LOG);
		#system "/usr/local/bin/mafft --auto $muscleinfile > $muscleoutfile";
	}
	
	elsif ($pauser eq "prank") {
		print STDERR "$path"."prank -d=$clustalinfile -o=$clustalinfile -F -noxml -notree\n";
		open(LOG, ">$out_file.$pauser.log") || die "Cannot create file outfile.pauser.log: $out_file.$pauser.log\n";
		my $pid = open(PH, "$path"."prank -d=$clustalinfile -o=$clustalinfile -F -noxml -notree 2>&1 |");              # or with an open pipe
		while(<PH>) {
			print STDERR $_;
			print LOG $_;
		}
		close(LOG);
		#system "$path/prank -d=$clustalinfile -o=$clustalinfile -F -noxml -notree";
	}
	
	else {	print STDERR "not a valid choice\n";
			print H      "<font color=\"red\">not a valid choice</font>\n";
			exit 1;
	}


# Open up the clustal alignment output file and put the ALIGNED amino acids into a hash keyed by the species name
	if ($pauser eq "clustalw" || $pauser eq "tcoffee") {
		$clustaloutfile = $out_file . '.aa_ali.pir';
	} elsif ($pauser eq "mafft") {
		$clustaloutfile = $out_file . '.aa_ali.fasta';
	} elsif ($pauser eq "muscle") {
		$clustaloutfile = $out_file . '.aa_ali.fasta';
	} elsif ($pauser eq "prank") {
		$clustaloutfile = $out_file . '.aaseqs.2.fas';
	}
	
} #POR HACER: creo que esto viene del primer "if($function == 1)"


# #########################################################################
# 
# 					Get the right file names for options 1 and 2 
# 								
# 
# #########################################################################


if ($function == 2) { #POR HACER: probar que esto tambin funciona
	print STDERR "We are at function=2.... (aa_file is $aa_file)\n";
	$clustaloutfile = $aa_file;
	system "perl -pi -e 's/\\r\\n?/\\n/g' $clustaloutfile";
}


# #########################################################################
# 
# Insert dashes into the unaligned nucleotide file according to the aligned 
# 								amino acids
# 
# #########################################################################


open (CLUSTALOUTFILE, $clustaloutfile) || print_in_H_and_exit(\*H, "Unable to open clustaloutfile: $clustaloutfile\n"); # put nbrf file into hash

$/ = '>'; # alter record separator for fasta/nbrf file
my $toFasta = 0;
if(!-e $out_file.".aa_ali.fasta") {
	$toFasta = 1;
	open(ToF, ">$out_file".".aa_ali.fasta") || print_in_H_and_exit(\*H, "Unable to open the file $out_file.aa_ali.fasta\n"); 
}
my $count = 1;
while (<CLUSTALOUTFILE>) {
	#print STDERR $count++, "   -> $_\n";
	chomp;
	if ($_ =~ m/\w/i) {
		if ($_ =~ m/\w\w;/ and $_ =~ /\*/i) { # This means it is an nbrf/pir file
			($aa_name, $discardline, @aa_seq) = split "\n"; # name into $name, several seq lines into @Sequence
		}
		else { # its a FASTA file
			($aa_name, @aa_seq) = split "\n"; # name into $name, several seq lines into @Sequence
		}
		$aa_seq = join ('',@aa_seq); #now join up elements in @n_seq to create $n_seq.
		$aa_name =~s/ [0-9]+ bp ?//g;
		$aa_name =~s/[ \t\r]+$//g;
		$aa_name =~ s/\s/_/g;
		$aa_name =~ s/P1;//;
		$aa_name =~ s/[\(\)]//g; # not getting rid of underscores
		$aa_seq =~ s/\s//g;
		$aa_seq =~ s/\*//g; # get rid of * at end of NBRF file
		$aa_nameseq{$aa_name} = $aa_seq;
		push @aa_spnames, $aa_name;
		#print STDERR "$aa_name in aaspnames: @aa_spnames\n";
		if($toFasta) {	print ToF ">$aa_name\n$aa_seq\n";	}
		$discardline='';
	}
	else {
		next;
	}
}
if($toFasta)	{	close(ToF);	}

# #########################################################################
# 
# Check each nt seq is 3x length of equivalent aa seq
# 							
# 
# #########################################################################

foreach (@aa_spnames) {
	#print STDERR "aa_spname $_\n";
	$tempnt = $nt_nameseq{$_};
	$tempaa = $aa_nameseq{$_};
	$tempnt =~ s/[\-\.]//g;
	$tempaa =~ s/[\-\.]//g;
	$tempntlength = length $tempnt;
	$tempntlength = int ($tempntlength/3);
	$tempaalength = length $tempaa;
	#print STDERR "$tempntlength - $tempaalength\n";
	if ($tempntlength != $tempaalength) {
		print STDERR "Length ($tempaalength) of amino acid seq from $_ does not correspond to 1/3 of nucelotide seq ($tempntlength)\n";
		print H      "<font color=\"red\">Error: length ($tempaalength) of amino acid seq from $_ does not correspond to 1/3 of nucelotide seq ($tempntlength)</font><br>There might be some problem with the format of your sequences.<br><b>Exiting...</b><br>\n";
		exit;
	}
	else {
		next;
	}
	$tempnt =();
	$tempaa =();
	$tempntlength = ();
	$tempaalength = ();
}




# Now compare original NUC file with newly aligned AMINO ACID file and put dashes in as appropriate.

foreach (@nt_spnames) {
	if (! $aa_nameseq{$_}) {
		print STDERR "Warning:	Sequence $_ not found in amino acid file\n";
		print H      "<font color=\"red\">Warning: Sequence $_ not found in amino acid file</font>\n";
		my @tmp = keys %aa_nameseq;
		print STDERR "@tmp\n";
	}
}


my $codon_coloured_file = "$out_file.aa_based_codon_coloured.html";
open(CC, ">$codon_coloured_file") || print_in_H_and_exit(\*H, "Unable to create the codon_coloured_file: $codon_coloured_file\n"); 
print CC "<html>\n";
print CC <<HTML_OTRO;
	<head><style type='text/css'>
		TD, TR
		{
			font-size: 12px;
			font-family: Verdana, Arial, Helvetica, Helv, sans-serif;
		}
		body{
			font-size: 12px;
			font-family: Verdana, Arial, Helvetica, Helv, sans-serif;
			margin: 0;
			padding: 0;
			border: 0;
			overflow: hidden;
			height: 100%; 
			max-height: 100%; 
		}
		
		#framecontent{
			position: absolute;
			top: 0;
			bottom: 0; 
			left: 0;
			width: 100px; /*Width of frame div*/
			height: 100%;
			overflow: hidden; /*Disable scrollbars. Set to "scroll" to enable*/
			background: navy;
			color: white;
		}
		
		#maincontent{
			position: fixed;
			top: 0; 
			left: 100px; /*Set left value to WidthOfFrameDiv*/
			right: 0;
			bottom: 0;
			overflow: auto; 
			background: #fff;
		}
		
		.innertube{
			margin: 15px; /*Margins for inner DIV inside each DIV (to provide padding)*/
		}
		
		* html body{ /*IE6 hack*/
			padding: 0 0 0 100px; /*Set value to (0 0 0 WidthOfFrameDiv)*/
		}
		
		* html #maincontent{ /*IE6 hack*/
			height: 100%; 
			width: 100%; 
		}
	</style>
	<title>AA-coloured codon alignment</title>
	</head>
HTML_OTRO

print CC "<body>\n";
	$fcounter = 0;
		print CC "<div id=\"framecontent\"><div class=\"innertube\">\n";
		print CC "<table><tr>\n";
	foreach my $aa ( @aas ) { 
		print CC "<td bgcolor=\"$aa2colors{$aa}\">$aa</td>\n";
		$fcounter++;
		if($fcounter == 5) {
			print CC "</tr><tr>\n"; 
			$fcounter = 0;
		}
	}
print CC "</tr></table></div></div><div id=\"maincontent\"><div id=\"innertube\"><br>\n";
print CC "<table border=0><tr>\n";

$fastaoutfile     = $out_file.".nt_ali.fasta";
my $firstoutfile  = $out_file.".nt1_ali.fasta";
my $secondoutfile = $out_file.".nt2_ali.fasta";
my $thirdoutfile  = $out_file.".nt3_ali.fasta";
my $firstsecondoutfile  = $out_file.".nt12_ali.fasta";
open(F1,">$firstoutfile")  || die "Unable to write to firstoutfile: $firstoutfile\n";
open(F2,">$secondoutfile") || die "Unable to write to secondoutfile: $secondoutfile\n";
open(F3,">$thirdoutfile")  || die "Unable to write to thirdoutfile: $thirdoutfile\n";
open(F12,">$firstsecondoutfile")  || die "Unable to write to firstsecondoutfile: $firstsecondoutfile\n";
print STDERR "Ready to print results on $fastaoutfile!!!\n";
open (FASTAOUT, ">$fastaoutfile");
foreach (@aa_spnames) {
	print CC "</tr><tr><td>$_</td>\n";
	if ($nt_nameseq{$_}) {
		$aa_posn = 0;
		@aa_array = split(//, $aa_nameseq{$_});
		@nt_array = $nt_nameseq{$_} =~ /(.{3})/g;
		$aa_length = scalar(@aa_array);
		@finalnt_array = ();
		my(@firstnt_array,@secondnt_array,@thirdnt_array,@firstsecondnt_array);
			while  ($aa_posn < $aa_length) {
				if ($aa_array[$aa_posn] =~ /\-/) {
					my $current_aa = uc($aa_array[$aa_posn]);
					print CC "<td bgcolor=\"$aa2colors{$current_aa}\">---</td>\n";
					push @finalnt_array, "---";
					push @firstnt_array,  "-";
					push @secondnt_array, "-";
					push @thirdnt_array,  "-";
					push @firstsecondnt_array, "--";
					$aa_posn++;
				} elsif ($aa_array[$aa_posn] =~ /[A-Za-z]/) {
					my $current_aa = uc($aa_array[$aa_posn]);
					$shifted_nts = shift @nt_array;
					print CC "<td bgcolor=\"$aa2colors{$current_aa}\">$shifted_nts</td>\n";
					$temptriplet = (uc $shifted_nts);
					my @tmp = split(//,$temptriplet);
					$gccontent_sp_1st{$_}->{$tmp[0]}++;
					$gccontent_sp_2nd{$_}->{$tmp[1]}++;
					$gccontent_sp_3rd{$_}->{$tmp[2]}++;
					$gccontent_sp{$_}    ->{$tmp[0]}++;
					$gccontent_sp{$_}    ->{$tmp[1]}++;
					$gccontent_sp{$_}    ->{$tmp[2]}++;
					$gccontent_1st{$tmp[0]}         ++;
					$gccontent_2nd{$tmp[1]}         ++;
					$gccontent_3rd{$tmp[2]}         ++;
					$gccontent{$tmp[0]}             ++;
					$gccontent{$tmp[1]}             ++;
					$gccontent{$tmp[2]}             ++;
					
					push @finalnt_array, $shifted_nts;
					push @firstnt_array,  $tmp[0];
					push @secondnt_array, $tmp[1];
					push @thirdnt_array,  $tmp[2];
					push @firstsecondnt_array, (@tmp)[0,1];
					$aa_posn++;
				} else {
					$aa_posn++;
				}
			}
			
		print FASTAOUT ">$_\n",@finalnt_array,"\n";
		print F1  ">$_\n",@firstnt_array,"\n";
		print F2  ">$_\n",@secondnt_array,"\n";
		print F3  ">$_\n",@thirdnt_array,"\n";
		print F12 ">$_\n",@firstsecondnt_array, "\n";
		#delete $nt_nameseq{$_};
	}

	else {
		print STDERR "Warning:	Sequence $_ not found in nucleotide file\n";
		print H      "<font color=\"red\">Warning: Sequence $_ not found in nucleotide file</font>\n";
	}
	
}
close FASTAOUT;
close(F1);	close(F2);	close(F3);	close(F12);

print CC "</tr></table></div></div></body></html>\n";
close(CC);



my $nt_ali_file         = "$out_file.nt_ali.fasta";
my $aa_ali_file         = "$out_file.aa_ali.fasta";

#print H "<div id=\"rounded-box2\" style=\"z-index: 1; padding: 5px; padding-left:40px; padding-bottom:14px; padding-top:14px; background-color:#CCAFAF; layer-background-color:#003366; \">\n";
print H "<table border=0 cellspacing=10 width=\"95%\"><tr valign=\"top\">\n";
print H "<td  width=\"50%\" style=\"border-radius: 7px; -ms-border-radius: 7px; -moz-border-radius:7px; -webkit-border-radius: 7px; -khtml-border-radius: 7px;z-index: 1; padding: 5px;  padding-left:40px; padding-right:15px;  padding-bottom:14px; padding-top:14px; height:90%; background-color:#CCAFAF;; layer-background-color:#CCCC88; \">\n";

print H "<H2>Nucleotide Alignment</H2>\n";
if($from_web_server) {	print H "<div style=\"display:none;\">\n";	}
print H "<textarea id=\"nt_textarea\" style=\"font-size=10px;\" ROWS=12 COLS=50 wrap=\"off\">\n";
	open(I, $nt_ali_file) || print_in_H_and_exit(\*H, "Unable to open the nt_ali_file: $nt_ali_file\n"); 
	while(<I>) {	print H "$_";	}
	close(I);
print H "</textarea><br>\n";
if($from_web_server) {	print H "</div>\n";	}


if($from_web_server) {
	my $tmp_file = $nt_ali_file;
	$tmp_file = (split(/\//,$tmp_file))[-1];
		print H	"<applet code=\"jalview.bin.JalviewLite\" archive=\"../jalviewApplet.jar\" width=\"90%\" height=\"200\"\n";
		print H "                   <param name=\"embedded\" value=\"true\">\n";
		print H	"                   <param name=\"file\" value=\"$tmp_file\">\n";
		print H	"                   <param name=\"defaultColour\" value=\"Clustal\">\n";
		print H	"                   <param name=\"label\"  value=\"View larger...\">\n";
		print H	"                   <param name=\"RGB\"  value=\"8899AA\">\n";
		print H	"                   <param name=\"showQuality\" value=\"false\">\n";
		print H	"                   <param name=\"showConsensus\" value=\"false\">\n";
		print H	"                   <param name=\"showConservation\" value=\"false\">\n";
		print H	"</applet><br>\n";

	print H "<table width=\"300\" padding=10><tr>\n";
	print H "<td align=\"left\"><input type=\"submit\" style=\"font-size:13px;\" value=\"Save fasta...\" onclick=\"send2file_webserver('$tmp_file');\">\n";
	print H "<img align=\"Absbottom\" src=\"/question.png\" width=\"28\" alt=\"Use Jalview &gt; File &gt; Output to Textbox to save in an alternative format \" title=\"Use Jalview &gt; File &gt; Output to Textbox to save in an alternative format\" />\n";
	print H "</td>\n";
	print H	"<td align=\"right\"><applet code=\"jalview.bin.JalviewLite\" archive=\"../jalviewApplet.jar\" width=\"140\" height=\"35\"\n";
	#print H "                   <param name=\"embedded\" value=\"true\">\n";
    print H	"                   <param name=\"file\" value=\"$tmp_file\">\n";
    print H	"                   <param name=\"defaultColour\" value=\"Clustal\">\n";
    print H	"                   <param name=\"label\"  value=\"View larger...\">\n";
    print H	"                   <param name=\"RGB\"  value=\"CCAFAF\">\n";
    print H	"                   <param name=\"showQuality\" value=\"false\">\n";
    print H	"                   <param name=\"showConsensus\" value=\"false\">\n";
	print H	"                   <param name=\"showConservation\" value=\"false\">\n";
    print H	"</applet></td>\n";
    print H "</tr></table>\n";
} else {
	print H "<input type=\"submit\" value=\"Save fasta...\" onclick=\"send2file('nt_textarea');\"> \n";
}
my $tmp_file = $out_file;
$tmp_file = (split(/\//,$tmp_file))[-1];
#print H "<ul><li><a href=\"$tmp_file.aa_based_codon_coloured.html\" onclick=\"window.open('$tmp_file.aa_based_codon_coloured.html','aa-based codon-coloured alignment','width=640,height=300,scrollbars=1')\">View aa-based codon-coloured alignment [html]</a></li></ul>\n";
print H "<ul><li><a href=\"\" onclick=\"window.open('$tmp_file.aa_based_codon_coloured.html','aa-based codon-coloured alignment','width=640,height=300,scrollbars=1')\">View aa-based codon-coloured alignment [html]</a></li></ul>\n";



$/ = $old_input_rec_sep;
my $jorl = $program;
$jorl = "user" if(!$jorl);
print H "</td><td  width=\"50%\" style=\"border-radius: 7px; -ms-border-radius: 7px; -moz-border-radius:7px; -webkit-border-radius: 7px; -khtml-border-radius: 7px;z-index: 1; padding: 5px;  padding-left:40px; padding-right:15px;  padding-bottom:14px; padding-top:14px; height:90%; background-color:#86AAD4;; layer-background-color:#CCCC88; \">\n";
print H "<H2>Amino acid alignment ($jorl)</H2>\n";
if($from_web_server) {	print H "<div style=\"display:none;\">\n";	}
print H "<textarea id=\"aa_textarea\" style=\"font-size=10px;\" width=\"90%\" ROWS=12 COLS=50 wrap=\"off\">\n";
	open(I, $aa_ali_file) || print_in_H_and_exit(\*H, "Unable to open the aa_ali_file: $aa_ali_file\n"); 
	$_ = <I>; print H "$_";
	while(<I>) {	if(/^>/){ print H "\n$_";	} else { chomp; print H "$_";} }
	close(I);
print H "</textarea><br>\n";
if($from_web_server) {	print H "</div>\n";	}


if($from_web_server) {
	my $tmp_file = $aa_ali_file;
	$tmp_file = (split(/\//,$tmp_file))[-1];
		print H	"<applet code=\"jalview.bin.JalviewLite\" archive=\"../jalviewApplet.jar\" width=\"90%\" height=\"200\"\n";
		print H "                   <param name=\"embedded\" value=\"true\">\n";
		print H	"                   <param name=\"file\" value=\"$tmp_file\">\n";
		print H	"                   <param name=\"defaultColour\" value=\"Clustal\">\n";
		print H	"                   <param name=\"label\"  value=\"View larger...\">\n";
		print H	"                   <param name=\"RGB\"  value=\"8899AA\">\n";
		print H	"                   <param name=\"showQuality\" value=\"false\">\n";
		print H	"                   <param name=\"showConsensus\" value=\"false\">\n";
		print H	"                   <param name=\"showConservation\" value=\"false\">\n";
		print H	"</applet><br>\n";

	print H "<table width=\"300\" padding=10><tr valign=\"bottom\">\n";
	print H "<td align=\"left\"><input type=\"submit\" style=\"font-size:13px;\" value=\"Save fasta...\" onclick=\"send2file_webserver('$tmp_file');\">\n";
	print H "<img align=\"Absbottom\" src=\"/question.png\" width=\"28\" alt=\"Use Jalview &gt; File &gt; Output to Textbox to save in an alternative format \" title=\"Use Jalview &gt; File &gt; Output to Textbox to save in an alternative format\" />\n";
	print H "</td>\n";
	print H	"<td align=\"right\"><applet code=\"jalview.bin.JalviewLite\" archive=\"../jalviewApplet.jar\" width=\"140\" height=\"35\"\n";
	#print H "                   <param name=\"embedded\" value=\"true\">\n";
    print H	"                   <param name=\"file\" value=\"$tmp_file\">\n";
    print H	"                   <param name=\"defaultColour\" value=\"Clustal\">\n";
    print H	"                   <param name=\"label\"  value=\"View larger...\">\n";
    print H	"                   <param name=\"RGB\"  value=\"86AAD4\">\n";
    print H	"                   <param name=\"showQuality\" value=\"false\">\n";
    print H	"                   <param name=\"showConsensus\" value=\"false\">\n";
	print H	"                   <param name=\"showConservation\" value=\"false\">\n";
    print H	"</applet></td>\n";
    print H "</tr></table>\n";
} else {
	print H "<input type=\"submit\" value=\"Save fasta...\" onclick=\"send2file('aa_textarea');\"> \n";
}
print H "</td></tr></table>\n";
#print H "</div>\n";




if($gblocks_params ne "") {
	#GBlocks results...
	#print H "<div id=\"rounded-box2\" style=\"z-index: 1; padding: 5px; padding-left:40px; padding-bottom:14px; padding-top:14px; background-color:#8899AA; layer-background-color:#003366; \">\n";
	
	#Ahora ejecutamos lo de gblocks....
	#Hay que pasarle el $aa_ali_file y los parmetros $gblocks_params
	system "$path"."Gblocks $aa_ali_file $gblocks_params";
	system "$path"."Gblocks $aa_ali_file $gblocks_params -p=s";
	#system "$path"."Gblocks $aa_ali_file $gblocks_params -p=y";
	
	my $tmp_file = $aa_ali_file;
	$tmp_file = (split(/\//,$tmp_file))[-1];

	my $gb_aa_ali_file = $aa_ali_file."-gb";
	print STDERR "Going to rename $gb_aa_ali_file to $out_file.aa_cleanali.fasta\n\n";
	rename($gb_aa_ali_file, $out_file.".aa_cleanali.fasta");
	$gb_aa_ali_file = $out_file.".aa_cleanali.fasta";
	print STDERR "Ready to remove spaces....";
	system "perl -pi -e 's/ //g' $gb_aa_ali_file";
	print STDERR "...done\n";
	
	my $gb_coordinates_file = $aa_ali_file."-gb.txts";
	my $ali_length;
	
	$/ = $old_input_rec_sep;
	open(I, $aa_ali_file) || print_in_H_and_exit(\*H, "Unable unable\n"); 
		my $seq = "";
		while(<I>) {
			if(/^>/) {
				if($seq ne "") {
					$ali_length = length($seq);
					last;
				}
			} else {
				s/[ \t\n]//g;
				$seq .= $_;
			}
		}
	close(I);
	
	$/ = $old_input_rec_sep;
	open(I, $gb_coordinates_file) || print_in_H_and_exit(\*H, "Unable to open the gb_coordinates file: $gb_coordinates_file\n");
	my @tmp;
	my $count = 0;
	while(<I>) {
		if(/Flank positions/) {
			$_ = <I>;
			s/^Flanks:[ \t]+//;
			s/[ \t\n\[\]]+/ /g;
			s/^ +//;
			@tmp = split(/ /,$_);
			last;
			last;
		}
	}
	close(I);
	my(@discarded);
	for(my $i=0; $i<$ali_length; $i++) {
		$discarded[$i] = 1;
	}
	for(my $i=0; $i<scalar(@tmp); $i=$i+2) {
		my $start = $tmp[$i];
		my $end   = $tmp[$i+1];
		for(my $j=$start-1; $j<$end; $j++) {
			$discarded[$j] = 0;
		}
	}
	

	my $codon_coloured_file = "$out_file.aa_based_codon_coloured-gb.html";
	open(CC, ">$codon_coloured_file") || print_in_H_and_exit(\*H, "Unable to create the codon_coloured-gb_file: $codon_coloured_file\n");
	print CC "<html>\n";
	print CC <<HTML_OTRO;
	<head><style type='text/css'>
		TD, TR
		{
			font-size: 12px;
			font-family: Verdana, Arial, Helvetica, Helv, sans-serif;
		}
		body{
			font-size: 12px;
			font-family: Verdana, Arial, Helvetica, Helv, sans-serif;
			margin: 0;
			padding: 0;
			border: 0;
			overflow: hidden;
			height: 100%; 
			max-height: 100%; 
		}
		
		#framecontent{
			position: absolute;
			top: 0;
			bottom: 0; 
			left: 0;
			width: 100px; /*Width of frame div*/
			height: 100%;
			overflow: hidden; /*Disable scrollbars. Set to "scroll" to enable*/
			background: navy;
			color: white;
		}
		
		#maincontent{
			position: fixed;
			top: 0; 
			left: 100px; /*Set left value to WidthOfFrameDiv*/
			right: 0;
			bottom: 0;
			overflow: auto; 
			background: #fff;
		}
		
		.innertube{
			margin: 15px; /*Margins for inner DIV inside each DIV (to provide padding)*/
		}
		
		* html body{ /*IE6 hack*/
			padding: 0 0 0 100px; /*Set value to (0 0 0 WidthOfFrameDiv)*/
		}
		
		* html #maincontent{ /*IE6 hack*/
			height: 100%; 
			width: 100%; 
		}
	</style>
	<title>AA-coloured codon alignment (gblocks)</title>
	</head>
HTML_OTRO
	
	print CC "<body>\n";
		my $fcounter = 0;
		print CC "<div id=\"framecontent\"><div class=\"innertube\">\n";
		print CC "<table><tr>\n";
		foreach my $aa ( @aas ) { 
			print CC "<td bgcolor=\"$aa2colors{$aa}\">$aa</td>\n";
			$fcounter++;
			if($fcounter == 5) {
				print CC "</tr><tr>\n"; 
				$fcounter = 0;
			}
		}
	print CC "</tr></table></div></div><div id=\"maincontent\"><div id=\"innertube\"><br>\n";
	print CC "<table border=0><tr>\n";
	
	$fastaoutfile = $out_file.".nt_cleanali.fasta";
	$firstoutfile  = $out_file.".nt1_cleanali.fasta";
	$secondoutfile = $out_file.".nt2_cleanali.fasta";
	$thirdoutfile  = $out_file.".nt3_cleanali.fasta";
	$firstsecondoutfile = $out_file.".nt12_cleanali.fasta";
	open(F1, ">$firstoutfile")        || die "Unable to write to firstoutfile  (clean): $firstoutfile\n";
	open(F2, ">$secondoutfile")       || die "Unable to write to secondoutfile (clean): $secondoutfile\n";
	open(F3, ">$thirdoutfile")        || die "Unable to write to thirdoutfile  (clean): $thirdoutfile\n";
	open(F12,">$firstsecondoutfile")  || die "Unable to write to firstsecondoutfile  (clean): $firstsecondoutfile\n";
	print STDERR "Ready to print results on $fastaoutfile!!!\n";
	open (FASTAOUT, ">$fastaoutfile") || print_in_H_and_exit(\*H, "Unable to write to fastaoutfile: $fastaoutfile\n");
	my @clean_nt_array;
	foreach (@aa_spnames) {
		print CC "</tr><tr><td>$_</td>\n";
		if ($nt_nameseq{$_}) {
			$aa_posn = 0;
			@aa_array = split(//, $aa_nameseq{$_});
			@nt_array = $nt_nameseq{$_} =~ /(.{3})/g;
			$aa_length = scalar(@aa_array);
			@clean_nt_array = ();
			my(@firstnt_array,@secondnt_array,@thirdnt_array,@firstsecondnt_array);
				while  ($aa_posn < $aa_length) {
					#print STDERR "$aa_posn($_): $aa_array[$aa_posn]\n";
					if($discarded[$aa_posn] == 1) {
						if($aa_array[$aa_posn] =~ /[A-Za-z]/) {
							$shifted_nts = shift @nt_array;
						}
						$aa_posn++;
						next;
					}
					if ($aa_array[$aa_posn] =~ /\-/) {
						push @clean_nt_array, "---";
						push @firstnt_array,  "-";
						push @secondnt_array, "-";
						push @thirdnt_array,  "-";
						push @firstsecondnt_array, "--";
						my $current_aa = uc($aa_array[$aa_posn]);
						print CC "<td bgcolor=\"$aa2colors{$current_aa}\">---</td>\n";
						$aa_posn++;
					} elsif ($aa_array[$aa_posn] =~ /[A-Za-z]/) {
						$shifted_nts = shift @nt_array;
						my $current_aa = uc($aa_array[$aa_posn]);
						print CC "<td bgcolor=\"$aa2colors{$current_aa}\">$shifted_nts</td>\n";
						my $temptriplet = (uc $shifted_nts);
						my @tmp = split(//,$temptriplet);
						#my @tmp = split(//,$shifted_nts);
						$clean_gccontent_sp_1st{$_}->{$tmp[0]}++;
						$clean_gccontent_sp_2nd{$_}->{$tmp[1]}++;
						$clean_gccontent_sp_3rd{$_}->{$tmp[2]}++;
						$clean_gccontent_sp{$_}    ->{$tmp[0]}++;
						$clean_gccontent_sp{$_}    ->{$tmp[1]}++;
						$clean_gccontent_sp{$_}    ->{$tmp[2]}++;
						$clean_gccontent_1st{$tmp[0]}         ++;
						$clean_gccontent_2nd{$tmp[1]}         ++;
						$clean_gccontent_3rd{$tmp[2]}         ++;
						$clean_gccontent{$tmp[0]}             ++;
						$clean_gccontent{$tmp[1]}             ++;
						$clean_gccontent{$tmp[2]}             ++;
						push @clean_nt_array, $shifted_nts;
						push @firstnt_array,  $tmp[0];
						push @secondnt_array, $tmp[1];
						push @thirdnt_array,  $tmp[2];
						push @firstsecondnt_array, (@tmp)[0,1];
						$aa_posn++;
					} else {
						$aa_posn++;
					}
				}

			#delete $nt_nameseq{$_};
			print FASTAOUT ">$_\n",@clean_nt_array,"\n";
			print F1  ">$_\n",@firstnt_array,"\n";
			print F2  ">$_\n",@secondnt_array,"\n";
			print F3  ">$_\n",@thirdnt_array,"\n";
			print F12 ">$_\n",@firstsecondnt_array,"\n";
		}
		else {
			print STDERR "Warning:	Sequence $_ not found in nucleotide file\n";
			print H      "<font color=\"red\">Warning: Sequence $_ not found in nucleotide file</font>\n";
		}
	}
	close FASTAOUT;
	close(F1);	close(F2);	close(F3);	close(F12);
	
	print CC "</tr></table></div></div></body></html>\n";
	close(CC);
	
	

			my $nttmp = $out_file.".nt_cleanali.fasta";
			my $aatmp = $out_file.".aa_cleanali.fasta";

			print H "<table border=0 width=\"95%\" cellspacing=10 ><tr>\n";
			print H "<td width=\"50%\" valign=\"top\"  style=\"border-radius: 7px; -ms-border-radius: 7px; -moz-border-radius:7px; -webkit-border-radius: 7px; -khtml-border-radius: 7px;z-index: 1; padding: 5px;  padding-left:40px; padding-right:15px;  padding-bottom:14px; padding-top:14px; background-color:#8899AA;; layer-background-color:#CCCC88; \">\n";
			
			print H "<H2>Clean nucleotide alignment (Gblocks)</H2>\n";
			if($from_web_server) {	print H "<div style=\"display:none;\">\n";	}
			print H "<textarea id=\"clean_nt_textarea\" style=\"font-size=10px;\" ROWS=12 COLS=50 wrap=\"off\">\n";
				open(I, $nttmp) || print_in_H_and_exit(\*H, "Unable to open the nttmp: $nttmp\n");
				while(<I>) {	print H "$_";	}
				close(I);
			print H "</textarea><br>\n";
			if($from_web_server) {	print H "</div>\n";	}

			if($from_web_server) {
				my $tmp_file = $nttmp;
				$tmp_file = (split(/\//,$tmp_file))[-1];
					print H	"<applet code=\"jalview.bin.JalviewLite\" archive=\"../jalviewApplet.jar\" width=\"90%\" height=\"200\"\n";
					print H "                   <param name=\"embedded\" value=\"true\">\n";
					print H	"                   <param name=\"file\" value=\"$tmp_file\">\n";
					print H	"                   <param name=\"defaultColour\" value=\"Clustal\">\n";
					print H	"                   <param name=\"label\"  value=\"View larger...\">\n";
					print H	"                   <param name=\"RGB\"  value=\"8899AA\">\n";
					print H	"                   <param name=\"showQuality\" value=\"false\">\n";
					print H	"                   <param name=\"showConsensus\" value=\"false\">\n";
					print H	"                   <param name=\"showConservation\" value=\"false\">\n";
					print H	"</applet><br>\n";

				print H "<table width=\"300\" padding=10><tr valign=\"bottom\">\n";
				print H "<td align=\"left\"><input type=\"submit\" style=\"font-size:13px;\" value=\"Save fasta...\" onclick=\"send2file_webserver('$tmp_file');\">\n";
				print H "<img align=\"Absbottom\" src=\"/question.png\" width=\"28\" alt=\"Use Jalview &gt; File &gt; Output to Textbox to save in an alternative format \" title=\"Use Jalview &gt; File &gt; Output to Textbox to save in an alternative format\" />\n";
				print H "</td>\n";
				print H	"<td align=\"right\"><applet code=\"jalview.bin.JalviewLite\" archive=\"../jalviewApplet.jar\" width=\"140\" height=\"35\"\n";
				#print H "                   <param name=\"embedded\" value=\"true\">\n";
				print H	"                   <param name=\"file\" value=\"$tmp_file\">\n";
				print H	"                   <param name=\"defaultColour\" value=\"Clustal\">\n";
				print H	"                   <param name=\"label\"  value=\"View larger...\">\n";
				print H	"                   <param name=\"RGB\"  value=\"8899AA\">\n";
				print H	"                   <param name=\"showQuality\" value=\"false\">\n";
				print H	"                   <param name=\"showConsensus\" value=\"false\">\n";
				print H	"                   <param name=\"showConservation\" value=\"false\">\n";
				print H	"</applet></td>\n";
				print H "</tr></table>\n";
			} else {
				print H "<input type=\"submit\" value=\"Save fasta...\" onclick=\"send2file('clean_nt_textarea');\"> \n";
			}

			$tmp_file = $out_file;
			$tmp_file = (split(/\//,$tmp_file))[-1];
			#print H "<ul><li><a href=\"$tmp_file.aa_based_codon_coloured-gb.html\" onclick=\"window.open('$tmp_file.aa_based_codon_coloured-b.html','aa-based codon-coloured alignment','width=640,height=300,scrollbars=1')\">View aa-based codon-coloured alignment [html]</a></li></ul>\n";
			print H "<ul><li><a href=\"\" onclick=\"window.open('$tmp_file.aa_based_codon_coloured-gb.html','aa-based codon-coloured gblocks alignment','width=640,height=300,scrollbars=1')\">View aa-based codon-coloured alignment [html]</a></li></ul>\n";

			
			print H "</td><td width=\"50%\" valign=\"top\" style=\"border-radius: 7px; -ms-border-radius: 7px; -moz-border-radius:7px; -webkit-border-radius: 7px; -khtml-border-radius: 7px;z-index: 1; padding: 5px;  padding-left:40px; padding-right:15px;  padding-bottom:14px; padding-top:14px; height:90%; background-color:#83C2C0;; layer-background-color:#CCCC88; \">\n";
			print H "<H2>Clean amino acid alignment (Gblocks)</H2>\n";
			if($from_web_server) {	print H "<div style=\"display:none;\">\n";	}
			print H "<textarea id=\"clean_aa_textarea\" style=\"font-size=10px;\" width=\"90%\" ROWS=12 COLS=50 wrap=\"off\">\n";
				open(I, $aatmp) || print_in_H_and_exit(\*H, "Unable to open the aatmp: $aatmp\n");
				$_ = <I>; print H "$_";
				while(<I>) {	if(/^>/){ print H "\n$_";	} else { chomp; print H "$_";} }
				close(I);
			print H "</textarea><br>\n";
			if($from_web_server) {	print H "</div>\n";	}
			
			if($from_web_server) {
				my $tmp_file = $aatmp;
				$tmp_file = (split(/\//,$tmp_file))[-1];

					print H	"<applet code=\"jalview.bin.JalviewLite\" archive=\"../jalviewApplet.jar\" width=\"90%\" height=\"200\"\n";
					print H "                   <param name=\"embedded\" value=\"true\">\n";
					print H	"                   <param name=\"file\" value=\"$tmp_file\">\n";
					print H	"                   <param name=\"defaultColour\" value=\"Clustal\">\n";
					print H	"                   <param name=\"label\"  value=\"View larger...\">\n";
					print H	"                   <param name=\"RGB\"  value=\"8899AA\">\n";
					print H	"                   <param name=\"showQuality\" value=\"false\">\n";
					print H	"                   <param name=\"showConsensus\" value=\"false\">\n";
					print H	"                   <param name=\"showConservation\" value=\"false\">\n";
					print H	"</applet><br>\n";

				print H "<table width=\"300\" padding=10><tr valign=\"bottom\">\n";
				print H "<td align=\"left\"><input type=\"submit\" style=\"font-size:13px;\" value=\"Save fasta...\" onclick=\"send2file_webserver('$tmp_file');\">\n";
				print H "<img align=\"Absbottom\" src=\"/question.png\" width=\"28\" alt=\"Use Jalview &gt; File &gt; Output to Textbox to save in an alternative format \" title=\"Use Jalview &gt; File &gt; Output to Textbox to save in an alternative format\" />\n";
				print H "</td>\n";
				print H	"<td align=\"right\"><applet code=\"jalview.bin.JalviewLite\" archive=\"../jalviewApplet.jar\" width=\"140\" height=\"35\"\n";
				#print H "                   <param name=\"embedded\" value=\"true\">\n";
				print H	"                   <param name=\"file\" value=\"$tmp_file\">\n";
				print H	"                   <param name=\"defaultColour\" value=\"Clustal\">\n";
				print H	"                   <param name=\"label\"  value=\"View larger...\">\n";
				print H	"                   <param name=\"RGB\"  value=\"83C2C0\">\n";
				print H	"                   <param name=\"showQuality\" value=\"false\">\n";
				print H	"                   <param name=\"showConsensus\" value=\"false\">\n";
				print H	"                   <param name=\"showConservation\" value=\"false\">\n";
				print H	"</applet></td>\n";
				print H "</tr></table>\n";
			} else {
				print H "<input type=\"submit\" value=\"Save fasta...\" onclick=\"send2file('clean_aa_textarea');\"> \n";
			}
			$tmp_file = $out_file.".aa_ali.fasta";
			$tmp_file = (split(/\//,$tmp_file))[-1];
			print H "<ul><li><a href=\"$tmp_file-gb.htm\">GBlocks Results: graphical representation in html</a></li>\n";
			print H "<li><a href=\"$tmp_file-gb.txts\">Blocks coordinates</a></li></ul>\n";
			print H "</td></tr></table>\n";

	

	
	#print H "</div>\n";
}


#-----------------------------------------------
#-----  GC content     -------------------------
#-----------------------------------------------
#print H "<div id=\"rounded-box2\" style=\"z-index: 1; padding: 5px; padding-left:40px; padding-bottom:14px; padding-top:14px; background-color:#86AAD4; layer-background-color:#003366; \">\n";


my $sum = 0;
foreach my $nt ( keys %gccontent ) {	$sum += $gccontent{$nt};	}
$gccontent{'GC'} = sprintf("%.2f",($gccontent{'G'}+$gccontent{'C'})/$sum*100);

$sum = 0;
foreach my $nt ( keys %gccontent_1st ) {	$sum += $gccontent_1st{$nt};	}
$gccontent_1st{'GC'} = sprintf("%.2f",($gccontent_1st{'G'}+$gccontent_1st{'C'})/$sum*100);

$sum = 0;
foreach my $nt ( keys %gccontent_2nd ) {	$sum += $gccontent_2nd{$nt};	}
$gccontent_2nd{'GC'} = sprintf("%.2f",($gccontent_2nd{'G'}+$gccontent_2nd{'C'})/$sum*100);

$sum = 0;
foreach my $nt ( keys %gccontent_3rd ) {	$sum += $gccontent_3rd{$nt};	}
$gccontent_3rd{'GC'} = sprintf("%.2f",($gccontent_3rd{'G'}+$gccontent_3rd{'C'})/$sum*100);


foreach (@aa_spnames) {
	$sum = 0;
	foreach my $nt ( keys %{$gccontent_sp{$_}} ) {	$sum += $gccontent_sp{$_}->{$nt};	}
	$gccontent_sp{$_}->    {'GC'} = sprintf("%.2f",($gccontent_sp{$_}->    {'G'}+$gccontent_sp{$_}->    {'C'})/$sum*100);
	
	$sum = 0;
	foreach my $nt ( keys %{$gccontent_sp_1st{$_}} ) {	$sum += $gccontent_sp_1st{$_}->{$nt};	}
	$gccontent_sp_1st{$_}->{'GC'} = sprintf("%.2f",($gccontent_sp_1st{$_}->{'G'}+$gccontent_sp_1st{$_}->{'C'})/$sum*100);
	
	$sum = 0;
	foreach my $nt ( keys %{$gccontent_sp_2nd{$_}} ) {	$sum += $gccontent_sp_2nd{$_}->{$nt};	}
	$gccontent_sp_2nd{$_}->{'GC'} = sprintf("%.2f",($gccontent_sp_2nd{$_}->{'G'}+$gccontent_sp_2nd{$_}->{'C'})/$sum*100);
	
	$sum = 0;
	foreach my $nt ( keys %{$gccontent_sp_3rd{$_}} ) {	$sum += $gccontent_sp_3rd{$_}->{$nt};	}
	$gccontent_sp_3rd{$_}->{'GC'} = sprintf("%.2f",($gccontent_sp_3rd{$_}->{'G'}+$gccontent_sp_3rd{$_}->{'C'})/$sum*100);
}


if($gblocks_params ne "") {
	$sum = 0;
	foreach my $nt ( keys %clean_gccontent ) {	$sum += $clean_gccontent{$nt};	}
	$clean_gccontent{'GC'} = sprintf("%.2f",($clean_gccontent{'G'}+$clean_gccontent{'C'})/$sum*100);
	
	$sum = 0;
	foreach my $nt ( keys %clean_gccontent_1st ) {	$sum += $clean_gccontent_1st{$nt};	}
	$clean_gccontent_1st{'GC'} = sprintf("%.2f",($clean_gccontent_1st{'G'}+$clean_gccontent_1st{'C'})/$sum*100);
	
	$sum = 0;
	foreach my $nt ( keys %clean_gccontent_2nd ) {	$sum += $clean_gccontent_2nd{$nt};	}
	$clean_gccontent_2nd{'GC'} = sprintf("%.2f",($clean_gccontent_2nd{'G'}+$clean_gccontent_2nd{'C'})/$sum*100);
	
	$sum = 0;
	foreach my $nt ( keys %clean_gccontent_3rd ) {	$sum += $clean_gccontent_3rd{$nt};	}
	$clean_gccontent_3rd{'GC'} = sprintf("%.2f",($clean_gccontent_3rd{'G'}+$clean_gccontent_3rd{'C'})/$sum*100);
	
	
	foreach (@aa_spnames) {
		$sum = 0;
		foreach my $nt ( keys %{$clean_gccontent_sp{$_}} ) {	$sum += $clean_gccontent_sp{$_}->{$nt};	}
		$clean_gccontent_sp{$_}->    {'GC'} = sprintf("%.2f",($clean_gccontent_sp{$_}->    {'G'}+$clean_gccontent_sp{$_}->    {'C'})/$sum*100);
		
		$sum = 0;
		foreach my $nt ( keys %{$clean_gccontent_sp_1st{$_}} ) {	$sum += $clean_gccontent_sp_1st{$_}->{$nt};	}
		$clean_gccontent_sp_1st{$_}->{'GC'} = sprintf("%.2f",($clean_gccontent_sp_1st{$_}->{'G'}+$clean_gccontent_sp_1st{$_}->{'C'})/$sum*100);
		
		$sum = 0;
		foreach my $nt ( keys %{$clean_gccontent_sp_2nd{$_}} ) {	$sum += $clean_gccontent_sp_2nd{$_}->{$nt};	}
		$clean_gccontent_sp_2nd{$_}->{'GC'} = sprintf("%.2f",($clean_gccontent_sp_2nd{$_}->{'G'}+$clean_gccontent_sp_2nd{$_}->{'C'})/$sum*100);
		
		$sum = 0;
		foreach my $nt ( keys %{$clean_gccontent_sp_3rd{$_}} ) {	$sum += $clean_gccontent_sp_3rd{$_}->{$nt};	}
		$clean_gccontent_sp_3rd{$_}->{'GC'} = sprintf("%.2f",($clean_gccontent_sp_3rd{$_}->{'G'}+$clean_gccontent_sp_3rd{$_}->{'C'})/$sum*100);
	}

}
my $cols = 1;
$cols = 2 if($gblocks_params ne "");
print H "<table border=0 width=\"95%\" cellspacing=10 ><tr><td  colspan=\"$cols\" style=\"border-radius: 7px; -ms-border-radius: 7px; -moz-border-radius:7px; -webkit-border-radius: 7px; -khtml-border-radius: 7px;z-index: 1; padding: 5px;  padding-left:40px; padding-right:15px;  padding-bottom:14px; padding-top:14px; height:90%; background-color:#6CB381;; layer-background-color:#CCCC88; \">\n";
#print H "<table border=0 width=\"90%\"><tr><td colspan=\"$cols\">\n";

print H "<H2>Compositional bias (GC content) and 1st, 2nd and 3rd position nt-alignments</H2>\n";
print H "</td></tr><tr><td style=\"border-radius: 7px; -ms-border-radius: 7px; -moz-border-radius:7px; -webkit-border-radius: 7px; -khtml-border-radius: 7px;z-index: 1; padding: 5px;  padding-left:40px; padding-right:15px;  padding-bottom:14px; padding-top:14px; height:90%; background-color:#7CC391;; layer-background-color:#CCCC88; \"><h3>Original data</h3>\n" if($gblocks_params ne "");
print H "<table bgcolor=\"B2CFD4\" Border=0 CELLSPACING=\"6\" style=\"border : thin solid; \">\n";
print H "<tr><th>Taxa</th><th>GC</th><th>GC 1st</th><th>GC 2nd</th><th>GC 3rd</th></tr>\n";
foreach ( sort { $gccontent_sp{$a}->{'GC'} <=> $gccontent_sp{$b}->{'GC'} } keys %gccontent_sp ) {
	print H "<tr><td><i>$_</i></td><td>$gccontent_sp{$_}->{'GC'}</td><td>$gccontent_sp_1st{$_}->{'GC'}</td><td>$gccontent_sp_2nd{$_}->{'GC'}</td><td>$gccontent_sp_3rd{$_}->{'GC'}</td></tr>\n";
}
#foreach (@aa_spnames) {
#	print H "<tr><td><i>$_</i></td><td>$gccontent_sp{$_}->{'GC'}</td><td>$gccontent_sp_1st{$_}->{'GC'}</td><td>$gccontent_sp_2nd{$_}->{'GC'}</td><td>$gccontent_sp_3rd{$_}->{'GC'}</td></tr>\n";
#}
print H "<tr><td><b>Overall</b></td><td><b>$gccontent{'GC'}</b></td><td><b>$gccontent_1st{'GC'}</b></td><td><b>$gccontent_2nd{'GC'}</b></td><td><b>$gccontent_3rd{'GC'}</b></td></tr>\n";
print H "</table>\n";

$tmp_file = $out_file;
$tmp_file = (split(/\//,$tmp_file))[-1];
print H "<table border=0><tr><td>1st position nt alignment</td> \n";
print H "<td><applet code=\"jalview.bin.JalviewLite\" archive=\"../jalviewApplet.jar\" width=\"140\" height=\"35\"\n";
print H	"                   <param name=\"file\" value=\"$tmp_file.nt1_ali.fasta\">\n";
print H	"                   <param name=\"defaultColour\" value=\"% Identity\">\n";
print H	"                   <param name=\"label\"  value=\"View alignment\">\n";
print H	"                   <param name=\"RGB\"  value=\"7CC391\">\n";
print H	"                   <param name=\"showQuality\" value=\"false\">\n";
print H	"                   <param name=\"showConsensus\" value=\"false\">\n";
print H	"                   <param name=\"showConservation\" value=\"false\">\n";
print H	"</applet></td></tr>\n";
print H "<tr><td>2nd position nt alignment</td> \n";
print H "<td><applet code=\"jalview.bin.JalviewLite\" archive=\"../jalviewApplet.jar\" width=\"140\" height=\"35\"\n";
print H	"                   <param name=\"file\" value=\"$tmp_file.nt2_ali.fasta\">\n";
print H	"                   <param name=\"defaultColour\" value=\"% Identity\">\n";
print H	"                   <param name=\"label\"  value=\"View alignment\">\n";
print H	"                   <param name=\"RGB\"  value=\"7CC391\">\n";
print H	"                   <param name=\"showQuality\" value=\"false\">\n";
print H	"                   <param name=\"showConsensus\" value=\"false\">\n";
print H	"                   <param name=\"showConservation\" value=\"false\">\n";
print H	"</applet></td></tr>\n";
print H "<tr><td>3rd position nt alignment</td> \n";
print H "<td><applet code=\"jalview.bin.JalviewLite\" archive=\"../jalviewApplet.jar\" width=\"140\" height=\"35\"\n";
print H	"                   <param name=\"file\" value=\"$tmp_file.nt3_ali.fasta\">\n";
print H	"                   <param name=\"defaultColour\" value=\"% Identity\">\n";
print H	"                   <param name=\"label\"  value=\"View alignment\">\n";
print H	"                   <param name=\"RGB\"  value=\"7CC391\">\n";
print H	"                   <param name=\"showQuality\" value=\"false\">\n";
print H	"                   <param name=\"showConsensus\" value=\"false\">\n";
print H	"                   <param name=\"showConservation\" value=\"false\">\n";
print H	"</applet></td></tr>\n";
print H "<tr><td>1st+2nd positions nt alignment</td> \n";
print H "<td><applet code=\"jalview.bin.JalviewLite\" archive=\"../jalviewApplet.jar\" width=\"140\" height=\"35\"\n";
print H	"                   <param name=\"file\" value=\"$tmp_file.nt12_ali.fasta\">\n";
print H	"                   <param name=\"defaultColour\" value=\"% Identity\">\n";
print H	"                   <param name=\"label\"  value=\"View alignment\">\n";
print H	"                   <param name=\"RGB\"  value=\"7CC391\">\n";
print H	"                   <param name=\"showQuality\" value=\"false\">\n";
print H	"                   <param name=\"showConsensus\" value=\"false\">\n";
print H	"                   <param name=\"showConservation\" value=\"false\">\n";
print H	"</applet></td></tr></table>\n";



if($gblocks_params ne "") {
	print H "</td><td style=\"border-radius: 7px; -ms-border-radius: 7px; -moz-border-radius:7px; -webkit-border-radius: 7px; -khtml-border-radius: 7px;z-index: 1; padding: 5px;  padding-left:40px; padding-right:15px;  padding-bottom:14px; padding-top:14px; height:90%; background-color:#7CC391;; layer-background-color:#CCCC88; \">\n";
	print H "<h3>After Gblocks</h3>\n";
	print H "<table bgcolor=\"B2CFD4\" Border=0 CELLSPACING=\"6\" style=\"border : thin solid; \">\n";
	print H "<tr><th>Taxa</th><th>GC</th><th>GC 1st</th><th>GC 2nd</th><th>GC 3rd</th></tr>\n";
	foreach ( sort { $gccontent_sp{$a}->{'GC'} <=> $gccontent_sp{$b}->{'GC'} } keys %gccontent_sp ) {
		print H "<tr><td><i>$_</i></td><td>$clean_gccontent_sp{$_}->{'GC'}</td><td>$clean_gccontent_sp_1st{$_}->{'GC'}</td><td>$clean_gccontent_sp_2nd{$_}->{'GC'}</td><td>$clean_gccontent_sp_3rd{$_}->{'GC'}</td></tr>\n";
	}
	#foreach (@aa_spnames) {
	#	print H "<tr><td><i>$_</i></td><td>$clean_gccontent_sp{$_}->{'GC'}</td><td>$clean_gccontent_sp_1st{$_}->{'GC'}</td><td>$clean_gccontent_sp_2nd{$_}->{'GC'}</td><td>$clean_gccontent_sp_3rd{$_}->{'GC'}</td></tr>\n";
	#}
	print H "<tr><td><b>Overall</b></td><td><b>$clean_gccontent{'GC'}</b></td><td><b>$clean_gccontent_1st{'GC'}</b></td><td><b>$clean_gccontent_2nd{'GC'}</b></td><td><b>$clean_gccontent_3rd{'GC'}</b></td></tr>\n";
	print H "</table>\n";

	$tmp_file = $out_file;
	$tmp_file = (split(/\//,$tmp_file))[-1];
	print H "<table border=0><tr><td>1st position nt alignment</td> \n";
	print H "<td><applet code=\"jalview.bin.JalviewLite\" archive=\"../jalviewApplet.jar\" width=\"140\" height=\"35\"\n";
    print H	"                   <param name=\"file\" value=\"$tmp_file.nt1_cleanali.fasta\">\n";
    print H	"                   <param name=\"defaultColour\" value=\"% Identity\">\n";
    print H	"                   <param name=\"label\"  value=\"View alignment\">\n";
    print H	"                   <param name=\"RGB\"  value=\"7CC391\">\n";
    print H	"                   <param name=\"showQuality\" value=\"false\">\n";
    print H	"                   <param name=\"showConsensus\" value=\"false\">\n";
	print H	"                   <param name=\"showConservation\" value=\"false\">\n";
    print H	"</applet></td>\n";
	print H "<tr><td>2nd position nt alignment</td> \n";
	print H "<td><applet code=\"jalview.bin.JalviewLite\" archive=\"../jalviewApplet.jar\" width=\"140\" height=\"35\"\n";
    print H	"                   <param name=\"file\" value=\"$tmp_file.nt2_cleanali.fasta\">\n";
    print H	"                   <param name=\"defaultColour\" value=\"% Identity\">\n";
    print H	"                   <param name=\"label\"  value=\"View alignment\">\n";
    print H	"                   <param name=\"RGB\"  value=\"7CC391\">\n";
    print H	"                   <param name=\"showQuality\" value=\"false\">\n";
    print H	"                   <param name=\"showConsensus\" value=\"false\">\n";
	print H	"                   <param name=\"showConservation\" value=\"false\">\n";
    print H	"</applet></td>\n";
	print H "<tr><td>3rd position nt alignment</td> \n";
	print H "<td><applet code=\"jalview.bin.JalviewLite\" archive=\"../jalviewApplet.jar\" width=\"140\" height=\"35\"\n";
    print H	"                   <param name=\"file\" value=\"$tmp_file.nt3_cleanali.fasta\">\n";
    print H	"                   <param name=\"defaultColour\" value=\"% Identity\">\n";
    print H	"                   <param name=\"label\"  value=\"View alignment\">\n";
    print H	"                   <param name=\"RGB\"  value=\"7CC391\">\n";
    print H	"                   <param name=\"showQuality\" value=\"false\">\n";
    print H	"                   <param name=\"showConsensus\" value=\"false\">\n";
	print H	"                   <param name=\"showConservation\" value=\"false\">\n";
    print H	"</applet></td>\n";
	print H "<tr><td>1st+2nd positions nt alignment</td> \n";
	print H "<td><applet code=\"jalview.bin.JalviewLite\" archive=\"../jalviewApplet.jar\" width=\"140\" height=\"35\"\n";
    print H	"                   <param name=\"file\" value=\"$tmp_file.nt12_cleanali.fasta\">\n";
    print H	"                   <param name=\"defaultColour\" value=\"% Identity\">\n";
    print H	"                   <param name=\"label\"  value=\"View alignment\">\n";
    print H	"                   <param name=\"RGB\"  value=\"7CC391\">\n";
    print H	"                   <param name=\"showQuality\" value=\"false\">\n";
    print H	"                   <param name=\"showConsensus\" value=\"false\">\n";
	print H	"                   <param name=\"showConservation\" value=\"false\">\n";
    print H	"</applet></td></tr></table>\n";
}

print H "</td></tr></table>\n";



print H "</td></tr></table>\n";
#print H "</div>\n";







print H <<HTML1;
	</body></html>
HTML1


close(H);


# my(	%gccontent_sp_1st, 	%gccontent_sp_2nd, 	%gccontent_sp_3rd,
# 	%gccontent_1st, 	%gccontent_2nd, 	%gccontent_3rd,
# 	%gccontent_sp, 		%gccontent					);
# my(	%clean_gccontent_sp_1st, 	%clean_gccontent_sp_2nd, 	%clean_gccontent_sp_3rd,
# 	%clean_gccontent_1st, 		%clean_gccontent_2nd, 		%clean_gccontent_3rd,
# 	%clean_gccontent_sp, 		%clean_gccontent			);

#POR HACER: imprimir el HTML report!!!
#Imprimo nt_ali en caja de texto (con botn guardar) y con botn View alignment (jalview)
#idem para aa_ali
#Imprimo estadsticas de contenido GC

#			GC		GC1		GC2		GC3
#Overall
#Taxon1
#Taxon2
#Taxon3
#Taxon4
#...

