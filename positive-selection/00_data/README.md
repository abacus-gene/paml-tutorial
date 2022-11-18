# Obtaining and parsing molecular data

## 1. Downloading sequence data

We are using the "dataset 1" described in [Hou et al. 2007](https://pubmed.ncbi.nlm.nih.gov/17467195/). This study only lists the accession numbers of the myxovirus genes used in the "Materials and Methods" section, but does not provide the reader with the alignments they used. Consequently, we decided to download the sequences using the given accession numbers and parsed the data following the procedure described in [Hou et al. 2007](https://pubmed.ncbi.nlm.nih.gov/17467195/), whenever possible.

Below, you can find the link to the site where each myxovirus sequence was downloaded for the taxa present in the phylogeny studied by [Hou et al. 2007](https://pubmed.ncbi.nlm.nih.gov/17467195/). Some sequences have additional notes used to highlight issues we encountered when trying to downloda them:

* Human, *Homo sapiens*, Mx1 (GenBank accession no: NM_002462). Downloaded from [here](https://www.ncbi.nlm.nih.gov/nuccore/NM_002462.5?from=334&to=2322).
* Rhesus macaque, *Macaca mulata*, Mx1 (Ensembl:ENSMMUT00000021494). Downloaded from [here](http://www.ensembl.org/Macaca_mulatta/Gene/Summary?db=core;g=ENSMMUG00000015329;r=3:5184764-5224295).
* Chimpanzee, *Pan troglodites*, Mx1 (Ensembl:ENSPTRG00000013927). Downloaded from [here](http://www.ensembl.org/Pan_troglodytes/Gene/Summary?db=core;g=ENSPTRG00000013927;r=21:27867000-27965931).
* Orangutan, *Pongo abeli*, cDNA from clone DKFZp469O2020, "note: myxovirus resistance protein 1 (Homo sapiens)" (Gen-Bank accession no: CR860897). Downloded from [here](https://www.ncbi.nlm.nih.gov/nuccore/CR860897.1?from=351&to=2339).
* Sheep, *Ovis aries*, Mx homologue (GenBank accession no:X66093). Downloded from [here](https://www.ncbi.nlm.nih.gov/nuccore/X66093.1?from=72&to=2036).
* Cow, *Bos taurus*, Mx1 (GenBank accession no: NM_173940). Downladed from [here](https://www.ncbi.nlm.nih.gov/nuccore/NM_173940.2?from=99&to=2045).
* Dog, *Canis familiaris*, Mx1 (Gen-Bank accession no: AF239823). Downloaded from [here](https://www.ncbi.nlm.nih.gov/nuccore/AF239823.1?from=116&to=2089).
* Rat, *Rattus norvegicus*, Mx1 (GenBank accession no:BC099784). Downloaded from [here](https://www.ncbi.nlm.nih.gov/nuccore/BC099784.1?from=309&to=2267).
* Pig, *Sus scrofa*, Mx1 (GenBank accession no: AB164037). Downloaded from [here](https://www.ncbi.nlm.nih.gov/nuccore/AB164037.1?from=101&to=2089).
* Mouse, *Mus musculus* (GenBank accession no: NM_010846). Downloded from [here](https://www.ncbi.nlm.nih.gov/nuccore/NM_010846.1?from=214&to=2109).
* Chicken, *Gallus gallus*, Mx1 (GenBank accession no: ~NM_204069~ **NOTE: They had a typo. This sequence is a frog, *Xenopus tropicalis*. The accession no for the chicken is NM_204609.1**). Downloaded from [here](https://www.ncbi.nlm.nih.gov/nuccore/NM_204609.1?from=141&to=2258).
* Duck, *Anas platyrhynchos*, clone 13, Mx protein (GenBank accession no:Z21549). Downloaded from [here](https://www.ncbi.nlm.nih.gov/nuccore/Z21549.1?from=61&to=2226).

In [Hou et al. 2007](https://pubmed.ncbi.nlm.nih.gov/17467195/), the sequences of rhesus macaque and chimpanzee were extracted from the ENSEMBL predicted transcripts (version 39, 2006), whereas all other mRNA sequences or cDNA sequences were retrieved from GenBank.

>> **NOTE**: For ENSEMBL, we could not access ENSEMBL v39, so we downloaded the two sequences from the current ENSEMBL version (ENSEMBL v104, May 2021 -- links provided above).
>> Rhesus macaque: It said that the gene to export was `ENSMMUG00000015329.4 (MX1)`. Therefore this is the variant downloaded.
>> Chimpanzee: It said that the gene to export was `ENSPTRG00000013927.5 (MX1)`. Therefore this is the variant downloaded.

## 2. Generating gene alignments

According to [Hou et al. 2007](https://pubmed.ncbi.nlm.nih.gov/17467195/), the coding sequences were aligned based on translated protein sequences for multispecies data. The alignment of multiple DNA sequences based on predicted protein sequences was generated using the `ClustalW` program.

Given that the alignments or the code to generate the alignments were not provided as supplementary material, we proceeded to generate the gene alignments based on the details they wrote in their paper, which are summarised in the paragraph above.

### Processing raw input data

First, we had to concatenate all the sequences we had previously downloaded (see links provided above) in a unique FASTA file, which we had saved [here](raw_data/data1). In order to do this, we ran the following code from directory [`raw_data`](raw_data):

```sh
# Run from `raw_data` directory 
for i in data1/*fasta
do
dir=$( echo $i | sed 's/\/..*//' )
name=$( echo $i | sed 's/..*\///' | sed 's/\.fasta//' )
sed 's/^>..*/\>'${name}'/' $i >> $dir/$dir"_unaln_raw.fasta"
printf "\n" >> $dir/$dir"_unaln_raw.fasta"
done 
```

Then, we generated an alignment file in which each sequence was written in one line, one after the other:

```sh
# Run from `raw_data` directory
for i in data1/*raw.fasta
do
name=$( echo $i | sed 's/\.fasta//' )
name2=$( echo $name | sed 's/\_raw//' )
../scripts/one_line_fasta.pl $i 
mv $name"_one_line.fa" $name2".fasta"
done
```

### Aligning sequence data

#### First procedure followed to align dataset 1

Once we had our analigned sequences in a unique FASTA file, we used the pipeline [`TranslatorX`](http://www.translatorx.co.uk/) to get the protein alignment. The settings we chose are the following:

* Build the alignment with `ClustalW`.
* Set the same genetic code for all taxa and select the option to guess the most likely reading frame.
* We did not choose to clean the alignment in blocks.

Files were saved as `data1_nuc_aln.fasta` and `data1_prot_aln.fasta`. We also saved the log generated by this tool ([`translatorx_3423.clustalw.log.txt`](alignments/data1/translatorx_3423.clustalw.log.txt)) when `ClustalW` was used to generate the alignment as well as the log file with the compositional bias ([`GCcontent_log_translatorx.txt`](alignments/data1/GCcontent_log_translatorx.txt)). All these output files are saved in [`alignments/data1`](alignments/data1).

Then, we also ran `PAL2NAL` to convert the alignments in codon data (just to double check the nucleotide alignment is the same). This is a perl script (see [here](scripts/pal2nal_v14)) that uses a protein alignment and a file with the DNA sequence to match, and then outputs a codon-based DNA alignment. The commands we used are the following:

```sh
# Run from the directory where this
# README.md file is
cd alignments/data1
mkdir pal2nal_checks
../../scripts/pal2nal_v14/pal2nal.pl data1_prot_aln.fasta data1_nuc_aln.fasta -output fasta > pal2nal_checks/data1pal2nal_out.fasta
# Convert output file to one line fasta file
../../scripts/one_line_fasta.pl pal2nal_checks/data1pal2nal_out.fasta
```

The last step is to convert the alignment from FASTA to PHYLIP format. For that purpose, we use the perl script `FASTAtoPHYL.pl` (saved [here](scripts)) as it follows:

```sh
# Run from this directory where
# this README.md file is
cd alignments/data1
num=$( grep '>' data1_nuc_aln.fasta | wc -l )
len=$( sed -n '2,2p' data1_nuc_aln.fasta | sed 's/\r//' | sed 's/\n//' | wc --m )
perl ../../scripts/FASTAtoPHYL.pl data1_nuc_aln.fasta $num $len 
mv data1_nuc_aln.phy data1.phy
```

#### Second procedure followed to align dataset 1

We detected that the alignment generated when using `ClustalW` resulted in a section of ~189 gaps for most of the mammal taxa except for the outgroup taxa and the mouse sequence. Therefore, we decided to use `MAFFT` as an alternative approach to generate the alignment for dataset 1. First, we downloaded the `translatorX.pl` perl script and saved it in the directory [`translatorX_perl`](translatorX_perl). We copied the unaligned fasta file we had already generated for dataset 1 (i.e., [`data1_unaln.fasta`](translatorX_perl/data1_unaln.fasta)) and ran the following:

```sh
# Run from `00_data/translatorX_perl` directory
./translatorX.pl -i data1_unaln.fasta -p F -o mafft_translatorx
mkdir ../alignments_mafft
cp mafft_translatorx.nt_ali.fasta ../alignments_mafft/data1_nuc_mafft_aln.fasta
cp mafft_translatorx.aa_ali.fasta ../alignments_mafft/data1_prot_mafft_aln.fasta
```

>> **NOTE**: We saved the screen output in file [`log.out_translatorx.txt`](translatorX_perl/log.out_translatorx.txt).

Then, we ran `PAL2NAL` to convert the alignments in codon data (just to double check the nucleotide alignment is the same). This is a perl script (see [here](scripts/pal2nal_v14)) that uses a protein alignment and a file with the DNA sequence to match, and then outputs a codon-based DNA alignment. The commands we used are the following:

```sh
# Run from this directory where this
# README.md file is
cd alignments_mafft
mkdir pal2nal_checks
../scripts/pal2nal_v14/pal2nal.pl data1_prot_mafft_aln.fasta data1_nuc_mafft_aln.fasta -output fasta > pal2nal_checks/data1_pal2nal_mafft_out.fasta
# Convert output file to one line fasta file
../scripts/one_line_fasta.pl pal2nal_checks/data1_pal2nal_mafft_out.fasta
```

The last step is to convert the FASTA file in PHYLIP format. For that purpose, we use the perl script `FASTAtoPHYL.pl` (saved [here](scripts)) as it follows:

```sh
# Run from this directory where
# this README.md file is
cd alignments_mafft
num=$( grep '>' data1_nuc_mafft_aln.fasta | wc -l )
len=$( sed -n '2,2p' data1_nuc_mafft_aln.fasta | sed 's/\r//' | sed 's/\n//' | wc --m )
perl ../scripts/FASTAtoPHYL.pl data1_nuc_mafft_aln.fasta $num $len 
mv data1_nuc_mafft_aln.phy data1_mafft.phy
```

#### Third procedure followed to align dataset 1

We detected that both alignments generated with `ClustalW` and `MAFFT` resulted in a section with gaps (more gaps with `ClustalW` than with `MAFFT`). Therefore, we took the alignment generated with `MAFFT`, `data1_mafft.phy`, and manually got rid of the first 213 nucleotides (71 codons). The resulting alignment (1992 bp) can be found in directory [`alignments_mafft_nogaps`](alignments_mafft_nogaps) as `data1_mafft_nogaps.phy`.

## 3. Generating gene tree

### Procedure described in Hou et al. 2007

The phylogenetic tree was constructed by the neighbor-joining (NJ) method (Saitou and Nei, 1987) based on predicted protein sequences. The p-distances (proportion of differences) with a high resolution of branching pattern were calculated (Nei and Kumar, 2000). One thousand bootstrap replicates were carried out to test the support for each node in the tree.

### Procedure followed here

Here, we followed a maximum-likelihood approach to estimate the best-scoring maximum likelihood tree with `RAxML v8.2.10` ([Stamatakis 2014](https://academic.oup.com/bioinformatics/article/30/9/1312/238053?login=true)). We used the following code:

```sh
# NOTE: This command was included in a script 
# (`scripts/run_raxml.sh`) as part of a job 
# array ran in an HPC.
# Directory `1` corresponds to `RAxML` directory,
# `2` to `RAxML_mafft`, and `3` to
# `RAxML_mafft_nogaps`.
# The first directory contains the alignment
# generated with `ClustalW`, the second the one 
# generated with `MAFFT`, and the third has the
# same alignment as in directory `MAFFT` but 
# without the gap section.
# The input file, variable `$seq`, is `data1.phy`, 
# `data1_mafft.phy`, or `data1_mafft_nogaps.phy`; respectively.
# Variable `$name` is `data1`, `data1_mafft`, or 
# `data1_mafft_nogaps`; respectively. 
raxmlHPC -f a -m GTRGAMMA -p 12345 -# 100 -x 12345 -# 500 -s $seq -o Duck,Chicken -n $name
```

For each analysis with each gene alignment, the best-scoring ML trees are `RAxML/1/RAxML_bestTree.data1`, `RAxML_mafft/1/RAxML_bestTree.data1_mafft`, and `RAxML_mafft_nogaps/1/RAxML_bestTree.data1_mafft_nogaps`; respectively. We viewed these trees in `FigTree` and made sure that the branches for duck and chicken were used to root the tree. We also oredered the nodes in increasing order so it matched the way the tree had been displayed in the original paper.

The rooted trees were then saved as `data1.tree` (Newick format) in each corresponding directory. All of the inferred gene trees agree with the tree topology.

>> **NOTE**: The best-scoring ML tree for dataset 1 in [Hou et al. 2007](https://pubmed.ncbi.nlm.nih.gov/17467195/) (Fig. S1) has *C. familiaris* as an outgroup to the clade with pig, cow, and sheep. We did not obtain this tree topology. We have kept the tree topology estimated by `RAxML` as the same tree topology was inferred regardless of the alignment used.

Then, we just used the next command to get rid of the branch lengths and obtain the tree topologies to be used in `CODEML`:

```sh
# Run this code from the 
# directory where this
# `README.md` file is 
for i in RAxML* 
do
echo Parsing tree in directory $i
printf "12  1\n" > $i/data1_nobl.tree
sed 's/\:[0-9]*\.[0-9]*//g' $i/data1.tree >> $i/data1_nobl.tree
done
```

## 4. Summary

For the analyses with `CODEML`, we decided to use the alignment generated with `MAFFT` and in which the gap section at the beginning of the alignment was removed. Before getting started, however, we decided to reorder the position of the sequences in the alignment so they matched the order in which the taxa are found in the phylogeny. In that way, we expect the output data by `CODEML` to be read in a more intuitive manner. For that purpose, we did the following:

```sh
# Create an alignment file with the PHYLIP 
# header and the next line, then grep the 
# sequences according to phylogeny
# Run this code from the directory where this 
# `README.md` file is found
sed -n '1,2p' alignments_mafft_nogaps/data1_mafft_nogaps.phy > Mx_aln.phy
grep 'Rhesus' alignments_mafft_nogaps/data1_mafft_nogaps.phy >> Mx_aln.phy
grep 'Orangutan' alignments_mafft_nogaps/data1_mafft_nogaps.phy >> Mx_aln.phy
grep 'Chimpanzee' alignments_mafft_nogaps/data1_mafft_nogaps.phy >> Mx_aln.phy
grep 'Human' alignments_mafft_nogaps/data1_mafft_nogaps.phy >> Mx_aln.phy
grep 'Dog' alignments_mafft_nogaps/data1_mafft_nogaps.phy >> Mx_aln.phy
grep 'Pig' alignments_mafft_nogaps/data1_mafft_nogaps.phy >> Mx_aln.phy
grep 'Cow' alignments_mafft_nogaps/data1_mafft_nogaps.phy >> Mx_aln.phy
grep 'Sheep' alignments_mafft_nogaps/data1_mafft_nogaps.phy >> Mx_aln.phy
grep 'Rat' alignments_mafft_nogaps/data1_mafft_nogaps.phy >> Mx_aln.phy
grep 'Mouse' alignments_mafft_nogaps/data1_mafft_nogaps.phy >> Mx_aln.phy
grep 'Chicken' alignments_mafft_nogaps/data1_mafft_nogaps.phy >> Mx_aln.phy
grep 'Duck' alignments_mafft_nogaps/data1_mafft_nogaps.phy >> Mx_aln.phy

# Add an identifier in the tag 
# to show they are myxovirus 
# sequences
sed -i 's/      /\_Mx      /g' Mx_aln.phy
```

Last, as we modified the tag names in the alignment, we need to do the same for the tree file. Therefore:

```sh
# Run this code from the directory where this 
# `README.md` file is found
cp RAxML_mafft_nogaps/data1_nobl.tree Mx_root.tree
sed -i 's/Chimpanzee/Chimpanzee\_Mx/' Mx_root.tree
sed -i 's/Human/Human\_Mx/' Mx_root.tree
sed -i 's/Orangutan/Orangutan\_Mx/' Mx_root.tree
sed -i 's/Rhesus\_macaque/Rhesus\_macaque\_Mx/' Mx_root.tree
sed -i 's/Sheep/Sheep\_Mx/' Mx_root.tree
sed -i 's/Cow/Cow\_Mx/' Mx_root.tree
sed -i 's/Pig/Pig\_Mx/' Mx_root.tree
sed -i 's/Dog/Dog\_Mx/' Mx_root.tree
sed -i 's/Mouse/Mouse\_Mx/' Mx_root.tree
sed -i 's/Rat/Rat\_Mx/' Mx_root.tree
sed -i 's/Duck/Duck\_Mx/' Mx_root.tree
sed -i 's/Chicken/Chicken\_Mx/' Mx_root.tree
```

After that, we manually unroot the tree before starting the analyses in `CODEML` (i.e., we remove the parentheses that include the duck and chicken tags) and save it as `Mx_unroot.tree`. You can also generate the tree if you run the next command:

```sh
# Run this code from the directory where this 
# `README.md` file is found
cp Mx_root.tree Mx_unroot.tree
sed -i 's/(Duck/Duck/' Mx_unroot.tree
sed -i 's/Chicken\_Mx)/Chicken\_Mx/' Mx_unroot.tree
printf "\n" >> Mx_unroot.tree
# Resulting unrooted tree you should see in the
# output file:
#   ((((((Chimpanzee_Mx,Human_Mx),Orangutan_Mx),Rhesus_macaque_Mx),(((Sheep_Mx,Cow_Mx),Pig_Mx),Dog_Mx)),(Mouse_Mx,Rat_Mx)),Duck_Mx,Chicken_Mx);
```

Once the input files for `CODEML` have been generated (i.e., `Mx_aln.phy`, `Mx_root.tree`, and `Mx_unroot.tree`), we can copy them into the [`00_protocol_CODEML`](../01_protocol_analyses) directory in this GitHub repository. Please note that the `Mx_aln.phy` alignment file you will find here still has STOP codons, and hence it has 1992 nucleotides. Nevertheless, the `Mx_aln.phy` you will find in the [`00_protocol_CODEML`](../01_protocol_analyses) directory does not have STOP codons (i.e., the alignment length is reduced 3 nucleotides, and so the PHYLIP header is changed). We manually removed the last three nucleotides and changed the header so the alignment is now ready to be used with `CODEML`.
