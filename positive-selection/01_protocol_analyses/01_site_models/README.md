# Brief introduction: site models 
Site models treat the $\omega$ ratio for any site (codon) in the gene as a random variable 
from a statistical distribution, thus allowing $\omega$ to vary among codons 
([Nielsen and Yang 1998](https://pubmed.ncbi.nlm.nih.gov/9539414/), 
[Yang et al. 2005](https://academic.oup.com/mbe/article/22/4/1107/1083468)).
There are different site models implemented
in `CODEML`, which you can use to analyse the same data in a single run.

For instance, we can compare the simplest model (`M0`, homogenous model as a null model) against an alternative
site model within which it is nested: `M1a`. 
If `M1a` fits the data better than `M0`, then we can compare two other pair of 
nested models: `M1a` (Nearly Neutral) vs. 
`M2a` (Positive Selection) 
([Nielsen and Yang 1998](https://pubmed.ncbi.nlm.nih.gov/9539414/), 
[Wong et al. 2004](https://academic.oup.com/genetics/article/168/2/1041/6059620), 
[Yang et al. 2005](https://academic.oup.com/mbe/article/22/4/1107/1083468)). 
If there are conflicting results, an additional comparison can be performed: 
`M7` (beta) vs. `M8` (beta&$\omega$) 
([Yang et al. 2000](https://pubmed.ncbi.nlm.nih.gov/10790415/)). 
Both `M7` and `M8` models assume that $\omega$ 
is beta-distributed among sites, in which the beta distribution in `M7` is approximated with 10 site 
classes and in `M8` there is an extra site class category to account for positively selected sites 
(i.e., two free parameters under `M7`and four under `M8`).

# Inference with `CODEML`
Before starting our analyses, we create the main working directory:

```sh
# Run from `01_site_models`
mkdir Site_models
cd Site_models 
```

## Setting the control file 
We use the [template control file](../../templates/template_CODEML.ctl) provided in this GitHub 
repository ([here](../../templates/). Then, we use the command `sed` to find the 
variable names defined in the template file so we can replace them
with the correct value for each option:

```sh
# Run from `01_site_models/Site_models`

# 1. Copy the template control file that 
# we will be later modifying
cp ../../../templates/template_CODEML.ctl codeml-sites.ctl 

# 2. Replace variable names with the 
# values needed to run the analysis 

# 2.1. Set path to input files
sed -i 's/ALN/\.\.\/\.\.\/myxovirus\.aln/' codeml-sites.ctl
sed -i 's/TREE/\.\.\/\.\.\/myxovirus\.tree/' codeml-sites.ctl
sed -i 's/OUT/out\_sites\.txt/' codeml-sites.ctl 

# 2.2. Set data to 1 (only one loci)
sed -i 's/NDAT/1/' codeml-sites.ctl 

# 2.3. Specify all the sites model under which
# we want our data set to be analysed
sed -i 's/CODMOD/0/' codeml-sites.ctl 
sed -i 's/NSSIT/0\ 1\ 2\ 7\ 8/' codeml-sites.ctl 
sed -i 's/CODFREQ/7/' codeml-sites.ctl
sed -i 's/ESTFREQ/0/' codeml-sites.ctl
sed -i 's/CLOCK/0/' codeml-sites.ctl 

# 2.3. Set the starting values that are to  
# be used when the model parameters are 
# estimated
sed -i 's/FIXOME/0/' codeml-sites.ctl 
sed -i 's/INITOME/0\.5/' codeml-sites.ctl 
```  

**NOTE:** In this example, we have not copied the alignment or the tree files 
in the working directory `Model_M0`. Instead, we have specified the path to 
these files as you can see in step 2.1 (see code snippet above). We have decided
to do this so we do not keep several copies of the same input files to carry out the different
tests for positive selection in different directories. If you were to run this analysis while having 
the two input files in the main directory, then you would not need to type the relative path to 
these files in the control file as shown above (e.g., `../../myxovirus.aln` or `../../myxovirus.tree` in this example) but
only the file names (e.g., `myxovirus.aln` or `myxovirus.tree`).

## Running `CODEML`
Now that we have the control file, we only need to run `CODEML`. 
If you want to use the compiled version of `CODEML` provided in 
this repository, you can use the code provided in the snippet below.
Otherwise, please modify the code so you can execute 
`CODEML` on your PC:

```sh
# Execute `CODEML` v4.10.5
../../../src/CODEML/codeml4.10.5 codeml-sites.ctl | tee logfile_codeml-sites.txt

# Remove unnecessary files 
rm 2N*
```

# Likelihood Ratio Test 
First, we need to extract those lines that include the term `lnL`. Then, we need 
to remove unnecessary information so that only the likelihood values are kept, which we 
save in a variable called `lnL_vals`. We also extract those lines that include the term `Model` 
and the term `np` so that we can match the likelihood 
values in `lnL_vals` with the corresponding models
under which they were calculated and with the number of parameters in said models. We save the 
collected values in `header` and `np_vals`, respectively.
Last, we output the content of these variables in a file (`lnL_sites.txt`) that we can later process
in `R`:

```sh
# Run from `01_site_models/Site_models`
lnL_vals=$( grep 'lnL' out_sites.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' )
np_vals=$( grep 'lnL' out_sites.txt | sed 's/..*np\:\ //' | sed 's/)..*//' )
header=$( grep 'Model ' out_sites.txt | sed 's/\:..*//' | sed 's/\ /\_/' )
echo $header > lnL_sites.txt
echo $lnL_vals >> lnL_sites.txt
echo $np_vals >> lnL_sites.txt
# Update header names with model names accordingly 
sed -i 's/NSsites\_//g' lnL_sites.txt 
sed -i 's/Model\ 1/Model\_1a/' lnL_sites.txt
sed -i 's/Model\ 2/Model\_2a/' lnL_sites.txt
sed -i 's/Model\ /Model\_/g' lnL_sites.txt
```

The `R` script to compute LRT statistics for each model comparison 
can be found [here](Find_bestmodel.R).
The results are already included as part of commented lines.

If you run this script, you will see that the preferred model is `M8`. 
If you want to extract the estimated parameters in this model, you can use the following command:

```sh 
p0=$( grep 'Parameters in M8' -A1 out_sites.txt | sed -n '2,2p' | sed 's/  p =..*//' | sed 's/p0..* //' )
p=$( grep 'Parameters in M8' -A1 out_sites.txt | sed -n '2,2p' | sed 's/..* p = *//' | sed 's/ ..*//' )
q=$( grep 'Parameters in M8' -A1 out_sites.txt | sed -n '2,2p' | sed 's/..* q =..* //' )
p1=$( grep 'Parameters in M8' -A2 out_sites.txt | sed -n '3,3p' | sed 's/..*p1 = *//' | sed 's/)..*//' )
w=$( grep 'Parameters in M8' -A2 out_sites.txt | sed -n '3,3p' | sed 's/..*w = *//' )
printf "p0\tp\tq\tp1\tw\n" > estpars_M8.tsv
printf $p0"\t"$p"\t"$q"\t"$p1"\t"$w"\n" >> estpars_M8.tsv 
```