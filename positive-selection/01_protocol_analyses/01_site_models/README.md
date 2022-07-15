# Brief introduction: site models 
We will run `CODEML` under different site models and then compute the LRT 
statistic among nested pairs to identify those that best fit the data. We 
run the software under the homogenous `M0` model and the following site models:
`M0`, `M1a`, `M2a`, `M7`, and `M8`. 

The first pairwise comparison to test is `M1a` vs `M0`. If `M1a` fits the data better than 
`M0`, then it is worth running two additional tests for positive selection:
(i) `M1a` (Nearly Neutral) vs. `M2a` (Positive Selection) models 
and (ii) `M7` (Beta) vs. `M8` (Beta&$\omega$) models. In the first test, the `M1a`
model (null model) has 2 free parameters and two site classes to which codon sites can
be assigned to. In one of these site classes, the value of $\omega$ is always fixed
to 1 ($\omega_{1}=1$), while in the other site class it is allowed to be $\omega_{0}<1$,
and thus under purifying selection. On the other hand, the `M2a` model (alternative
model) has 4 free parameters and three site classes.
As with the `M1a` model, one of these site classes is fixed to 1 ($\omega_{1}=1$) and another
has a value of $\omega$ under purifying selection ($\omega_{0}<1$). Sites assigned to the
third site class, however, are under positive selection with $\omega_{2}>1$. 
In the second test, both `M7` and `M8` models assume that $\omega$ is beta-distributed among 
sites; the former is approximated with 10 site classes while the latter allows for an extra 
site class category to account for positively selected sites (i.e., two free parameters under 
the `M7` model and four under the `M8` model).

# Inference with `CODEML`
We will run our analyses under the 5 site models described above as we aim 
to use these results to find the best-fitting model for each pairwise comparison. 

First, we will create the main working directory:

```sh
# Run from `01_site_models`
mkdir Site_models
cd Site_models 
```

## Setting the control file 
We will use the template provided in this GitHub repository to learn how to 
specify the different options. We will use the command `sed` to find the 
variable names we defined in the template file so we can replace them
with the correct value for each option: 

```sh
# Relative paths from `CODEML` directory 
# created above

# 1. Copy the template control file that 
# we will be later modifying
cp ../../../../templates/template_CODEML.ctl codeml-sites.ctl 

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

## Running `CODEML`
Now that you have the control file, you only need to run `CODEML`. 
If you want to use the compiled version of `CODEML` provided in 
this repository, you can use the code provided in the snippet below.
Otherwise, please modify the code so you can execute 
`CODEML` on your PC:

```sh
# Execute `CODEML`
../../../../src/CODEML/codeml4.9j codeml-sites.ctl | tee logfile_codeml-sites.txt

# Remove unnecessary files 
rm 2N*
```

# Likelihood Ratio Test 
First, we need to extract those lines that include the term `lnL` and then 
remove unnecessary information to keep only the likelihood value. We can 
also match this value with the model under which it was calculated.
Then, we will save the output in a variabe called `header`:
 
```sh
# Keep running from `Site_models`
lnL_vals=$( grep 'lnL' out_sites.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' )
header=$( grep 'Model ' out_sites.txt | sed 's/\:..*//' | sed 's/\ /\_/' )
```

Last, we can output the content of both variables in a file that we can later process
in `R`:

```sh
lnL_vals=$( grep 'lnL' out_sites.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' )
np_vals=$( grep 'lnL' out_sites.txt | sed 's/..*np\:\ //' | sed 's/)..*//' )
header=$( grep 'Model ' out_sites.txt | sed 's/\:..*//' | sed 's/\ /\_/' )
echo $header > lnL_sites.txt
echo $lnL_vals >> lnL_sites.txt
echo $np_vals >> lnL_sites.txt
# Update header names with model names 
# accordingly 
sed -i 's/Model\_1/Model\_1a/' lnL_sites.txt
sed -i 's/Model\_2/Model\_2a/' lnL_sites.txt
```

The `R` script to compute LRT statistics for each model comparison 
can be found [here](Find_bestmodel.R).
The results are already included as part of commented lines.

Once you run this script, you will see that the preferred models for the two positive 
selection tests were `M2a` and `M8`. If we want to extract the 
$\omega$ estimates, we can use the next commands:

```sh 
grep 'MLEs..*(K=3)' -A3 out_sites.txt  > omegas_M2a_M8.txt
printf '\n\n' >> omegas_M2a_M8.txt
grep 'MLEs..*(K=11)' -A3 out_sites.txt >> omegas_M2a_M8.txt
```