# Brief introduction: branch models 
Branch models assume that positive selection may occur
at specific sites in specific branches of the phylogeny. 
Consequently, we need to specify which branches in our phylogeny 
are to have different $\omega$ values. As we do not allow
$\omega$ to vary across sites, we will select two branches that correspond to the two lineages 
for which we want to test for positive selection: the duck and the chicken lineages. We will 
use the snippet code below to generate two copies of the main tree file, 
[`myxovirus.tree`](../myxovirus.tree), and add the tag `#1` to identify these two lineages,
one tag per tree:

```sh 
# Run from `02_branch_models`
cp ../myxovirus.tree ../myxovirus_branch_chicken.tree 
cp ../myxovirus.tree ../myxovirus_branch_duck.tree 

# Add tags 
sed -i 's/Chicken\_Mx/Chicken\_Mx\ \#1/' ../myxovirus_branch_chicken.tree
sed -i 's/Duck\_Mx/Duck\_Mx\ \#1/' ../myxovirus_branch_duck.tree
```

# Inference with `CODEML`
We will run two analyses: one having the branch leading to the duck lineage selected 
and another with the chicken lineage.

First, we will create the corresponding main directories:

```sh
# Run from `02_branch_models`
mkdir -p Branch_model_chicken/CODEML Branch_model_duck/CODEML
```

## Setting the control file 
We will use the template provided in this GitHub repository to learn how to 
specify the different options. We will use the command `sed` to find the 
variable names we defined in the template file so we can replace them
with the correct value for each option: 

```sh
# Relative paths from `02_branch_models` directory

# 1. Copy the template control file that 
# we will be later modifying
cp ../../../templates/template_CODEML.ctl Branch_model_chicken/CODEML/codeml-branch.ctl 
cp ../../../templates/template_CODEML.ctl Branch_model_duck/CODEML/codeml-branch.ctl 

# 2. Replace variable names with the 
# values needed to run the analysis 

# 2.1. Set path to input files
sed -i 's/ALN/\.\.\/\.\.\/\.\.\/myxovirus\.aln/' Branch_model_chicken/CODEML/codeml-branch.ctl 
sed -i 's/TREE/\.\.\/\.\.\/\.\.\/myxovirus\_branch\_chicken\.tree/' Branch_model_chicken/CODEML/codeml-branch.ctl 
sed -i 's/OUT/out\_chicken\_branch\.txt/' Branch_model_chicken/CODEML/codeml-branch.ctl 

sed -i 's/ALN/\.\.\/\.\.\/\.\.\/myxovirus\.aln/' Branch_model_duck/CODEML/codeml-branch.ctl 
sed -i 's/TREE/\.\.\/\.\.\/\.\.\/myxovirus\_branch\_duck\.tree/' Branch_model_duck/CODEML/codeml-branch.ctl 
sed -i 's/OUT/out\_duck\_branch\.txt/' Branch_model_duck/CODEML/codeml-branch.ctl 

# 2.2. Go over all control files
for i in Branch_model_duck Branch_model_chicken
do
# 2.2.1 Set data to 1 (only one loci)
sed -i 's/NDAT/1/' $i/CODEML/codeml-branch.ctl 

# 2.2.2 Specify all the sites model under which
# we want our data set to be analysed
sed -i 's/CODMOD/2/' $i/CODEML/codeml-branch.ctl 
sed -i 's/NSSIT/0/' $i/CODEML/codeml-branch.ctl 
sed -i 's/CODFREQ/7/' $i/CODEML/codeml-branch.ctl 
sed -i 's/ESTFREQ/0/' $i/CODEML/codeml-branch.ctl 
sed -i 's/CLOCK/0/' $i/CODEML/codeml-branch.ctl 

# 2.2.3 Set the starting values that are to  
# be used when the model parameters are 
# estimated
sed -i 's/FIXOME/0/' $i/CODEML/codeml-branch.ctl 
sed -i 's/INITOME/0\.5/' $i/CODEML/codeml-branch.ctl  

done
```  

## Running `CODEML`
Now that you have the control file, you only need to run `CODEML`. 
If you want to use the compiled version of `CODEML` provided in 
this repository, you can use the code provided in the snippet below.
Otherwise, please modify the code so you can execute 
`CODEML` on your PC:

```sh
# Run from `02_branch_models`
home_dir=$( pwd )
for i in Branch_model_duck Branch_model_chicken
do
# Move to dir 
cd $i/CODEML/
# Execute `CODEML`
name=$( echo $i | sed 's/..*\_//' )
../../../../../src/CODEML/codeml4.9j codeml-branch.ctl | tee logfile_codeml-branch_$name.txt
# Remove unnecessary files 
rm 2N*
# Go back to home_dir 
cd $home_dir
done
```

# Likelihood Ratio Test
Now, we can test if the branch model fits the data better than the null model, the `M0` model, which 
we previously ran. The results for the `M0` model can be found in the `../00_homogenous_model/M0_model` 
directory.

First, we need to extract those lines that have `lnL` and then 
remove unnecessary information to keep only the `lnL` value. We will do this for both 
analyses:

```sh
# Run from `02_branch_models`
grep 'lnL' Branch_model_chicken/CODEML/out_chicken_branch.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' > lnL_branch_mods.txt
grep 'lnL' Branch_model_duck/CODEML/out_duck_branch.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> lnL_branch_mods.txt
grep 'lnL' ../00_homogenous_model/M0_model/out_M0.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> lnL_branch_mods.txt
```

The R script `Find_bestmodel.R` to compute LRT can be found in the `02_branch_models` directory.
The results are already printed out there in commented lines.

In addition, we can extract the $\omega$ ratios when running the branch model and 
having the chicken lineage as a foreground lineage and when the foreground lineage 
is the duck: 

```sh
# Run from `02_branch_models`
grep 'w ratios' -A2 Branch_model_chicken/CODEML/out_chicken_branch.txt | sed -n '2,2p' > Branch_model_chicken/tree_chicken_wratios.tree
grep 'w ratios' -A2 Branch_model_duck/CODEML/out_duck_branch.txt | sed -n '2,2p' > Branch_model_duck/tree_duck_wratios.tree
```