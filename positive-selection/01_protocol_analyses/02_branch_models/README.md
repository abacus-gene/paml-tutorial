# Brief introduction: branch models 
Branch models assume that $\omega$ varies across the branches of the tree
([Yang 1998](https://pubmed.ncbi.nlm.nih.gov/9580986/),
[Yang and Nielsen 1998](https://link.springer.com/article/10.1007/PL00006320)). 
Consequently, we need to specify which branches in our phylogeny 
are to have different $\omega$ values. In this example, we will run two 
analyses with `CODEML`: one in which the chicken lineage is labelled as the foreground branch (i.e., branch 
that will have a different value of $\omega) and another with the 
duck lineage instead. 

Before we get started with the tutorial, we need to generate two copies of the main tree file, 
[`myxovirus.tree`](../myxovirus.tree), and add the tag `#1` to identify the foreground 
branch in each case:   

```sh 
# Run from `02_branch_models`
cp ../myxovirus.tree ../myxovirus_branch_chicken.tree 
cp ../myxovirus.tree ../myxovirus_branch_duck.tree 

# Add tags 
sed -i 's/Chicken\_Mx/Chicken\_Mx\ \#1/' ../myxovirus_branch_chicken.tree
sed -i 's/Duck\_Mx/Duck\_Mx\ \#1/' ../myxovirus_branch_duck.tree
```

# Inference with `CODEML`
Given that we want to test two hypotheses (i.e., one in which the foreground lineage is the 
duck lineage and another with the chicken lineage), we create two different directories 
to run the branch model under the two hypotheses:

```sh
# Run from `02_branch_models`
mkdir -p Branch_model_chicken/CODEML Branch_model_duck/CODEML
```

## Setting the control file 
We use the [template control file](../../templates/template_CODEML.ctl) provided in this GitHub 
repository ([here](../../templates/). Then, we use the command `sed` to find the 
variable names defined in the template file so we can replace them
with the correct value for each option:

```sh
# Run from `02_branch_models`

# 1. Copy the template control file that 
# we will be later modifying
cp ../../templates/template_CODEML.ctl Branch_model_chicken/CODEML/codeml-branch.ctl 
cp ../../templates/template_CODEML.ctl Branch_model_duck/CODEML/codeml-branch.ctl 

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

**NOTE:** In this example, we have not copied the alignment or the tree files 
in the working directory `Model_M0`. Instead, we have specified the path to 
these files as you can see in step 2.1 (see code snippet above). We have decided
to do this so we do not keep several copies of the same input files to carry out the different
tests for positive selection in different directories. If you were to run this analysis while having 
the two input files in the main directory, then you would not need to type the relative path to 
these files in the control file as shown above (e.g., `../../../myxovirus.aln` or `../../../myxovirus.tree` in this example) but
only the file names (e.g., `myxovirus.aln` or `myxovirus.tree`).

## Running `CODEML`
Now that we have the control file for each analysis, we only need to run `CODEML`. 
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
../../../../src/CODEML/codeml4.10.5 codeml-branch.ctl | tee logfile_codeml-branch_$name.txt
# Remove unnecessary files 
rm 2N*
# Go back to home_dir 
cd $home_dir
done
```

# Likelihood Ratio Test
Now, we can test if the branch model fits the data better than the null model (`M0` model), which 
we previously ran (i.e., the results for the `M0` model can be found in the `../00_homogenous_model/Model_M0` 
directory).

First, we extract those lines that have the `lnL` term and then 
remove unnecessary information to keep only the `lnL` value. We will do this for the three analyses 
which likelihood values we want to compare:

```sh
# Run from `02_branch_models`
grep 'lnL' Branch_model_chicken/CODEML/out_chicken_branch.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' > lnL_branch_mods.txt
grep 'lnL' Branch_model_duck/CODEML/out_duck_branch.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> lnL_branch_mods.txt
grep 'lnL' ../00_homogenous_model/Model_M0/out_M0.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> lnL_branch_mods.txt
```

The R script [`Find_bestmodel.R`](Find_bestmodel.R) is used to compute the LRT. You can also find 
the results for each test written there as commented lines.

In addition, we can extract the $\omega$ ratios when the chicken lineage is set as the 
foreground lineage and then when the chicken is the foreground lineage:

```sh
# Run from `02_branch_models`
w_c_back=$( grep 'w (dN/dS) for branches'  Branch_model_chicken/CODEML/out_chicken_branch.txt | sed 's/..*: *//' | sed 's/ ..*//' )
w_c_for=$( grep 'w (dN/dS) for branches'  Branch_model_chicken/CODEML/out_chicken_branch.txt | sed 's/..* //')
w_d_back=$( grep 'w (dN/dS) for branches'  Branch_model_duck/CODEML/out_duck_branch.txt | sed 's/..*: *//' | sed 's/ ..*//' )
w_d_for=$( grep 'w (dN/dS) for branches'  Branch_model_duck/CODEML/out_duck_branch.txt | sed 's/..* //')
printf "w_back_chicken\tw_fore_chicken\tw_back_duck\tw_fore_duck\n" > w_est_branches.tsv 
printf $w_c_back"\t"$w_c_for"\t"$w_d_back"\t"$w_d_for"\n" >> w_est_branches.tsv 
```