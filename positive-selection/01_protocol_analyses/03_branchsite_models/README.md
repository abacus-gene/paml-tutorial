# Brief introduction: branch-site models 
Branch-site models assume that $\omega\$ can vary both across sites and across lineages
([Yang and Nielsen 2002](https://link.springer.com/article/10.1007/PL00006320)). 
In these analyses, the aim is to detect positive selection at specific sites along the 
so-called _foreground_ branches (i.e., user-specified branches with tags along which 
positive selection at specific sites will be tested).

As we did with the tutorial to run the branch model, here we will run two analyses: one in which the foreground branch is the chicken lineage and 
another with the duck lineage instead. If you have run the tutorial for branch models 
[here](../02_branch_models/README.md), you will already have the tree files with the tags selecting 
the foreground branches. Otherwise, please run the following code snippet before getting started 
with the tutorial:

```sh 
# Run from `02_branch_models`
cp ../myxovirus.tree ../myxovirus_branch_chicken.tree 
cp ../myxovirus.tree ../myxovirus_branch_duck.tree 

# Add tags 
sed -i 's/Chicken\_Mx/Chicken\_Mx\ \#1/' ../myxovirus_branch_chicken.tree
sed -i 's/Duck\_Mx/Duck\_Mx\ \#1/' ../myxovirus_branch_duck.tree
```

# Inference with `CODEML`

## Part 1: alternative models 
Given that we want to test two hypotheses (i.e., one in which the foreground lineage will be the 
duck lineage and another with the chicken lineage), we create two different directories 
to run the corresponding analyses under the branch model:

```sh
# Run from `03_branchsite_models`
mkdir -p Branchsite_model_chicken/CODEML Branchsite_model_duck/CODEML
```

### Setting the control file 
We use the [template control file](../../templates/template_CODEML.ctl) provided in this GitHub 
repository ([here](../../templates/). Then, we use the command `sed` to find the 
variable names defined in the template file so we can replace them
with the correct value for each option:

```sh
# Relative paths from `03_branchsite_models` directory 
# created above

# 1. Copy the template control file that 
# we will be later modifying
cp ../../templates/template_CODEML.ctl Branchsite_model_chicken/CODEML/codeml-branchsite.ctl 
cp ../../templates/template_CODEML.ctl Branchsite_model_duck/CODEML/codeml-branchsite.ctl 

# 2. Replace variable names with the 
# values needed to run the analysis 

# 2.1. Set path to input files
sed -i 's/ALN/\.\.\/\.\.\/\.\.\/myxovirus\.aln/' Branchsite_model_chicken/CODEML/codeml-branchsite.ctl 
sed -i 's/TREE/\.\.\/\.\.\/\.\.\/myxovirus\_branch\_chicken\.tree/' Branchsite_model_chicken/CODEML/codeml-branchsite.ctl 
sed -i 's/OUT/out\_chicken\_branchsite\.txt/' Branchsite_model_chicken/CODEML/codeml-branchsite.ctl 

sed -i 's/ALN/\.\.\/\.\.\/\.\.\/myxovirus\.aln/' Branchsite_model_duck/CODEML/codeml-branchsite.ctl 
sed -i 's/TREE/\.\.\/\.\.\/\.\.\/myxovirus\_branch\_duck\.tree/' Branchsite_model_duck/CODEML/codeml-branchsite.ctl 
sed -i 's/OUT/out\_duck\_branchsite\.txt/' Branchsite_model_duck/CODEML/codeml-branchsite.ctl 

# 2.2. Go over all control files
for i in Branchsite_model_duck Branchsite_model_chicken
do
# 2.2.1 Set data to 1 (only one loci)
sed -i 's/NDAT/1/' $i/CODEML/codeml-branchsite.ctl 

# 2.2.2 Specify all the sites model under which
# we want our data set to be analysed
sed -i 's/CODMOD/2/' $i/CODEML/codeml-branchsite.ctl 
sed -i 's/NSSIT/2/' $i/CODEML/codeml-branchsite.ctl 
sed -i 's/CODFREQ/7/' $i/CODEML/codeml-branchsite.ctl 
sed -i 's/ESTFREQ/0/' $i/CODEML/codeml-branchsite.ctl 
sed -i 's/CLOCK/0/' $i/CODEML/codeml-branchsite.ctl 

# 2.2.3 Set the starting values that are to  
# be used when the model parameters are 
# estimated
sed -i 's/FIXOME/0/' $i/CODEML/codeml-branchsite.ctl 
sed -i 's/INITOME/0\.5/' $i/CODEML/codeml-branchsite.ctl  

done
```

**NOTE:** In this example, we have not copied the alignment or the tree files 
in the working directories (`Branchsite_model_chicken` or `Branchsite_model_duck`). Instead, we have specified the path to 
these files as you can see in step 2.1 in the code snippet above. We have decided
to do this so we do not keep several copies of the same input files to carry out the different
tests for positive selection in different directories. Nevertheless, 
you may copy these two input files in the working directory of this analysis or other analysis 
you may run for other projects, 
in which case you do not need to type the relative path to these files in the control
file as shown above (e.g., `../../../myxovirus.aln` or `../../../myxovirus.tree` in this example) but
only the file names (e.g., `myxovirus.aln` or `myxovirus.tree`).

### Running `CODEML`
Now that we have the control file, we only need to run `CODEML`. 
If you want to use the compiled version of `CODEML` provided in 
this repository, you can use the code provided in the snippet below.
Otherwise, please modify the code so you can execute 
`CODEML` on your PC:

```sh
# Run from `03_branchsite_models`
home_dir=$( pwd )
for i in Branchsite_model_duck Branchsite_model_chicken
do
# Move to dir 
cd $i/CODEML/
# Execute `CODEML`
name=$( echo $i | sed 's/..*\_//' )
../../../../src/CODEML/codeml4.10.5 codeml-branchsite.ctl | tee logfile_codeml-branchsite_$name.txt
# Remove unnecessary files 
rm 2N*
# Go back to home_dir 
cd $home_dir
done
```  

## Part 2: null models 
We now run `CODEML` under a modified version of the `M2a` model where 
the value of $\omega$ for site classes "2a" and "2b" is fixed to $\omega=1$. Consequently, 
as this parameter is not estimated, this model has one parameter less than the 
so-called branch-site model A 
([Zhang et al. 2005](https://academic.oup.com/mbe/article/22/12/2472/1009544)). 
We then create one directory for each analysis within 
the corresponding main directories previously created:

```sh
# Run from `03_branchsite_models`
mkdir -p Branchsite_model_chicken/CODEML_2 Branchsite_model_duck/CODEML_2
```

### Setting the control file 
We use the [template control file](../../templates/template_CODEML.ctl) provided in this GitHub 
repository ([here](../../templates/). Then, we use the command `sed` to find the 
variable names defined in the template file so we can replace them
with the correct value for each option:

```sh
# Relative paths from `03_branchsite_models` directory 
# created above

# 1. Copy the template control file that 
# we will be later modifying
cp ../../templates/template_CODEML.ctl Branchsite_model_chicken/CODEML_2/codeml-branchsite_null.ctl 
cp ../../templates/template_CODEML.ctl Branchsite_model_duck/CODEML_2/codeml-branchsite_null.ctl 

# 2. Replace variable names with the 
# values needed to run the analysis 

# 2.1. Set path to input files
sed -i 's/ALN/\.\.\/\.\.\/\.\.\/myxovirus\.aln/' Branchsite_model_chicken/CODEML_2/codeml-branchsite_null.ctl 
sed -i 's/TREE/\.\.\/\.\.\/\.\.\/myxovirus\_branch\_chicken\.tree/' Branchsite_model_chicken/CODEML_2/codeml-branchsite_null.ctl 
sed -i 's/OUT/out\_chicken\_branchsite\.txt/' Branchsite_model_chicken/CODEML_2/codeml-branchsite_null.ctl 

sed -i 's/ALN/\.\.\/\.\.\/\.\.\/myxovirus\.aln/' Branchsite_model_duck/CODEML_2/codeml-branchsite_null.ctl 
sed -i 's/TREE/\.\.\/\.\.\/\.\.\/myxovirus\_branch\_duck\.tree/' Branchsite_model_duck/CODEML_2/codeml-branchsite_null.ctl 
sed -i 's/OUT/out\_duck\_branchsite\.txt/' Branchsite_model_duck/CODEML_2/codeml-branchsite_null.ctl 

# 2.2. Go over all control files
for i in Branchsite_model_duck Branchsite_model_chicken
do
# 2.2.1 Set data to 1 (only one loci)
sed -i 's/NDAT/1/' $i/CODEML_2/codeml-branchsite_null.ctl 

# 2.2.2 Specify all the sites model under which
# we want our data set to be analysed
sed -i 's/CODMOD/2/' $i/CODEML_2/codeml-branchsite_null.ctl 
sed -i 's/NSSIT/2/' $i/CODEML_2/codeml-branchsite_null.ctl 
sed -i 's/CODFREQ/7/' $i/CODEML_2/codeml-branchsite_null.ctl 
sed -i 's/ESTFREQ/0/' $i/CODEML_2/codeml-branchsite_null.ctl 
sed -i 's/CLOCK/0/' $i/CODEML_2/codeml-branchsite_null.ctl 

# 2.2.3 Fix value of omega to 1
sed -i 's/FIXOME/1/' $i/CODEML_2/codeml-branchsite_null.ctl 
sed -i 's/INITOME/1/' $i/CODEML_2/codeml-branchsite_null.ctl 

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


### Running `CODEML`
Now that we have the control file, we only need to run `CODEML`. 
If you want to use the compiled version of `CODEML` provided in 
this repository, you can use the code provided in the snippet below.
Otherwise, please modify the code so you can execute 
`CODEML` on your PC:

```sh
# Run from `03_branchsite_models`
home_dir=$( pwd )
for i in Branchsite_model_duck Branchsite_model_chicken
do
# Move to dir 
cd $i/CODEML_2/
# Execute `CODEML`
name=$( echo $i | sed 's/..*\_//' )
../../../../src/CODEML/codeml4.10.5 codeml-branchsite_null.ctl  | tee logfile_codeml-branchsite_$name"_null.txt"
# Remove unnecessary files 
rm 2N*
# Go back to home_dir 
cd $home_dir
done
```  

# Likelihood Ratio Test 
Now, we can test whether the branch-site model A fits the data better than its modified version 
with one less parameter. 

First, we extract those lines that have the `lnL` term and then 
remove unnecessary information to keep only the `lnL` value. We will do this for the three analyses 
which likelihood values we want to compare:

```sh
# Run from `03_branchsite_models`
grep 'lnL' Branchsite_model_chicken/CODEML/out_chicken_branchsite.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' > lnL_branchsite_mods.txt
grep 'lnL' Branchsite_model_duck/CODEML/out_duck_branchsite.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> lnL_branchsite_mods.txt
grep 'lnL' Branchsite_model_chicken/CODEML_2/out_chicken_branchsite.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> lnL_branchsite_mods.txt
grep 'lnL' Branchsite_model_duck/CODEML_2/out_duck_branchsite.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> lnL_branchsite_mods.txt
```

The R script [`Find_bestmodel.R`](Find_bestmodel.R) is then used to compute LRT.
The results are already written in the script as commented lines.

In addition, you can also extract the $\omega$ estimates for each analysis:

```sh 
# Run from `03_branchsite_models`
grep 'MLEs of' -A5 Branchsite_model_chicken/CODEML/out_chicken_branchsite.txt > Branchsite_model_chicken/branchsite_chicken_MLEs.txt
grep 'MLEs of' -A5 Branchsite_model_chicken/CODEML_2/out_chicken_branchsite.txt > Branchsite_model_chicken/branchsite_chicken_MLEs_2.txt
grep 'MLEs of' -A5 Branchsite_model_duck/CODEML/out_duck_branchsite.txt > Branchsite_model_duck/branchsite_duck_MLEs.txt
grep 'MLEs of' -A5 Branchsite_model_duck/CODEML_2/out_duck_branchsite.txt > Branchsite_model_duck/branchsite_duck_MLEs_2.txt
```