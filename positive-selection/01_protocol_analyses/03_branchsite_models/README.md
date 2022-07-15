# Brief introduction: branch-site models 
Branch-site models assume that positive selection may happen at 
specific sites and also across specific branches of the given phylogeny. 
The main difference between branch-site and branch models is that the individual models 
in the former model that are assigned to a specific branch model will have more than
one site class. While $\omega$ is not allowed to change across sites 
in branch models, you can find more than one category for $\omega$ in branch-site 
models (i.e., as many categories as site classes are to model).

# Inference with `CODEML`
We will run two analyses: one in which the foreground branch is the chicken lineage and 
another in which the duck lineage is. In case you have not run 
the analyses under the branch models
(see [here](../02_branch_models/README.md)), you will need to generate the input tree files 
needed for this analysis as the foreground lineages need to be identified with a tag. You can 
run the next commands if you need to do this: 

```sh 
# Run the code below from `03_branchsite_models` 
# only if you have not generated before the input 
# tree files identifying the foreground branches 
cp ../myxovirus.tree ../myxovirus_branch_chicken.tree 
cp ../myxovirus.tree ../myxovirus_branch_duck.tree 

# Add tags 
sed -i 's/Chicken\_Mx/Chicken\_Mx\ \#1/' ../myxovirus_branch_chicken.tree
sed -i 's/Duck\_Mx/Duck\_Mx\ \#1/' ../myxovirus_branch_duck.tree
```

## Part 1: alternative models 
First, we need to create the corresponding working directories where the analyses 
will take place for each scenario described above:

```sh
# Run from `03_branchsite_models`
mkdir -p Branchsite_model_chicken/CODEML Branchsite_model_duck/CODEML
```

### Setting the control file 
We will use the template provided in this GitHub repository to learn how to 
specify the different options. We will use the command `sed` to find the 
variable names we defined in the template file so we can replace them
with the correct value for each option: 

```sh
# Relative paths from `03_branchsite_models` directory 
# created above

# 1. Copy the template control file that 
# we will be later modifying
cp ../../../templates/template_CODEML.ctl Branchsite_model_chicken/CODEML/codeml-branchsite.ctl 
cp ../../../templates/template_CODEML.ctl Branchsite_model_duck/CODEML/codeml-branchsite.ctl 

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

### Running `CODEML`
Now that you have the control file, you only need to run `CODEML`. 
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
../../../../../src/CODEML/codeml4.9j codeml-branchsite.ctl | tee logfile_codeml-branchsite_$name.txt
# Remove unnecessary files 
rm 2N*
# Go back to home_dir 
cd $home_dir
done
```  

## Part 2: null models 
We will now run `CODEML` under a modified version of the `M2a` model where 
the value of $\omega$ for site classes "2a" and "2b" is fixed to $\omega=1$. Consequently, 
as this parameter is not estimated, this model has one parameter less than the 
so-called branch-site model A.

```sh
# Run from `03_branchsite_models`
mkdir -p Branchsite_model_chicken/CODEML_2 Branchsite_model_duck/CODEML_2
```

### Setting the control file 
We will use the template provided in this GitHub repository to learn how to 
specify the different options. We will use the command `sed` to find the 
variable names we defined in the template file so we can replace them
with the correct value for each option: 

```sh
# Relative paths from `03_branchsite_models` directory 
# created above

# 1. Copy the template control file that 
# we will be later modifying
cp ../../../templates/template_CODEML.ctl Branchsite_model_chicken/CODEML_2/codeml-branchsite_null.ctl 
cp ../../../templates/template_CODEML.ctl Branchsite_model_duck/CODEML_2/codeml-branchsite_null.ctl 

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

### Running `CODEML`
Now that you have the control file, you only need to run `CODEML`. 
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
../../../../../src/CODEML/codeml4.9j codeml-branchsite_null.ctl  | tee logfile_codeml-branchsite_$name"_null.txt"
# Remove unnecessary files 
rm 2N*
# Go back to home_dir 
cd $home_dir
done
```  

# Likelihood Ratio Test 
First, we need to extract those lines that have `lnL` and then 
remove unnecessary information to keep only the `lnL` value. We will do this for both 
analyses. The files will be saved in `03_branchsite_models`:

```sh
# Run from `03_branchsite_models`
grep 'lnL' Branchsite_model_chicken/CODEML/out_chicken_branchsite.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' > lnL_branchsite_mods.txt
grep 'lnL' Branchsite_model_duck/CODEML/out_duck_branchsite.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> lnL_branchsite_mods.txt
grep 'lnL' Branchsite_model_chicken/CODEML_2/out_chicken_branchsite.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> lnL_branchsite_mods.txt
grep 'lnL' Branchsite_model_duck/CODEML_2/out_duck_branchsite.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> lnL_branchsite_mods.txt
```

The R script `Find_bestmodel.R` to compute LRT can be found in the `03_branchsite_models` directory.
The results are already printed out there in commented lines.

In addition, you can also extract the $\omega$ estimates for each analysis:

```sh 
# Run from `03_branchsite_models`
grep 'MLEs of' -A5 Branchsite_model_chicken/CODEML/out_chicken_branchsite.txt > Branchsite_model_chicken/branchsite_chicken_MLEs.tree
grep 'MLEs of' -A5 Branchsite_model_chicken/CODEML_2/out_chicken_branchsite.txt > Branchsite_model_chicken/branchsite_chicken_MLEs_2.tree
grep 'MLEs of' -A5 Branchsite_model_duck/CODEML/out_duck_branchsite.txt > Branchsite_model_duck/branchsite_duck_MLEs.tree
grep 'MLEs of' -A5 Branchsite_model_duck/CODEML_2/out_duck_branchsite.txt > Branchsite_model_duck/branchsite_duck_MLEs_2.tree
```