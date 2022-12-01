# Brief introduction: branch-site models

Branch-site models assume that $\omega\$ can vary both across sites and across lineages ([Yang and Nielsen 2002](https://link.springer.com/article/10.1007/PL00006320)). In these analyses, the aim is to detect positive selection at specific sites along the so-called _foreground_ branches (i.e., user-specified branches with tags along which positive selection at specific sites will be tested).

As we did with the tutorial to run the branch model, here we will run four analyses: one in which the foreground branch is the chicken lineage, another with the duck lineage, a third one in which both the chicken and the duck lineages are labelled simultaneously as foreground branches, and a last one in which all the branches part of the bird clade are labelled as foreground branches (i.e., duck and chicken branches and the branch from the root to this clade).

> **NOTE**: If you are a Mac user (i.e., UNIX-based system), you will see that the code snippets below include the command `sed`. By default, this command is different from Linux-based systems, and hence you will not be able to execute them properly. There are two approaches that you can follow, being one easier than the other:
>
>1. (EASY): Instead of running the commands below using the format `sed -i 's/PATTERN/REPLACEMENT/'`, you should include `''` between `-i` and `'s/PATTERN/REPLACEMENT/'`: `sed -i '' 's/PATTERN/REPLACEMENT/'`. Remember to modify the commands in this tutorial accordingly before you paste them on your terminal!
>2. (MORE DIFFICULT): You should install `GNU sed` and establish it as the "standard" `sed` command instead of the one you will have by default. At the time of writing, [this post](https://medium.com/@bramblexu/install-gnu-sed-on-mac-os-and-set-it-as-default-7c17ef1b8f64) is available and has a detailed explanation that you could follow for this purpose. Nevertheless, there are many other tutorials out there that you could follow to achieve the same goal.

If you have run the tutorial for branch models [here](../02_branch_models/README.md), you will already know that we will use an unrooted tree when testing the first three hypotheses and a rooted tree when testing the fourth one. In addition, you will already have the four tree files with the tags selecting the foreground branches generated. If not, please run the following code snippet before getting started with the tutorial:

```sh
# Run from `03_branchsite_models`
cp ../Mx_unroot.tree ../Mx_branch_chicken.tree 
cp ../Mx_unroot.tree ../Mx_branch_duck.tree 
cp ../Mx_unroot.tree ../Mx_branch_duckchicken.tree 
cp ../Mx_root.tree ../Mx_branch_bird.tree 

# Add tags 
sed -i 's/Chicken\_Mx/Chicken\_Mx\ \#1/' ../Mx_branch_chicken.tree
sed -i 's/Duck\_Mx/Duck\_Mx\ \#1/' ../Mx_branch_duck.tree
sed -i 's/Chicken\_Mx/Chicken\_Mx\ \#1/' ../Mx_branch_duckchicken.tree
sed -i 's/Duck\_Mx/Duck\_Mx\ \#1/' ../Mx_branch_duckchicken.tree
sed -i 's/Chicken\_Mx/Chicken\_Mx\ \#1/' ../Mx_branch_bird.tree
sed -i 's/Duck\_Mx/Duck\_Mx\ \#1/' ../Mx_branch_bird.tree
sed -i 's/Chicken\_Mx\ \#1)/Chicken\_Mx\ \#1)\ \#1/' ../Mx_branch_bird.tree
```

## 1. Inference with `CODEML`

### Part 1: alternative models

Given that we want to test four hypotheses, we create four different directories to run the branch-site model under the corresponding hypothesis:

```sh
# Run from `03_branchsite_models`
mkdir -p Branchsite_model_chicken/CODEML Branchsite_model_duck/CODEML Branchsite_model_duckchicken/CODEML Branchsite_model_bird/CODEML
```

#### 1.1. Setting the control file

We use the [template control file](../../templates/template_CODEML.ctl) provided in this GitHub repository ([here](../../templates/). Then, we use the command `sed` to find the variable names defined in the template file so we can replace them with the correct value for each option:

```sh
# Relative paths from `03_branchsite_models` directory 
# created above

# 1. Copy the template control file that 
# we will be later modifying
cp ../../templates/template_CODEML.ctl Branchsite_model_chicken/CODEML/codeml-branchsite.ctl 
cp ../../templates/template_CODEML.ctl Branchsite_model_duck/CODEML/codeml-branchsite.ctl 
cp ../../templates/template_CODEML.ctl Branchsite_model_duckchicken/CODEML/codeml-branchsite.ctl 
cp ../../templates/template_CODEML.ctl Branchsite_model_bird/CODEML/codeml-branchsite.ctl 

# 2. Replace variable names with the 
# values needed to run the analysis 

# 2.1. Set path to input files
sed -i 's/ALN/\.\.\/\.\.\/\.\.\/Mx\_aln\.phy/' Branchsite_model_chicken/CODEML/codeml-branchsite.ctl 
sed -i 's/TREE/\.\.\/\.\.\/\.\.\/Mx\_branch\_chicken\.tree/' Branchsite_model_chicken/CODEML/codeml-branchsite.ctl 
sed -i 's/OUT/out\_chicken\_branchsite\.txt/' Branchsite_model_chicken/CODEML/codeml-branchsite.ctl 

sed -i 's/ALN/\.\.\/\.\.\/\.\.\/Mx\_aln\.phy/' Branchsite_model_duck/CODEML/codeml-branchsite.ctl 
sed -i 's/TREE/\.\.\/\.\.\/\.\.\/Mx\_branch\_duck\.tree/' Branchsite_model_duck/CODEML/codeml-branchsite.ctl 
sed -i 's/OUT/out\_duck\_branchsite\.txt/' Branchsite_model_duck/CODEML/codeml-branchsite.ctl 

sed -i 's/ALN/\.\.\/\.\.\/\.\.\/Mx\_aln\.phy/' Branchsite_model_duckchicken/CODEML/codeml-branchsite.ctl 
sed -i 's/TREE/\.\.\/\.\.\/\.\.\/Mx\_branch\_duckchicken\.tree/' Branchsite_model_duckchicken/CODEML/codeml-branchsite.ctl 
sed -i 's/OUT/out\_duckchicken\_branchsite\.txt/' Branchsite_model_duckchicken/CODEML/codeml-branchsite.ctl 

sed -i 's/ALN/\.\.\/\.\.\/\.\.\/Mx\_aln\.phy/' Branchsite_model_bird/CODEML/codeml-branchsite.ctl 
sed -i 's/TREE/\.\.\/\.\.\/\.\.\/Mx\_branch\_bird\.tree/' Branchsite_model_bird/CODEML/codeml-branchsite.ctl 
sed -i 's/OUT/out\_bird\_branchsite\.txt/' Branchsite_model_bird/CODEML/codeml-branchsite.ctl 

# 2.2. Go over all control files
for i in Branchsite_model_duck Branchsite_model_chicken Branchsite_model_duckchicken Branchsite_model_bird
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

**NOTE:** In this example, we have not copied the alignment or the tree files in the working directories (e.g., `Branchsite_model_chicken`). Instead, we have specified the path to these files as you can see in step 2.1 in the code snippet above. We have decided to do this so we do not keep several copies of the same input files to carry out the different tests for positive selection in different directories. Nevertheless, you may copy these two input files in the working directory of this analysis or other analysis you may run for other projects, in which case you do not need to type the relative path to these files in the control file as shown above (e.g., `../../../Mx_aln.phy`) but only the file names (e.g., `Mx_aln.phy`).

#### 1.2. Running `CODEML`

Now that we have the control file, we only need to run `CODEML`. If you want to use the compiled version of `CODEML` provided in this repository, you can use the code provided in the snippet below. Otherwise, please modify the code so you can execute `CODEML` on your PC:

```sh
# Run from `03_branchsite_models`
home_dir=$( pwd )
for i in Branchsite_model_duck Branchsite_model_chicken Branchsite_model_duckchicken Branchsite_model_bird
do
# Move to dir 
cd $i/CODEML/
# Execute `CODEML` v4.10.6
# 1. If you have the latest version of `CODEML` installed on your PC and 
#    you have exported its location to your PATH, please modify the command
#    below accordingly so it is this version of `CODEML` the one you run.
# 2. The link below redirects to an executable file that was compiled with 
#    the Windows Linux Subsystem. If you are using another operating system,
#    the command below will not work. Please modify the command below so you
#    run the version of `CODEML` compatible with your operating system.
#
# Once you know which `CODEML` executable you are running on your PC, copy 
# and paste the command below on your terminal either as it is or after
# having modified the path to the `CODEML` executable you want to run
name=$( echo $i | sed 's/..*\_//' )
../../../../src/CODEML/codeml4.10.6 codeml-branchsite.ctl | tee logfile_codeml-branchsite_$name.txt &
# Go back to home_dir 
cd $home_dir
done
# Remove unnecessary files 
rm */CODEML/2N*
```  

### Part 2: null models

We now run `CODEML` under a modified version of the `M2a` model where the value of $\omega$ for site classes "2a" and "2b" is fixed to $\omega=1$. Consequently, as this parameter is not estimated, this model has one parameter less than the so-called branch-site model A ([Zhang et al. 2005](https://academic.oup.com/mbe/article/22/12/2472/1009544)). We then create one directory for each analysis within the corresponding main directories previously created:

```sh
# Run from `03_branchsite_models`
mkdir -p Branchsite_model_chicken/CODEML_2 Branchsite_model_duck/CODEML_2 Branchsite_model_duckchicken/CODEML_2 Branchsite_model_bird/CODEML_2
```

#### 2.1. Setting the control file

We use the [template control file](../../templates/template_CODEML.ctl) provided in this GitHub repository ([here](../../templates/). Then, we use the command `sed` to find the variable names defined in the template file so we can replace them with the correct value for each option:

```sh
# Relative paths from `03_branchsite_models` directory 
# created above

# 1. Copy the template control file that 
# we will be later modifying
cp ../../templates/template_CODEML.ctl Branchsite_model_chicken/CODEML_2/codeml-branchsite_null.ctl 
cp ../../templates/template_CODEML.ctl Branchsite_model_duck/CODEML_2/codeml-branchsite_null.ctl 
cp ../../templates/template_CODEML.ctl Branchsite_model_duckchicken/CODEML_2/codeml-branchsite_null.ctl 
cp ../../templates/template_CODEML.ctl Branchsite_model_bird/CODEML_2/codeml-branchsite_null.ctl 

# 2. Replace variable names with the 
# values needed to run the analysis 

# 2.1. Set path to input files
sed -i 's/ALN/\.\.\/\.\.\/\.\.\/Mx\_aln\.phy/' Branchsite_model_chicken/CODEML_2/codeml-branchsite_null.ctl 
sed -i 's/TREE/\.\.\/\.\.\/\.\.\/Mx\_branch\_chicken\.tree/' Branchsite_model_chicken/CODEML_2/codeml-branchsite_null.ctl 
sed -i 's/OUT/out\_chicken\_branchsite\.txt/' Branchsite_model_chicken/CODEML_2/codeml-branchsite_null.ctl 

sed -i 's/ALN/\.\.\/\.\.\/\.\.\/Mx\_aln\.phy/' Branchsite_model_duck/CODEML_2/codeml-branchsite_null.ctl 
sed -i 's/TREE/\.\.\/\.\.\/\.\.\/Mx\_branch\_duck\.tree/' Branchsite_model_duck/CODEML_2/codeml-branchsite_null.ctl 
sed -i 's/OUT/out\_duck\_branchsite\.txt/' Branchsite_model_duck/CODEML_2/codeml-branchsite_null.ctl 

sed -i 's/ALN/\.\.\/\.\.\/\.\.\/Mx\_aln\.phy/' Branchsite_model_duckchicken/CODEML_2/codeml-branchsite_null.ctl 
sed -i 's/TREE/\.\.\/\.\.\/\.\.\/Mx\_branch\_duckchicken\.tree/' Branchsite_model_duckchicken/CODEML_2/codeml-branchsite_null.ctl 
sed -i 's/OUT/out\_duckchicken\_branchsite\.txt/' Branchsite_model_duckchicken/CODEML_2/codeml-branchsite_null.ctl

sed -i 's/ALN/\.\.\/\.\.\/\.\.\/Mx\_aln\.phy/' Branchsite_model_bird/CODEML_2/codeml-branchsite_null.ctl 
sed -i 's/TREE/\.\.\/\.\.\/\.\.\/Mx\_branch\_bird\.tree/' Branchsite_model_bird/CODEML_2/codeml-branchsite_null.ctl 
sed -i 's/OUT/out\_bird\_branchsite\.txt/' Branchsite_model_bird/CODEML_2/codeml-branchsite_null.ctl

# 2.2. Go over all control files
for i in Branchsite_model_duck Branchsite_model_chicken Branchsite_model_duckchicken Branchsite_model_bird
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

**NOTE:** In this example, we have not copied the alignment or the tree files in the working directories. Instead, we have specified the path to these files as you can see in `step 2.1` as detailed in the code snippet above. We have decided to do this so we do not keep several copies of the same input files to carry out the different tests for positive selection in different directories. If you were to run this analysis while having the two input files in the main directory, then you would not need to type the relative path to these files in the control file as shown above (e.g., `../../../Mx_aln.phy`) but only the file names (e.g., `Mx_aln.phy`).

#### 2.2. Running `CODEML`

Now that we have the control file, we only need to run `CODEML`. If you want to use the compiled version of `CODEML` provided in this repository, you can use the code provided in the snippet below. Otherwise, please modify the code so you can execute `CODEML` on your PC:

```sh
# Run from `03_branchsite_models`
home_dir=$( pwd )
for i in Branchsite_model_duck Branchsite_model_chicken Branchsite_model_duckchicken Branchsite_model_bird
do
# Move to dir 
cd $i/CODEML_2/
# Execute `CODEML`
name=$( echo $i | sed 's/..*\_//' )
../../../../src/CODEML/codeml4.10.6 codeml-branchsite_null.ctl  | tee logfile_codeml-branchsite_$name"_null.txt" &
# Go back to home_dir 
cd $home_dir
done
# Remove unnecessary files 
rm */CODEML_2/2N*
```  

## 2. Likelihood Ratio Test (LRT)

Now, we can test whether the branch-site model A fits the data better than its modified version with one less parameter.

First, we extract those lines that have the `lnL` term and then remove unnecessary information to keep only the `lnL` value. We will do this for the analyses which likelihood values we want to compare:

```sh
# Run from `03_branchsite_models`
grep 'lnL' Branchsite_model_chicken/CODEML/out_chicken_branchsite.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' > lnL_branchsite_mods.txt
grep 'lnL' Branchsite_model_duck/CODEML/out_duck_branchsite.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> lnL_branchsite_mods.txt
grep 'lnL' Branchsite_model_duckchicken/CODEML/out_duckchicken_branchsite.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> lnL_branchsite_mods.txt
grep 'lnL' Branchsite_model_bird/CODEML/out_bird_branchsite.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> lnL_branchsite_mods.txt
grep 'lnL' Branchsite_model_chicken/CODEML_2/out_chicken_branchsite.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> lnL_branchsite_mods.txt
grep 'lnL' Branchsite_model_duck/CODEML_2/out_duck_branchsite.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> lnL_branchsite_mods.txt
grep 'lnL' Branchsite_model_duckchicken/CODEML_2/out_duckchicken_branchsite.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> lnL_branchsite_mods.txt
grep 'lnL' Branchsite_model_bird/CODEML_2/out_bird_branchsite.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> lnL_branchsite_mods.txt
```

The R script [`Find_bestmodel.R`](Find_bestmodel.R) is then used to compute LRT. The results are already written in the script as commented lines.

In addition, you can also extract the $\omega$ estimates for each analysis:

```sh
# Run from `03_branchsite_models`
grep 'MLEs of' -A5 Branchsite_model_chicken/CODEML/out_chicken_branchsite.txt > Branchsite_model_chicken/branchsite_chicken_MLEs.txt
grep 'MLEs of' -A5 Branchsite_model_chicken/CODEML_2/out_chicken_branchsite.txt > Branchsite_model_chicken/branchsite_chicken_MLEs_2.txt
grep 'MLEs of' -A5 Branchsite_model_duck/CODEML/out_duck_branchsite.txt > Branchsite_model_duck/branchsite_duck_MLEs.txt
grep 'MLEs of' -A5 Branchsite_model_duck/CODEML_2/out_duck_branchsite.txt > Branchsite_model_duck/branchsite_duck_MLEs_2.txt
grep 'MLEs of' -A5 Branchsite_model_duckchicken/CODEML/out_duckchicken_branchsite.txt > Branchsite_model_duckchicken/branchsite_duckchicken_MLEs.txt
grep 'MLEs of' -A5 Branchsite_model_duckchicken/CODEML_2/out_duckchicken_branchsite.txt > Branchsite_model_duckchicken/branchsite_duckchicken_MLEs_2.txt
grep 'MLEs of' -A5 Branchsite_model_bird/CODEML/out_bird_branchsite.txt > Branchsite_model_bird/branchsite_outgroup_MLEs.txt
grep 'MLEs of' -A5 Branchsite_model_bird/CODEML_2/out_bird_branchsite.txt > Branchsite_model_bird/branchsite_outgroup_MLEs_2.txt
```
