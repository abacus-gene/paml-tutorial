# Brief introduction: branch models

Branch models assume that $\omega$ varies across the branches of the tree ([Yang 1998](https://pubmed.ncbi.nlm.nih.gov/9580986/), [Yang and Nielsen 1998](https://link.springer.com/article/10.1007/PL00006320)). Consequently, we need to specify which branches in our phylogeny are to have different $\omega$ values. In this example, we will run three analyses with `CODEML`: one in which the chicken lineage is labelled as the foreground branch (i.e., branch that will have a different value of $\omega), another in which the foreground branch is the duck lineage instead, another in which these two branches are simultaneously labeled as foreground branches, and a last one in which the branch leading to the outgroup together with the duck and chicken branches are treated as a foreground branches (i.e., we refer to this group of branches as the bird clade).

Note that, for the first the three hypotheses mentioned above (chicken or duck lineages as foreground, both separately and simultaneously), we will need to use an unrooted tree. If we were using a rooted tree, we would have an identifiability issue as the two branch lengths leading to the root cannot be properly estimated, only their sum can. Nevertheless, the hypothesis testing the outgroup branch as a foreground can only be tested using a rooted tree. In other words, one of the two branch lengths leading to the root (the one from the outgroup, duck and chicken) is clearly placed in a rooted tree, but we do not know where exactly in the unrooted tree this branch finishes (i.e., we do not know how close this branch can be from the outgroup or from the inner clade). Consequently, we will be using a rooted tree when testing this hypothesis.

> **NOTE**: If you are a Mac user (i.e., UNIX-based system), you will see that the code snippets below include the command `sed`. By default, this command is different from Linux-based systems, and hence you will not be able to execute them properly. There are two approaches that you can follow, being one easier than the other:
>
>1. (EASY): Instead of running the commands below using the format `sed -i 's/PATTERN/REPLACEMENT/'`, you should include `''` between `-i` and `'s/PATTERN/REPLACEMENT/'`: `sed -i '' 's/PATTERN/REPLACEMENT/'`. Remember to modify the commands in this tutorial accordingly before you paste them on your terminal!
>2. (MORE DIFFICULT): You should install `GNU sed` and establish it as the "standard" `sed` command instead of the one you will have by default. At the time of writing, [this post](https://medium.com/@bramblexu/install-gnu-sed-on-mac-os-and-set-it-as-default-7c17ef1b8f64) is available and has a detailed explanation that you could follow for this purpose. Nevertheless, there are many other tutorials out there that you could follow to achieve the same goal.

Before we get started with the tutorial, we need to generate three copies of the unrooted tree file, [`Mx_unroot.tree`](../Mx_unroot.tree), and one copy of the rooted tree file, [`Mx_root.tree`](../Mx_root.tree), so we can add the tag `#1` to identify the foreground branches in each case:

```sh
# Run from `02_branch_models`
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

## Inference with `CODEML`

Given that we want to test four hypotheses, we create four different directories to run the branch model under the corresponding hypothesis:

```sh
# Run from `02_branch_models`
mkdir -p Branch_model_chicken/CODEML Branch_model_duck/CODEML Branch_model_duckchicken/CODEML Branch_model_bird/CODEML
```

### 1. Setting the control file

We use the [template control file](../../templates/template_CODEML.ctl) provided in this GitHub repository ([here](../../templates/). Then, we use the command `sed` to find the variable names defined in the template file so we can replace them with the correct value for each option:

```sh
# Run from `02_branch_models`

# 1. Copy the template control file that 
# we will be later modifying
cp ../../templates/template_CODEML.ctl Branch_model_chicken/CODEML/codeml-branch.ctl 
cp ../../templates/template_CODEML.ctl Branch_model_duck/CODEML/codeml-branch.ctl 
cp ../../templates/template_CODEML.ctl Branch_model_duckchicken/CODEML/codeml-branch.ctl 
cp ../../templates/template_CODEML.ctl Branch_model_bird/CODEML/codeml-branch.ctl 

# 2. Replace variable names with the 
# values needed to run the analysis 

# 2.1. Set path to input files
sed -i 's/ALN/\.\.\/\.\.\/\.\.\/Mx\_aln\.phy/' Branch_model_chicken/CODEML/codeml-branch.ctl 
sed -i 's/TREE/\.\.\/\.\.\/\.\.\/Mx\_branch\_chicken\.tree/' Branch_model_chicken/CODEML/codeml-branch.ctl 
sed -i 's/OUT/out\_chicken\_branch\.txt/' Branch_model_chicken/CODEML/codeml-branch.ctl 

sed -i 's/ALN/\.\.\/\.\.\/\.\.\/Mx\_aln\.phy/' Branch_model_duck/CODEML/codeml-branch.ctl 
sed -i 's/TREE/\.\.\/\.\.\/\.\.\/Mx\_branch\_duck\.tree/' Branch_model_duck/CODEML/codeml-branch.ctl 
sed -i 's/OUT/out\_duck\_branch\.txt/' Branch_model_duck/CODEML/codeml-branch.ctl 

sed -i 's/ALN/\.\.\/\.\.\/\.\.\/Mx\_aln\.phy/' Branch_model_duckchicken/CODEML/codeml-branch.ctl 
sed -i 's/TREE/\.\.\/\.\.\/\.\.\/Mx\_branch\_duckchicken\.tree/' Branch_model_duckchicken/CODEML/codeml-branch.ctl 
sed -i 's/OUT/out\_duckchicken\_branch\.txt/' Branch_model_duckchicken/CODEML/codeml-branch.ctl 

sed -i 's/ALN/\.\.\/\.\.\/\.\.\/Mx\_aln\.phy/' Branch_model_bird/CODEML/codeml-branch.ctl 
sed -i 's/TREE/\.\.\/\.\.\/\.\.\/Mx\_branch\_bird\.tree/' Branch_model_bird/CODEML/codeml-branch.ctl 
sed -i 's/OUT/out\_bird\_branch\.txt/' Branch_model_bird/CODEML/codeml-branch.ctl 

# 2.2. Go over all control files
for i in Branch_model_duck Branch_model_chicken Branch_model_duckchicken Branch_model_bird
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

**NOTE:** In this example, we have not copied the alignment or the tree files in the working directories (e.g., `Branch_model_chicken`). Instead, we have specified the path to these files as you can see in `step 2.1` described in the code snippet above. We have decided to do this so we do not keep several copies of the same input files to carry out the different tests for positive selection in different directories. If you were to run this analysis while having the two input files in the main directory, then you would not need to type the relative path to these files in the control file as shown above (e.g., `../../../Mx_aln.phy`) but only the file names (e.g., `Mx_aln.phy`).

### 2. Running `CODEML`

Now that we have the control file for each analysis, we only need to run `CODEML`. If you want to use the compiled version of `CODEML` provided in this repository, you can use the code provided in the snippet below. Otherwise, please modify the code so you can execute `CODEML` on your PC:

```sh
# Run from `02_branch_models`
home_dir=$( pwd )
for i in Branch_model_duck Branch_model_chicken Branch_model_duckchicken Branch_model_bird
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
../../../../src/CODEML/codeml4.10.6 codeml-branch.ctl | tee logfile_codeml-branch_$name.txt &
# Go back to home_dir 
cd $home_dir
done
# Remove unnecessary files 
rm */CODEML/2N*
```

## Likelihood Ratio Test (LRT)

Now, we can test if the branch model fits the data better than the null model (`M0` model), which we previously ran (i.e., the results for the `M0` model can be found in the `../00_homogenous_model/Model_M0` directory) for each of our hypotheses.

First, we extract those lines that have the `lnL` term and then remove unnecessary information to keep only the `lnL` value. We will do this for the three analyses which likelihood values we want to compare:

```sh
# Run from `02_branch_models`
grep 'lnL' Branch_model_chicken/CODEML/out_chicken_branch.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' > lnL_branch_mods.txt
grep 'lnL' Branch_model_duck/CODEML/out_duck_branch.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> lnL_branch_mods.txt
grep 'lnL' Branch_model_duckchicken/CODEML/out_duckchicken_branch.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> lnL_branch_mods.txt
grep 'lnL' Branch_model_bird/CODEML/out_bird_branch.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> lnL_branch_mods.txt
grep 'lnL' ../00_homogenous_model/Model_M0/out_M0.txt | sed 's/..*\:\ *//' | sed 's/\ ..*//' >> lnL_branch_mods.txt
```

The R script [`Find_bestmodel.R`](Find_bestmodel.R) is used to compute the LRT. You can also find the results for each test written there as commented lines.

In addition, we can extract the $\omega$ ratios when the chicken lineage is set as the foreground lineage and then when the chicken is the foreground lineage:

```sh
# Run from `02_branch_models`
w_c_back=$( grep 'w (dN/dS) for branches'  Branch_model_chicken/CODEML/out_chicken_branch.txt | sed 's/..*: *//' | sed 's/ ..*//' )
w_c_for=$( grep 'w (dN/dS) for branches'  Branch_model_chicken/CODEML/out_chicken_branch.txt | sed 's/..* //')
w_d_back=$( grep 'w (dN/dS) for branches'  Branch_model_duck/CODEML/out_duck_branch.txt | sed 's/..*: *//' | sed 's/ ..*//' )
w_d_for=$( grep 'w (dN/dS) for branches'  Branch_model_duck/CODEML/out_duck_branch.txt | sed 's/..* //')
w_dc_back=$( grep 'w (dN/dS) for branches'  Branch_model_duckchicken/CODEML/out_duckchicken_branch.txt | sed 's/..*: *//' | sed 's/ ..*//' )
w_dc_for=$( grep 'w (dN/dS) for branches'  Branch_model_duckchicken/CODEML/out_duckchicken_branch.txt | sed 's/..* //')
w_b_back=$( grep 'w (dN/dS) for branches'  Branch_model_bird/CODEML/out_bird_branch.txt | sed 's/..*: *//' | sed 's/ ..*//' )
w_b_for=$( grep 'w (dN/dS) for branches'  Branch_model_bird/CODEML/out_bird_branch.txt | sed 's/..* //')
printf "w_back_chicken\tw_fore_chicken\tw_back_duck\tw_fore_duck\tw_back_duckchicken\tw_fore_duckchicken\tw_back_bird\tw_fore_bird\n" > w_est_branches.tsv 
printf $w_c_back"\t"$w_c_for"\t"$w_d_back"\t"$w_d_for"\t"$w_dc_back"\t"$w_dc_for"\t"$w_b_back"\t"$w_b_for"\n" >> w_est_branches.tsv 
```
