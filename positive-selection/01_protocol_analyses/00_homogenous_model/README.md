# Brief introduction: the homogenous model

The simplest codon model implemented in `CODEML` is the so-called `M0`, which assumes that  $\omega$ does not vary across sites or across lineages. All alignment sites of a gene have evolved under the same evolutionary pressure in all taxa.

In this tutorial, we will learn how to run `CODEML` under the `M0` model.

## Inference with `CODEML`

First, we create the main working directory:

```sh
# Run from `00_homogenous_model`
mkdir Model_M0
```

> **NOTE**: If you are a Mac user (i.e., UNIX-based system), you will see that the code snippets below include the command `sed`. By default, this command is different from Linux-based systems, and hence you will not be able to execute them properly. There are two approaches that you can follow, being one easier than the other:
>
>1. (EASY): Instead of running the commands below using the format `sed -i 's/PATTERN/REPLACEMENT/'`, you should include `''` between `-i` and `'s/PATTERN/REPLACEMENT/'`: `sed -i '' 's/PATTERN/REPLACEMENT/'`. Remember to modify the commands in this tutorial accordingly before you paste them on your terminal!
>2. (MORE DIFFICULT): You should install `GNU sed` and establish it as the "standard" `sed` command instead of the one you will have by default. At the time of writing, [this post](https://medium.com/@bramblexu/install-gnu-sed-on-mac-os-and-set-it-as-default-7c17ef1b8f64) is available and has a detailed explanation that you could follow for this purpose. Nevertheless, there are many other tutorials out there that you could follow to achieve the same goal.

### 1. Setting the control file

We use the [template control file](../../templates/template_CODEML.ctl) provided in this GitHub repository ([here](../../templates/). Note that, if we were using a rooted tree, the two branch lengths leading to the root are unidentifiable and only their sum can be properly estimated. Consequently, to avoid this identifiability issue, we will be using an unrooted tree (which also results in reducing the number of parameters to be estimated to one). Once we have located our input files (alignment and tree files) and know which model we want to run, we use the command `sed` to find the variable names defined in the template file so we can replace them with the correct value for each option:

```sh
# Run from `Model_M0` directory 
cd Model_M0 

# 1. Copy the template control file that 
# we will be later modifying
cp ../../../templates/template_CODEML.ctl codeml-M0.ctl 

# 2. Replace variable names with the 
# values needed to run the analysis 

# 2.1. Set path to input files
sed -i 's/ALN/\.\.\/\.\.\/Mx\_aln\.phy/' codeml-M0.ctl 
sed -i 's/TREE/\.\.\/\.\.\/Mx\_unroot\.tree/' codeml-M0.ctl 
sed -i 's/OUT/out\_M0\.txt/' codeml-M0.ctl 

# 2.2. Set data to 1 (only one loci)
sed -i 's/NDAT/1/' codeml-M0.ctl 

# 2.3. Specify all the sites model under which
# we want our data set to be analysed
sed -i 's/CODMOD/0/' codeml-M0.ctl 
sed -i 's/NSSIT/0/' codeml-M0.ctl 
sed -i 's/CODFREQ/7/' codeml-M0.ctl 
sed -i 's/ESTFREQ/0/' codeml-M0.ctl 
sed -i 's/CLOCK/0/' codeml-M0.ctl 

# 2.3. Set the starting values that are to  
# be used when the model parameters are 
# estimated
sed -i 's/FIXOME/0/' codeml-M0.ctl 
sed -i 's/INITOME/0\.5/' codeml-M0.ctl 
```

**NOTE:** In this example, we have not copied the alignment or the tree files in the working directory `Model_M0`. Instead, we have specified the path to these files as you can see in `step 2.1` described in the code snippet above. We have decided to do this so we do not keep several copies of the same input files to carry out the different tests for positive selection in different directories. If you were to run this analysis while having the two input files in the main directory, then you would not need to type the relative path to these files in the control file as shown above (e.g., `../../Mx_aln.phy`) but only the file names (e.g., `Mx_aln.phy`).

### 2. Running `CODEML`

Now that we have the control file, we only need to run `CODEML`. If you want to use the compiled version of `CODEML` provided in this repository, you can use the code provided in the snippet below. Otherwise, please modify the code so you can execute `CODEML` on your PC:

```sh
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
../../../src/CODEML/codeml4.10.6 codeml-M0.ctl | tee logfile_codemlM0.txt

# Remove unnecessary files 
rm 2N*
```

If you want to quickly extract the estimated $\omega$ in an output file, you can do the following:

```sh
printf "omega\n" > omega_est.txt
grep 'omega ' out_M0.txt | sed 's/..*= *//' >> omega_est.txt 
```
