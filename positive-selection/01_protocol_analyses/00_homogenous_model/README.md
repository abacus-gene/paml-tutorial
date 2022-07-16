# Brief introduction: the homogenous model
The simplest codon model implemented in `CODEML` is the so-called `M0`, which assumes that 
$\omega$ does not vary across sites or across lineages. All alignment sites of a gene have 
evolved under the same evolutionary pressure in all taxa. 

In this tutorial, we will learn how to run `CODEML` under the `M0` model.

# Inference with `CODEML`
First, we create the main working directory:

```sh
# Run from `00_homogenous_model`
mkdir Model_M0
```

## Setting the control file 
We use the [template control file](../../templates/template_CODEML.ctl) provided in this GitHub 
repository ([here](../../templates/). Then, we use the command `sed` to find the 
variable names defined in the template file so we can replace them
with the correct value for each option:

```sh
# Run from `Model_M0` directory 
cd Model_M0 

# 1. Copy the template control file that 
# we will be later modifying
cp ../../../templates/template_CODEML.ctl codeml-M0.ctl 

# 2. Replace variable names with the 
# values needed to run the analysis 

# 2.1. Set path to input files
sed -i 's/ALN/\.\.\/\.\.\/myxovirus\.aln/' codeml-M0.ctl 
sed -i 's/TREE/\.\.\/\.\.\/myxovirus\.tree/' codeml-M0.ctl 
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

**NOTE:** In this example, we have not copied the alignment or the tree files 
in the working directory `Model_M0`. Instead, we have specified the path to 
these files as you can see in step 2.1 (see code snippet above). We have decided
to do this so we do not keep several copies of the same input files to carry out the different
tests for positive selection in different directories. If you were to run this analysis while having 
the two input files in the main directory, then you would not need to type the relative path to 
these files in the control file as shown above (e.g., `../../myxovirus.aln` or 
`../../myxovirus.tree` in this example) but only the file names (e.g., `myxovirus.aln` or `myxovirus.tree`).

## Running `CODEML`
Now that we have the control file, we only need to run `CODEML`. 
If you want to use the compiled version of `CODEML` provided in 
this repository, you can use the code provided in the snippet below.
Otherwise, please modify the code so you can execute 
`CODEML` on your PC:

```sh
# Execute `CODEML` v4.10.5
../../../src/CODEML/codeml4.10.5 codeml-M0.ctl | tee logfile_codemlM0.txt

# Remove unnecessary files 
rm 2N*
```

If you want to quickly extract the estimated $\omega$ in an output file, you can do the following:

```sh 
printf "omega\n" > omega_est.txt
grep 'omega ' out_M0.txt | sed 's/..*= *//' >> omega_est.txt 
```
