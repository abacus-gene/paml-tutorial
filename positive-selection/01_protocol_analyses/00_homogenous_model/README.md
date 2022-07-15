# Inference with `CODEML`
We will run our analyses under the homogenous, `M0`, model. First, we will create the
main working directory:

```sh
# Run from `00_homogenous_model`
mkdir Model_M0
```

## Setting the control file 
We will use the template provided in this GitHub repository to learn how to 
specify the different options. We will use the command `sed` to find the 
variable names we defined in the template file so we can replace them
with the correct value for each option: 

```sh
# Relative paths from `Model_M0` directory 
# created above
cd Model_M0 

# 1. Copy the template control file that 
# we will be later modifying
cp ../../../../templates/template_CODEML.ctl codeml-M0.ctl 

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

## Running `CODEML`
Now that you have the control file, you only need to run `CODEML`. 
If you want to use the compiled version of `CODEML` provided in 
this repository, you can use the code provided in the snippet below.
Otherwise, please modify the code so you can execute 
`CODEML` on your PC:

```sh
# Execute `CODEML`
../../../../src/CODEML/codeml4.9j codeml-M0.ctl | tee logfile_codemlM0.txt

# Remove unnecessary files 
rm 2N*
```


