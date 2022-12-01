# Analyses with more than one gene

The different analyses that we have discussed in the protocol (see the step-by-step tutorial with the data and code [here](../01_protocol_analyses/)) are based on detecting positive selection when there is only one gene. Note that, when using one gene, there is only one tree topology associated with this gene (e.g., either the corresponding gene tree or the species tree). In other words, all species that are part of the tree topology have available sequence data in the corresponding gene alignments (i.e., there are no taxa for which the sequence is all gaps). Sometimes, however, you may have more than one gene alignment for the taxa that you are interested in to test for positive selection. In that case, some of these gene alignments might even have some missing taxa (i.e., that gene is not available for one or more species in your gene alignment). In that case, there is not a unique tree topology that can be used for all the gene alignments as it happened when having only one gene alignment. In addition, the more gene alignments you have, the more time it takes to generate individual pruned tree topologies for each gene alignment. Now, the new version of `CODEML` (PAML v4.10.6, you can download the latest version from the PAML GitHub repository, [here](https://github.com/abacus-gene/paml)) has a new implemented feature that can help you with this type of analysis!

Below, we will describe the four scenarios that you may encounter when running analyses with `CODEML`.

## Case **a**

This scenario is the implementation that has been already available in `CODEML` in previous releases:

* You have more than one gene alignment and no missing taxa.
* As there are no missing taxa, the tree topology is the same for all the gene alignments.

[Here](case_a), you have an example dataset with the following files:

1. [Control file to execute `CODEML`](case_a/case_a.ctl): Note that we are using the homogenous model to speed up the analyses. Nevertheless, when using your own data, you should change the options according to the positive selection test that you want to run. The most important option in the control file that you need to remember to enable this scenario is **`ndata`**. In **case a**, you only need to type the amount of blocks of gene alignments included in the alignment file (e.g., three in this example, `ndata = 3`).
2. [Alignment file](case_a/case_a.phy): As many alignment blocks as available gene alignments. In this case, there are 3 blocks (i.e., 3 gene alignments) and each block has sequence data for all taxa (no missing taxa). Note that the gene alignments provided as part of this alignment file have no biological meaning: this file is used as an example to show how this case runs quick in `CODEML`.
3. [Tree file](case_a/case_a.tree): There is only one tree topology as there are no missing taxa in the gene alignments. We refer to the PHYLIP header and the tree topology in Newick file (i.e., the content of this tree file) as a unique "tree block".

Now, you can execute `CODEML` using the command below:

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
cd case_a
../../src/CODEML/codeml4.10.6 case_a.ctl | tee logfile_caseA.txt
```

Once `CODEML` finishes, you will see that, for each gene, you will have the corresponding results in the [output file](case_a/out_caseA_M0.txt) separated by a header that starts with `Dataset`. In this example, you will have three headers: `Dataset 1`, `Dataset 2`, and `Dataset 3`. Nevertheless, you would have as many sections with the corresponding headers as alignment blocks in the alignment file.

Note that, depending on the test for positive selection you run and the options you enable for this purpose, the output file will change. Please read the protocol for more information about what to expect under each type of analysis and run the examples using the data available in the corresponding directories in this GitHub repository:

* [Homogenous model](../01_protocol_analyses/00_homogenous_model/README.md)
* [Site models](../01_protocol_analyses/01_site_models/README.md)
* [Branch model](../01_protocol_analyses/02_branch_models/README.md)
* [Branch-site model](../01_protocol_analyses/03_branchsite_models/README.md)

## Case **b**

This scenario is part of the new features implemented in `CODEML`. You should enable this case if your dataset consists of the following:

* You have more than one gene alignment and some of these genes have missing taxa.
* You have already generated one individual gene tree for each of your gene alignments (e.g., in this example, for each of the three gene alignments there is a matching tree topology).

[Here](case_b), you have an example dataset with the following files:

1. [Control file to execute `CODEML`](case_b/case_b.ctl): Note that we are using the homogenous model to speed up the analyses. Nevertheless, when using your own data, you should change the options according to the positive selection test that you want to run. The most important option in the control file that you need to remember to enable this scenario is **`ndata`**. In **case b**, you need to type the amount of blocks of gene alignments included in the alignment file followed by `separate_trees` (e.g., there are three gene alignments in this example, `ndata = 3 separate_trees`). The second argument, `separate_trees`, indicates that the tree file will have one tree topology (i.e., tree block) for each of the gene alignments in your alignment file (i.e., alignment block).
2. [Alignment file](case_b/case_b.phy): As many alignment blocks as available gene alignments. In this case, there are 3 blocks (i.e, 3 gene alignments) and the last two blocks have missing data (i.e., the sequence data for `S1` are missing in the second gene alignment and the sequence data for `S1` and `S3` are missing in the third gene alignment). Note that the gene alignments provided as part of this alignment file have no biological meaning: this file is used as an example to show how this case runs quick in `CODEML`.
3. [Tree file](case_b/case_b.tree): There is **one tree topology for each gene tree** and they are written in the same order as the alignment blocks are found in the alignment file (i.e., the first tree topology will correspond to the first gene alignment, and so on). In total, there are three tree blocks (i.e., PHYLIP header + Tree in Newick format) in this file as there are three gene alignments.

Now, you can execute `CODEML` using the command below:

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
cd case_b
../../src/CODEML/codeml4.10.6 case_b.ctl | tee logfile_caseB.txt
```

Once `CODEML` finishes, you will see that, for each gene, you will have the corresponding results in the [output file](case_b/out_caseB_M0.txt) separated by a header that starts with `Dataset`. In this example, you will have three headers: `Dataset 1`, `Dataset 2`, and `Dataset 3`. Nevertheless, you would have as many sections with the corresponding headers as alignment blocks in the alignment file.

Note that, depending on the test for positive selection you run and the options you enable for this purpose, the output file will change. Please read the protocol for more information about what to expect under each type of analysis and run the examples using the data available in the corresponding directories in this GitHub repository:

* [Homogenous model](../01_protocol_analyses/00_homogenous_model/README.md)
* [Site models](../01_protocol_analyses/01_site_models/README.md)
* [Branch model](../01_protocol_analyses/02_branch_models/README.md)
* [Branch-site model](../01_protocol_analyses/03_branchsite_models/README.md)

> **IMPORTANT**: Note that, if you run the branch or the branch-site model (i.e., you select foreground branches with the tag `#1`), you should not include gene alignments with missing taxa which tree topology is incompatible with these model. For instance, you should always check the following:
>
> 1. The order in which you write the tree topologies in the tree file corresponds to the order in which the gene alignments appear in the alignment file. In other words, the position of the gene alignment matches the position of the corresponding gene tree. Otherwise, `CODEML` will not run.
> 2. Branch and branch-site models need at least two branch types in the tree to be executed. For instance, let's imagine that the tree with all taxa is `((S1,S2 #1),S3,S4);`. If one gene alignment did not have sequence data for taxa `S1` and `S3`, the resulting tree would be `(S2, S4);`. In that case, the pruned tree does not have the foreground anymore because it is pruned when taxa `S1` and `S3` are not included anymore. If a tree like that were found in the tree file when running a branch or a branch-site model, `CODEML` will not run. If you have your gene trees ready but are unsure whether the `#1` is lost or kept after pruning, please read sections `Case c` and `Case d` in this `README.md` file to learn what you should do to avoid any mistakes!

## Case **c**

This scenario is part of the new features implemented in `CODEML`. You should enable this case if your dataset consists of the following:

* You have more than one gene alignment and some of these genes have missing taxa.
* You **do not have** an individual gene tree for each of your gene alignments, you only have the species tree (e.g., in this example, the species tree has four taxa but, some of the gene trees, have missing taxa).

[Here](case_c), you have an example dataset with the following files:

1. [Control file to execute `CODEML`](case_c/case_c.ctl): Note that we are using the homogenous model to speed up the analyses. Nevertheless, when using your own data, you should change the options according to the positive selection test that you want to run. The most important option in the control file that you need to remember to enable this scenario is **`ndata`**. In **case c**, you need to type the amount of blocks of gene alignments included in the alignment file followed by `maintree` or `maintree 1` (e.g., there are three gene alignments in this example, `ndata = 3 maintree` or `ndata = 3 maintree 1`). The second argument, `maintree`, indicates that the tree file will have a unique tree block (i.e., PHYLIP header + Tree in Newick format) that would correspond to the species tree. This tree topology is the one that `CODEML` will use as a "main" tree and, when a gene alignment that has missing taxa is read, `CODEML` will prune those missing taxa and generate the corresponding tree topology (i.e., you do not need to have individual gene trees, `CODEML` will do this for you). Note that this control file includes `ndata = 3 maintree` as the third argument, `1`, is the default behaviour and hence not needed. Nevertheless, you can also include the third argument if you want to.
2. [Alignment file](case_c/case_c.phy): As many alignment blocks as available gene alignments. In this case, there are 3 blocks (i.e., 3 gene alignments) and the last two blocks have missing data (i.e., the sequence data for `S1` are missing in the second gene alignment and the sequence data for `S1` and `S3` are missing in the third gene alignment). Note that the gene alignments provided as part of this alignment file have no biological meaning: this file is used as an example to show how this case runs quick in `CODEML`.
3. [Tree file](case_c/case_c.tree): There is only one tree block with one tree topology (the "main" tree topology), which would correspond to the species tree. Whenever a gene alignment with missing taxa is read with `CODEML`, (i) the "main" tree will be pruned accordingly, (ii) the resulting pruned tree will be used by `CODEML` to run the analysis for this gene alignment.

Now, you can execute `CODEML` using the command below:

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
cd case_c
../../src/CODEML/codeml4.10.6 case_c.ctl | tee logfile_caseC.txt
```

Once `CODEML` finishes, you will see that, for each gene, you will have the corresponding results in the [output file](case_c/out_caseC_M0.txt) separated by a header that starts with `Dataset`. In this example, you will have three headers: `Dataset 1`, `Dataset 2`, and `Dataset 3`. Nevertheless, you would have as many sections with the corresponding headers as alignment blocks in the alignment file.

Note that, depending on the test for positive selection you run and the options you enable for this purpose, the output file will change. Please read the protocol for more information about what to expect under each type of analysis and run the examples using the data available in the corresponding directories in this GitHub repository:

* [Homogenous model](../01_protocol_analyses/00_homogenous_model/README.md)
* [Site models](../01_protocol_analyses/01_site_models/README.md)
* [Branch model](../01_protocol_analyses/02_branch_models/README.md)
* [Branch-site model](../01_protocol_analyses/03_branchsite_models/README.md)

> **IMPORTANT**: Note that, if you run the branch or the branch-site model (i.e., you select foreground branches with the tag `#1`), you should not include gene alignments with missing taxa which tree topology is incompatible with these model. While sometimes this situation can be obvious and it is easy to manually check which tree topologies will be generated for each gene alignment, it gets very time consuming with more taxa and more genes. Under `case c`, `CODEML` will automatically prune the "main" tree you have included in the tree file when a gene tree has missing taxa. If you have enabled a branch or a branch-site model and the pruned tree is not compatible with the branch or the branch-site model (i.e., the `#1` that selects the foreground branch/es has been lost after pruning and there is no `#1` left), `CODEML` will stop. If you want to avoid this from happening, you need to read the details about `case d`.
>
> There is another important outcome of the pruning that you should consider. If you have more than one label (e.g., several branches have been labelled with `#1` or you have different `#X` tags), some labels could disappear after pruning but others may remain. In this case, as there is at least one branch marked as foreground, `CODEML` will run. You **must always check** the results you obtain with gene alignments with missing taxa: while some labels may remain (and hence your hypothesis), others may disappear (and affect your main hypothesis).

## Case **d**

This scenario is part of the new features implemented in `CODEML`. You should enable this case if your dataset consists of the following:

* You have more than one gene alignment and some of these genes have missing taxa.
* You **do not have** an individual gene tree for each of your gene alignments, you only have the species tree (e.g., in this example, the species tree has four taxa but, some of the gene trees, have missing taxa).
* You want to get a tree file with as many tree blocks that have been pruned according to the missing taxa in the corresponding gene alignment in the alignment file. You may want to do this because you want to run a branch or a branch-site model and you are not sure whether the foreground branches are kept after pruning.

[Here](case_d), you have an example dataset with the following files:

1. [Control file to execute `CODEML`](case_d/case_d.ctl): Note that, despite we are using the homogenous model, this option does not matter under `case d`. The most important option in the control file that you need to remember to enable this scenario is **`ndata`**. In `case d`, you need to type the amount of blocks of gene alignments included in the alignment file followed by `maintree 0` (e.g., there are three gene alignments in this example, `ndata = 3 maintree 0`). The second argument, `maintree`, indicates that the tree file will have a unique tree block (i.e., PHYLIP header + Tree in Newick format) that would correspond to the species tree. The second argument, `0`, tells `CODEML` that only a tree file with the individual pruned tree topologies for all the gene alignments in the alignment file should be generated (i.e., `CODEML` will not run, only an output tree file compatible with `CODEML` will be generated). Please note that tree topology included in the tree file is the one that `CODEML` will use as a "main" tree and, when a gene alignment that has missing taxa is read, `CODEML` will prune those missing taxa and generate the corresponding tree topology.
2. [Alignment file](case_d/case_d.phy): As many alignment blocks as available gene alignments. In this case, there are 3 blocks (i.e., 3 gene alignments) and the last two blocks have missing data (i.e., the sequence data for `S1` are missing in the second gene alignment and the sequence data for `S1` and `S3` are missing in the third gene alignment). Note that the gene alignments provided as part of this alignment file have no biological meaning: this file is used as an example to show how this case runs quick in `CODEML`.
3. [Tree file](case_d/case_d.tree): There is only one tree block with one tree topology (the "main" tree topology), which would correspond to the species tree. Whenever a gene alignment with missing taxa is read with `CODEML`, (i) the "main" tree will be pruned accordingly, (ii) the tree topology will be included in the tree file that will be output by `CODEML` (i.e., `genetree.trees`).

Now, you can execute `CODEML` using the command below:

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
cd case_d
../../src/CODEML/codeml4.10.6 case_d.ctl | tee logfile_caseD.txt
```

Once `CODEML` finishes, you will see that an output file called [`genetrees.trees`](case_d/genetrees.trees) will be generated. This is the tree file that will include as many tree blocks as alignment blocks, with the corresponding pruned branches according to the missing taxa. Note that, if branches had been labelled as foreground but they are gone after pruning, the foreground branches will not be included. If that occurs, you should check whether the pruning has affected your hypothesis for the branch or the branch-site model you wanted to carry out before running `CODEML`. If there are no incompatibilities with any of the models you want to run with `CODEML`, you should then run `CODEML` under `case b` (i.e., you now have a tree file with the matching tree topologies for each gene alignment).
