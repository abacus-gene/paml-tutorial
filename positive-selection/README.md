# Beginner's guide on the use of PAML to detect positive selection

## Introduction

The ratio between nonsynonymous and synonymous substitution rates, $\omega=dN/dS$, has been widely used to measure the effect of natural selection on protein-coding genes ([Kimura 1977](https://www.nature.com/articles/267275a0), [Miyaga & Yasunaga 1980](https://link.springer.com/article/10.1007/BF01732067)). Synonymous mutations (also known as "silent" mutations) do not change the amino acid, and synonymous substitution rate may be expected to equal the neutral mutation rate. Nonsynonymous mutations (also known as replacement mutations) change the amino acid, and may often be deleterious and purged from the population by natural selection, resulting in reduced nonsynonymous substitution rate. The $\omega$ ratio is a measure of such selection on nonsynonymous mutations. Specifically, $\omega$ can be used to identify genes ...

* ... under positive (or diversifying) selection (or adaptive evolution), if $\omega>1$;
* ... under neutral evolution, if $\omega=1$;
* ... under negative (or purifying) selection, if $\omega<1$.

As every functional protein has some structural constraints, the $\omega$ ratio average over the whole protein sequence is often less than 1, even if positive selection operates at some sites or over certain time intervals. Thus, statistical tests have been designed to detect positive selection that targets only a small subset of the amino acid residues or affect a limited time interval, leading to the so-called sites test and branch test.

## What can I find in this repository?

In this repository, you will find the data, the code, and step-by-step guidelines for reproducing the results in the `CODEML` protocol ([Álvarez-Carretero, et al. 2023](https://academic.oup.com/mbe/article/40/4/msad041/7140562?searchresult=1&login=true)). We performed all positive selection analyses with `CODEML`, in the [PAML v4.10.6 package](https://github.com/abacus-gene/paml) ([Yang 2007](https://doi.org/10.1093/molbev/msm088)).

We use the alignment and tree files for the myxovirus gene sequences from ten mammal species and two birds (outgroup) analysed by Huo et al. ([2007](https://pubmed.ncbi.nlm.nih.gov/17467195/)).
[In direcyory `00_data`](00_data), we explain how we downloaded and parsed these sequences before we generated the alignment and the gene tree. Then, we carried out different tests for positive selection under the following models:

* **Homogenous model**: all alignment sites and taxa have evolved under the same evolutionary pressure. This model, also known as `M0` model, assumes that $\omega$ is constant across all sites and lineages.
* **Site models** assume that different (amino acid or codon) sites are under different selective pressures and have different $\omega$ values. Positive selection is detected when a subset of sites in the protein-coding gene have $\omega>1$.
* **Branch models** assume that $\omega$ varies among branches of the phylogeny and positive selection is detected along specific lineages if $\omega$ for the branches is $>1$.
* **Branch-site models** assume that $\omega$ varies among branches of the phylogeny and across sites of the gene, and positive selection is detected if a subset of sites for specific branches of the phylogeny have $\omega>1$.

[In directory `01_protocol_analyses`](01_protocol_analyses), you can find one directory for each of the analyses mentioned above, with the corresponding `README.md` file. All code snippets and explanations needed to run `CODEML` under each scenario according to our protocol ([Álvarez-Carretero et al. 2023](https://academic.oup.com/mbe/article/40/4/msad041/7140562?searchresult=1&login=true)) are provided. [In directory `02_extra_analyses`](02_extra_analyses), you can find one directory for each of the new analyses that can be carried out with `CODEML` when (i) there are several genes that need to be analysed at once and (ii) some of these genes may have missing taxa.

We hope that the protocol will be useful for illustrating the control-file settings and interpretations of the program outputs, enabling you to apply similar analyses to your own data.

## References

If you use any of the code we provide in this GitHub repository or consult the protocol for your own analyses, please cite:

* [Álvarez-Carretero S, Kapli P, Yang Z. 2023. Beginner's guide on the use of PAML to detect positive selection, _Mol Biol Evol_, 40(4):msad041](https://academic.oup.com/mbe/article/40/4/msad041/7140562?searchresult=1&login=true).
* [Yang Z. 2007. PAML 4: Phylogenetic analysis by maximum likelihood. _Mol Biol Evol_, 24(8):1586-1591](https://doi.org/10.1093/molbev/msm088).
