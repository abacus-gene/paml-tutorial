                               * [PATHS TO INPUT/OUTPUT FILES]
      seqfile = ALN            * Sequence data filename
     treefile = TREE           * Tree file name
      outfile = out_codeml.txt * Main result file name
   
                               * [TYPE OF MESSAGES IN LOG FILE / SCREEN]
        noisy = 3              * 0,1,2,3,9: how much rubbish on the screen
      verbose = 1              * 0: concise | 1: detailed | 2: too much

                               * [DATA TYPE]
      seqtype = 1              * 1:codons | 2:AAs | 3:codons-->AAs
        ndata = NDAT           * Number of data sets or loci in the alignment
        icode = 0              * 0: universal code | 1: mammalian mt | 2-10: see below
	  
                               * [DEFINE THE MODEL]
        model = CODMOD         * <<OMEGA VARIES ACROSS LINEAGES>>
		                       * -- MODELS FOR CODON DATA --
		                       * 0: one ω ratio for all branches
					           * 1: one ratio for each branch
					           * 2: two or more dN/dS ratios for branches
                               * -- MODELS FOR AAs or CODON-TRANSLATED AAs --
                               * 0: poisson   | 1: proportional | 2: Empirical | 3: Empirical+F
                               * 6: FromCodon | 7: AAClasses    | 8: REVaa_0   | 9: REVaa(nr=189)
	  NSsites = NSSIT          * <<OMEGA VARIES ACROSS SITES>>
	                           * -- "SITE CODON" MODELS --
	                           * 0: one w  (sites.m)        | 1: neutral (M1a, sites.m) | 2: selection (M2a, sites.m) 
				               * 3: discrete (sites.m)      | 4: freqs                  | 5: gamma         
				               * 6: 2gamma                  | 7: beta (sites.m)         | 8: beta&w (sites.m)
				               * 9: beta&gamma              | 10: beta&gamma+1          | 11: beta&normal>1
				               * 12: 0&2normal>1            | 13: 3normal>0             | 22: like M2a but omega>2 (sites.m) 
                               * 23: Tgamma | 24: Tinvgamma | 25: Tgamma+1 | 26: Tinvgamma+1 
        ncatG = NCATG          * -- NUMBER OF CATEGORIES --
		                       * Number of categories in the omega distribution under the site models.
		                       * It works together with `NSsites` option.
					           * E.g., for M3 model, ncatG could be 2 or 3 (2 or 3 site classes).
		                       * NOTE: For branch-site models, the variable ncatG is ignored by the program
					           *	   since the number of site classes is fixed under both models A and B.
							   
    CodonFreq = CODFREQ        * -- SET EQUILIBRIUM CODON FREQENCIES --
	                           * 0: 1/61 each (std genetic code)
						       * 1: F1X4 (calculated from avg nucleotide freq.
						       * 2: F3X4 (calculated from avg nucleotide freq. at the three codon pos.
						       * 3: codon table (free parameters)
                               * 4: F1x4MG | 5: F3x4MG | 6: FMutSel0 | 7: FMutSel
       aaDist = AADIST         * -- SET WHICH AA DISTANCES ARE TO BE USED --
                               * 0: equal | +: geometric   | -: linear
					           * 1-6: G1974,Miyata,c,p,v,a | 7: AAClasses
        clock = CLOCK          * -- CLOCK -- 
                               * 0: no clock | 1: global clock | 2: local clock
      
   aaRatefile = DATF           * -- SET WHICH AA. DIST. MATRIX IS TO BE USED --
                               * Only used for AA seqs under empirical models (i.e., model set to 2 or 3)
					           * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, 
					           * path to a file with your own distance matrix

                               * [SET MODEL PARAMETERS]
    fix_kappa = 0              * 1: kappa fixed, 0: kappa to be estimated
        kappa = INITKAP        * Initial or fixed kappa
    fix_omega = 0              * 1: omega or omega_1 fixed, 0: estimate
        omega = INITOME        * Initial or fixed omega for codons or codon-based AAs
		                       * any value > 1 if it is to be estimated for posSel tests
    fix_alpha = FIXA           * 0: estimate gamma shape parameter | 1: fix it at alpha
        alpha = ALPHA          * Initial or fixed alpha, 0:infinity (constant rate)
		                       * NOTE: For the sites model, it is best to run the NSsites models 
					           *       instead, which is specified by setting fix_alpha to 1
					           *       setting alpha to 0. If alpha equals to 0 (where 0 means infinity),
					           *       it assumes the model has one rate for all sites, so the 
					           *       gamma-rates model is not used.
       Malpha = 0              * 0: one gamma distribution will be applied across all sites 
	                           *    when subst. rates are assumed to vary from site to site.
				               *    Change this to 1 if a different gamma dist. is to be used 
				   	           *    for each codon position ot fot each gene .
	  fix_rho = FIXR           * 0: estimate rho | 1: fix it at rho
          rho = RHO            * Initial or fixed rho, 0: no correlation

                               * [SET PARTITION MODELS] 
        Mgene = MGENE          * The sequence file will need to have option `G` appended to header
				               * Setup of partition models of codon substituion:
				               *
				               * SEQ.FILE | CTL.FILE   | PARAMETERS.ACROSS.GENES
				               * --------   ----------   ---------------------------------------
				               * No G       Mgene = 0    Everything equal
                               * Option G   Mgene = 0    The same (κ, ω) and π, but different cs
                               *                         (proportional branch lengths)
                               * Option G   Mgene = 2    The same (κ,ω), but different πs and cs
                               * Option G   Mgene = 3    The same π, but different (κ, ω) and cs
                               * Option G   Mgene = 4    Different (κ, ω), πs, and cs
                               * Option G   Mgene = 1    Separate analysis

                               * [ESTIMATE dN AND dS] 
	runmode = RUNMOD           * 0: Uses discrete gamma for ML estimation
                               * -2: pairwise comparisons of prot-coding sequences.
				               * -3: Bayesian method to estimate distance t and dN/dS 
				               *     ratio omega in pairwise comparisons.
				               *     Default gamma priors: t~Γ(1.1,1.1) ω~ Γ(1.1,2.2)
			                   *       runmode = -3 1.1 1.1 1.1 2.2

                               * [MISCELLANEOUS SETTINGS]
        getSE = 0              * 0: don't want them | 1: want S.E.s of estimates
 RateAncestor = 0              * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)
   Small_Diff = 1e-8           * Small value used in the difference approximation of derivatives
                               * in the iteration algorithm
       method = 1              * Optimization method:
	                           * 0: simultaneous | 1: one branch a time
  fix_blength = 0              * What to do if the tree has branch lengths:
                               * 0: ignore b.l.         | -1: random search
					           * 1: use them as initial | 2: fix b.l.	
    cleandata = 0              * Remove sites with ambiguity data (1:yes, 0:no)?
				  


* [GENETIC CODES IMPLEMENTED IN CODEML]
* 0: universal        | 1:mammalian mt.          | 2:yeast mt.,
* 3: mold mt.         | 4: invertebrate mt.      | 5: ciliate nuclear 
* 6: echinoderm mt.   | 7: euplotid mt.          | 8: alternative yeast nu. 
* 9: ascidian mt.     | 10: blepharisma nu.      | 11: Yang's regularized code 
* These codes correspond to transl_table 1 to 11 of GenBank.
