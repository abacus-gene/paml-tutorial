      seqfile = ALN            * Path to the alignment file
     treefile = TREE           * Path to the tree file
      outfile = OUT            * Path to the output file
   
        noisy = 3              * How much rubbish on the screen
      verbose = 1              * More or less detailed report

      seqtype = 1              * Data type
        ndata = NDAT           * Number of data sets or loci
        icode = 0              * Genetic code 
    cleandata = 0              * Remove sites with ambiguity data?
		
        model = CODMOD         * Models for ω varying across lineages
	  NSsites = NSSIT          * Models for ω varying across sites
    CodonFreq = CODFREQ        * Codon frequencies
	  estFreq = ESTFREQ        * Use observed freqs or estimate freqs by ML
        clock = CLOCK          * Clock model
    fix_omega = FIXOME         * Estimate or fix omega
        omega = INITOME        * Initial or fixed omega
