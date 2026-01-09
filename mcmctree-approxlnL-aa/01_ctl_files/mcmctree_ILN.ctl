* [[ INPUT/OUTPUT FILES AND SEED NUMBER ]]

          seed = -1                                           * Seed number, type -1 to generate a random seed number
       seqfile = ../../../../../00_inp_data/mtCDNApri.phy     * Path to the input sequence file
      treefile = ../../../../../00_inp_data/mtCDNApri.trees   * Path to the input tree file
      mcmcfile = mcmc.txt                                     * Path to the output file with MCMC samples
       outfile = out.txt                                      * Path to the log output file

* [[ DATA FORMAT AND LIKELIHOOD CALCULATION ]]

         ndata = 1                              * Number of partitions (alignment blocks) in sequence file
       seqtype = 2                              * 0: nucleotides; 1: codons; 2: AAs
       usedata = 2 ../../../../00_CODEML/in.BV  * 0: no data (prior);
                                                * 1: exact likelihood;
                                                * 2 <path_to_inBV>: enable approximate likelihood, add
                                                *                   absolute/relative path to "in.BV"
                                                *                   file as a second argument;
                                                * 3: out.BV (generate "in.BV" file)
     cleandata = 0                              * remove sites with ambiguity data (1:yes, 0:no)?

* [[ EVOLUTIONARY MODEL ]]

         model = 3                                   * [ NUC DATA ]
                                                     * 0: JC69; 1: K80; 2: F81; 3: F84; 4: HKY85
                                                     * [ AA DATA ]
                                                     * 0: poisson; 1: proportional; 2: Empirical; 3:Empirical+F
                                                     * 6: FromCodon; 8: REVaa_0; 9:REVaa(nr=189)
         alpha = 0.5                                 * Value for parameter alpha for gamma rates at sites
         ncatG = 5                                   * Number of categories in discrete gamma
         clock = 2                                   * 1: global clock; 2: independent rates; 3: correlated rates
    aaRatefile = ../../../../../00_inp_data/lg.dat   * Path to the file with LG substitution matrix
                                                     * Find available matrices under "dat" folder or
                                                     * visit the GitHub repository:
                                                     * https://github.com/abacus-gene/paml/tree/master/dat

* [[ PRIORS ]]

       BDparas = 1 1 0.1 m  * [BIRTH-DEATH PROCESS WITH SPECIES SAMPLING]
                            * <lambda> (birth rate);
                            * <mu> (death rate);
                            * <rho> (sampling fraction for extant taxa);
                            * <psi> (sampling fraction for extinct taxa, only add if enabling `TipDating`) 
                            * The last argument enables the conditional or the multiplicative construction:
                            *   c: conditional construction (default until PAML v4.10.7)
                            *   m: multiplicative construction (availble from PAML v4.10.8)
                            * There is no default, users must specify either "c" (or "C") to enable
                            * the conditional construction or "m" (or "M") to enable the multiplicative
                            * construction
   rgene_gamma = 2 20       * [RATE PRIOR]
                            * <shape> <scale> for gammaDir prior for locus rates mu_{i} (rates for genes)
  sigma2_gamma = 1 10       * [PRIOR ON THE VARIANCE IN THE LOGARITHM OF THE RATE]
                            * <shape> <scale> for gamma prior for sigma_{i}^2, enabled only
                            *  if "clock" equals 2 or 3

* [[ MCMC SETTINGS ]]

         print = 1        * 0: no mcmc sample
                          * 1: everything except branch rates
                          * 2: everything
                          * -1: do not run MCMC, summarise samples from file
                          *     specified in "mcmcfile" and map them onto tree file
                          *     specified in "treefile"; an alignment file, while 
                          *     not used, is required (e.g., dummy alignment)
        burnin = 1000     * Number of iterations that will be discarding as part of
                          * the burn-in phase; these samples are not saved in the
                          * specified output file with variable "mcmcfile"
      sampfreq = 100      * Sampling frequency
       nsample = 20000    * Number of samples that you want to collect in the specified
                          * output file with variable "mcmcfile"
                          * The total number of iterations is calculated as follows: 
                          *   num_iter = burnin + sampfreq x nsample