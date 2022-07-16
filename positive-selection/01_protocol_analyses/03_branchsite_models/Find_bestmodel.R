#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------#
# SET ENVIRONMENT #
#-----------------#
library( rstudioapi ) 
# Get the path to current open R script and find main dir "00_Gene_filtering"
path_to_file <- getActiveDocumentContext()$path
wd <- paste( dirname( path_to_file ), "/", sep = "" )
setwd( wd )

#-----------#
# LOAD DATA #
#-----------#
# 1. Load our text file with the lnL values
# Chicken | Duck | Chicken_2 | Duck_2
lnL_vals <- read.table( file = "lnL_branchsite_mods.txt", sep= " ", stringsAsFactors = FALSE, 
                        header = FALSE )
rownames( lnL_vals ) <- c( "BS_chicken", "BS_duck", "BS-w1-chicken", "BS-w1-duck" )

# 2. We can now compute the LRT statistic.

## BS_chicken vs BS-w1-chicken ##
diff_BSvsBSw1_chick <- 2*( lnL_vals[1,] - lnL_vals[3,] )
# diff =  4.108244
pchisq( diff_BSvsBSw1_chick, df = 1, lower.tail=F )
# p-val = 0.04267465 < 0.05
# Barely significantly different

## BS-duck vs BS-w1-duck ##
diff_BSvsBSw1_duck <- 2*( lnL_vals[2,] - lnL_vals[4,] )
# diff = 13.2206
pchisq( diff_BSvsBSw1_duck, df = 1, lower.tail=F )
# p-val = 0.0002768896 < 0.05
# NOTE: The LRT statistic should be compared with the 50:50 mixture of
# point mass 0 and 1,5%2=2.71 and 1,1%2=5.41(Self and Liang 1987). 
# Nevertheless, we use critical values \chi_{1,5%}^2=3.84 and \chi_{1,1%}^2=5.99
# to avoid to guide against violations of model assumption instead of the mixture
# distribution mentioned above
Chisq.crit1 <- 3.84
Chisq.crit2 <- 5.99

# 3. Plot results
par( mfrow = c( 1, 2 ) )
# Chicken
curve( dchisq( x, df = 1 ), from = 0, to =  7 )
abline( v = c( Chisq.crit1, Chisq.crit2, diff_BSvsBSw1_chick ), col = c( "darkgray", "brown", "red" ) )
coords_dev     <- c( 1.35, 1.3 )
coords_pval    <- c( 1.25, 1.2 )
coords_alphac  <- c( 1.28, 1.0 )
coords_alphac2 <- c( 1.28, 0.9 )
text( x = coords_dev[1], y = coords_dev[2],
      labels = expression( atop( paste( '2', Delta, 'l = 4.11', sep = "" ) ) ),
      cex = 1.2, col = "red" )
text( x = coords_alphac[1], y = coords_alphac[2],
      labels = expression( atop( paste( chi["1,0.05"]^"2", "= 3.84", sep = " " ) ) ),
      cex = 1.2, col = "gray" )
text( x = coords_alphac2[1], y = coords_alphac2[2],
      labels = expression( atop( paste( chi["1,0.01"]^"2", "= 5.99", sep = " " ) ) ),
      cex = 1.2, col = "brown" )
text( x = coords_pval[1], y = coords_pval[2],
      labels = expression( atop( ''*italic( pval )*' = 0.04' ) ),
      cex = 1.2, col = "black" )
title( expression( 'A) '*italic(Chicken)*': branch-site model VS branch-site model '*omega*'=1' ) )
# Duck
curve( dchisq( x, df = 1 ), from = 0, to =  15 )
abline( v = c( Chisq.crit1, Chisq.crit2, diff_BSvsBSw1_duck ), col = c( "darkgray", "brown", "red" ) )
coords_dev     <- c( 10.4, 0.85 )
coords_pval    <- c( 10.15, 0.78 )
coords_alphac  <- c( 10, 0.67 )
coords_alphac2 <- c( 10, 0.60 )
text( x = coords_dev[1], y = coords_dev[2],
      labels = expression( atop( paste( '2', Delta, 'l = 13.22', sep = "" ) ) ),
      cex = 1.2, col = "red" )
text( x = coords_alphac[1], y = coords_alphac[2],
      labels = expression( atop( paste( chi["1,0.05"]^"2", "= 3.84", sep = " " ) ) ),
      cex = 1.2, col = "gray" )
text( x = coords_alphac2[1], y = coords_alphac2[2],
      labels = expression( atop( paste( chi["1,0.01"]^"2", "= 5.99", sep = " " ) ) ),
      cex = 1.2, col = "brown" )
text( x = coords_pval[1], y = coords_pval[2],
      labels = expression( atop( ''*italic( pval )*' = 2.7e-4' ) ),
      cex = 1.2, col = "black" )
title( expression( 'B) '*italic(Duck)*': branch-site model VS branch-site model '*omega*'=1' ) )
