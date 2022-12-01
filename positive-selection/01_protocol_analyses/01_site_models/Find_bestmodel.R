#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------#
# SET ENVIRONMENT #
#-----------------#
library( rstudioapi ) 
# Get the path to current open R script and find main dir
path_to_file <- getActiveDocumentContext()$path
wd <- paste( dirname( path_to_file ), "/", sep = "" )
setwd( wd )

#-----------#
# LOAD DATA #
#-----------#
# 1. Load our text file with the lnL values
lnL_vals <- read.table( file = "Site_models/lnL_sites.txt", sep= " ", stringsAsFactors = FALSE, 
                        header = TRUE )

# 2. We can now compute the LRT statistic, but this can 
# only be done between models that are nested:
#     M0 < M1a < M2a 
#     M7 < M8 

# First, we will compare M1a (alternative) to the M0 (null). 
# If the null is rejected, then we can compare M1a (null) to 
# M2a (alternative).
## M0 vs M1a ##
diff_M0vsM1a <- 2*( lnL_vals$Model_1a[1] - lnL_vals$Model_0[1] )
# diff = 559.2598
pchisq( diff_M0vsM1a, df = 1, lower.tail=F )
# p-val = 1.217966e-123 < 0.05
Chisq.crit.M0vsM1a <- qchisq( p = 0.95, df = 1 )
# alpha critical value at 5% = 3.841459
Chisq.crit.M0vsM1a_2 <- qchisq( p = 0.99, df = 1 )
# alpha critical value at 1% = 6.634897

# As M1a is a better fit to the data than M0, we can compare M1a 
# (Nearly Neutral) against M2a (Positive Selection). 
## M1 vs M2a ##
diff_M1avsM2a <- 2*( lnL_vals$Model_2a[1] - lnL_vals$Model_1a[1] )
# diff = 0
pchisq( diff_M1avsM2a, df = 2, lower.tail=F )
# p-val = 1 > 0
Chisq.crit.M1vsM2a <- qchisq( p = 0.95, df = 2 )
# alpha critical value at 5% level = 5.991465
Chisq.crit.M1vsM2a_2 <- qchisq( p = 0.99, df = 2 )
# alpha critical value at 1% level = 9.21034

# In addition, we can run an additional comparison between 
# M7 (beta) and M8 (beta&omega).
## M7 vs M8 ##
diff_M7vsM8 <- 2*( lnL_vals$Model_8[1] - lnL_vals$Model_7[1] )
# 12.5435
pchisq( diff_M7vsM8, df = 2, lower.tail=F )
# p-val = 0.001888922 < 0.05
Chisq.crit.M7vsM8 <- qchisq( p = 0.95, df = 2 )
# alpha critical value at 5% level = 5.991465
Chisq.crit.M7vsM8_2 <- qchisq( p = 0.99, df = 2 )
# alpha critical value at 1% level = 9.21034

# 3. Plot results 
par( mfrow = c(1, 3) )

# M0 vs M1a
curve( dchisq( x, df = 1 ), from = 0, to =  595 )
abline( v = c( Chisq.crit.M0vsM1a, Chisq.crit.M0vsM1a_2, diff_M0vsM1a ), col = c( "darkgray", "brown", "red" ) )
coords_dev    <- c( 407, 0.008 )
coords_pval   <- c( 410.4, 0.0075 )
coords_alphac <- c( 390, 0.0068 )
coords_alphac2 <- c( 390, 0.0064 )
text( x = coords_dev[1], y = coords_dev[2],
      labels = expression( atop( paste( '2', Delta, 'l = 559.26', sep = "" ) ) ),
      cex = 1.2, col = "red" )
text( x = coords_alphac[1], y = coords_alphac[2],
      labels = expression( atop( paste( chi["1,0.05"]^"2", "= 3.84", sep = " " ) ) ),
      cex = 1.2, col = "gray" )
text( x = coords_alphac2[1], y = coords_alphac2[2],
      labels = expression( atop( paste( chi["1,0.01"]^"2", "= 6.63", sep = " " ) ) ),
      cex = 1.2, col = "brown" )
text( x = coords_pval[1], y = coords_pval[2],
      labels = expression( atop( ''*italic( pval )*' = 1.22e-123' ) ),
      cex = 1.2, col = "black" )
title( "A) M0 vs M1a" )

#M1a vs M2a 
curve( dchisq( x, df = 2 ), from = 0, to =  10 )
abline( v = c( Chisq.crit.M1vsM2a, Chisq.crit.M1vsM2a_2, diff_M1avsM2a ), col = c( "darkgray", "brown", "red" ) )
coords_dev     <- c( 7.85, 0.48 )
coords_pval    <- c( 7.8, 0.45 )
coords_alphac  <- c( 8, 0.41 )
coords_alphac2 <- c( 8, 0.38 )
text( x = coords_dev[1], y = coords_dev[2],
      labels = expression( atop( paste( '2', Delta, 'l = 0', sep = "" ) ) ),
      cex = 1.2, col = "red" )
text( x = coords_alphac[1], y = coords_alphac[2],
      labels = expression( atop( paste( chi["2,0.05"]^"2", "= 5.99", sep = " " ) ) ),
      cex = 1.2, col = "gray" )
text( x = coords_alphac2[1], y = coords_alphac2[2],
      labels = expression( atop( paste( chi["2,0.01"]^"2", "= 9.21", sep = " " ) ) ),
      cex = 1.2, col = "brown" )
text( x = coords_pval[1], y = coords_pval[2],
      labels = expression( atop( ''*italic( pval )*' = 1' ) ),
      cex = 1.2, col = "black" )
title( "B) M1a vs M2a" )

# M7 vs M8
curve( dchisq( x, df = 2 ), from = 0, to =  15 )
abline( v = c( Chisq.crit.M7vsM8, Chisq.crit.M7vsM8_2, diff_M7vsM8 ), col = c( "darkgray", "brown", "red" ) )
coords_dev     <- c( 3.5, 0.48 )
coords_pval    <- c( 3.8, 0.45 )
coords_alphac  <- c( 3.4, 0.41 )
coords_alphac2 <- c( 3.4, 0.38 )
text( x = coords_dev[1], y = coords_dev[2],
      labels = expression( atop( paste( '2', Delta, 'l = 12.54', sep = "" ) ) ),
      cex = 1.2, col = "red" )
text( x = coords_alphac[1], y = coords_alphac[2],
      labels = expression( atop( paste( chi["2,0.05"]^"2", "= 5.99", sep = " " ) ) ),
      cex = 1.2, col = "gray" )
text( x = coords_alphac2[1], y = coords_alphac2[2],
      labels = expression( atop( paste( chi["2,0.01"]^"2", "= 9.21", sep = " " ) ) ),
      cex = 1.2, col = "brown" )
text( x = coords_pval[1], y = coords_pval[2],
      labels = expression( atop( ''*italic( pval )*' = 0.0019' ) ),
      cex = 1.2, col = "black" )
title( "C) M7 vs M8" )

