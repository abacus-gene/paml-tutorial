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
#    Order in columns: Chicken | Duck | M0
lnL_vals <- read.table( file = "lnL_branch_mods.txt", sep= " ", stringsAsFactors = FALSE, 
                        header = FALSE )
row.names( lnL_vals ) <- c( "Chicken-branch", "Duck-branch", "DuckChicken-branch", "Outgroup-branch", "M0-nobranch" )

# 2. We can now compute the LRT statistic
## 2.1. Chicken-branch vs M0-nobranch ##
diff_chickVSM0nob <- 2*( lnL_vals[1,] - lnL_vals[5,] )
# diff = 43.13729
pchisq( diff_chickVSM0nob, df = 1, lower.tail=F )
# p-val = 5.103027e-11 < 0.05
Chisq.crit.chick <- qchisq( p = 0.95, df = 1 )
# alpha critical value = 3.841459
## 2.2. Chicken-branch vs M0-nobranch ##
diff_duckVSM0nob <- 2*( lnL_vals[2,] - lnL_vals[5,] )
# diff = 36.96113
pchisq( diff_duckVSM0nob, df = 1, lower.tail=F )
# p-val = 1.205079e-09 < 0.05
Chisq.crit.duck <- qchisq( p = 0.95, df = 1 )
# alpha critical value = 3.841459
## 2.3. DuckChicken-branch vs M0-nobranch ##
diff_duckchickVSM0nob <- 2*( lnL_vals[3,] - lnL_vals[5,] )
# diff = 51.03423
pchisq( diff_duckchickVSM0nob, df = 1, lower.tail=F )
# p-val = 9.076912-05 < 0.05
Chisq.crit.duckchick <- qchisq( p = 0.95, df = 1 )
# alpha critical value = 3.841459
## 2.4. Outgroup-branch vs M0-nobranch ##
diff_outgroupVSM0nob <- 2*( lnL_vals[4,] - lnL_vals[5,] )
# diff = 17.15947
pchisq( diff_outgroupVSM0nob, df = 2, lower.tail=F )
# p-val = 0.0001878751 < 0.05
Chisq.crit.outgroup <- qchisq( p = 0.95, df = 2 )
# alpha critical value = 5.991465

# 3. Plot results
par( mfrow = c( 2, 2 ) )
# Chicken branch
curve( dchisq( x, df = 1 ), from = 0, to =  45 )
abline( v = c( Chisq.crit.chick, diff_chickVSM0nob ), col = c( "darkgray", "red" ) )
coords_dev    <- c( 33.8, 0.4 )
coords_pval   <- c( 34.4, 0.35 )
coords_alphac <- c( 33, 0.28 )
text( x = coords_dev[1], y = coords_dev[2],
      labels = expression( atop( paste( '2', Delta, 'l = 43.14', sep = "" ) ) ),
      cex = 1.2, col = "red" )
text( x = coords_alphac[1], y = coords_alphac[2],
      labels = expression( atop( paste( chi["1,0.05"]^"2", "= 3.84", sep = " " ) ) ),
      cex = 1.2, col = "gray" )
text( x = coords_pval[1], y = coords_pval[2],
      labels = expression( atop( ''*italic( pval )*' = 5.10e-11' ) ),
      cex = 1.2, col = "black" )
#title( "A) Branch model (chicken) VS M0 model" )
title( expression( 'A) '*italic(Chicken)*': Branch model VS M0 model ' ) )

# Duck branch
curve( dchisq( x, df = 1 ), from = 0, to =  40 )
abline( v = c( Chisq.crit.duck, diff_duckVSM0nob ), col = c( "darkgray", "red" ) )
coords_dev    <- c( 28.6, 0.4 )
coords_pval   <- c( 29.4, 0.35 )
coords_alphac <- c( 28, 0.28 )
text( x = coords_dev[1], y = coords_dev[2],
      labels = expression( atop( paste( '2', Delta, 'l = 36.96', sep = "" ) ) ),
      cex = 1.2, col = "red" )
text( x = coords_alphac[1], y = coords_alphac[2],
      labels = expression( atop( paste( chi["1,0.05"]^"2", "= 3.84", sep = " " ) ) ),
      cex = 1.2, col = "gray" )
text( x = coords_pval[1], y = coords_pval[2],
      labels = expression( atop( ''*italic( pval )*' = 1.21e-09' ) ),
      cex = 1.2, col = "black" )
#title( "B) Branch model (duck) VS M0 model" )
title( expression( 'B) '*italic(Duck)*': Branch model VS M0 model ' ) )

# DuckChicken branch
curve( dchisq( x, df = 1 ), from = 0, to =  60 )
abline( v = c( Chisq.crit.duckchick, diff_duckchickVSM0nob ), col = c( "darkgray", "red" ) )
coords_dev    <- c( 28.6, 0.3 )
coords_pval   <- c( 29.4, 0.25 )
coords_alphac <- c( 28, 0.20 )
text( x = coords_dev[1], y = coords_dev[2],
      labels = expression( atop( paste( '2', Delta, 'l = 51.03', sep = "" ) ) ),
      cex = 1.2, col = "red" )
text( x = coords_alphac[1], y = coords_alphac[2],
      labels = expression( atop( paste( chi["1,0.05"]^"2", "= 3.84", sep = " " ) ) ),
      cex = 1.2, col = "gray" )
text( x = coords_pval[1], y = coords_pval[2],
      labels = expression( atop( ''*italic( pval )*' = 9.08-05' ) ),
      cex = 1.2, col = "black" )
#title( "C) Branch model (duck-chicken) VS M0 model" )
title( expression( 'C) '*italic(Duck-Chicken)*': Branch model VS M0 model ' ) )

# Outgroup branch
curve( dchisq( x, df = 1 ), from = 0, to =  60 )
abline( v = c( Chisq.crit.outgroup, diff_outgroupVSM0nob ), col = c( "darkgray", "red" ) )
coords_dev    <- c( 28.6, 0.3 )
coords_pval   <- c( 29.4, 0.25 )
coords_alphac <- c( 28, 0.20 )
text( x = coords_dev[1], y = coords_dev[2],
      labels = expression( atop( paste( '2', Delta, 'l = 17.15', sep = "" ) ) ),
      cex = 1.2, col = "red" )
text( x = coords_alphac[1], y = coords_alphac[2],
      labels = expression( atop( paste( chi["2,0.05"]^"2", "= 5.99", sep = " " ) ) ),
      cex = 1.2, col = "gray" )
text( x = coords_pval[1], y = coords_pval[2],
      labels = expression( atop( ''*italic( pval )*' = 1.87e-03' ) ),
      cex = 1.2, col = "black" )
#title( "C) Branch model (duck-chicken) VS M0 model" )
title( expression( 'D) '*italic(Outgroup)*': Branch model VS M0 model ' ) )
