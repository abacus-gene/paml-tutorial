#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------------------------------------#
# LOAD PACKAGES, FUNCTIONS, AND SET ENVIRONMENT #
#-----------------------------------------------#
# This package lets you find automatically the path to a specific location
# in your file structure
# If you have not installed this package, you will need to install it. 
# You can uncomment the following line to do this:
#install.packages( "rstudioapi" )
library( rstudioapi ) 
# Get the path of current open file
path_to_file <- getActiveDocumentContext()$path 
# Get working directory path
wd      <- paste( dirname( path_to_file ), "/", sep = "" )
wd.name <- dirname( path_to_file )
# Set working directory
setwd( wd )

#-------------------------------------------------------------#
# DEFINE GLOBAL VARIABLES -- modify according to your dataset #
#-------------------------------------------------------------#
# Name that you want the output calibrated tree file to have.
# E.g. the file name will have the following format
# "<out_name>_calib_MCMCtree.tree".
out_name <- c( "mtcdnapri" )
num_dat  <- length( out_name )

# Path to your input text file that allows you to match the flags you have 
# used to constrain node ages with the calibration you 
# want to use in `MCMCtree` notation. The format you need to follow is given
# below:
#
#   - Header.
#   - One row per calibration.
#   - No spaces at all, semi-colon separated.
#   - There are 4 columns:
#       - Name you want to give to the calibrated node (no spaces!).
#       - Name of one of the tips (e.g., tip 1) that leads to MRCA (no spaces!).
#       - Name of the other tip (e.g., tip 2) that leads to MRCA (no spaces!).
#       - Calibration in `MCMCtree` notation (no spaces!). More details on the 
#         `MCMCtree` notation you need to use in the fourth column in the PAML
#         documentation:
#         https://github.com/abacus-gene/paml/blob/master/doc/pamlDOC.pdf
# 
# E.g.: Header and one row in a calibration text file:
#
# ```
# name;tip1;tip2;MCMCtree
# root;sp1;sp2;'B(0.256,1.34,0.025,1e-300)'
# ```
#
# NOTE: Always check that there is at least one blank line at the 
# end of the this text file! Otherwise, you will get an error telling you that 
# there is an incomplete final line in these files. This file needs to be
# already in PHYLIP format. Please follow the same format as used in the 
# example tree file provided.
path_textconv <- c( "../raw_calibs/calibrations.txt" )
calibrations_all <- read.table( file = path_textconv,
                                stringsAsFactors = FALSE, sep = ";",
                                blank.lines.skip = TRUE, header = TRUE,
                                colClasses = rep( "character", 4 ) )

calibrations <- keep_indexes <- ind_dup <- nodes_dup <- tt_all <-
  vector( mode = "list", num_dat )
calibrations[[ 1 ]] <- calibrations_all
names( calibrations ) <- names( keep_indexes )  <- names( ind_dup ) <- 
  names( nodes_dup ) <- names( tt_all ) <- c( "mtcdnapri" )

# Path to tree
path_tree <- c( "../../mtCDNApri.trees" )
for( c in 1:length( calibrations ) ){
  cat( "\n[[ ANALYSING CALIBRATION FILE ", names(calibrations)[c], " ]]\n" )
  tt_ape <- ape::read.tree( file = path_tree[c] )
  keep_indexes[[c]] <- matrix( 0, nrow = length(rownames(calibrations[[c]])),
                               ncol = 3 )
  # Generate empty vector with as many entries as nodes in the tree
  tt_ape$node.label <- rep( NA, tt_ape$Nnode )
  for( i in 1:length(rownames(calibrations[[c]])) ){
    # Get calibration in the same format input by the user
    node_lab <- calibrations[[c]][i,4]
    # Get MRCA for these two tips
    mrca <- ape::getMRCA( phy = tt_ape, tip = c( calibrations[[c]][i,2],
                                                 calibrations[[c]][i,3]) )
    keep_indexes[[c]][i,1] <- mrca-ape::Ntip(tt_ape)
    keep_indexes[[c]][i,2] <- calibrations[[c]][i,1]
    keep_indexes[[c]][i,3] <- paste( calibrations[[c]][i,2], "-",
                                calibrations[[c]][i,3], "-",
                                node_lab, sep = "" )
    print( mrca-ape::Ntip( tt_ape ) )
    # Replace node label accordingly
    tt_ape$node.label[mrca-ape::Ntip(tt_ape)] <- paste0( "[",
                                                         calibrations[[c]][i,1],
                                                         "]", collapse = "" )
  }
  # Find duplicates
  ind_dup[[c]]   <- which( duplicated(keep_indexes[[c]][,1]) == TRUE )
  nodes_dup[[c]] <- which( keep_indexes[[c]][,1] %in% as.numeric( keep_indexes[[c]][ind_dup[[c]],1] ) )
  keep_indexes[[c]][nodes_dup[[c]],]
  # Save tree with node labels
  tt_all[[c]] <- tt_ape
} 

##>> ---
## CHECK for duplicates
ind_dup
## NO duplicates, success!
##>> ---

# Remove "NA" from the labs
for( c in 1:length( calibrations ) ){
  ind_na_bools <- is.na( x = tt_all[[c]]$node.label )
  ind_na       <- which( ind_na_bools == TRUE )
  tt_all[[c]]$node.label[ind_na] <- ""
  # Write PHYLIP header, then the calibrated tree
  writeLines( text = paste( length(tt_all[[c]]$tip.label ), " 1", sep = "" ), 
              con = paste( "../raw_calibs/cals_only_",
              names( calibrations )[c], ".tree", sep = "" ) )
  ape::write.tree( phy = tt_all[[c]],
                   file = paste( "../raw_calibs/cals_only_",
                   names( calibrations )[c], ".tree", sep = "" ),
                   append = TRUE )
}

#>> TEST
tt_all[[1]]$node.label[c(1,2,4)]
# plot.phylo(tt_all[[c]])
# nodelabels(text=1:tt_all[[c]]$Nnode,node=1:tt_all[[c]]$Nnode+Ntip(tt_all[[c]]),
#            cex = 0.7, frame = "none" )
#>> SUCCESS!
# PLAN: Load the tree later to then replace with the corresponding
# calibrations following my old script

#---------------------------------#
# READ TREE AND CALIBRATIONS FILE #
#---------------------------------#
# Read tree and get phylip header
# NOTE: Make always sure that there is at least one blank line at the 
# end of the tree file! Otherwise, you will get an error telling you that 
# there is an incomplete final line in these files.
tt_name <- c( "../raw_calibs/cals_only_mtcdnapri.tree" )
for( t in 1:length( tt_name ) ){
  
  #-------------------#
  # Get PHYLIP header #
  #-------------------#
  tt            <- readLines( tt_name[t] )
  phylip.header <- tt[1]
  tt            <- tt2 <- tt3 <- tt[2]
  
  #--------------------------------#
  # REPLACE TAGS WITH CALIBRATIONS #
  #--------------------------------#
  # Replace calibration names with corresponding calibration
  for( j in 1:length( rownames( calibrations[[t]] ) ) ){
    # Get node label as input by user
    node_lab <- calibrations[[t]][j,4]
    # Get rid of unnecessary notation for susbequent formatting
    tmp_calib <- gsub( x = node_lab, pattern = "\\(..*",
                       replacement = "" )
    tmp_calib <- gsub( x = tmp_calib, pattern = "[0-9]..*",
                       replacement = "" )
    # Conditional is used so that the single quotation marks are only kept 
    # in the upper-bound calibration for the root. Inequality calibrations
    # do not require single quotation marks
    if( tmp_calib == 'B' || tmp_calib == 'U' || tmp_calib == 'L' ){
      tt <- gsub( pattern = paste0( "\\[", calibrations[[t]][j,1], "\\]" ),
                  x = tt,
                  replacement = paste( "'", node_lab, "'", sep = "" ) )
    }else{ # For cross-braced nodes
      tt <- gsub( pattern = paste0( "\\[", calibrations[[t]][j,1], "\\]" ),
                  x = tt,
                  replacement = paste( node_lab, sep = "" ) )
    }
    # Copy to visualise in `FigTree`
    reps <- gsub( x = gsub( x = gsub( x = gsub( x = gsub( x = node_lab,
                                                          pattern = "\\{",
                                                          replacement = "(" ),
                                                pattern = "\\}",
                                                replacement = ")" ), 
                                      pattern = "\\[|\\]", replacement = "" ),
                            pattern = "\\#", replacement = "flag" ),
                  pattern = " ", replacement = "-" )
    # For cross-braced calibrations without fossil
    if( tmp_calib == '#' ){
      reps <- gsub( x = gsub( x = reps, pattern = "\\#", replacement = "flag" ),
                    pattern = "\\]", replacement = "" )
      tt2 <- gsub( pattern = paste0( "\\[", calibrations[[t]][j,1], "\\]" ),
                   x = tt2,
                   replacement = paste0( "'", reps, "-", calibrations[[t]][j,1],
                                         "'", collapse = "" ) )
    }else{ # For the rest of calibrations
      tt2 <- gsub( pattern = paste0( "\\[", calibrations[[t]][j,1], "\\]" ),
                   x = tt2,
                   replacement = paste0( "'", reps, "-", calibrations[[t]][j,1],
                                         "'", collapse = "" ) )
    }
    # Generate an uncalibrated tree for `BASEML`/`CODEML`!
    tt3 <- gsub( pattern = paste0( "\\[", calibrations[[t]][j,1], "\\]" ),
                 x = tt3,
                 replacement = "" )
  }
  
  #-------------------------------#
  # WRITE CALIBRATED TREE IN FILE #
  #-------------------------------#
  out_dir <- "../../tree_display/"
  if( ! dir.exists( "../../tree_display/" ) ){
    dir.create( "../../tree_display/" )
  }
  # Check if `out_name` != `names(calibrations)[t]`
  if( out_name == names( calibrations )[t] ){
    out_fullname <- paste( out_name )
  }else{
    out_fullname <- paste( out_name,  names( calibrations )[t], sep = "_" )
  }
  # Write calibrated tree file and file to visualise in `FigTree`
  write( x = phylip.header, file = paste( out_dir, out_fullname,
                                          "_calib_MCMCtree.tree", sep = "" ) )
  write( x = tt, file = paste( out_dir, out_fullname,
                               "_calib_MCMCtree.tree", sep = "" ),
         append = TRUE )
  write( x = phylip.header, file = paste( out_dir, out_fullname,
                                          "_fordisplay_calib_MCMCtree.tree",
                                          sep = "" ) )
  write( x = tt2, file = paste( out_dir, out_fullname,
                                "_fordisplay_calib_MCMCtree.tree", sep = "" ),
         append = TRUE )
  # Write an uncalibrated tree file, only once (first iteration)!
  if( t == 1 ){
    write( x = phylip.header, file = paste( out_dir, out_fullname,
                                            "_uncalib.tree", sep = "" ) )
    write( x = tt3, file = paste( out_dir, out_fullname,
                                  "_uncalib.tree", sep = "" ),
           append = TRUE )
  }
  
}




