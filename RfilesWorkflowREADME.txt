# R files work flow:

# 1. Start with the R file HMP_metadata_updated_taxonomy_StR_CSTs_020619.R
#  and run it without running all the lines saving RData files or text files.
#  The final product from that file is the list "integrated.data2" that is read in the next R file:

# 2. continue with the file "AllKalmanPreds.R" that feeds on the integrated.data2 list.  In this file, the ML estimates of the noise-filtered (log) abundances are calculated.  These will be used next in the file "LactobacillusVsRest.R"

# 3. Open the LactobacillusVsRest.R file and run it to do the stability calculations and the PCA in the space of the stability metrics.  Skip the lines writing files 

