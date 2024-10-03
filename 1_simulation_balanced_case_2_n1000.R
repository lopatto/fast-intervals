home_dir = "/home/hc654/BinaryOutcomePermutationTest"
setwd(home_dir)

#the script below has the compare_ATE function to run simulations
source("2_sim_raw.R")  #compare_ATE



#command line, take input from command line
args = commandArgs(trailingOnly=TRUE)
index = as.integer(args[1])

#set random seed
set.seed(index)

#output file name
save_path_1000=  paste0("/home/hc654/palmer_scratch/binary_outcome_sim/balanced_s2/n1000/",index,'_1000.csv')

#This script contains simulations for the balanced case with SATE=0 and v11_share=0.5, v01_share=0, v10_share=0.
#The sample sizes is n=1000
#The number of treated units is nt=n/2.
#This test case is designed such that the Wald test is expected to perform poorly even with relatively large sample sizes.

############################################################
#########Balanced Case and SATE=0###########################
############################################################
v11_share=0.08
v10_share=0
v01_share=0
test_type=1
nsim=10

#############################
#####n=1000###################
#############################
n=1000
nt=500
print(n)
result_1000=compare_ATE(nsim,n,nt,v11_share,v10_share,v01_share,1,display_progress=FALSE,warning_msg=FALSE)
write.csv(result_1000,file = save_path_1000)


