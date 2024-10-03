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
save_path_50= paste0("/home/hc654/palmer_scratch/binary_outcome_sim/unbalanced_s1/n50/",index,'_50.csv')
save_path_100= paste0("/home/hc654/palmer_scratch/binary_outcome_sim/unbalanced_s1/n100/",index,'_100.csv')
save_path_200=  paste0("/home/hc654/palmer_scratch/binary_outcome_sim/unbalanced_s1/n200/",index,'_200.csv')

#This script contains simulations for the balanced case with SATE=0 and v11_share=0.5, v01_share=0, v10_share=0.
#The sample sizes is n=50,100,200
#The number of treated units is nt=0.6*n.
#This test case is designed such that the Wald test can perform poorly even with relatively large sample sizes.

############################################################
#########Balanced Case and SATE=0###########################
############################################################
print('Unbalanced')
v11_share=0.5
v10_share=0
v01_share=0
test_type=1
nsim=10

#############################
#####n=100###################
#############################
n=50
nt=30
print(n)
result_50=compare_ATE_unbalanced(nsim,n,nt,v11_share,v10_share,v01_share,1,display_progress=FALSE,warning_msg=FALSE)
write.csv(result_50,file = save_path_50)


#############################
#####n=100###################
#############################
n=100
nt=60
print(n)
result_100=compare_ATE_unbalanced(nsim,n,nt,v11_share,v10_share,v01_share,1,display_progress=FALSE,warning_msg=FALSE)
write.csv(result_100,file = save_path_100)

#############################
#####n=200###################
#############################
n=200
nt=120
print(n)
start.time=proc.time()
result_200=compare_ATE_unbalanced(nsim,n,nt,v11_share,v10_share,v01_share,1,display_progress=FALSE,warning_msg=FALSE)
end.time=proc.time()
print(end.time-start.time)
write.csv(result_200,file = save_path_200)





