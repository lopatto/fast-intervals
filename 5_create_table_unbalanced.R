rm(list=ls())
library(xtable)
###########################################
##########Scenario 1: ununbalanced#############
###########################################
######
#n=50#
######
table_unbalanced_s1=c()
setwd('/home/hc654/palmer_scratch/binary_outcome_sim/unbalanced_s1/n50')
files_s1_50=list.files()
output_50_s1=c()
for (i in 1:length(files_s1_50)){
  output_50_s1=rbind(output_50_s1, read.csv(files_s1_50[i])[,2:15])
}
result_50_s1=c(mean(output_50_s1[,1]),median(output_50_s1[,2]),mean(output_50_s1[,3]),median(output_50_s1[,4]),mean(output_50_s1[,5]),median(output_50_s1[,6]))
mean(output_50_s1[,13])
mean(output_50_s1[,14])

########
#n=100##
########
setwd('/home/hc654/palmer_scratch/binary_outcome_sim/unbalanced_s1/n100')
files_s1_100=list.files()
output_100_s1=c()
for (i in 1:length(files_s1_100)){
  output_100_s1=rbind(output_100_s1, read.csv(files_s1_100[i])[,2:15])
}
result_100_s1=c(mean(output_100_s1[,1]),median(output_100_s1[,2]),mean(output_100_s1[,3]),median(output_100_s1[,4]),mean(output_100_s1[,5]),median(output_100_s1[,6]))
mean(output_100_s1[,13])
mean(output_100_s1[,14])

########
#n=200##
########
setwd('/home/hc654/palmer_scratch/binary_outcome_sim/unbalanced_s1/n200')
files_s1_200=list.files()
output_200_s1=c()
for (i in 1:length(files_s1_200)){
  output_200_s1=rbind(output_200_s1, read.csv(files_s1_200[i])[,2:15])
}
result_200_s1=c(mean(output_200_s1[,1]),median(output_200_s1[,2]),mean(output_200_s1[,3]),median(output_200_s1[,4]),mean(output_200_s1[,5]),median(output_200_s1[,6]))
mean(output_200_s1[,13])
mean(output_200_s1[,14])

table_unbalanced_s1=rbind(result_50_s1,result_100_s1,result_200_s1)


xtable(round(table_unbalanced_s1,3))
table_unbalanced_s1






###########################################
##########Scenario 2: unbalanced#############
###########################################
######
#n=50#
######
setwd('/home/hc654/palmer_scratch/binary_outcome_sim/unbalanced_s2/n50')
files_s2_50=list.files()
output_50_s2=c()
for (i in 1:length(files_s2_50)){
  output_50_s2=rbind(output_50_s2, read.csv(files_s2_50[i])[,2:15])
}
result_50_s2=c(mean(output_50_s2[,1]),median(output_50_s2[,2]),mean(output_50_s2[,3]),median(output_50_s2[,4]),mean(output_50_s2[,5]),median(output_50_s2[,6]))
mean(output_50_s2[,13])
mean(output_50_s2[,14])
########
#n=100##
########
setwd('/home/hc654/palmer_scratch/binary_outcome_sim/unbalanced_s2/n100')
files_s2_100=list.files()
output_100_s2=c()
for (i in 1:length(files_s2_100)){
  output_100_s2=rbind(output_100_s2, read.csv(files_s2_100[i])[,2:15])
}
result_100_s2=c(mean(output_100_s2[,1]),median(output_100_s2[,2]),mean(output_100_s2[,3]),median(output_100_s2[,4]),mean(output_100_s2[,5]),median(output_100_s2[,6]))
mean(output_100_s2[,13])
mean(output_100_s2[,14])

########
#n=200##
########
setwd('/home/hc654/palmer_scratch/binary_outcome_sim/unbalanced_s2/n200')
files_s2_200=list.files()
output_200_s2=c()
for (i in 1:length(files_s2_200)){
  output_200_s2=rbind(output_200_s2, read.csv(files_s2_200[i])[,2:15])
}
result_200_s2=c(mean(output_200_s2[,1]),median(output_200_s2[,2]),mean(output_200_s2[,3]),median(output_200_s2[,4]),mean(output_200_s2[,5]),median(output_200_s2[,6]))
mean(output_200_s2[,13])
mean(output_200_s2[,14])
table_unbalanced_s2=rbind(result_50_s2,result_100_s2,result_200_s2)

xtable(round(table_unbalanced_s2,3))






