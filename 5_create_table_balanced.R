rm(list=ls())
library(xtable)
###########################################
##########Scenario 1: Balanced#############
###########################################
######
#n=50#
######
table_balanced_s1=c()
setwd('/home/hc654/palmer_scratch/binary_outcome_sim/balanced_s1/n50')
files_s1_50=list.files()
output_50_s1=c()
for (i in 1:length(files_s1_50)){
  output_50_s1=rbind(output_50_s1, read.csv(files_s1_50[i])[,2:9])
}
result_50_s1=c(mean(output_50_s1[,1]),median(output_50_s1[,2]),mean(output_50_s1[,3]),median(output_50_s1[,4]))
result_50_s1
########
#n=100##
########
setwd('/home/hc654/palmer_scratch/binary_outcome_sim/balanced_s1/n100')
files_s1_100=list.files()
output_100_s1=c()
for (i in 1:length(files_s1_100)){
  output_100_s1=rbind(output_100_s1, read.csv(files_s1_100[i])[,2:9])
}
result_100_s1=c(mean(output_100_s1[,1]),median(output_100_s1[,2]),mean(output_100_s1[,3]),median(output_100_s1[,4]))
result_100_s1


########
#n=200##
########
setwd('/home/hc654/palmer_scratch/binary_outcome_sim/balanced_s1/n200')
files_s1_200=list.files()
output_200_s1=c()
for (i in 1:length(files_s1_200)){
  output_200_s1=rbind(output_200_s1, read.csv(files_s1_200[i])[,2:9])
}
result_200_s1=c(mean(output_200_s1[,1]),median(output_200_s1[,2]),mean(output_200_s1[,3]),median(output_200_s1[,4]))
result_200_s1



#########
#n=1000##
#########
setwd('/home/hc654/palmer_scratch/binary_outcome_sim/balanced_s1/n1000')
files_s1_1000=list.files()
output_1000_s1=c()
for (i in 1:length(files_s1_1000)){
  output_1000_s1=rbind(output_1000_s1, read.csv(files_s1_1000[i])[,2:9])
}
result_1000_s1=c(mean(output_1000_s1[,1]),median(output_1000_s1[,2]),mean(output_1000_s1[,3]),median(output_1000_s1[,4]))
result_1000_s1


table_s1=rbind(result_50_s1,result_100_s1,result_200_s1,result_1000_s1)

xtable( round(table_s1,2))

###########################################
##########Scenario 2: Balanced#############
###########################################
######
#n=50#
######
setwd('/home/hc654/palmer_scratch/binary_outcome_sim/balanced_s2/n50')
files_s2_50=list.files()
output_50_s2=c()
for (i in 1:length(files_s2_50)){
  output_50_s2=rbind(output_50_s2, read.csv(files_s2_50[i])[,2:9])
}
result_50_s2=c(mean(output_50_s2[,1]),median(output_50_s2[,2]),mean(output_50_s2[,3]),median(output_50_s2[,4]))
result_50_s2

########
#n=100##
########
setwd('/home/hc654/palmer_scratch/binary_outcome_sim/balanced_s2/n100')
files_s2_100=list.files()
output_100_s2=c()
for (i in 1:length(files_s2_100)){
  output_100_s2=rbind(output_100_s2, read.csv(files_s2_100[i])[,2:9])
}
result_100_s2=c(mean(output_100_s2[,1]),median(output_100_s2[,2]),mean(output_100_s2[,3]),median(output_100_s2[,4]))
result_100_s2


########
#n=200##
########
setwd('/home/hc654/palmer_scratch/binary_outcome_sim/balanced_s2/n200')
files_s2_200=list.files()
output_200_s2=c()
for (i in 1:length(files_s2_200)){
  output_200_s2=rbind(output_200_s2, read.csv(files_s2_200[i])[,2:9])
}

result_200_s2=c(mean(output_200_s2[,1]),median(output_200_s2[,2]),mean(output_200_s2[,3]),median(output_200_s2[,4]))
result_200_s2


table_200_s2=rbind(result_50_s2,result_100_s2,result_200_s2)

xtable( round(table_200_s2,2))

########
#n=1000#
########

setwd('/home/hc654/palmer_scratch/binary_outcome_sim/balanced_s2/n1000')
files_s2_1000=list.files()
output_1000_s2=c()
for (i in 1:length(files_s12_1000)){
  output_1000_s2=rbind(output_1000_s2, read.csv(files_s2_1000[i])[,2:9])
}
result_1000_s2=c(mean(output_1000_s2[,1]),median(output_1000_s2[,2]),mean(output_1000_s2[,3]),median(output_1000_s2[,4]))
result_1000_s2


table_s2=rbind(result_50_s2,result_100_s2,result_200_s2,result_1000_s2)

xtable( round(table_s2,2))
xtable( round(table_s1,2))




