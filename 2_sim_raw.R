#This script contains function used to compare the coverage of the SATE CIs with the permutation test and the Wald test

library(pbapply) #pbsapply
library(sandwich) #vcovHC
source("2_permutation_test.R")


compare_ATE_inner = function(i,n,nt,v11_share,v10_share,v01_share,test_type,display_progress=TRUE,warning_msg=FALSE){
   
    print(i) 
    #This function generates a test case for the permutation test
    #We are testing two different method for calculating the SATE CI

    nc=n-nt #number of control units
    index=1:n #index of experimental units
    v11 = round(v11_share*n,0) #alway positve
    v10 = round(v10_share*n,0) #positive only if treated
    v01 = round(v01_share*n,0) #positive only if control
    v00 = n-v11-v10-v01 #always negative

    po_vector = c(rep('11',v11),rep('10',v10),rep('01',v01),rep('00',v00))
    y1_vector = c(rep(1,v11),rep(1,v10),rep(0,v01),rep(0,v00))
    y0_vector = c(rep(1,v11),rep(0,v10),rep(1,v01),rep(0,v00))

    SATE = (v10-v01)/n

    treatment_status=rep(0,n)
    treatment_status[sample(index,nt)]=1
    outcome = treatment_status*y1_vector+(1-treatment_status)*y0_vector

    #observed data
     data = as.data.frame(cbind(index,treatment_status,outcome))
     colnames(data)=c('index','treatment','outcome')

    #permutation test
     start_time=proc.time()
     result_test_case = permutation_test(data, test_type , alpha=0.05, epsilon=0.005, nperm=0,return_interval = TRUE ,display_progress = display_progress,warning_msg=warning_msg  )
     end_time=proc.time()
     dim_time=(end_time-start_time)[3]
     CS= result_test_case[['CS']] 

     cov_perm = (SATE<= CS[2]) && (SATE>= CS[1]) #coverage of the permutation test
     perm_width =CS[2]-CS[1] #width of the permutation test

     #Wald test
     temp=lm(outcome~treatment,data=data)
     ate_est = temp$coefficients[2] #estimated treatment effect
     std=sqrt(vcovHC(temp, type = "HC1")[2,2]) #standard error of the treatment effect  

     cov_wald = (SATE <= ate_est + qnorm(0.975)*std) & (SATE >= ate_est - qnorm(0.975)*std) #coverage of the Wald test
     wald_width = qnorm(0.975)*std*2 #width of the permutation test

     output = c(cov_perm=cov_perm,perm_width=perm_width, cov_wald=cov_wald, wald_width=wald_width, CS_wald_lb =ate_est - qnorm(0.975)*std, CS_wald_ub =ate_est + qnorm(0.975)*std,CS_lb=CS[1],CS_ub=CS[2],dim_time=dim_time) #output as a vector
     return(output)

}

compare_ATE = function(nsim,n,nt,v11_share,v10_share,v01_share,test_type,display_progress=TRUE,warning_msg=FALSE){

    result_mat = pbsapply(1:nsim, function(i) compare_ATE_inner(i,n,nt,v11_share,v10_share,v01_share,test_type,display_progress=display_progress,warning_msg=warning_msg))

    result_mat = t(result_mat)

    colnames(result_mat) = c('cov_perm','perm_width', 'cov_wald','wald_width','CS_wald_lb','CS_wald_ub', 'CS_lb','CS_ub','dim_time')
    
    return(result_mat)
}


compare_ATE_inner_unbalanced = function(i,n,nt,v11_share,v10_share,v01_share,test_type,display_progress=TRUE,warning_msg=FALSE){
  
  print(i) 
  #This function generates a test case for the permutation test
  #We are testing two different method for calculating the SATE CI
  
  nc=n-nt #number of control units
  index=1:n #index of experimental units
  v11 = round(v11_share*n,0) #alway positve
  v10 = round(v10_share*n,0) #positive only if treated
  v01 = round(v01_share*n,0) #positive only if control
  v00 = n-v11-v10-v01 #always negative
  
  po_vector = c(rep('11',v11),rep('10',v10),rep('01',v01),rep('00',v00))
  y1_vector = c(rep(1,v11),rep(1,v10),rep(0,v01),rep(0,v00))
  y0_vector = c(rep(1,v11),rep(0,v10),rep(1,v01),rep(0,v00))
  
  SATE = (v10-v01)/n
  
  treatment_status=rep(0,n)
  treatment_status[sample(index,nt)]=1
  outcome = treatment_status*y1_vector+(1-treatment_status)*y0_vector
  
  #observed data
  data = as.data.frame(cbind(index,treatment_status,outcome))
  colnames(data)=c('index','treatment','outcome')
  
  #permutation test
  start_time=proc.time() #time the code
  result_test_case = permutation_test(data, test_type , alpha=0.05, epsilon=0, nperm=20000,return_interval = TRUE ,display_progress = display_progress,warning_msg=warning_msg  )
  end_time=proc.time()
  dim_time=(end_time-start_time)[3]
  
  CS=rep(NA,2)
  CS[1]=min(result_test_case$CS[which(result_test_case$CS[,2]==1),1])
  CS[2]=max(result_test_case$CS[which(result_test_case$CS[,2]==1),1])
  
  cov_perm = (SATE<= CS[2]) && (SATE>= CS[1]) #coverage of the permutation test
  perm_width =CS[2]-CS[1] #width of the permutation test
  
  
  #studentized statistics
  start_time_std=proc.time() #time the code
  result_test_case_std = permutation_test(data, test_type=2 , alpha=0.05, epsilon=0, nperm=20000,return_interval = TRUE ,display_progress = display_progress,warning_msg=warning_msg  )
  end_time_std=proc.time()
  std_time=(end_time_std-start_time_std)[3]
  
  CS_std=rep(NA,2)
  CS_std[1]=min(result_test_case_std$CS[which(result_test_case_std$CS[,2]==1),1])
  CS_std[2]=max(result_test_case_std$CS[which(result_test_case_std$CS[,2]==1),1])
  
  cov_perm_std = (SATE<= CS_std[2]) && (SATE>= CS_std[1]) #coverage of the permutation test
  perm_width_std =CS_std[2]-CS_std[1] #width of the permutation test  
  
  #Wald test
  temp=lm(outcome~treatment,data=data)
  ate_est = temp$coefficients[2] #estimated treatment effect
  std=sqrt(vcovHC(temp, type = "HC1")[2,2]) #standard error of the treatment effect  
    
  cov_wald = (SATE <= ate_est + qnorm(0.975)*std) & (SATE >= ate_est - qnorm(0.975)*std) #coverage of the Wald test
  wald_width = qnorm(0.975)*std*2 #width of the permutation test

  output = c(cov_perm=cov_perm,perm_width=perm_width, cov_perm_std=cov_perm_std,perm_width_std=perm_width_std,cov_wald=cov_wald, wald_width=wald_width, CS_wald_lb =ate_est - qnorm(0.975)*std, CS_wald_ub =ate_est + qnorm(0.975)*std,CS_lb=CS[1],CS_ub=CS[2],CS_std_lb=CS_std[1],CS_std_ub=CS_std[2],dim_time=dim_time,std_time=std_time) #output as a vector
  return(output)
  
}

compare_ATE_unbalanced = function(nsim,n,nt,v11_share,v10_share,v01_share,test_type,display_progress=TRUE,warning_msg=FALSE){
  
  
  result_mat = pbsapply(1:nsim, function(i) compare_ATE_inner_unbalanced(i,n,nt,v11_share,v10_share,v01_share,test_type,display_progress=display_progress,warning_msg=warning_msg))
  
  result_mat = t(result_mat)
  
  colnames(result_mat) = c('cov_perm','perm_width','cov_perm_std', 'perm_width_std',   'cov_wald','wald_width','CS_wald_lb','CS_wald_ub', 'CS_lb','CS_ub','CS_std_lb','CS_std_ub','dim_time','std_time')
  return(result_mat)
}
