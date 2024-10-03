library(Rcpp) #load Rcpp package
library(RcppProgress) #load RcppProgress package for tracking progress

print('Compling the Rcpp code')
sourceCpp('permutation_test.cpp')

permutation_test = function(data, test_type, alpha=0.05, epsilon=0.01,nperm=0,return_interval=TRUE,display_progress=TRUE,warning_msg=FALSE){

#This function calculates the confidence interval for the general estimand using the permutation test
# Inputs:
# data: a data frame contains columns 'treatment', and 'outcome'. 
#       'treatment' is a binary variable indicating the treatment status of the unit.
#       'outcome' is a binary variable indicating the outcome of the unit.
# test_type: an integer indicating the type of test to be conducted.
#            1: difference-in-means test for the SATE
# alpha: the significance level of the test. 
# epsilon: the slack probability to calculate the number of simulated draws.
#        The implied confidence level has  1-alpha-epsilon.
# nperm: the number of permutation draws. If nperm=0, the function will calculate the number of permutation draws based on epsilon.



     
     n =nrow(data) #number of experimental units
     nt = sum(data[,'treatment']) #number of treated


    if (nperm==0){

        if (warning_msg){
            print(paste0('Using the default number of permutation draws with epsilon ', epsilon))
            print(nperm)

        }
        nperm = round(1/(epsilon^2)* log(4*n*log(n,2)/epsilon))+1 #number of simulated draws
        print(nperm)
        
         }   
    n11= sum(data[data$treatment==1,'outcome'])  #number of treated units with outcome 1
    n10 = sum(1-data[data$treatment==1,'outcome']) #number of treated units with outcome 0
    n01 = sum(data[data$treatment==0,'outcome']) #number of control units with outcome 1
    n00 = sum(1-data[data$treatment==0,'outcome']) #number of control units with outcome 0

    Cn_max =  (sum(data[,'outcome'] * (2*data[,'treatment']-1))-nt+n) #generating largest possible tau0 (up to normalization by sample size)
    Cn_min =  (sum(data[,'outcome'] * (2*data[,'treatment']-1))-nt) #generating smallest possible tau0 (up to normalization sample size)
    Cn= (Cn_min:Cn_max)/n #generating all possible tau0 


    if (test_type==1){

        #if test_type is 1, we are using the difference-in-means statistic to construct a CI for tau0

        test_statistic = mean(data[data$treatment==1,'outcome'])-mean(data[data$treatment==0,'outcome'])

        #pruning the grid by tail inequality
        Cn_preprocessed = c()

        var_max = n^2/(n-1) * 1/(nt*(n-nt)) #largest possible variance
        for (i in 1:length(Cn)){

            diff_temp= abs(test_statistic-Cn[i])

            if (diff_temp<sqrt(var_max/alpha)){
                Cn_preprocessed=c(Cn_preprocessed,Cn[i])
            }
        }
    
        result_table=perm_test_interface(test_statistic, test_type , alpha, epsilon,n11,n10,n01,n00,n,nt,nperm,Cn_preprocessed,length(Cn_preprocessed),display_progress=display_progress,return_interval=return_interval,warning_msg=warning_msg)

        return(list(CS=result_table,grid=Cn_preprocessed))
        
    }else if (test_type==2){

        #if test_type is 2, we are using the studentized difference-in-means statistic to construct a CI for tau0

        numerator = mean(data[data$treatment==1,'outcome'])-mean(data[data$treatment==0,'outcome']) #difference in means

        denominator = sqrt(var(data[data$treatment==1,'outcome'])/nt + var(data[data$treatment==0,'outcome'])/(n-nt)) #estimated std

        if (warning_msg && abs(denominator)<1e-10){
            print('Warning: The estimated standard deviation is too small. The studentized difference-in-means statistic may not be well-defined.')
        }
        #pruning the grid by tail inequality
        Cn_preprocessed = c()

        var_max = n^2/(n-1) * 1/(nt*(n-nt)) #largest possible variance
        for (i in 1:length(Cn)){

            diff_temp= abs(numerator-Cn[i])

            if (diff_temp<sqrt(var_max/alpha)){
                Cn_preprocessed=c(Cn_preprocessed,Cn[i])
            }
        }
    
        result_table=perm_test_STD_interface(numerator,denominator , test_type , alpha, epsilon,n11,n10,n01,n00,n,nt,nperm,Cn_preprocessed,length(Cn_preprocessed),display_progress=display_progress,return_interval=return_interval,warning_msg=warning_msg)

        CS = result_table
        #CS = unique(result_table[which(result_table[,2]==1),1])

    
        return(list(CS=CS,grid=Cn_preprocessed))
        #return(list(CS_SATE=CS,grid=Cn_preprocessed))
        #return(list(CS_general=CS2,grid=Cn_preprocessed))

    }else if(test_type==3){


        test_statistic = (n11/nt) / (n01/(n-nt))
        
        result_table=perm_test_interface(test_statistic , test_type , alpha, epsilon,n11,n10,n01,n00,n,nt,nperm,Cn,length(Cn),display_progress=display_progress,return_interval=return_interval,warning_msg=warning_msg)

        CS = result_table

        return(list(CS=CS,grid=Cn))

    }else if (test_type==4){


        test_statistic = n11 * n00 / (n10*n01)

        result_table=perm_test_interface(test_statistic , test_type , alpha, epsilon,n11,n10,n01,n00,n,nt,nperm,Cn,length(Cn),display_progress=display_progress,return_interval=return_interval,warning_msg=warning_msg)

        CS = result_table  

        return(list(CS=CS,grid=Cn))
     

    }






} 