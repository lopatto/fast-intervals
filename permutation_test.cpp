#include <iostream>
#include <stdexcept>
#include <algorithm> //max_element, min_element and sample
#include <chrono>
#include <thread>


// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include <progress_bar.hpp>

////[[Rcpp::depends(RcppClock)]]
//#include <RcppClock.h>
#include <Rcpp.h>
#include <R.h>

using namespace Rcpp;
using namespace R;
using namespace std;
using namespace this_thread; // sleep_for, sleep_until



double compute_parameter_value(double v11, double v10, double v01, double v00,int test_type, bool warning_msg=false){

// This function computes the parameter value given the po table and the test type
// When the test type is 1 and 2, the parameter value is SATE;
// When the test type is 3, the parameter value is the risk ratio;
// When the test type is 4, the parameter value is the odds ratio.

   double parameter_value=0;

    if (test_type==1 || test_type==2){

        parameter_value = (v10-v01)/(v11+v10+v01+v00);

    }else if (test_type==3){
        
        if ((v11+v01)<1e-10){

            if (warning_msg){
                cout << "zero in the denominator: imputing infinity" << endl;
            }
            parameter_value =  std::numeric_limits<double>::infinity();

        }else{
            parameter_value = (v11+v10)/(v11+v01);
        }

    }else if(test_type==4){

        if (abs(v01+v00)<1e-10 || abs(v11+v01)<1e-10){

            if (warning_msg){
                cout << "zero in the denominator: imputing infinity" << endl;
            }
            parameter_value =  std::numeric_limits<double>::infinity();

        }else{
            parameter_value = (v11+v10)/(v01+v00) * (v10+v00)/(v11+v01);
        }
    }

    return parameter_value;
}


double compute_parameter_value_perm(double v11, double v10, double v01, double v00,int test_type, bool warning_msg=false){

// This function computes the parameter value given the po table and the test type
// When the test type is 1 and 2, the parameter value is SATE;
// When the test type is 3, the parameter value is the risk ratio;
// When the test type is 4, the parameter value is the odds ratio.

   double parameter_value=0;

    if (test_type==1 || test_type==2){

        parameter_value = (v10-v01)/(v11+v10+v01+v00);

    }else if (test_type==3){
        
        if ((v11+v01)<1e-10){

            if (warning_msg){
                cout << "zero in the denominator: imputing 0" << endl;
            }
            parameter_value =  0;

        }else{
            parameter_value = (v11+v10)/(v11+v01);
        }

    }else if(test_type==4){

        if ( abs(v01+v00)<1e-10 || abs(v11+v01)<1e-10){

            if (warning_msg){
                cout << "zero in the denominator: imputing 0" << endl;
            }

            parameter_value =   0;

        }else{
            parameter_value = (v11+v10)/(v01+v00) * (v10+v00)/(v11+v01);
        }
    }

    return parameter_value;
}

double compute_test_statistics(double n11, double n10, double n01, double n00, int test_type, bool warning_msg=false){

    /*
    This function computes the test statistics for the observed data.
    test_type: 1 for two-sided test using a difference-in-mean statistics
               2 for two-sided test using a studentized difference-in-mean statistics
    */
    double test_statistics=0;

    if (test_type==1){
        
        test_statistics = n11/(n11+n10)-n01/(n01+n00);

    }else if(test_type==3){

        if ( n01*(n11+n10)<1e-10){

            if (warning_msg){
                cout << "zero in the denominator: imputing infinity" << endl;
            }
            test_statistics =   std::numeric_limits<double>::infinity();
        }else{
            test_statistics = n11*(n01+n00)/(n01*(n11+n10));
        }

    }else if (test_type==4){

        if (n10*n01<1e-10){

            if (warning_msg){
                cout << "zero in the denominator: imputing infinity" << endl;
            }

            test_statistics =   std::numeric_limits<double>::infinity();
        }else{
            test_statistics = n11*n00/(n10*n01);
        }


    }

    return test_statistics;

}

double compute_test_statistics_STD(double n11, double n10, double n01, double n00, double parameter_value, int test_type, bool warning_msg=false){

    /*
    This function computes the test statistics for the observed data.
    test_type: 1 for two-sided test using a difference-in-mean statistics
               2 for two-sided test using a studentized difference-in-mean statistics
    */
    double test_statistics=0;

    if(test_type==2){

        double T= n11/(n11+n10)-n01/(n01+n00);
        double Q1= n11 * (n11/(n11+n10) - 1)*(n11/(n11+n10) - 1) + n10 *(n11/(n11+n10))*(n11/(n11+n10));
        double Q0= n01 * (n01/(n01+n00) - 1)*(n01/(n01+n00) - 1) + n00 *(n01/(n01+n00))*(n01/(n01+n00));
        double sigma = sqrt(Q1/((n11+n10)*(n11+n10)) + Q0/((n01+n00)*(n01+n00)));
        
        if (sigma<1e-10){

            if (warning_msg){
                cout << "zero in the denominator: imputing infinity" << endl;
            }

            test_statistics =   std::numeric_limits<double>::infinity();
        }else{
            test_statistics = abs((T-parameter_value)/sigma);
        }   

    }

    return test_statistics;

}

bool check_j(int j, double tau0, double n, double n11, double n10, double n01, double n00){

/* 
This function checks the four inequalities given in (B.3) in the paper.
j: an iterator defined in Algorithm 5.3
tau0: SATE (v10-v01)/n
n: total number of units
n11: number of treated units with outcome 1
n10: number of treated units with outcome 0
n01: number of control units with outcome 1
n00: number of control units with outcome 0
*/
    bool check1 = (j>= n*tau0+n01);
    bool check2 = (j>= n11);
    bool check3 = (j<= n-n10);
    bool check4 = (j<= n11+n*tau0+n10+n01);

    if (check1 && check2 && check3 && check4){
      return true;
    }
    else{
      return false;
    }

}

double generate_v10(double j, double tau0,double n11, double n10, double n01, double n00, double n){

/*

j: an iterator defined in Algorithm 5.3
tau0: SATE (v10-v01)/n
n11: number of treated units with outcome 1
n10: number of treated units with outcome 0
n01: number of control units with outcome 1
n00: number of control units with outcome 0
*/
  //Find smallest v10
  double lb [4] ={0,n*tau0,j-n11-n01,n11+n01+n*tau0-j}; // This follows from (B.3).
  double ub [4] ={j, n11+n00, n10+n01+n*tau0,n+n*tau0-j}; // This follows from (B.3).
   
   /*
    cout << "lb: " << lb[0] << endl;
    cout << "ub: " << ub[0] << endl;
    cout << "lb2: " << lb[1] << endl;
    cout << "ub2: " << ub[1] << endl;
    cout << "lb3: " << lb[2] << endl;
    cout << "ub3: " << ub[2] << endl;
    cout << "lb4: " << lb[3] << endl;
    cout << "ub4: " << ub[3] << endl;
*/
    
   double max_lb = *max_element(lb, lb+4);
   double min_ub = *min_element(ub, ub+4);
   //cout << max_lb << endl;
   // cout << min_ub << endl;

    if (min_ub< max_lb){
      //cout << "The intersection is empty. Return n+1." << endl;

      return n+1;  // return n+1 if the intersection of the intervals is empty 
    }
    else{

      return max_lb; // return the maximum of the lower bounds of the intervals if nonempty
    }
  


}

bool LD_check(double n, double n11, double n10, double n01, double n00, double v10, double v01, double v11){

/*
This function checks the compatability of the po tables with data. The inequalities are given in Lemma 2.4 of the paper.
*/
      double lb [4] ={0,n11-v10,v11-n01,v11+v01-n10-n01}; // This follows from (B.3).
      double ub [4] ={v11, n11, v11+v01-n01,n-v10-n01-n10}; // This follows from (B.3).


      double max_lb = *max_element(lb, lb+4);
      double min_ub = *min_element(ub, ub+4);


    if (min_ub< max_lb-1e-10){
      //cout << "The intersection is empty. Return n+1." << endl;
     //cout <<"False" <<endl;
     //cout << max_lb << endl;
     //cout << min_ub << endl;
      return false;  // return n+1 if the intersection of the intervals is empty 
    }
    else{
     //cout <<"True" <<endl;
     //cout << max_lb << endl;
     //cout << min_ub << endl;
      return true; // return the maximum of the lower bounds of the intervals if nonempty
    }
  
}

NumericMatrix sim_q(double K, double v11, double v10, double v01, double v00, double n, double nt){

    

/*
This function creates simulated data for the permutation test.

 q: 8* (num of simulations) matrice of integers. Each row is a Monte Carlo simulation. 
 q[0]: number of always positive units in treated group  
 q[1]: number of always positive units in control group  
 q[2]: number of positive-only-if-treated units in the treated group
 q[3]: number of positive-only-if-treated units in the control group
 q[4]: number of positive-only-if-control units in the treated group
 q[5]: number of positive-only-if-control units in the control group
 q[6]: number of always negative units in the treated group
 q[7]: number of always negative units in the control group 
*/

//initialize the vector of potential outcomes
    NumericVector po_vector(n);
    NumericMatrix q(K, 8); // matrix to store the simulated data

    for (int i = 0; i < n; i++){

        if (i<v11){
            po_vector[i] = 1;
        }
        else if (i<v11+v10){
            po_vector[i] = 2;
        }
        else if (i<v11+v10+v01){
            po_vector[i] = 3;
        }
        else{
            po_vector[i] = 4;
        }
    }

    //cout<<po_vector<<endl;

    for (int i = 0; i < K; i++){
        
        //note this is 1 based, not 0 based
        //sample units that are assigned to treatment
        IntegerVector sample_temp=sample(n, nt, false,  R_NilValue); 

        for (int j=0; j<sample_temp.length(); j++){
            
            if (po_vector[sample_temp[j]-1]==1){
                q(i,0) += 1;
            }
            else if (po_vector[sample_temp[j]-1]==2){
                q(i,2) += 1;
            }
            else if (po_vector[sample_temp[j]-1]==3){
                q(i,4) += 1;
            }
            else{
                q(i,6) += 1;
            }
        }
        q(i,1) = v11-q(i,0);
        q(i,3) = v10-q(i,2);
        q(i,5) = v01-q(i,4);
        q(i,7) = v00-q(i,6);


    }



    return q;

    

}


NumericMatrix sim_q2(double K, double v11, double v10, double v01, double v00, double n, double nt){

    

/*
This function creates simulated data for the permutation test generated using multivariate hypergeometric distribution.

 q: 8* (num of simulations) matrice of integers. Each row is a Monte Carlo simulation. 
 q[0]: number of always positive units in treated group  
 q[1]: number of always positive units in control group  
 q[2]: number of positive-only-if-treated units in the treated group
 q[3]: number of positive-only-if-treated units in the control group
 q[4]: number of positive-only-if-control units in the treated group
 q[5]: number of positive-only-if-control units in the control group
 q[6]: number of always negative units in the treated group
 q[7]: number of always negative units in the control group 
*/

//initialize the vector of potential outcomes
    NumericMatrix q(K, 8); // matrix to store the simulated data

    for (int i = 0; i < K; i++){

     
        q(i,0)=round(R::rhyper(v11, n-v11, nt));

        q(i,1)=round(v11-q(i,0));

        q(i,2)=round(R::rhyper(v10, n-v10-v11, nt-q(i,0)));

        q(i,3)=round(v10-q(i,2));

        q(i,4)=round(R::rhyper(v01, n-v11-v10-v01, nt-q(i,0)-q(i,2)));

        q(i,5)=round(v01-q(i,4));

        q(i,6)=round(nt-q(i,0)-q(i,2)-q(i,4));

        q(i,7)=round(v00-q(i,6));


    }

    return q;

    

}


double compute_v10_max(double v11, double v10, double v01, double v00, double n11, double n10, double n01, double n00, double n){

    /*
    This function finds the largest v10 that is compatible with the data, along the line given in (5.18)
    */
       
    double v10_max = v10; //prescreen largest compatible v10

    double v11_new, v10_new, v01_new;
    for (int k=1;k<n+1;k++){

        //update po table
        v11_new = v11-k; 
        v10_new = v10+k;
        v01_new = v01+k;

        if (LD_check(n, n11, n10, n01, n00, v10_new, v01_new, v11_new)){
           
            v10_max = v10_new;
           
            }

        }

    return v10_max;

    
}

int pre_calculate_size(double n11, double n10, double n01, double n00, double n,  NumericVector cn, double cn_length ){

    
   int size=0;

   double tau0, v10, v01, v11,v00, v10_max;

   for (int i = 0; i < cn_length; i++){
      
      tau0 = cn[i]; // SATE
    for (int j=0; j<n+1; j++){

        if (check_j(j, tau0, n, n11, n10, n01, n00)){ // if j is comptable

            v10 = generate_v10(j, tau0, n11, n10, n01, n00, n); // find the smallest possible v10

            //po tables 
            v01 =v10 - n*tau0;
            v11 = j-v10;
            v00 = n-v11-v10-v01;

            if (LD_check(n, n11, n10, n01, n00, v10, v01, v11)){
                size = size + 1;
            }
                
            //line search
            v10_max = compute_v10_max(v11,v10,v01,v00,n11,n10,n01,n00,n); //compute the largest v10 that is comptable

            for (int k=1; k<(v10_max-v10+1); k++){
                    
                        //update po table
                        v11 = v11-1; 
                        v10 = v10+1;
                        v01 = v01+1;
                        v00 = v00-1;
                    
                        if (LD_check(n, n11, n10, n01, n00, v10, v01, v11)){      
                            size = size+1;         
                        }

                        
            }
        } 
            
    }
    }


    return size;
}

//[[Rcpp::export()]]
NumericMatrix update_q(double K, NumericMatrix q, bool warning_msg=false){

    /*
This function creates simulated data for the permutation test.

 q: 8* (num of simulations) matrice of integers. Each row is a Monte Carlo simulation. 
 q[0]: number of always positive units in treated group  
 q[1]: number of always positive units in control group  
 q[2]: number of positive-only-if-treated units in the treated group
 q[3]: number of positive-only-if-treated units in the control group
 q[4]: number of positive-only-if-control units in the treated group
 q[5]: number of positive-only-if-control units in the control group
 q[6]: number of always negative units in the treated group
 q[7]: number of always negative units in the control group 
*/


    NumericMatrix qnew(K, 8); // matrix to store the simulated data
    qnew = Rcpp::clone(q);
    double p00;
    double p11;
    double B1;
    double B2;
    //cout << "update: " << endl;
    //cout << "q: " << q << endl;
    for (int i = 0; i < K; i++){

        if ( (q(i,6)+q(i,7))>1e-10){

            p00 = q(i,7)/(q(i,6)+q(i,7));

            if (p00<0){
                p00=0;
            }
            
            B1 =  R::rbinom( 1, p00) ;


            if ( abs(B1-1)<1e-10){
                qnew(i,7) = round(qnew(i,7)- 1);
                qnew(i,5) = round(qnew(i,5)+1);
             }else{
                qnew(i,6) = round(qnew(i,6)- 1);
                qnew(i,4) = round(qnew(i,4)+1);
            }

        }

        if ( (q(i,0)+q(i,1))>1e-10){

            p11 = q(i,1)/(q(i,0)+q(i,1));
 
            if (p11<0){
                p11=0;
            }
            B2 =  R::rbinom(  1, p11);
  
            if (abs(B2-1)<1e-10){
                qnew(i,1) = round(qnew(i,1)- 1);
                qnew(i,3) = round(qnew(i,3)+1);
            } else{
                qnew(i,0) = round(qnew(i,0)- 1);
                qnew(i,2) = round(qnew(i,2)+1);
        }
  
        }

    
       if (warning_msg){
 
        if (qnew(i,0)<-1e-10 || qnew(i,1)<-1e-10 || qnew(i,2)<-1e-10 || qnew(i,3)<-1e-10 || qnew(i,4)<-1e-10 || qnew(i,5)<-1e-10 || qnew(i,6)<-1e-10 || qnew(i,7)<-1e-10){
            cout << "Negative values in the po table when updating." << endl;
            cout << q(i,0) << endl;
            cout << q(i,1) << endl;
            cout << q(i,2) << endl;
            cout << q(i,3) << endl;
            cout << q(i,4) << endl;
            cout << q(i,5) << endl;
            cout << q(i,6) << endl;
            cout << q(i,7) << endl;
            }
            
       }
    }

    return qnew;

}

// Algorithm 5.1 in the paper
double perm_test_inner(NumericMatrix q, double n_mc, double test_statistics, double alpha, double epsilon, double v11 , double v10, double v01, double v00,int test_type,bool warning_msg=false){
 

/*
 q: 8* (num of simulations) matrice of integers. Each row is a Monte Carlo simulation. 
 q[0]: number of always positive units in treated group  
 q[1]: number of always positive units in control group  
 q[2]: number of positive-only-if-treated units in the treated group
 q[3]: number of positive-only-if-treated units in the control group
 q[4]: number of positive-only-if-control units in the treated group
 q[5]: number of positive-only-if-control units in the control group
 q[6]: number of always negative units in the treated group
 q[7]: number of always negative units in the control group 

 n_mc: number of Monte Carlo simulations
 test_statistics: test statistics computed from the observed data
 alpha: significance level

 test_type: 1 for two-sided test using a difference-in-mean statistics
            2 for two-sided test using a studentized difference-in-mean statistics
            3 for two-sided test using a risk ratio
            4 for two-sided test using an odds ratio

*/
    // Initialize variables
    double n11,n10,n01,n00 ; // number of treated units with outcome 1
    double test_statistics_mc; // test statistics for each Monte Carlo simulation
    double term_obs, term_mc; // term in the observed test statistics
    double count = 0; // number of Monte Carlo simulations that have test statistics greater than the observed test statistics
    double p_value=0; // p-value



    double parameter_value= compute_parameter_value_perm(v11,v10,v01,v00,test_type,warning_msg); // parameter value


    term_obs = abs(test_statistics - parameter_value); // term in the observed test statistics
    //main loop: calculate test statistics for each Monte Carlo simulation

    for (int i = 0; i < q.nrow(); i++){

      n11 = q(i,0) + q(i,2); 
      n10 = q(i,4) + q(i,6);
      n01 = q(i,1) + q(i,5);
      n00 = q(i,3) + q(i,7);
      
      test_statistics_mc = compute_test_statistics(n11,n10,n01,n00, test_type,warning_msg); // test statistics for each Monte Carlo simulation

      if (test_statistics_mc == std::numeric_limits<double>::infinity()){
        count += 1; // if test statistics for each Monte Carlo simulation is infinity, increment count by 1
      }else{

        term_mc = abs(test_statistics_mc - parameter_value); // term in the Monte Carlo test statistics
        
        if (term_mc >= term_obs){
          count += 1; // if test statistics for each Monte Carlo simulation is weakly greater than the observed test statistics, increment count by 1
        }
      }



   }
   //Calculate p_value
   p_value = count/n_mc; // p-value
   
   //cout << "p_value: " << p_value << endl;
   if ( (p_value + epsilon )< alpha){

     return 0; // reject null hypothesis
   }
   else{

     return 1; //fail to reject null hypothesis
   }
    
}

// Algorithm 5.1 in the paper
double perm_test_inner_STD(NumericMatrix q, double n_mc, double test_statistics_num, double test_statistics_den,double alpha, double epsilon, double v11 , double v10, double v01, double v00,int test_type, bool warning_msg=false){
 

/*
 q: 8* (num of simulations) matrice of integers. Each row is a Monte Carlo simulation. 
 q[0]: number of always positive units in treated group  
 q[1]: number of always positive units in control group  
 q[2]: number of positive-only-if-treated units in the treated group
 q[3]: number of positive-only-if-treated units in the control group
 q[4]: number of positive-only-if-control units in the treated group
 q[5]: number of positive-only-if-control units in the control group
 q[6]: number of always negative units in the treated group
 q[7]: number of always negative units in the control group 

 n_mc: number of Monte Carlo simulations
 test_statistics: test statistics computed from the observed data
 alpha: significance level

 test_type: 1 for two-sided test using a difference-in-mean statistics
            2 for two-sided test using a studentized difference-in-mean statistics
            3 for two-sided test using a risk ratio
            4 for two-sided test using an odds ratio

*/
    // Initialize variables
    double n11,n10,n01,n00 ; // number of treated units with outcome 1
    double test_statistics_mc; // test statistics for each Monte Carlo simulation
    double term_obs, term_mc; // term in the observed test statistics
    double count = 0; // number of Monte Carlo simulations that have test statistics greater than the observed test statistics
    double p_value=0; // p-value

    double parameter_value= compute_parameter_value(v11,v10,v01,v00,test_type,warning_msg); // parameter value


    term_obs = abs( (test_statistics_num - parameter_value)/test_statistics_den); // term in the observed test statistics
    //main loop: calculate test statistics for each Monte Carlo simulation

    for (int i = 0; i < q.nrow(); i++){

      n11 = q(i,0) + q(i,2); 
      n10 = q(i,4) + q(i,6);
      n01 = q(i,1) + q(i,5);
      n00 = q(i,3) + q(i,7);
      
      test_statistics_mc = compute_test_statistics_STD(n11,n10,n01,n00, parameter_value, test_type,warning_msg); // test statistics for each Monte Carlo simulation

      if (test_statistics_mc == std::numeric_limits<double>::infinity()){
        count += 1; // if test statistics for each Monte Carlo simulation is infinity, increment count by 1
      }else{

        term_mc = test_statistics_mc; // term in the Monte Carlo test statistics
        
        if (term_mc >= term_obs){
          count += 1; // if test statistics for each Monte Carlo simulation is weakly greater than the observed test statistics, increment count by 1
        }
      }



   }
   //Calculate p_value
   p_value = count/n_mc; // p-value
   
   //cout << "p_value: " << p_value << endl;
   if ( (p_value + epsilon) < alpha){
     return 0; // reject null hypothesis
   }
   else{
     return 1; //fail to reject null hypothesis
   }
    
}

//[[Rcpp::export()]]
NumericMatrix perm_test_general_interval( double test_statistics, double test_type, double alpha, double epsilon,double n11, double n10, double n01, double n00, double n,double nt, double K, NumericVector cn, double cn_length,bool display_progress=true, bool warning_msg=false){


  Progress p(cn_length, display_progress);
  NumericMatrix q;
  //Initialize output intervals
  double interval_upper_bound, interval_lower_bound;
  interval_upper_bound = -std::numeric_limits<double>::infinity(); 
  interval_lower_bound = std::numeric_limits<double>::infinity();
  
  double parameter_value;

  for (int i = 0; i < cn_length; i++){

    
    double tau0 = cn[i]; // SATE
    //cout << "tau0: " << tau0 << endl;

    for (int j=0; j<n+1; j++){ //iterate over j


        if (check_j(j, tau0, n, n11, n10, n01, n00)){ // if j is comptable

            double v10 = generate_v10(j, tau0, n11, n10, n01, n00, n); // find the smallest v10

            //po tables 
            double v01 =v10 - n*tau0;
            double v11 = j-v10;
            double v00 = n-v11-v10-v01;


            NumericMatrix q = sim_q2(K, v11, v10, v01, v00, n, nt);


            if (LD_check(n, n11, n10, n01, n00, v10, v01, v11)){ //if no v10 compatiable, v10=n+1 and this loop wont be executed 
                

                parameter_value = compute_parameter_value(v11,v10,v01,v00,test_type,warning_msg);

                //cout << "interval_upper_bound:  " << interval_upper_bound << endl;
                //cout << "interval_lower_bound:  " << interval_lower_bound << endl; 

                if ( (parameter_value > interval_upper_bound) || parameter_value < interval_lower_bound){

                    double result = perm_test_inner(q, K, test_statistics, alpha, epsilon,v11,v10,v01,v00, test_type,warning_msg);

                    if (abs(result-1)<1e-10){
                        if (parameter_value > interval_upper_bound){
                            interval_upper_bound = parameter_value;
                        }
                        if (parameter_value < interval_lower_bound){
                            interval_lower_bound = parameter_value;
                        }
                    }
                  
                }

                }
          
                double v10_max = compute_v10_max(v11,v10,v01,v00,n11,n10,n01,n00,n); //compute the largest v10 that is comptable
 
 
            for (int k=1; k<(v10_max-v10+1); k++){ //iterate from v10 to v10_max
                    
                //update po table 
                v11 = v11-1; 
                v10 = v10+1;
                v01 = v01+1;
                v00 = v00-1;
                    
                //q = sim_q(K, v11, v10, v01, v00, n, nt);
                q = update_q(K,q,warning_msg); //update the permutation draws 

                if (LD_check(n, n11, n10, n01, n00, v10, v01, v11)){ // if the po table is comptable

                    //main permutation test

                    parameter_value = compute_parameter_value(v11,v10,v01,v00,test_type,warning_msg);

                    //cout << "interval_upper_bound:  " << interval_upper_bound << endl;
                    //cout << "interval_lower_bound:  " << interval_lower_bound << endl; 
                    if ( (parameter_value > interval_upper_bound) || parameter_value < interval_lower_bound){
                            
                        double result = perm_test_inner(q, K, test_statistics, alpha, epsilon,v11,v10,v01,v00, test_type,warning_msg);
                        if (abs(result-1)<1e-10){
                            if (parameter_value > interval_upper_bound){
                                    interval_upper_bound = parameter_value;
                            }
                            if (parameter_value < interval_lower_bound){
                                    interval_lower_bound = parameter_value;
                            }
                        }

                }
            }   
            
        }

}
    p.increment(); //progress bar
    
    } 

  }
    NumericMatrix  output_vector(2,1);
    output_vector(0,0) = interval_lower_bound;
    output_vector(1,0) = interval_upper_bound;

    
    return output_vector;

}

//[[Rcpp::export()]]
NumericMatrix perm_test_general_set( double test_statistics, double test_type, double alpha, double epsilon, double n11, double n10, double n01, double n00, double n,double nt, double K, NumericVector cn, double cn_length,bool display_progress=true, bool warning_msg=false){



  //Initialize output matrices
  int size = pre_calculate_size(n11, n10, n01, n00, n, cn, cn_length); //calculate the size of the confidence set
  NumericMatrix output_matrix(size,6); // vector to store the result
  //cout << "Initalizing a matrix of  " << size << " times 2" << endl; Don't forget to change me back
  int count = 0;
  
  Progress p(size, display_progress); //tracking progress

  for (int i = 0; i < cn_length; i++){

    double tau0 = cn[i]; // SATE

    for (int j=0; j<n+1; j++){ //iterate over j


        if (check_j(j, tau0, n, n11, n10, n01, n00)){ // if j is comptable

            double v10 = generate_v10(j, tau0, n11, n10, n01, n00, n); // find the smallest v10


            //po tables 
            double v01 =v10 - n*tau0;
            double v11 = j-v10;
            double v00 = n-v11-v10-v01;


            //generate draws
            NumericMatrix q = sim_q2(K, v11, v10, v01, v00, n, nt);

            if (LD_check(n, n11, n10, n01, n00, v10, v01, v11)){ //if no v10 compatiable, v10=n+1 and this loop wont be executed 
                
                double result = perm_test_inner(q, K, test_statistics, alpha, epsilon ,v11,v10,v01,v00,test_type,warning_msg);
                p.increment(); //progress bar

                output_matrix(count,0)=compute_parameter_value(v11,v10,v01,v00,test_type);
                output_matrix(count,1)=result;
                output_matrix(count,2)=v11;
                output_matrix(count,3)=v10;
                output_matrix(count,4)=v01;
                output_matrix(count,5)=v00;
                count += 1;

                }
          
                double v10_max = compute_v10_max(v11,v10,v01,v00,n11,n10,n01,n00,n); //compute the largest v10 that is comptable
 
 
            for (int k=1; k<(v10_max-v10+1); k++){ //iterate from v10 to v10_max
                    
                //update po table 
                v11 = v11-1; 
                v10 = v10+1;
                v01 = v01+1;
                v00 = v00-1;
                    
                q = update_q(K,q,warning_msg); //update the permutation draws

                if (LD_check(n, n11, n10, n01, n00, v10, v01, v11)){ // if the po table is comptable

                    //main permutation test
                    double result = perm_test_inner(q, K, test_statistics, alpha, epsilon,v11,v10,v01,v00, test_type,warning_msg);
                    p.increment(); //progress bar

                    output_matrix(count,0)=compute_parameter_value(v11,v10,v01,v00,test_type,warning_msg);
                    output_matrix(count,1)=result;
                    output_matrix(count,2)=v11;
                    output_matrix(count,3)=v10;
                    output_matrix(count,4)=v01;
                    output_matrix(count,5)=v00;
                    count += 1;
                }
            }   
            
        }

}
    p.increment(); //progress bar
    
}
    //cout << "Initalizing a matrix of  " << size << " times 2" << endl; Don't forget to change me back
    //cout << "count"<< count<<endl;
    //cout << "size" << size<<endl;
    return output_matrix;

}

//[[Rcpp::export()]]
NumericMatrix perm_test_SATE( double test_statistics, double test_type, double alpha, double epsilon,double n11, double n10, double n01, double n00, double n,double nt, double K, NumericVector cn, double cn_length,bool display_progress=true,bool warning_msg=false){


  Progress p(cn_length, display_progress);
   
  NumericMatrix my_matrix(cn_length,4);


  for (int i = 0; i < cn_length; i++){

    
    double tau0 = cn[i]; // SATE
    my_matrix(i,0) = tau0; 
    int result = 0; //result = 0 means it is currently being rejected
    int j = 0;

    while( (j<n+1) && (abs(result)<1e-10) ){ //as soon as we find a j that is not rejected, we stop

        if (check_j(j, tau0, n, n11, n10, n01, n00)){ // if j is comptable

            double v10 = generate_v10(j, tau0, n11, n10, n01, n00, n); // find the smallest v10

            //po tables 
            double v01 =v10 - n*tau0;
            double v11 = j-v10;
            double v00 = n-v11-v10-v01;


            if (v10 < n+1){ // if v10 is comptable

                NumericMatrix q = sim_q2(K, v11, v10, v01, v00, n, nt);

                if (LD_check(n, n11, n10, n01, n00, v10, v01, v11)){ //if no v10 compatiable, v10=n+1 and this loop wont be executed 
                
                    double result = perm_test_inner(q, K, test_statistics, alpha, epsilon,v11,v10,v01,v00,test_type,warning_msg);

                    my_matrix(i,1)=result;

                   // if (result==1){
                   //    my_matrix(i,2) = j;
                   //     my_matrix(i,3) = 0;
                   // }
                }
          
                double v10_max = compute_v10_max(v11,v10,v01,v00,n11,n10,n01,n00,n);
 
            if (v10_max!=v10 && abs(result)<1e-10){ //it is being rejected at the end point so continue

                int k = 1;
                while ( (k<(v10_max-v10+1)) && (abs(result)<1e-10)){
                    
                    //update po table
                    v11 = v11-1; 
                    v10 = v10+1;
                    v01 = v01+1;
                    v00 = v00-1;
                    
                    //update the permutation draws
                    q = update_q(K,q,warning_msg);
                    
                    if (LD_check(n, n11, n10, n01, n00, v10, v01, v11)){ // if the po table is comptable

                        //main permutation test
                        result = perm_test_inner(q, K, test_statistics, alpha, epsilon,v11,v10,v01,v00,test_type,warning_msg);

                        my_matrix(i,1) = result;

                    //if (result==1){
                    //    my_matrix(i,2) = j;
                    //    my_matrix(i,3) = k;
                    //s}

                        }
                    k = k +1;    
                    }
                }
            
            } 
            
        }

    j = j + 1;

}
    p.increment(); //progress bar
}


    return my_matrix;

}

//[[Rcpp::export()]]
double test_tau0_balanced(double tau0, double test_statistics, double test_type, double alpha, double epsilon, double n11, double n10, double n01, double n00, double n,double nt, double K, bool warning_msg=false){


    double result = 0; //result = 0 means it is currently being rejected
    int j = 0;
    

    while( (j<n+1) && (abs(result)<1e-10) ){ //as soon as we find a j that is not rejected, we stop


        if (check_j(j, tau0, n, n11, n10, n01, n00)){ // if j is comptable

            double v10 = generate_v10(j, tau0, n11, n10, n01, n00, n); // find the smallest v10

            //po tables 
            double v01 =v10 - n*tau0;
            double v11 = j-v10;
            double v00 = n-v11-v10-v01;

            if (v10 < n+1){ // if v10 is comptable

                NumericMatrix q = sim_q2(K, v11, v10, v01, v00, n, nt);

                if (LD_check(n, n11, n10, n01, n00, v10, v01, v11)){ //if no v10 compatiable, v10=n+1 and this loop wont be executed 
                
                    result = perm_test_inner(q, K, test_statistics, alpha, epsilon,v11,v10,v01,v00,test_type,warning_msg);

                }

                if ( abs(result)<1e-10 && abs(v10)<1e-10 && abs(v01)<1e-10){ //if the result is rejected and v10=0 and v01=0, then we do one more permutation test

                    v10 = 1;
                    double v01 =v10 - n*tau0;
                    double v11 = j-v10;
                    double v00 = n-v11-v10-v01;

                    NumericMatrix q = sim_q2(K, v11, v10, v01, v00, n, nt);
                    result = perm_test_inner(q, K, test_statistics, alpha, epsilon,v11,v10,v01,v00,test_type,warning_msg);
                }
                       
        }

}
    j = j + 1;

}
    return result;

}

//[[Rcpp::export()]]
NumericMatrix perm_test_SATE_balanced( double test_statistics, double test_type, double alpha, double epsilon,double n11, double n10, double n01, double n00, double n,double nt, double K, NumericVector cn, double cn_length,bool display_progress=true,bool warning_msg=false){


   //Asserts floating point compatibility at compile time
   //static_assert(std::numeric_limits<float>::is_iec559, "IEEE 754 required"); //This line checks our way of defining negative infinity is ok

   double upper_bound_lb=test_statistics;
   double lower_bound_ub=test_statistics;
   double upper_bound_ub = max(cn);
   double lower_bound_lb = min(cn);
   double upper_bound_mid, lower_bound_mid, upper_bound_result,lower_bound_result;

   while ( abs(upper_bound_lb-upper_bound_ub) > (1/n+ 1e-10)){

        upper_bound_mid = round((upper_bound_ub + upper_bound_lb)/2 * n)/n ;



        upper_bound_result = test_tau0_balanced(upper_bound_mid, test_statistics, test_type, alpha, epsilon, n11, n10, n01, n00, n, nt, K, warning_msg);
        if (abs(upper_bound_result-1)<1e-10){ 
            upper_bound_lb = upper_bound_mid;  //if the mid point is accepted, update the lower end point.
        }else{
            upper_bound_ub= upper_bound_mid; //if the mid point is rejected, update the upper end point.
        }
   }

   while (abs(lower_bound_lb-lower_bound_ub) > (1/n+ 1e-10)){

        lower_bound_mid = round((lower_bound_ub + lower_bound_lb)/2 * n)/n ;


        lower_bound_result = test_tau0_balanced(lower_bound_mid, test_statistics, test_type, alpha, epsilon, n11, n10, n01, n00, n, nt, K, warning_msg);
        if (abs(lower_bound_result-1)<1e-10){ 
            lower_bound_ub = lower_bound_mid;  //if the mid point is accepted, update the lower end point.
        }else{
            lower_bound_lb= lower_bound_mid; //if the mid point is rejected, update the upper end point.
        }
   }

    NumericMatrix output_matrix(2,1);

    upper_bound_result = test_tau0_balanced(upper_bound_ub, test_statistics, test_type, alpha, epsilon, n11, n10, n01, n00, n, nt, K, warning_msg);
    if ( abs(upper_bound_result-1)<1e-10){
       output_matrix(1,0) = upper_bound_ub; //if upper bound is accepted, return the upper bound
    }else{
       output_matrix(1,0) = upper_bound_lb; //if upper bound is rejected, return the lower bound
    }


  
    lower_bound_result = test_tau0_balanced(lower_bound_lb, test_statistics, test_type, alpha, epsilon, n11, n10, n01, n00, n, nt, K, warning_msg);
    if ( abs(lower_bound_result-1)<1e-10){
       output_matrix(0,0) = lower_bound_lb; //if lower bound is accepted, return the lower bound
    }else{
       output_matrix(0,0) = lower_bound_ub; //if lower bound is rejected, return the upper bound
     }

    return output_matrix;
     
}

//[[Rcpp::export()]]
NumericMatrix perm_test_SATE_STD( double test_statistics_num, double test_statistics_den, double test_type, double alpha, double epsilon,double n11, double n10, double n01, double n00, double n,double nt, double K, NumericVector cn, double cn_length,bool display_progress=true,bool warning_msg=false){

  
  //This function is used for the studentized difference-in-means statistics. 

  Progress p(cn_length, display_progress);
   
  NumericMatrix my_matrix(cn_length,4);


  for (int i = 0; i < cn_length; i++){

    
    double tau0 = cn[i]; // SATE
    my_matrix(i,0) = tau0; 
    int result = 0; //result = 0 means it is currently being rejected
    int j = 0;

    while( (j<n+1) && (abs(result)<1e-10) ){ //as soon as we find a j that is not rejected, we stop

        if (check_j(j, tau0, n, n11, n10, n01, n00)){ // if j is comptable

            double v10 = generate_v10(j, tau0, n11, n10, n01, n00, n); // find the smallest v10

            //po tables 
            double v01 =v10 - n*tau0;
            double v11 = j-v10;
            double v00 = n-v11-v10-v01;


            if (v10 < n+1){ // if v10 is comptable

                NumericMatrix q = sim_q2(K, v11, v10, v01, v00, n, nt);

                if (LD_check(n, n11, n10, n01, n00, v10, v01, v11)){ //if no v10 compatiable, v10=n+1 and this loop wont be executed 
                
                    double result = perm_test_inner_STD(q, K, test_statistics_num, test_statistics_den, alpha, epsilon,v11,v10,v01,v00,test_type,warning_msg);

                    my_matrix(i,1)=result;

                    //if (result==1){
                    //    my_matrix(i,2) = j;
                    //   my_matrix(i,3) = 0;
                    //}

                }
          
                double v10_max = compute_v10_max(v11,v10,v01,v00,n11,n10,n01,n00,n);
 
            if (v10_max!=v10 && abs(result)<1e-10){ //it is being rejected at the end point so continue

                int k = 1;
                while ( (k<(v10_max-v10+1)) && (abs(result)<1e-10)){
                    
                    //update po table
                    v11 = v11-1; 
                    v10 = v10+1;
                    v01 = v01+1;
                    v00 = v00-1;
                    
                    //update the permutation draws
                    q = update_q(K,q,warning_msg);
                    
                    if (LD_check(n, n11, n10, n01, n00, v10, v01, v11)){ // if the po table is comptable

                        //main permutation test
                        result = perm_test_inner_STD(q, K, test_statistics_num, test_statistics_den, alpha, epsilon,v11,v10,v01,v00,test_type,warning_msg);

                        my_matrix(i,1) = result;

                       // if (result==1){
                       // my_matrix(i,2) = j;
                       // my_matrix(i,3) = k;
                        //}

                        }
                    k = k +1;    
                    }
                }
            
            } 
            
        }

    j = j + 1;

}
    p.increment(); //progress bar
}


    return my_matrix;

}

//[[Rcpp::export()]]
NumericMatrix perm_test_interface(double test_statistics, double test_type, double alpha, double epsilon,double n11, double n10, double n01, double n00, double n, double nt, double K,  NumericVector cn, double cn_length,bool display_progress=true, bool return_interval=true, bool warning_msg=false){

/*
This function performs the permutation test for the general set of parameter values. 
test_statistics: test statistics computed from the observed data
test_type: 1 for two-sided test using a difference-in-mean statistics
           2 for two-sided test using a studentized difference-in-mean statistics
           3 for two-sided test using a risk ratio
           4 for two-sided test using an odds ratio
n11: number of treated units with outcome 1
n10: number of treated units with outcome 0
n01: number of control units with outcome 1
n00: number of control units with outcome 0
n: total number of units
nt: number of treated units
K: number of Monte Carlo simulations
cn: vector of parameter values
cn_length: length of the vector of parameter values
display_progress: boolean to display progress bar
return_interval: boolean to return confidence intervals or confidence sets
*/

    if (test_type==1){

        if ( abs(nt-n/2)<1e-10){
            if (warning_msg==true){
                cout<<"Balanced Sample Size. The program will return a confidence interval. "<<endl;
            }
            NumericMatrix output_matrix = perm_test_SATE_balanced(test_statistics, test_type, alpha, epsilon,n11, n10, n01, n00, n, nt, K, cn, cn_length, display_progress,warning_msg);
            return output_matrix;
        }

        NumericMatrix output_matrix = perm_test_SATE(test_statistics, test_type, alpha, epsilon,n11, n10, n01, n00, n, nt, K, cn, cn_length, display_progress,warning_msg);
        return output_matrix;
    }
    if (return_interval==false){

        if (n>1000){
            cout<<"Large Sample Size (n>1000). The program will return a confidence set. This may require a big memory."<<endl;
        }
        
        //cout<<"The program will return a confidence set. "<<endl;
        NumericMatrix output_matrix = perm_test_general_set(test_statistics, test_type, alpha, epsilon, n11, n10, n01, n00, n, nt, K, cn, cn_length, display_progress,warning_msg);
        return output_matrix;

    }else if (return_interval==true){
         
        if (warning_msg){
            cout<<"The program will return a confidence interval. If a confidence set is desired, include the option return_interval=false "<<endl;
        }
        NumericMatrix output_matrix = perm_test_general_interval(test_statistics, test_type, alpha, epsilon, n11, n10, n01, n00, n, nt, K, cn, cn_length, display_progress,warning_msg);
        return output_matrix;

    }else{

        cout<<"Invalid input"<<endl;
        NumericMatrix output_matrix(1,1);
        return output_matrix;
    }
    

}

//[[Rcpp::export()]]
NumericMatrix perm_test_STD_interface(double test_statistics_num, double test_statistics_den ,double test_type, double alpha, double epsilon,double n11, double n10, double n01, double n00, double n, double nt, double K,  NumericVector cn, double cn_length,bool display_progress=true, bool return_interval=true, bool warning_msg=false){

/*
This function performs the permutation test for the general set of parameter values. 
test_statistics: test statistics computed from the observed data
test_type: 1 for two-sided test using a difference-in-mean statistics
           2 for two-sided test using a studentized difference-in-mean statistics
           3 for two-sided test using a risk ratio
           4 for two-sided test using an odds ratio
n11: number of treated units with outcome 1
n10: number of treated units with outcome 0
n01: number of control units with outcome 1
n00: number of control units with outcome 0
n: total number of units
nt: number of treated units
K: number of Monte Carlo simulations
cn: vector of parameter values
cn_length: length of the vector of parameter values
display_progress: boolean to display progress bar
return_interval: boolean to return confidence intervals or confidence sets
*/
    if (test_type==2){
        NumericMatrix output_matrix = perm_test_SATE_STD(test_statistics_num, test_statistics_den, test_type, alpha, epsilon,n11, n10, n01, n00, n, nt, K, cn, cn_length, display_progress, warning_msg);
        return output_matrix;        
    }else{

        cout<<"Invalid input"<<endl;
        NumericMatrix output_matrix(1,1);
        return output_matrix;
    }

}
