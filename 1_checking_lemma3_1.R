###########################################################################
LD = function(v11,v10,v01,v00,n11,n10,n01,n00){
  
  max_lb = max(c(0,n11-v10,v11-n01,v11+v01-n10-n01))
  min_ub = min(c(v11,n11,v11+v01-n01,n-v10-n01-n10))
  
  if (min_ub<max_lb-1e-10){
    
    return(FALSE)
  }else{
    
    return(TRUE)
  }
}

Test_SATE = function(SATE,n11,n10,n01,n00,dim){
  
  max_count=0
  possible_output=FALSE
  
  n = n11+n10+n01+n00
  
  for (j in 0:n){
    
    for (v10 in 0:n){
      
      
      v11 = j - v10
      v01 = v10-n*SATE
      v00 = n-j-v10+n*SATE
      
      po_table = c(v11,v10,v01,v00)
      possible = LD(v11,v10,v01,v00,n11,n10,n01,n00)

      count=0
      if (possible==TRUE ){
        
        possible_output=TRUE
        y1=c(rep(1,v11),rep(1,v10),rep(0,v01),rep(0,v00))
        y0=c(rep(1,v11),rep(0,v10),rep(1,v01),rep(0,v00))
        
        
        for (i in 1:nrow(assignment)){
            
            dim_temp = sum(y1*assignment[i,])/7 - sum(y0*assignment[i,])/2
            if (abs(dim_temp-SATE)>=abs(dim-SATE)){
              
              count = count + 1
            }
          
        max_count = max(c(max_count,count))
        print(SATE)
        print(po_table)
        print(count)
        }
        
        
      }
      
      
    }
    
  }
  
  return(c(max_count,possible_output))
  
}

#data
n=9
nt=7
nc=2
v11=0
v10=1
v00=6
v01=0
true_SATE=v10/n-v01/n

n11=1
n10=6
n01=0
n00=2
dim=n11/7-n01/2

#generate all possible assignment vectors
x=1:n
assignment_to_control=t(combn(x, 2)) #36 possible assignments
assignment = matrix(1,nrow(assignment_to_control),9)
for (i in 1:nrow(assignment_to_control)){
  assignment[i,assignment_to_control[i,]]=0
}

#sate grid
sate_grid = ((n10-n01-nt):(n10-n01-nt+n))
count=matrix(0,length(sate_grid),3)
count[,1]=sate_grid
for (i in 1:length(sate_grid)){
  
  SATE=sate_grid[i]/n
  
  output=Test_SATE(SATE,n11,n10,n01,n00,dim)
  count[i,2]=output[1]
  count[i,3]=output[2]
  }





