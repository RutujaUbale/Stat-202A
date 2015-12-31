#Sweep Operator
mySweep <- function(A, m)
{
  n <- dim(A)[1]
  c <- dim(A)[2]
  for (k in 1:m) 
  {
    for (i in 1:n)     
      for (j in 1:n)   
        if (i!=k  & j!=k)     
          A[i,j] <- A[i,j] - A[i,k]*A[k,j]/A[k,k]    
        
        for (i in 1:n) 
          if (i!=k) 
            A[i,k] <- A[i,k]/A[k,k]  
          
          for (j in 1:n) 
            if (j!=k) 
              A[k,j] <- A[k,j]/A[k,k]  
            
            A[k,k] <- - 1/A[k,k] 
  }
  return(A)
}

A <- matrix(c(1,2,3,7,11,13,17,21,23), 3,3)
mySweep(A,3)
solve(A)

#Optional 1: Sweep Operator by checking whether the input matrix is a squared matrix and if it is invertible
mySweep <- function(A, m)
{
  n <- dim(A)[1]
  c <- dim(A)[2]
  try (if (n != c) stop ("The matrix is not a square matrix"))
  for (k in 1:m) 
  {
    try (if (abs(A[k,k]) < 0.001) stop ("The matrix is not invertible"))
    
    for (i in 1:n)     
      for (j in 1:n)   
        if (i!=k  & j!=k)     
          A[i,j] <- A[i,j] - A[i,k]*A[k,j]/A[k,k]    
        
        for (i in 1:n) 
          if (i!=k) 
            A[i,k] <- A[i,k]/A[k,k]  
          
          for (j in 1:n) 
            if (j!=k) 
              A[k,j] <- A[k,j]/A[k,k]  
            
            A[k,k] <- - 1/A[k,k] 
  }
  return(A)
}

A <- matrix(c(1,2,3,7,11,13,17,21,23), 3,3)
mySweep(A,3)
solve(A)

#Optional 2: To find the Determinant of the Matrix
myDet <- function(A)
{
  n <- dim(A)[1]
  c <- dim(A)[2]
  Det <- 1
  try (if (n != c) stop ("The matrix is not a square matrix"))
  for (k in 1:n) 
  {
    try (if (abs(A[k,k]) < 0.001) stop ("The matrix is not invertible"))
    Det <- Det * A[k,k]
    for (i in 1:n)     
      for (j in 1:n)   
        if (i!=k  & j!=k)     
          A[i,j] <- A[i,j] - A[i,k]*A[k,j]/A[k,k]    
    
    for (i in 1:n) 
      if (i!=k) 
        A[i,k] <- A[i,k]/A[k,k]  
    
    for (j in 1:n) 
      if (j!=k) 
        A[k,j] <- A[k,j]/A[k,k]  
    
    A[k,k] <- - 1/A[k,k] 
  }
  return(Det)
}

A = matrix(c(1,2,3,7,11,13,17,21,23), 3,3)
myDet(A)

#Optional 3: Reverse Sweep operator
mySweep <- function(A, m)
{
  n <- dim(A)[1]
  c <- dim(A)[2]
  try (if (n != c) stop ("The matrix is not a square matrix"))
  for (k in 1:m) 
  {
    try (if (abs(A[k,k]) < 0.001) stop ("The matrix is not invertible"))
    for (i in 1:n)     
      for (j in 1:n)   
        if (i!=k  & j!=k)     
          A[i,j] <- A[i,j] - A[i,k]*A[k,j]/A[k,k]    
        
        for (i in 1:n) 
          if (i!=k) 
            A[i,k] <- A[i,k]/A[k,k]  
          
          for (j in 1:n) 
            if (j!=k) 
              A[k,j] <- A[k,j]/A[k,k]  
            
            A[k,k] <- - 1/A[k,k] 
  }
  return(A)
}

A = matrix(c(1,2,3,7,11,13,17,21,23), 3,3)
B <- mySweep(A,3)
solve(A)

revSweep <- function(A, m)
{
  n <- dim(A)[1]
  
  for (k in 1:m) 
  {
    
    for (i in 1:n)     
      for (j in 1:n)   
        if (i!=k  & j!=k)     
          A[i,j] <- A[i,j] - A[i,k]*A[k,j]/A[k,k]    
        
        for (i in 1:n) 
          if (i!=k) 
            A[i,k] <- A[i,k]/A[k,k]  
          
          for (j in 1:n) 
            if (j!=k) 
              A[k,j] <- A[k,j]/A[k,k]  
            
            A[k,k] <- - 1/A[k,k] 
  }
  return(A)
}
revSweep(B,3)

#Gauss-Jordan elimination plain version
myGaussJordan <- function(A,m)
{
  n<- dim(A)[1]
  B<- cbind(A, diag(rep(1,n)))
  for (k in 1:m)
  {
    a <- B[k,k]
    for (j in 1:(n*2))
      B[k,j] = B[k,j]/a
    for (i in 1:n)
      if (i != k)
      {
        a <- B[i,k]
        for (j in 1:(n*2))
          B[k,j] <- B[i,j] - B[k,j]*a
      }
  }
  return (B)
}
A = matrix(c(1,2,3,7,11,13,17,21,23), 3,3)
solve(A)
myGaussJordan(A,3)

#Gauss-Jordan elimination vectorized form
myGaussJordan_vec <- function(A,m)
{
  n <- dim(A)[1]
  B<- cbind(A, diag(rep(1,n)))
  for(k in 1:m)
  {
    B[k, ] <- B[k, ]/B[k,k]
    for (i in 1:n)
      if (i!=k)
      {
        B[i, ] <- B[1, ] - B[k, ]*B[i,k]
      }
  }
  return(B)
}
A = matrix(c(1,2,3,7,11,13,17,21,23), 3,3)
solve(A)
myGaussJordan_vec(A,3)
  


