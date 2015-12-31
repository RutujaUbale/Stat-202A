# coding: utf-8

# In[1]:

#Sweep Operator 
import numpy as np
def mySweep(B, m):
    A = np.copy(B)
    n, c = A.shape  
    for k in range(m):
        for i in range(n):
            for j in range(n):
                if i!=k and j!=k:
                    A[i, j] = A[i, j] - A[i, k]*A[k,j]/ A[k,k]
        for i in range(n):
            if i!=k:
                A[i, k] = A[i,k]/A[k,k]
        for j in range(n):
            if j!=k:
                A[k, j]= A[k,j]/A[k,k]
        A[k,k] = -1/A[k,k]
    return A

A = np.array([[1,2,3],[7,11,13],[17,21,23]], dtype=float)
print mySweep(A, 3)


#Optional 1: Sweep Operator by checking whether the input matrix is a squared matrix and if it is invertible
import numpy as np
import sys
def mySweep(B, m):
    A = np.copy(B)
    n, c = A.shape
    if n != c:
        sys.exit("The matrix is not a square matrix")   
    for k in range(m):
        if abs(A[k,k]) < 0.001:
            sys.exit("The matrix is not invertible")
        for i in range(n):
            for j in range(n):
                if i!=k and j!=k:
                    A[i, j] = A[i, j] - A[i, k]*A[k,j]/ A[k,k]
        for i in range(n):
            if i!=k:
                A[i, k] = A[i,k]/A[k,k]
        for j in range(n):
            if j!=k:
                A[k, j]= A[k,j]/A[k,k]
        A[k,k] = -1/A[k,k]
    return A

A = np.array([[1,2,3],[7,11,13],[17,21,23]], dtype=float)
print mySweep(A, 3)


#Optional 2: To find the Determinant of the Matrix
def myDet(B):
    A = np.copy(B)
    n, c = A.shape
    Det = 1
    if n != c:
        sys.exit("The matrix is not a square matrix") 
    for k in range(n):
        Det = Det * A[k,k]
        if abs(A[k,k]) < 0.001:
            sys.exit("The matrix is not invertible")
        for i in range(n):
            for j in range(n):
                if i!=k and j!=k:
                    A[i, j] = A[i, j] - A[i, k]*A[k,j]/ A[k,k]
        for i in range(n):
            if i!=k:
                A[i, k] = A[i,k]/A[k,k]
        for j in range(n):
            if j!=k:
                A[k, j]= A[k,j]/A[k,k]
        A[k,k] = -1/A[k,k]
    return Det

A = np.array([[1,2,3],[7,11,13],[17,21,23]], dtype=float)
print myDet(A)


#Optional 3: Reverse Sweep operator
def mySweep(B, m):
    A = np.copy(B)
    n, c = A.shape
    if n != c:
        sys.exit("The matrix is not a square matrix")   
    for k in range(m):
        if abs(A[k,k]) < 0.001:
            sys.exit("The matrix is not invertible")
        for i in range(n):
            for j in range(n):
                if i!=k and j!=k:
                    A[i, j] = A[i, j] - A[i, k]*A[k,j]/ A[k,k]
        for i in range(n):
            if i!=k:
                A[i, k] = A[i,k]/A[k,k]
        for j in range(n):
            if j!=k:
                A[k, j]= A[k,j]/A[k,k]
        A[k,k] = -1/A[k,k]
    return A

A = np.array([[1,2,3],[7,11,13],[17,21,23]], dtype=float)
B = mySweep(A, 3)

def revSweep(B, m):
    A = np.copy(B)
    n, c = A.shape
    if n != c:
        sys.exit("The matrix is not a square matrix")   
    for k in range(m):
        if abs(A[k,k]) < 0.001:
            sys.exit("The matrix is not invertible")
        for i in range(n):
            for j in range(n):
                if i!=k and j!=k:
                    A[i, j] = A[i, j] - A[i, k]*A[k,j]/ A[k,k]
        for i in range(n):
            if i!=k:
                A[i, k] = A[i,k]/A[k,k]
        for j in range(n):
            if j!=k:
                A[k, j]= A[k,j]/A[k,k]
        A[k,k] = -1/A[k,k]
    return A
revSweep(B, 3)


#Gauss-Jordan elimination plain version
import numpy as np
def myGaussJordan(A, m):
    n = A.shape[0]
    B  = np.hstack((A, np.identity(n)))
    
    for k in range(m):
        a = B[k, k]
        for j in range(n*2):
            B[k, j] = B[k, j] / a
        for i in range(n):
            if i != k:
                a = B[i, k]
                for j in range(n*2):
                    B[i, j] = B[i, j] - B[k, j]*a;
    return B


#Gauss-Jordan elimination vectorized form
def myGaussJordanVec(A, m):
    n = A.shape[0]
    B  = np.hstack((A, np.identity(n)))
    
    for k in range(m):
        B[k, :] = B[k, ] / B[k, k]
        for i in range(n):
            if i != k:
                B[i, ] = B[i, ] - B[k, ]*B[i, k];
    
    return B

A = np.array([[1,2,3],[7,11,13],[17,21,23]], dtype=float).T
print myGaussJordan(A, 3)
print myGaussJordanVec(A, 3)




