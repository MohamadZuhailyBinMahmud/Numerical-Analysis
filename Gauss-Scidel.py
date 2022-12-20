### import packages ###
import numpy as np

### define values ###
size = int(5) # matrix size including the border
top=100 # top edge value
right=100 # right edge value
bot=50 # bottom edge value
left=75 # left edge value
err_test=0.01 # % error in decimal

### define function ###

def fx(t,b,l,r): # numerical solution to the laplace equation
    return (1/4)*(r+l+t+b)

def init(n,t,b,l,r): # create the problem matrix
    MatA = np.zeros((n,n)) # create a matrix filled with zeros
    for i in range(0,n): # fill the edges of the matrix with values
        MatA[0][i] = t # assign top value 
        MatA[n-1][i] = b # assign bottom value
        MatA[i][0] = l #assign left value
        MatA[i][n-1] = r # assign right value
    MatA[0][0]=(t+l)/2 # top-left corner value
    MatA[0][n-1]=(t+r)/2 # top-right corner value
    MatA[n-1][0]=(l+b)/2 # bottom-left corner value
    MatA[n-1][n-1]=(r+b)/2 # bottom-right corner value
    return MatA
    
def gauss_sidel(MatA,err_test,n): # solve the matrix using gauss sidel method
    err = np.ones((n-2,n-2)) # create error matrix to be filled up with % error
    
    for i in range(1,n-1): # loop over the rows of the matrix
        for j in range(1,n-1): #loop over the column of the matrix
            x = fx(MatA[i+1][j],MatA[i-1][j],MatA[i][j+1],MatA[i][j-1]) # solve the value for particular elements in the matrix
            if MatA[i][j]!=0: # to avoid zero division
                err[i-1][j-1]=abs((MatA[i][j]-x)/MatA[i][j]) # calculate the error as convergence indicator 
            MatA[i][j]=x  # save the new value back 
    
                             
    return MatA,err

### initialize the matrix ###

Mat=init(size,top,bot,left,right) # create the matrix that need to be solved
print(Mat,'\n\n')

### run the numerical calculation ###

final , error = gauss_sidel(Mat,0.01,size) # first run to initialize the final and error array
count = 1 # initialize the iteration count starting from 1

while (any(x>=err_test for x in error[0])) or (any(x>=err_test for x in error[1])) or (any(x>=err_test for x in error[2])):
    # iterate as long as the condition is true (which is > 1% error)
    final , error = gauss_sidel(Mat,0.01,size)
    count+=1
    
print('Results\n',final,'\n\n Error for each elements in the matrix \n',error,'\n\n number of iteration :',count)
