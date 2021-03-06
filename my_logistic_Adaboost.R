library(data.table)# allows us to use function fread,
# which quickly reads data from csv files 

myAdaboost <- function(X, Y, Xtest, Ytest, m=100, num_iterations=4000,
                  learning_rate=1e-1 , lambda=0.00000001) 
{
  n = dim(X)[1] 
  p = dim(X)[2]+1
  X1 = cbind(rep(1, n), X)
  th = 0.8
  X1 = X1 > th
  X1 = 2*X1-1
  Y = 2*Y - 1
  
  X1test = cbind(rep(1, n), Xtest)
  X1test = X1test > th
  X1test = 2*X1test - 1
  Ytest = 2*Ytest - 1
  
  beta = matrix(rep(0,p), nrow=p)
  
  for (i in 1:num_iterations)
  {
    S = X1%*%beta
    Stest = X1test%*%beta
    W = exp(-Y*S)
    a=matrix(rep(1,n), nrow=1)%*%(matrix(W*Y,n,p)*X1)/n
    e=(1-a)/2
    j=which.min(e)
    beta[j]=beta[j]+log((1-e[j])/e[j])/2
  }
  acc_train = mean(sign(S)*Y)
  acc_test = mean(sign(Stest)*Ytest)
  
  #Training Accuracies
  #accuracies = append(accuracies,acc)
  
  #predict_test = 1/(1+exp(-(Xtest1%*%beta)))
  #Testing Accuracies
  #accuracies2 = append(accuracies2,accuracy(predict_test, Ytest))
  
  return(list(beta, acc_train, acc_test))
}


# load data
load_digits <- function(subset=NULL, normalize=TRUE) {
  
  #Load digits and labels from digits.csv.
  
  #Args:
  #subset: A subset of digit from 0 to 9 to return.
  #If not specified, all digits will be returned.
  #normalize: Whether to normalize data values to between 0 and 1.
  
  #Returns:
  #digits: Digits data matrix of the subset specified.
  #The shape is (n, p), where
  #n is the number of examples,
  #p is the dimension of features.
  #labels: Labels of the digits in an (n, ) array.
  #Each of label[i] is the label for data[i, :]
  
  # load digits.csv, adopted from sklearn.
  
  df <- fread("Downloads/digits.csv") 
  df <- as.matrix(df)
  
  ## only keep the numbers we want.
  if (length(subset)>0) {
    
    c <- dim(df)[2]
    l_col <- df[,c]
    index = NULL
    
    for (i in 1:length(subset)){
      
      number = subset[i]
      index = c(index,which(l_col == number))
    }
    sort(index)
    df = df[index,]
  }
  
  # convert to arrays.
  digits = df[,-1]
  labels = df[,c]
  
  # Normalize digit values to 0 and 1.
  if (normalize == TRUE) {
    digits = digits - min(digits)
    digits = digits/max(digits)}
  
  
  # Change the labels to 0 and 1.
  for (i in 1:length(subset)) {
    labels[labels == subset[i]] = i-1
  }
  
  return(list(digits, labels))
  
}

split_samples <- function(digits,labels) {
  
  # Split the data into a training set (70%) and a testing set (30%).
  
  num_samples <- dim(digits)[1]
  num_training <- round(num_samples*0.7)
  indices = sample(1:num_samples, size = num_samples)
  training_idx <- indices[1:num_training]
  testing_idx <- indices[-(1:num_training)]
  
  return (list(digits[training_idx,], labels[training_idx],
               digits[testing_idx,], labels[testing_idx]))
}


#====================================
# Load digits and labels.
result = load_digits(subset=c(1, 7), normalize=TRUE)
digits = result[[1]]
labels = result[[2]]

result = split_samples(digits,labels)
training_digits = result[[1]]
training_labels = result[[2]]
testing_digits = result[[3]]
testing_labels = result[[4]]

# print dimensions
length(training_digits)
length(testing_digits)

# Train a net and display training accuracy.
final_results=myAdaboost(training_digits, training_labels, testing_digits, testing_labels)
beta= final_results[[1]]
beta
train_accuracy= final_results[[2]]
train_accuracy
test_accuracy= final_results[[3]]
test_accuracy
