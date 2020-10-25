readTsmInstance <- function(file, unitaryCost=FALSE) {
  chr<-scan(file,sep=";",what=character())
  index <- 1
  testCosts <- numeric()
  coverageMatrix <- numeric()
  
  # Read number of tests
  testsAndElements<-strsplit(chr[index]," ")[[1]]
  numTests <- as.numeric(testsAndElements[1])
  numElements <- as.numeric(testsAndElements[2])
  for (test in 1:numTests) {
    mask <- chr[index+2*test]
    if (length(mask)!=numElements) {
      simpleCondition("Error in the number of columns")
    }
    coverageMatrix<-c(coverageMatrix,as.numeric(strsplit(mask,"")[[1]]))
    testCosts <- c(testCosts, as.numeric(chr[index+2*test+1]))
  }
  
  if (unitaryCost) {
    testCosts <- rep(1,numTests)
  }
  
  coverageMatrix<-matrix(coverageMatrix,ncol=numElements,byrow=TRUE)
  result <- list(testCosts, coverageMatrix)
  names(result) <- c("testCosts", "coverageMatrix")
  return(result)
}

ilpModel4Tsm <- function(tsmInstance, costUpperBound=NULL,covLowerBound=NULL) {
  n <- dim(tsmInstance$coverageMatrix)[1]
  m <- dim(tsmInstance$coverageMatrix)[2]
  if (is.null(costUpperBound) && is.null(covLowerBound)) {
    stop("One of cost or coverage bounds must be specified")
  }
  if (!is.null(covLowerBound) && !is.null(costUpperBound)) {
    stop("Only one of cost or coverage bound must be specified")
  }
  
  constraints <- numeric()
  for (elem in 1:m) {
    row <- c(-tsmInstance$coverageMatrix[,elem],rep(0,m))
    row[n+elem] <- 1
    constraints <- c(constraints,row)
  }
  
  direction <- rep("<=",m)
  rightHandSide <- rep(0,m)
  
  coverageVector <- c(rep(0,n),rep(1,m))
  costVector <- c(tsmInstance$testCosts,rep(0,m))
  
  if (is.null(covLowerBound)) {
    max <- TRUE
    objectiveVector <- coverageVector
    constraints <- c(constraints,costVector)
    direction <- c(direction,"<=")
    rightHandSide <- c(rightHandSide,costUpperBound)
  } else {
    max <- FALSE
    objectiveVector <- costVector
    constraints <- c(constraints,coverageVector)
    direction <- c(direction, ">=")
    rightHandSide <- c(rightHandSide, covLowerBound)
  }
  
  # Prepare the output
  constraintMatrix <- matrix(constraints,ncol=n+m,byrow=TRUE)
  types <- rep("B", times=n+m)
  
  result <- list(objectiveVector, constraintMatrix, direction, rightHandSide, types, max)
  names(result) <- c("obj", "mat", "dir", "rhs", "types", "max")
  return(result)
}

solveModel <- function(tsmModel) {
  result<- Rsymphony_solve_LP(tsmModel$obj,tsmModel$mat,tsmModel$dir,tsmModel$rhs,types=tsmModel$types,max=tsmModel$max)
  return(result)
}

computeParetoFront <- function (tsmInstance) {
  # TODO
}

reduceInstance <- function(tsmInstance) {
  testCosts <- tsmInstance$testCosts
  coverageMatrix <- tsmInstance$coverageMatrix
  n <- dim(coverageMatrix)[1]
  m <- dim(coverageMatrix)[2]
  
  i<-1
  while (i <= length(testCosts)) {
    removeI <- FALSE
    j <- i+1
    while (j <= length(testCosts) && !removeI) {
      dom <- dominance(testCosts,coverageMatrix,i,j)
      if (dom < 0) {
        testCosts <- testCosts[-i]
        coverageMatrix <- coverageMatrix[-i,]
        removeI <- TRUE
      } else if (dom > 0) {
        testCosts <- testCosts[-j]
        coverageMatrix <- coverageMatrix[-j,]
      } else {
        j <- j+1
      }
    }
    if (!removeI) {
      i <- i+1
    }
  }
  result <- list(testCosts, coverageMatrix)
  names(result) <- c("testCosts", "coverageMatrix")
  return(result)
}

dominance <- function(testCosts, coverageMatrix, i, j) {
  if (dominates(testCosts, coverageMatrix, i, j)) {
    return(1)
  } else if (dominates(testCosts, coverageMatrix, j, i)) {
    return(-1)
  }
  return(0)
}

dominates <- function(testCosts, coverageMatrix,i,j) {
  return (all(coverageMatrix[i,] >= coverageMatrix[j,]) && (testCosts[i] <= testCosts[j]))
}