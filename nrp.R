readNrpInstance <- function(file) {
  numbers <- scan(file)
  index <- 1
  requirementCosts <- numeric()
  
  # Read requirement costs
  levels <- numbers[index]
  index <- index+1
  for (level in 1:levels) {
    numberOfRequirementsInLevel <- numbers[index]
    requirementCosts <- c(requirementCosts,numbers[(1:numberOfRequirementsInLevel)+index])
    index <- index+numberOfRequirementsInLevel+1
  }
  
  # Read dependencies
  numberOfDependencies <- numbers[index]
  if (numberOfDependencies > 0) {
    dependencies <- matrix(numbers[index+1:(2*numberOfDependencies)],ncol = 2,byrow=TRUE)
  } else {
    dependencies <- matrix(0,nrow=0,ncol=2)
  }
  index <- index+2*numberOfDependencies+1;
  
  # Read stakeholders
  numberOfStakeholders <- numbers[index]
  index <- index+1;
  stakeholdersWeights <- numeric()
  stakeholdersRequirements <- list()
  
  for (stakeholder in 1:numberOfStakeholders) {
    stakeholdersWeights[stakeholder] <- numbers[index]
    index <- index+1
    requirementsOfStakeholder <- numbers[index]
    stakeholdersRequirements[[stakeholder]] <- numbers[(1:requirementsOfStakeholder)+index]
    index <- index + requirementsOfStakeholder + 1;
  }
  result <- list(requirementCosts, dependencies, stakeholdersWeights, stakeholdersRequirements)
  names(result) <- c("requirementCosts", "dependencies", "stakeholdersWeights", "stakeholdersRequirements")
  return(result)
}

ilpModel <- function(nrpInstance, budgetLimitFraction) {
  n <- length(nrpInstance$requirementCosts)
  m <- length(nrpInstance$stakeholdersWeights)
  budget <- budgetLimitFraction * sum(nrpInstance$requirementCosts)
  # Set objective expression
  objectiveVector <- c(rep(0,n),nrpInstance$stakeholdersWeights)
  # Set cost constraint
  constraints <- c(nrpInstance$requirementCosts,rep(0,times=m))
  nrows <- 1
  # Set dependencies constraints
  if (dim(nrpInstance$dependencies)[1] > 0) {
    for (dep in 1:dim(nrpInstance$dependencies)[1]) {
      row <- rep(0,n+m)
      row[nrpInstance$dependencies[dep,1]] <- -1
      row[nrpInstance$dependencies[dep,2]] <- 1
      constraints <- c(constraints,row)
      nrows <- nrows+1;
    }
  }
  # Set stakeholders constraints
  for (stakeholder in 1:length(nrpInstance$stakeholdersRequirements)) {
    for (requirement in nrpInstance$stakeholdersRequirements[[stakeholder]]) {
      row <- rep(0,n+m)
      row[requirement] <- -1
      row[stakeholder+n] <- 1
      constraints <- c(constraints,row)
      nrows <- nrows+1;
    }
  }
  # Prepare the output
  constraintMatrix <- matrix(constraints,nrow=nrows,byrow=TRUE)
  direction <- rep("<=",nrows);
  rightHandSide <- c(budget,rep(0,nrows-1))
  types <- rep("B", times=n+m)
  
  result <- list(objectiveVector, constraintMatrix, direction, rightHandSide, types, TRUE)
  names(result) <- c("obj", "mat", "dir", "rhs", "types", "max")
  return(result)
}


ilpModelWithComponents <- function(nrpInstance, budgetLimitFraction, stakeholders=NULL) {
  if (is.null(stakeholders)) {
    stakeholders = 1:length(nrpInstance$stakeholdersWeights)
    #print (stakeholders)
  }
  n <- length(nrpInstance$requirementCosts)
  m <- length(stakeholders)
  budget <- budgetLimitFraction * sum(nrpInstance$requirementCosts)
  # Set objective expression
  objectiveVector <- c(rep(0,n),nrpInstance$stakeholdersWeights[stakeholders])
  # Set cost constraint
  constraints <- c(nrpInstance$requirementCosts,rep(0,times=m))
  nrows <- 1
  # Set dependencies constraints
  if (dim(nrpInstance$dependencies)[1] > 0) {
    for (dep in 1:dim(nrpInstance$dependencies)[1]) {
      row <- rep(0,n+m)
      row[nrpInstance$dependencies[dep,1]] <- -1
      row[nrpInstance$dependencies[dep,2]] <- 1
      constraints <- c(constraints,row)
      nrows <- nrows+1;
    }
  }
  # Set stakeholders constraints
  for (stakeholder in 1:m) {
    for (requirement in nrpInstance$stakeholdersRequirements[[stakeholders[stakeholder]]]) {
      row <- rep(0,n+m)
      row[requirement] <- -1
      row[stakeholder+n] <- 1
      constraints <- c(constraints,row)
      nrows <- nrows+1;
    }
  }
  # Prepare the output
  constraintMatrix <- matrix(constraints,nrow=nrows,byrow=TRUE)
  direction <- rep("<=",nrows);
  rightHandSide <- c(budget,rep(0,nrows-1))
  types <- rep("B", times=n+m)
  
  result <- list(objectiveVector, constraintMatrix, direction, rightHandSide, types, TRUE, stakeholders)
  names(result) <- c("obj", "mat", "dir", "rhs", "types", "max", "stakeholders")
  return(result)
}

computeAllRequirementsForStakeholders <- function(nrpInstance) {
  # Initialize data structures
  n <- length(nrpInstance$requirementCosts)
  m <- length(nrpInstance$stakeholdersWeights)
  unprocessedDirectDependencies <- numeric(0)
  toProcess <- numeric(0)
  dependents <- matrix(FALSE, nrow=n, ncol=n)
  dependencies <- matrix(FALSE, nrow=n, ncol=n)
  
  dependents[nrpInstance$dependencies] <- TRUE
  for (d in 1:n) {
    unprocessedDirectDependencies[d] = sum(dependents[,d])
    if (unprocessedDirectDependencies[d] == 0) {
      toProcess <- c(toProcess, d)
    }
  }
  # Start loop
  while (length (toProcess) > 0) {
    req <- toProcess[1]
    toProcess <- toProcess[-1]
    
    for (dep in (1:n)[dependents[req,]]) {
      dependencies[dep,] <- dependencies[dep,] | dependencies[req,]
      dependencies[dep,req] <- TRUE
      unprocessedDirectDependencies[dep] <- unprocessedDirectDependencies[dep]-1
      if (unprocessedDirectDependencies[dep] == 0) {
        toProcess <- c(toProcess, dep)
      }
    }
  }
  if (sum(unprocessedDirectDependencies) > 0) {
    print("Error!! Not all nodes have been processed")
  }
  
  stakeholdersDependencies <- list()
  for (stakeholder in 1:m) {
    requirements <- rep(FALSE, times=n)
    for (requirement in nrpInstance$stakeholdersRequirements[[stakeholder]]) {
      requirements <- requirements | dependencies[requirement,]
      requirements[requirement] <- TRUE
    }
    stakeholdersDependencies[[stakeholder]] <- (1:n)[requirements]
  }
  
  return(stakeholdersDependencies)
}

generateLocalOptima <- function(nrpInstance, budgetLimitFraction) {
  n <- length(nrpInstance$requirementCosts)
  m <- length(nrpInstance$stakeholdersWeights)
  budget <- budgetLimitFraction * sum(nrpInstance$requirementCosts)
  
  stakeholderRequirements <- computeAllRequirementsForStakeholders(nrpInstance)
  
  solution <- rep(FALSE,times=n+m)
  cost <- 0
  value <- 0
  
  for (stakeholder in sample(m)) {
    requirements <- stakeholderRequirements[[stakeholder]]
    requirementsToAdd <- requirements[!solution[requirements]]
    addedCost <- sum(nrpInstance$requirementCosts[requirementsToAdd])
    if (addedCost+cost <= budget) {
      value <- value + nrpInstance$stakeholdersWeights[stakeholder]
      cost <- cost + addedCost
      solution[n+stakeholder] <- TRUE
      solution[requirementsToAdd] <- TRUE
    }
  }
  result <- list((1:n)[solution[1:n]], (1:m)[solution[n+(1:m)]], cost, value)
  names(result) <- c("requirements", "stakeholders", "cost", "value")
  return(result)
}

cmsa <- function(nrpInstance, budgetLimitFraction, na, maxAge, ilpTimeLimit, totalTimeLimit) {
  n <- length(nrpInstance$requirementCosts)
  m <- length(nrpInstance$stakeholdersWeights)
  age <- rep(maxAge, times=m)
  
  startTime = proc.time()[3]
  while (proc.time()[3]-startTime < totalTimeLimit) {
    # Lines 5 to 10 of CMSA
    generatedStakeholders <- numeric()
    for (i in 1:na) {
      solution <- generateLocalOptima(nrpInstance, budgetLimitFraction)
      generatedStakeholders <- union(generatedStakeholders, solution$stakeholders)
    }
    for (stakeholder in generatedStakeholders) {
      if (age[stakeholder] >= maxAge) {
        age[stakeholder] <- 0
      }
    }
    
    # Apply exact solver (Line 11 of CMSA)
    stakeholders <- (1:m)[age < maxAge]
    mprime <- length(stakeholders)
    model <- ilpModelWithComponents(nrpInstance, budgetLimitFraction, stakeholders)
    ilpSolution <- Rsymphony_solve_LP(obj=model$obj, 
                                   mat=model$mat, 
                                   dir=model$dir, 
                                   rhs=model$rhs, 
                                   types=model$types, 
                                   max=model$max,
                                   time_limit = ilpTimeLimit)
    print(ilpSolution$objval)
    # Adapt (line 13 of CMSA)
    age <- age+1
    selectedStakeholders <- stakeholders[ilpSolution$solution[n+(1:mprime)] > 0.5]
    age[selectedStakeholders] <- 0
  }
  
  
}