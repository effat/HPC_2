source("GeneticAlg_helpers.R")


my_GeneticAlg.int <- function(genomeLen, codonMin, codonMax,
                              genomeMin=rep.int(codonMin, genomeLen), genomeMax=rep.int(codonMax, genomeLen),
                              suggestions=NULL,
                              popSize=50, 
                              iterations=100, terminationCost=NA,
                              mutationChance= 1/(genomeLen+1),
                              elitism=floor(popSize/10),
                              geneCrossoverPoints=NULL,
                              monitorFunc=NULL, 
                              pCross = 0.8, #newly added on 1/17
                              evalFunc,
                              allowrepeat = TRUE,
                              showSettings=FALSE, verbose=FALSE,
                              plapply = lapply) {
  # Optimizes an Integer chromosome using a genetic algorithm.
  #
  # popSize          = the population size
  # iterations       = number of generations
  # terminationCost  = The cost (error) that if reached, the GA should termiante
  # mutationChance   = chance that a var in the string gets mutated
  # geneCrossoverPoints = An array determining the genes to be swapped in crossover
  #
  # Partially based on "R Based Genetic Algorithm (genalg package)""
  # http://cran.r-project.org/web/packages/genalg/
  
  is.verbose = verbose
  verbose = function(...) { if (is.verbose) cat(...)}
  
  if (is.null(evalFunc)) {
    stop("A evaluation function must be provided. See the evalFunc parameter.");
  }
  
  stopifnot(genomeLen > 1)
  
  # do a variaty of sanity checks first
  verbose("Testing the sanity of parameters...\n");
  if (length(genomeMin) != length(genomeMax)) {
    stop("The vectors genomeMin and genomeMax must be of equal length.");
  }
  if (popSize < 5) {
    stop("The population size must be at least 5.");
  }
  if (iterations < 1) {
    stop("The number of iterations must be at least 1.");
  }
  if (!(elitism < popSize)) {
    stop("The population size must be greater than the elitism.");
  }
  if (elitism < 0) {
    stop("elitism must be at least 0.");
  }
  if ((mutationChance < 0) | (mutationChance  > 1)) {
    stop("mutationChance must be between 0 and 1.");
  }
  if (!is.null(geneCrossoverPoints)) {
    if (!is.numeric(geneCrossoverPoints) | length(geneCrossoverPoints) != genomeLen) {
      stop("Invalid geneCrossoverPoints.");
    }
  }
  
  if (showSettings) {
    verbose("The start conditions:\n");
    result = list(genomeMin=genomeMin, genomeMax=genomeMax, suggestions=suggestions,
                  popSize=popSize, iterations=iterations,
                  elitism=elitism, mutationChance=mutationChance);
    class(result) = "rbga";
    
    cat(summary(result));
  } else {
    verbose("Not showing GA settings...\n");
  }
  
  ##########
  # Creation
  population = matrix(nrow=popSize, ncol=genomeLen);
  
  if (!is.null(suggestions)) {
    verbose("Adding suggestions to first population...\n");
    suggestionCount = nrow(suggestions)
    population[1:suggestionCount,] = population[i,]
    verbose("Filling others with random values in the given domains...\n");
  } else {
    verbose("Starting with random values in the given domains...\n");
    suggestionCount = 0
  }
  
  for (i in (suggestionCount+1):popSize) {
    population[i,] = ga.new.chromosome(genomeLen, genomeMin, genomeMax, allowrepeat)
  }
  
  ############################################################################
  # do iterations
  bestEvals = rep(NA, iterations);
  meanEvals = rep(NA, iterations);
  evalVals = rep(NA, popSize);
  
  ###added on 12/18
  stagnantCount<-0
  stagnantList<-rep(NA, ((iterations)/20)*4);
  sLen<-length(stagnantList)
  
  for (iter in 1:iterations) {
    verbose(paste("Starting iteration", iter, "\n"));
    cat("Starting iteration", iter, "\n")
    
    ##########
    # Evaluation
    
    verbose("Calucating evaluation values... ");
    
    to.eval.Ids = which(is.na(evalVals))
    evalVals[to.eval.Ids] = unlist(plapply(to.eval.Ids, 
                                           function(i, population, evalFunc) evalFunc(population[i, ]),
                                           population, evalFunc))
    
    # check for invalid items
    if ((!all(is.numeric(evalVals))) |
        any(is.na(evalVals)) |
        any(is.nan(evalVals))) {
      stop("Invalid cost function return value (NA or NaN).")
    }
    
    # extract statistics about generation
    bestEvals[iter] = min(evalVals);
    meanEvals[iter] = mean(evalVals);
    bestInd = which.min(evalVals)
    verbose(" done.\n");
    
    collect.results <- function() {
      settings = list(genomeMin=genomeMin, genomeMax=genomeMax,
                      popSize=popSize, elitism=elitism, geneCrossoverPoints = geneCrossoverPoints,
                      iterations=iterations, suggestions=suggestions,
                      mutationChance=mutationChance)
      
      pop.info = list(population=population, evaluations=evalVals, best=bestEvals, mean=meanEvals, currentIteration=iter)
      
      best = list(genome=population[bestInd,], cost = evalVals[bestInd]);
      
      ret = list(settings = settings, population = pop.info, best = best)
      
      ##modifier by Effat, 11/8
      # class(ret) = "GeneticAlg.int";
      
      return (ret)
    }
    
    if (!is.null(monitorFunc)) {
      verbose("Sending current state to rgba.monitor()...\n");
      # report on GA results
      monitorFunc(collect.results());
    }
    
    ##########
    # check termination conditions
    if (iter == iterations) {
      verbose("End of generations iteration reached.\n");
      break
    }
    
    if (!is.na(terminationCost)) {
      if (bestEvals[iter] <= terminationCost) {
        verbose("Cost better than termination cost reached.\n");
        break
      }
    }
    
    ##########
    # Selection
    
    verbose("Creating next generation...\n");
    newPopulation = matrix(nrow=popSize, ncol=genomeLen);
    newEvalVals = rep(NA, popSize);
    
    verbose("  sorting results...\n");
    sortedEvaluations = sort(evalVals, index=TRUE);
    sortedPopulation  = matrix(population[sortedEvaluations$ix,], ncol=genomeLen);
    
    # save the best
    #  if (elitism > 0) {
    #   verbose("  applying elitism...\n");
    #   newPopulation[1:elitism,] = sortedPopulation[1:elitism,];
    #   newEvalVals[1:elitism] = sortedEvaluations$x[1:elitism]
    # } # ok, save nothing
    
    ###no elitism in steady-state GA, copy all pops
    
    newPopulation[1:popSize,] = sortedPopulation[1:popSize,]
    newEvalVals[1:popSize] = sortedEvaluations$x[1:popSize]
    
    
    
    ##########
    # Crossover
    # fill the rest by doing crossover
    
    parentsID<-rep(0, 2)
    
    
    for(p_i in 1:2){
      
      Rand_num <- runif(1) * abs(sum(newEvalVals)) #### fitness value is neg
      parent_select<- abs(newEvalVals[1])
      parent_index <- 1
      
      while((Rand_num > parent_select) && (parent_index < popSize)){
        parent_index <-parent_index + 1
        parent_select <- parent_select + abs(newEvalVals[parent_index])
        
      }##end while
      
      ###insert parent index 
      parentsID[p_i]<-parent_index
      
    }###end for selecting parents
    
    
    parents = sortedPopulation[parentsID,]
    ###apply crossover on 2 parents 2 get 2 children
    ###if crossoverpoint is 0, genomlen or no crossover is selected, then 2 children are just copy of parents
    child1<-parents[1,]
    child2<-parents[2,]
    
    if(runif(1) < pCross ){
      crossOverPoint = sample(0:genomeLen,1)
      if ((crossOverPoint != 0) && (crossOverPoint != genomeLen)) {
        
        child1<- c(parents[1, 1:crossOverPoint], parents[2, (crossOverPoint+1):genomeLen])
        child2<-c(parents[1, (crossOverPoint+1):genomeLen], parents[2, 1:crossOverPoint])
        
      }###end if
      
      
    }### end if crossover
    
    
    ##########
    # Mutation
    
    mut_ch1<-my_mutation(child1, mutationChance, genomeLen, codonMin, codonMax, allowrepeat, dempeningFactor)
    
    mut_ch2<-my_mutation(child2, mutationChance, genomeLen, codonMin, codonMax, allowrepeat, dempeningFactor)
    
    ###evaluate ne ch2
    eval_ch1<-evalFunc( mut_ch1)
    eval_ch2<-evalFunc( mut_ch2)
    
    ###by default child 1 is best
    best_ch_eval<-eval_ch1
    best_ch<-eval_ch1
    
    if(abs(eval_ch2) > abs (eval_ch1)){
      best_ch_eval<-eval_ch2
      best_ch<-eval_ch2
      
    }
    
    ###replace worst pop if best_eval is better
    
    newPopulation[1:popSize,] = sortedPopulation[1:popSize,]
    newEvalVals[1:popSize] = sortedEvaluations$x[1:popSize]
    
    if(abs(best_ch_eval) > abs(newEvalVals[popSize])){
      
      newPopulation[popSize,] = best_ch
      newEvalVals[popSize] = best_ch_eval
      
    }
    
    ########################### stagnation Count
    stagnantCount<-stagnantCount+1
    
    
    if(stagnantCount>sLen )
      stagnantCount<-1
    
    stagnantList[stagnantCount]<-evalVals[1]
    
    isStagnant<-FALSE
    sCount<-2
    while((sCount<=sLen) && (iter> sLen)){
      # print(paste0("mute2 ", sCount," ", sCount-1," ", stagnantList[sCount], " ",stagnantList[sCount-1]))
      if(stagnantList[sCount]!=stagnantList[1]){
        
        break;
        
      }
      
      sCount<-sCount+1
      
    }
    
    if(sCount>sLen)
      isStagnant<-TRUE
    
    if(isStagnant==TRUE)
      break;
    
  }###end for iteration
  
  # report on GA results
  result = collect.results()
  
  return(result);
}
