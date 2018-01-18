ptm<-proc.time()

library(cvTools)
library(TeachingSampling)

library(TunePareto)
library(Rcpp)
library(RWeka)
library(plyr)

### Cross validation init

currDataSet<-NULL;
defaultRule<-NULL
df<-NULL
classIndex<-NULL
not_inc<-c(2) #20, 32, 28, 34, 17
m<-c(1,4,5,7,10,6,3,9,2,8)
#c(1,10,11,7,4,8,14,9,6,5,12,15,13)#--thoracic
#c(4,5,9,10,11,12,13,21,28,36,44,45,46,47,48,49,51,52,54)#ger
#c(22,19,1,5,13,12,9,3,8,6,2,10,20,11,14,7,15,4,16,17,21,18)#park
#(1,4,5,7,10,6,3,9,2,8)#page



#c(13,21,8,22,16,17,18,20,7,3,2,12,11,10,6,15,4,1,14,5,9,19)#spect
#c(1,10,11,7,4,8,14,9,6,5,12,15,13)#--thoracic
#c(4,5,9,10,11,12,13,21,28,36,44,45,46,47,48,49,51,52,54)#sonar
#c(3:34,1)#ion#c(1,3,2,6,4,5,12,7,15,13,14,9,20,10,17,19)--ger


source("opt_cv_steady.R")


##updated in each cv
tuned_val_a<-0
tuned_val_b<-0


bbo.control <- function(pModify = 1, pMutate = 0.3, KEEP = 2, popSize = 50, maxGen = 100, numVar = 2, orderDep = TRUE, creepRate=0.5, creepFrac=0.01) {
  
  
  if (pModify > 1 | pModify < 0) {
    warning("'pModify' not in [0,1]; set to default value 1\n", immediate. = TRUE)
    pModify <- 1
  }
  
  if (pMutate < 0 | pMutate > 1) {
    warning("'pMutate' not in [0,1]; set to default value 0.3\n", immediate. = TRUE)
    pMutate <- 0.3
  }
  
  if (popSize <= 0 | popSize < KEEP) {
    warning("'popSize' must be > 0 and/or greater than 'KEEP'; 'popSize' set to default value 20\n", immediate. = TRUE)
    #popSize <- 20
    stop("Halted!")
  }
  
  if (KEEP < 0 | KEEP > popSize) {
    warning("'KEEP' cannot be negative or greater than 'popSize'; 'KEEP' default value is 5\n", immediate. = TRUE)
    stop("incorrect value passed!")
  }
  
  if (maxGen <= 0) {
    warning("'maxGen' must be > 0, set to default value 20\n", immediate. = TRUE)
    maxGen <- 20
  }
  
  
  if (numVar <= 0) {
    warning("'numVar' cannot be negative or O; numVar default value is 2\n", immediate. = TRUE)
    stop("incorrect value passed!")
  }
  
  if (missing(orderDep)){
    #warning("'orderDep' is set to TRUE\n", immediate. = TRUE)
    orderDep <- TRUE
  }
  
  
  list(pModify = pModify, pMutate = pMutate, KEEP = KEEP, popSize = popSize, maxGen = maxGen, numVar = numVar, orderDep = orderDep,creepRate=creepRate, creepFrac=creepFrac)
  
}






bbo2 <- function(fn, genomeLen,rMin,rMax, lower, upper, attType, sets,cLabel,DisplayFlag = TRUE, RandSeed, control = bbo.control(), ... ) {
  env <- new.env()
  ctrl <- do.call(bbo.control, as.list(control))
  selfSet <- FALSE
  if(missing(RandSeed)){
    RandSeed <- 1024
    selfSet <- TRUE
    
    cat("no randseed \n")
  }
  
  
  
  if (length(lower) != length(upper))
    stop("'lower' and 'upper' are not of same length")
  
  if (!is.vector(lower))
    lower <- as.vector(lower)
  
  if (!is.vector(upper))
    upper <- as.vector(upper)
  
  if (any(lower > upper))
    stop("'lower' > 'upper'")
  
  if (any(lower == "Inf"))
    warning("A component of 'lower' assigned 'Inf'. May imply 'NaN' results", immediate. = TRUE)
  
  if (any(lower == "-Inf"))
    warning("A component of 'lower' assigned '-Inf'. May imply 'NaN' results", immediate. = TRUE)
  
  if (any(upper == "Inf"))
    warning("A component of 'upper' assigned 'Inf'. May imply 'NaN' results", immediate. = TRUE)
  
  if (any(upper == "-Inf"))
    warning("A component of 'upper' assigned '-Inf'. May imply 'NaN' results", immediate. = TRUE)
  
  
  
  
  if(!missing(RandSeed)){
    if(!selfSet){
      set.seed(RandSeed)
    }
  }
  
  ## FeasibleFunction = #todo (Checks whether each candidate solution is a feasible solution; if not it bounds the individual parameter values by the minimum and maximum)
  
  ## Supporting functions::
  #TO DO# ClearDups : makes sure that the population does not have duplicates
  
  
  ## InitPopulation function code
  ## Next version:  we may allow the user to supply his own set of initial population, if that seems feasible considering all possible constraints
  
  # Creation
  #population = matrix(nrow=popSize, ncol=genomeLen);
  
  population<-vector("list",ctrl$popSize)
  confPop<-vector("list",ctrl$popSize)
  
  
  varRule<-function(seedData,val){
    #1:7 --
    
    max_att<-5
    attrN<-sample(1:max_att,1 , prob = rep(1/max_att, max_att )) # att number in chromosome
    ###do sampling attributes
    m1<-sample(m, size = length(m))#m[sample(length(m))]
    
    attrSet<-sample(m1, attrN, prob = rep(1/length(m1), length(m1))) # listing attrN atributes by index
    seedDataVal<-as.numeric(currDataSet[seedData,])
    
    
    
    varChrom=c()
    j<-1
    
    
    totalAtt<-length(attrSet)
    for(i in 1:totalAtt){
      
      attrIndex<-attrSet[i]
      
      if(attType[attrIndex]==0){
        
        spade<-(runif(1)*(.70 - .25) + .25)*(upper[attrIndex] - lower[attrIndex])
        geneL<-max(seedDataVal[attrIndex]-spade*.5, lower[attrIndex] )
        geneU<-min(seedDataVal[attrIndex]+spade*.5, upper[attrIndex] )
        gene<-runif(2)*(upper[attrIndex] - lower[attrIndex]) + lower[attrIndex]
        
        
        
        if(val==0)
        {
          #gene<-c()
          #mid<-.5*(upper[attrIndex] - lower[attrIndex])
          #gene[1]<-runif(1)*(mid - lower[attrIndex]) + lower[attrIndex]
          #gene[2]<-runif(1)*(upper[attrIndex] - mid) + mid
          
          ####changed 10/5
          temp   <-gene[1]
          gene[1]<- gene[2]
          gene[2]<-temp
          
          geneL=gene[1]
          geneU=gene[2]
          
          
        }
        
        varChrom[j]=attrIndex
        varChrom[j+1]=geneL
        varChrom[j+2]=geneU
        
        
      }
      
      else{
        
        ###sign: 0 means not =, 1 means eq
        sign<-ifelse(runif(1)<.05, 0,1)#sample(c(0, 1),1)
        attVal<-NULL
        
        if(val==1){
          ##generate value randomly excluding seed
          nomEnd<-setdiff(lowerRange[attrIndex]:upperRange[attrIndex], seedDataVal[attrIndex])
          ##shuffle
          nomEnd<-sample(nomEnd, size = length(nomEnd))   #nomEnd[sample(length(nomEnd))]
          
          
          if(length(nomEnd)>1)
            val_w_o_seed<-sample(nomEnd,1, prob = rep(1/length(nomEnd), length(nomEnd)))
          else
            val_w_o_seed<-nomEnd[1]
          
          if(sign==1)
            attVal<-seedDataVal[attrIndex]
          else
            attVal<-val_w_o_seed
          
        }### end if val==1 nom att
        else ## val!=1
          
        {
          all_nom <-c(lowerRange[attrIndex]:upperRange[attrIndex])
          all_nom <-sample(all_nom, size = length(all_nom))   #all_nom[sample(length(all_nom))] 
          
          if(length(all_nom)>1)
            attVal<-sample(all_nom,1, prob = rep(1/length(all_nom), length(all_nom)))
          else
            attVal<-all_nom[1]
          
        }
        
        varChrom[j]=attrIndex
        varChrom[j+1]=sign
        varChrom[j+2]=attVal
        
        
        
      }
      
      
      j<-j+3
      
      
    }
    
    varChrom[j]<-seedDataVal[classIndex]
    
    ret<-varChrom
    
  }
  
  
  totalPop<-ctrl$popSize
  for(i in 1: totalPop){
    
    possible_rules<-c(rMin:rMax)
    possible_rules<- possible_rules[sample(length(possible_rules))]
    
    ruleN<-sample(possible_rules, 1, prob = rep(1/length(possible_rules), length(possible_rules) ))
    
    # print(paste0("in totalpop ", i, "ruleN ", ruleN))
    trData<-rep(0,nrow(currDataSet))
    
    seeds<-NULL
    confCount<-rep(0, ruleN)
    classes<-allClass #c(1,2,3,4,5,7)
    len_class<-length(classes)
    
    seeded_rule<-floor(.8*ruleN)
    prop_class<-floor(seeded_rule/len_class)
    
    prop_cl_index<-rep(0, len_class)
    
    for(cl_in in 1:(len_class))
      prop_cl_index[cl_in]<-cl_in*prop_class
    
    
    ##rest of class upto ruleN is assigned last class
    # prop_cl_index[len_class]<-prop_cl_index[len_class]+(seeded_rule-prop_cl_index[len_class-1])---wrong
    prop_cl_index[len_class]<-prop_cl_index[len_class]+(ruleN-prop_cl_index[len_class])
    
    # cat("props ", prop_cl_index, "\n")
    
    ruleN_class<-rep(0,(len_class+1))
    
    # cat("prop =",prop_cl_index, "\n")
    j<-1
    while(j<=ruleN){
      
      ##get class index
      class_indices<- which(j<=prop_cl_index)
      class_indices<-class_indices[1] ## 
      
      currClassLabel<-classes[class_indices]  #sample(1:3, 1)
      currClass<-which(currDataSet[,classIndex]==currClassLabel)  
      
      unmark<-NULL 
      
      # for(k in 1:length(currClass))# -- why loop
      unmark<-setdiff(currClass,seeds)
      
      
      
      #seed tuple 
      if(!is.null(unmark)) {
        
        ###sample data indices
        unmark<-unmark[sample(length(unmark))]
        
        
        
        if(length(unmark)>1)
          seedData<-sample(unmark,1, prob = rep(1/length(unmark), length(unmark))  )
        else
          seedData<- unmark[1]#sample(unmark,1  )
      }
      
      else{
        currClass<-currClass[sample(length(currClass))]
        
        if(length(currClass)>1)
          seedData<-sample(currClass,1, prob = rep(1/length(currClass), length(currClass)))
        else
          seedData<-currClass[1]#sample(currClass,1)
      }
      
      
      # seedData<-sample(unmark,1)  
      
      
      val<-0
      if(j<=.8*ruleN)
        val<-1
      
      temp<-varRule(seedData,val)#varRule1(as.numeric(currDataSet[seedData,]), val, upper, lower,
      # attType, 1:7, 1:(genomeLen-1))#varRule(seedData,val)
      
      
      classCov<-c(0,0,0,0)
      covCount<-0
      confVal<-0
      index<- -1
      
      if(val==0)
      {
        evalRes <-confCountRand(temp,attType,df)#calcFitTst(population[[evalId]],attType,df)
        confVal<-evalRes[1]
        index<-evalRes[2]
        # cat(" ", index, " ")
      }    
      
      
      if(val==1){
        index=currClassLabel
        evalRes <-confCountRand(temp,attType,df)#calcFitTst(population[[evalId]],attType,df)
        confVal<-evalRes[1]
        # index<-evalRes[2] # done 8/16 to ensure equal prop rule in releSet
        #cat("=", index, " ")
      }
      
      # print( paste(c("returning varRule: ", val, "val ", temp,"index ",index,"len ", length(temp),"\n"), collapse = " "))
      
      
      if(index!=-1){
        
        seeds<-c(seeds,seedData)
        
        
        temp[length(temp)]<-index
        # if(covCount>0)
        confCount[j]<-confVal
        
        
        
        
        if(j==1)
          population[[i]]<-list(temp)
        
        else{
          list1<-population[[i]]
          list1[[j]]<-temp
          population[[i]]<-list1
        }
        
        
        ##update ruleNClass
        #  rule_j_cls<-temp[length(temp)]
        #  ruleN_class[rule_j_cls]<-ruleN_class[rule_j_cls]+1
        
        j<-j+1
        
      }
      
    }## end while for ruleN 
    
    # cat(" finished 1 ch \n")
    
    temp<-population[[i]]
    sortedEvaluations<- sort(confCount, decreasing=TRUE, index=TRUE)
    population[[i]]  = c(temp[sortedEvaluations$ix])
    confPop[[i]]<-c(sortedEvaluations$x)
    
    # cat(" ruleN= ", i, "  poplen ", length(population[[i]]), "conflen ",length(confPop[[i]]), "\n")
    
    
  } 
  
  ## Pass the additional arguments to the objective function
  costs = rep(NA, ctrl$popSize);
  
  #  to.eval.Ids = which(is.na(costs))
  # costs[to.eval.Ids] = unlist(lapply(to.eval.Ids, 
  #                                   function(i, population, fn) fn(population[[i]],attType),
  #                                  population, fn))
  
  
  
  
  for(evalId in 1:ctrl$popSize){
    
    #if(is.na(costs[evalId]))
    evalRes <-calcFitTst(population[[evalId]],attType,df)#calcFitTst(population[[evalId]],attType,df)
    costs[evalId]<-evalRes[1]
    # print(paste0(costs[evalId]))
  }
  ## 'members' is a matrix, 'costs' is a vector
  
  #  population <- list(members = members, costs = costs)
  
  
  #sort in descending order
  
  # sortedEvaluations<-sort(costs, index=TRUE)
  sortedEvaluations<-sort(costs, decreasing = TRUE, index=TRUE)
  population  = c(population[sortedEvaluations$ix])
  costs<-c(sortedEvaluations$x)
  confPop<-c(confPop[sortedEvaluations$ix])
  
  
  cat("entering in meta-GA \n")
  ###optimize migration here
  
 
  
  P<-ctrl$popSize
  
  get_data<<-migrate_data
  
  
  ###normalized ranking..will be used in evolving GA and original bbo
  norm_rank <- rep(Inf, P)
  x_cord_seq<-seq(from =0, to = 100, by = 10)
  
  
  for( i in 1:P ){
    
    norm_rank[i] <- (i - 1) / (P-1) # linear
  }
  
 
  ###corrected 12/5...population should not be altered, opt function should use a copy of it
  
  my_pop<-population
  
  
  
  ###tuned params by irace
  # tuned.confs <- irace(scenario = scenario, parameters = parameters)
  #best_ch <-Curveopt(norm_rank, x_cord_seq, my_pop, confPop, get_data)
  best_ch <-opt_partition(norm_rank, x_cord_seq, my_pop, confPop)#
  
  
  
  #best_ch_part<- length(best_ch)/2
  #tuned_a <-best_ch[1:best_ch_part]
  #tuned_b <-best_ch[(best_ch_part+1):length(best_ch)]
  
  ##change made on 1/18
  tuned_a <-best_ch
  tuned_b <-1 - best_ch
 
  
  cat(tuned_a," tuned a\n")
  cat(tuned_b, " tuned_b\n")
  #tuned_w1<-tuned.confs$w1[1]
  #tuned_mut<-tuned.confs$m[1]
  
  tuned_val_a<<-tuned_a ## global variable
  tuned_val_b<<-tuned_b
  
  cat("tuned res finished", "\n")
  
  eliteSolutions <- vector("list",ctrl$KEEP)
  eliteCosts <- rep(Inf, ctrl$KEEP)
  eliteConf<-vector("list",ctrl$KEEP) 
  stagnantList<-rep(NA, ((ctrl$maxGen)/20)*4);
  
  MinCostEachGen <- c()
  AvgCostEachGen <- c()
  BestMemberEachGen <- c()
  stagnantCount<-0
  
  
  
  ## Begin the optimization loop
  for( genIndex in 1:ctrl$maxGen ){
    
    if(costs[1]==-1)
      break;
    ## Save the best habitats in a temporary array/vector
    
    eliteSolutions[1:ctrl$KEEP] <- population[1:ctrl$KEEP]
    eliteCosts[1:ctrl$KEEP] <- costs[1:ctrl$KEEP]
    eliteConf[1:ctrl$KEEP]<-confPop[1:ctrl$KEEP]
    #elite <- list(eliteSolutions, eliteCosts)
    

    lambda <- rep(Inf, ctrl$popSize)
    mu <- rep(Inf, ctrl$popSize)
    
    ## mu(i) is the extinction rate for individual i.
    ## This routine assumes the population is sorted from most fit to least fit.
    ## getLambdaMu() is a separate function in Matlab impelementation, while here we do it in-place
    
    P <- ctrl$popSize
    
  
    lambda_curve <- tuned_a
    mu_curve    <- tuned_b
    ###calc lambda
    for(i in 1:P){
      
      cost_i<-norm_rank[i]*100#changed 11/15, must mult by 100
      discrete<-floor(cost_i/10)*10
      
      low_ind<-which(x_cord_seq == discrete)
      
      high_ind<-ifelse(discrete == 100, low_ind, low_ind+1)
      
      
      y_app<-( lambda_curve[high_ind])
      
      lambda[i]<-y_app*5
      
    }##for ends lambda
    lambda<-lambda/100
    
    
    ###calc mu
    for(i in 1:P){
      
      cost_i<-norm_rank[i]*100#changed 11/15, must mult by 100
      discrete<-floor(cost_i/10)*10
      
      low_ind<-which(x_cord_seq == discrete)
      
      high_ind<-ifelse(discrete == 100, low_ind, low_ind+1)
      
      
      y_app<-( mu_curve[high_ind])
      
      mu[i]<-y_app*5
    }##for ends mu
    
    mu<-mu/100
    
    
     
    ## getLambdaMu done
    ##plot evolved lambda mu
    
    if(genIndex == 1){
       pdf_name<-paste(RandSeed, ".pdf", sep = "")
       pdf(pdf_name)

      par( col="blue")
      plot(norm_rank, lambda, type="o", xlab="rank of pop", ylab = "lambda", main = "evolved norm ranking migration curve: lambda")
      plot(norm_rank, mu, type="o", xlab="rank of pop", ylab = "mu", main = "evolved norm ranking migration curve: mu")
      
      tuned_val_a<<- lambda
      tuned_val_b<<- mu
      
    }
    
    
    
    lambdaMin <- min(lambda)
    lambdaMax <- max(lambda)
    
    Island <- vector("list",ctrl$popSize)
    IslandConf <- vector("list",ctrl$popSize)
    
    for( k in 1:ctrl$popSize ){
      if(runif(1) > ctrl$pModify){
        ## do-nothing
        next
      }
      
      ## Normalize the immigration rate.
      denom<-ifelse((lambdaMax - lambdaMin)==0, 1, (lambdaMax - lambdaMin) )
      lambdaScale <- (lambda[k] - lambdaMin) / denom
      
      #lambdaScale <- (lambda[k] - lambdaMin) / (lambdaMax - lambdaMin)
      
      ## Probabilistically input new information into habitat i
      currentPop<-population[[k]]
      currPopLen<-length(currentPop)
      
      Island[[k]]<-NaN*seq(currPopLen)
      IslandConf[[k]]<-NaN*seq(currPopLen)
      
      
      for( j in 1:currPopLen ){
        
        
        if(runif(1) < lambdaScale){
          ## Pick a habitat from which to obtain a feature
          RandomNum <- runif(1) * sum(mu)
          Select <- mu[1]
          SelectIndex <- 1
          while((RandomNum > Select) & (SelectIndex < ctrl$popSize)){
            SelectIndex <- SelectIndex + 1
            Select <-Select + mu[SelectIndex]
          }## while-ends
          
          
          
          
          
          ##get target class of jth rule
          currPop_jth_rule<-currentPop[[j]]
          jth_rule_class<-currPop_jth_rule[length(currPop_jth_rule)]
          
          ####get same class best rule from donor
          donor_rule_index<-getTrgtClassRule(jth_rule_class, population[[SelectIndex]])
          
          
          Island[[k]][j] <- population[[SelectIndex]][donor_rule_index]
          IslandConf[[k]][j]<- confPop[[SelectIndex]][donor_rule_index]
          
          
          
          
          ##test
          migrated<-Island[[k]][[j]]
          if(is.nan(migrated))
            cat("problem in bbo same class", k, " ", j, "\n")
          
          #}###end runif mu[selectIndex]
          
          
        }## end if (runif(1)<lambda)
        else{
          Island[[k]][j] <- population[[k]][j]
          IslandConf[[k]][j]<- confPop[[k]][j]
          
          
          
        }## if-else-ends
        
        
      }## second-for-ends
      
      tempConf<-IslandConf[[k]]
      tempPop<-Island[[k]]
      sortedConf<- sort(tempConf, decreasing = TRUE, index=TRUE) # not -temponf
      
      Island[[k]]  = c(tempPop[sortedConf$ix])
      IslandConf[[k]]=c(sortedConf$x)
      
    }## first-for-ends
    
    
    
    creepOperation<-c(-1,1)
    creepBound<-c("U","L")
    for( k in 1:ctrl$popSize) {
      
      len<-length(Island[[k]])
      
      
      if(ctrl$pMutate> runif(1))
      {
        
        all_rules<-c(1:length(Island[[k]]))
        all_rules<-all_rules[sample(length(all_rules))]
        
        varRule<-sample(all_rules, 1)
        temp1  <-Island[[k]][[varRule]]
        len1   <- (length(temp1)-1)/3
        
        mutAtt     <-sample(1:len1,1) # att of rule to mutate
        posAtt     <-(mutAtt-1)*3+1 #pos of the mutation attribute
        attrIndex1 <-as.numeric(temp1[posAtt])
        lowBound   <-as.numeric(temp1[posAtt+1])
        upBound    <-as.numeric(temp1[posAtt+2])
        low        <-as.numeric(lower[attrIndex1])
        up         <-as.numeric(upper[attrIndex1])
        
        
        if(attType[attrIndex1]==0){
          
          creepFrac<-(runif(1)*(0.30 - 0.0) + 0.0)
          multiplier<-sample(creepOperation, 1)
          muttTerm<-(up-low)*creepFrac*multiplier#(up-low)*creepFrac*multiplier
          
          if(sample(creepBound,1)=="L")##low bound mutate
            temp1[posAtt+1]<-max(low, lowBound+muttTerm)
          else
            temp1[posAtt+2]<-min(up, upBound+muttTerm)
          
        }
        
        else{
          
          nomEnd<-c(upper[attrIndex1]:lower[attrIndex1],1) ## all values of nom att
          replaceAble<-setdiff(nomEnd, upBound) ## assign a new value other than previous
          replaceVal<-replaceAble[1]
          temp1[posAtt+2]=replaceVal # replace value
          
        }#### end of nom att mut
        
        
        Island[[k]][[varRule]] <- temp1
        evalRes <-confCount(temp1,attType,df)#calcFitTst(population[[evalId]],attType,df)
        confVal<-evalRes[1]
        
        
        IslandConf[[k]][[varRule]] <-confVal
        
        tempConf<-IslandConf[[k]]
        tempPop<-Island[[k]]
        
        
        sortedConf<- sort(tempConf, decreasing = TRUE, index=TRUE) # not -temponf
        Island[[k]]  = c(tempPop[sortedConf$ix])
        IslandConf[[k]]=c(sortedConf$x)
        
        
      }###end if of mutating 1 island
      
    }###end for mut of all pop
    
    
    
    
    ## Mutation ##
    
    
    
    ## Replace the habitats with their new versions.
    population = Island;
    confPop=IslandConf
    
    
    ## Make sure each individual is legal.
    #  population <- feasibilityCheck(population, lower, upper)
    
    ## Calculate cost
    newEvalVals = rep(NA, ctrl$popSize);
    #to.eval.Ids = which(is.na(newEvalVals))
    #costs[to.eval.Ids] = unlist(lapply(to.eval.Ids, 
    #                                  function(i, population, fn) fn(population[[i]],attType),
    #                                 population, fn))
    
    for(evalId in 1:ctrl$popSize){
      
      #if(is.na(costs[evalId]))
      evalRes <-calcFitTst(population[[evalId]],attType,df)#calcFitTst(population[[evalId]],attType,df)
      costs[evalId]<-evalRes[1]
      #print(paste0(costs[evalId]))
    }
    #population$costs <- apply(population$members, 1, fn)
    
    ## Sort from best to worst
    # sortedEvaluations = sort(costs, index=TRUE);
    sortedEvaluations = sort(costs, decreasing = TRUE, index=TRUE);
    sortedPopulation  = c(population[sortedEvaluations$ix])
    newConf<-c(confPop[sortedEvaluations$ix])
    
    costs=sortedEvaluations$x
    population=sortedPopulation
    confPop= newConf
    ## Replace the worst with the previous generation's elites.
    n <- ctrl$popSize
    
    if(ctrl$KEEP > 0){
      for( k in 1:ctrl$KEEP ){
        population[n-k+1] <- eliteSolutions[k]
        costs[n-k+1] <- eliteCosts[k]
        confPop[n-k+1]<-eliteConf[k]
      }
    }
    
    
    ## Make sure the population does not have duplicates. 
    ## Population <- clearDups(population, upper, lower, fn)
    
    ## Sort from best to worst
    #  sortedEvaluations = sort(costs, index=TRUE);
    sortedEvaluations = sort(costs, decreasing = TRUE, index=TRUE);
    sortedPopulation  = c(population[sortedEvaluations$ix])
    newConf<-c(confPop[sortedEvaluations$ix])
    
    costs=sortedEvaluations$x
    population=sortedPopulation
    confPop= newConf
    
    ## Compute the average cost
    ## nLegal: return value from computeAveCost is not used anywhere, so we dont collect it, for now.
    AverageCost <- mean(costs)
    ########################### stagnation Count
    stagnantCount<-stagnantCount+1
    sLen<-length(stagnantList)
    
    if(stagnantCount>sLen )
      stagnantCount<-1
    
    stagnantList[stagnantCount]<-costs[1]
    
    isStagnant<-FALSE
    sCount<-2
    while((sCount<=sLen) && (genIndex> sLen)){
      # print(paste0("mute2 ", sCount," ", sCount-1," ", stagnantList[sCount], " ",stagnantList[sCount-1]))
      if(stagnantList[sCount]!=stagnantList[1]){
        
        break;
        
      }
      
      sCount<-sCount+1
      
    }
    
    if(sCount>sLen)
      isStagnant<-TRUE
    
    
    ## Display info to screen
    MinCostEachGen <- c(MinCostEachGen, costs[1])
    AvgCostEachGen <- c(AvgCostEachGen, AverageCost)
    currBestPop <- population[[1]] 
    BestMemberEachGen <- rbind(BestMemberEachGen, currBestPop)
    rownames(BestMemberEachGen) <- NULL
    if( DisplayFlag ){
      cat('The best and mean of Generation # ', genIndex, ' are ', tail(MinCostEachGen, 1), ' and ', tail(AvgCostEachGen, 1), '\n')
    }
    
    if(isStagnant==TRUE)
      break;
    
    
    
  }## main-optimization-for-loop-upto-maxgen-ends
  
  bestMember <- population[[1]] 
  bestValue <-  tail(MinCostEachGen, 1)
  ## conclude function can be called here
  ## FINAL RETURN
  minCost <- list(bestMember = bestMember, bestValue = bestValue)
  bestHabitat <- list(itersBestValue = MinCostEachGen, itersBestMember = BestMemberEachGen, itersAvg = AvgCostEachGen)
  
  return(list(prop = ctrl, minCost = minCost, bestHabitat = bestHabitat))
  
}


getTrgtClassRule<-function(trgtClass, ruleSet)
{
  
  parent1_rules<-length(ruleSet)
  
  
  ###get all rule indices from 2 parents of trgtClass
  p1_trgtClass<-c()
  
  
  for(i in 1:parent1_rules){
    currRule<-ruleSet[[i]]
    endIndex<-length(currRule)
    
    if(currRule[endIndex]==trgtClass)
      p1_trgtClass<-c(p1_trgtClass, i)
  }##end for parent1_rules 
  
  #cat("trgtclas ", trgtClass, " ", p1_trgtClass, "\n")
  index<-p1_trgtClass[1]
  donor_rule<-ruleSet[[index]]
  
  # if(donor_rule[length(donor_rule)]==trgtClass)
  # cat("done selection successfully \n")
  # else
  #  cat("problem \n")
  
  trgtIn<-p1_trgtClass[1]
  donor_rule<-ruleSet[[trgtIn]]
  if(is.nan(donor_rule))
    cat("problem here \n")
  
  ##return 1st index of same class rule
  return (p1_trgtClass[1])
}


calcFit<-function(x,attType){
  
  count<-0
  totalRow<-nrow(currDataSet)
  
  for(i in 1:totalRow){
    
    
    
    stIndex<-1
    endIndex<-1
    classIndexTst<-0
    totalRules<-length(x)
    currRow<-currDataSet[i,]
    
    while(stIndex<=totalRules){
      
      retVal<-coverCpp(x[[stIndex]],as.numeric(currDataSet[i,]),attType)
      if(retVal[1]==1) {
        
        
        temp1<-x[[stIndex]]
        endIndex<-length(temp1)
        endIndex<-as.numeric(temp1[endIndex])
        
        classIndexTst<-endIndex
        break;
        
      }
      
      stIndex<-stIndex+1
    }
    
    
    
    #no rules matched
    if(stIndex>totalRules){
      classIndexTst<-defaultRule
      
    }
    
    
    if(classIndexTst==as.numeric(currRow[classIndex]))
      count<-count+1
    
  }
  
  
  
  fitValue<-(-1)*(count/totalRow)
  
}



