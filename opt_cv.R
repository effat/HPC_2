library(cvTools)
library(TeachingSampling)

library(TunePareto)
library(Rcpp)
library(RWeka)
library(plyr)

source("my_GeneticAlg.int.R")
#source("C:\\Drive E\\Fall 2017\\migration_optimize\\Pitt_migrate_norm_3.R")###for jth rule replacement



Curveopt<-function(norm_rank, x_cord_seq, my_pop, conf_pop, instancedf){
  
  norm_rank <- norm_rank
  x_cord_seq<-x_cord_seq
  P<-length(norm_rank)
  #pop<-my_pop
  #confPop<-conf_pop
  
  for(i in 1:P){
    cat("test P", "i = ",i," ",length(my_pop[[i]]), " conf len ", length(conf_pop[[i]]), "\n")
  }
  
  
  ###fitness function
  my_target.runner<-function(migration_vect, pop = my_pop, confPop = conf_pop){
    
    
    ch_partition<-length(migration_vect)/2
    lambda_curve <-migration_vect[1:ch_partition]
    mu_curve    <-migration_vect[(ch_partition + 1): length(migration_vect)]
    
    lambda_val<-rep(0, P)
    mu_val<-rep(0, P)
    
    ###calc lambda
    for(i in 1:P){
       cost_i<-norm_rank[i]*100#changed 11/15, must mult by 100
       discrete<-floor(cost_i/10)*10
      
      low_ind<-which(x_cord_seq == discrete)
      high_ind<-ifelse(discrete == 100, low_ind, low_ind+1)
  
      y_app<-(lambda_curve[low_ind]+ lambda_curve[high_ind])/2###this value [1...100], scale it down
      lambda_val[i]<-y_app*5 ##discritized y axis
     }##for ends lambda
    
    lambda_val<-lambda_val/100
    
    ###calc mu
    for(i in 1:P){
      cost_i<-norm_rank[i]*100#changed 11/15, must mult by 100
      discrete<-floor(cost_i/10)*10
      
      low_ind<-which(x_cord_seq == discrete)
      high_ind<-ifelse(discrete == 100, low_ind, low_ind+1)
  
      y_app<-(mu_curve[low_ind]+ mu_curve[high_ind])/2###this value [1...100], scale it down
      mu_val[i]<-y_app*5
    }##for ends mu
    mu_val<-mu_val/100
    
    ###combined lambda
    
    min_lambda <- min(lambda_val)
    max_lambda <-max(lambda_val)
    maxgen_val <-20
    best_cost  <-0
    
    best_member<-NULL
    
    
    ### apply BBO for 20 gen
    
    for(i in 1:maxgen_val){
      
      
      Island     <- vector("list",P)
      IslandConf <- vector("list",P)
      
      ##start migration for all pops
      for(k in 1:P){
        
        denom<-ifelse((max_lambda - min_lambda)==0, 1, (max_lambda - min_lambda) )
        lambdaScale <- (lambda_val[k] - min_lambda) / denom
        #cat(" test lambda ", lambda_val[k], " ", min_lambda, " ", max_lambda," ", k, "\n")
        
        ## Probabilistically input new information into habitat i
        currentPop<-pop[[k]]
        currPopLen<-length(currentPop)
        
        Island[[k]]    <-NaN*seq(currPopLen)
        IslandConf[[k]]<-NaN*seq(currPopLen)
        
        #cat("test P", "k = ",k," ",length(Island[[k]]), " conf len ", length(IslandConf[[k]]), "\n")
        
        
        for( j in 1:currPopLen ){
          
          if(runif(1) < lambdaScale){
            
            ## Pick a habitat from which to obtain a feature
            RandomNum <- runif(1) * sum(mu_val)
            Select <- mu_val[1]
            SelectIndex <- 1
            while((RandomNum > Select) & (SelectIndex < P)){
              SelectIndex <- SelectIndex + 1
              Select <-Select + mu_val[SelectIndex]
            }## while-ends
            
            
            ##get target class of jth rule
            currPop_jth_rule <-currentPop[[j]]
            jth_rule_class   <-currPop_jth_rule[length(currPop_jth_rule)]
            
            ####get same class best rule from donor
            donor_rule_index  <-getTrgtClassRule(jth_rule_class, pop[[SelectIndex]])
            
            Island[[k]][j]    <- pop[[SelectIndex]][donor_rule_index]
            IslandConf[[k]][j]<- confPop[[SelectIndex]][donor_rule_index]
            
          }###end if lamdaScale
          else{
            conf<-confPop[[k]]
          #  cat("conf test j",j, " conf ", length(conf), "\n")
            IslandConf[[k]][j]<- confPop[[k]][j]
            Island[[k]][j]   <- pop[[k]][j]
            
          }## if-else-ends
          
          
        }### for for j in currPopLen
        
        ###sort rules of one island based on confidence
        tempConf <-IslandConf[[k]]
        tempPop  <-Island[[k]]
        sortedConf<- sort(tempConf, decreasing = TRUE, index=TRUE) # not -temponf
        
        Island[[k]]    <- c(tempPop[sortedConf$ix])
        IslandConf[[k]]<-c(sortedConf$x)
        
        
      }##end migration for all pops
      
      
      
      ### replace pops by islands
      pop = Island;
      confPop=IslandConf
      
      ###added 12/5
      costs<-rep(0, P)
      
      ###calc costs of pops after migration
      for(evalId in 1:P){
        evalRes      <-calcFitTst(pop[[evalId]],attType1, instancedf)
        costs[evalId]<-evalRes[1]
      }## end for cost calc
      
      
      ## Sort from best to worst
      sortedEvaluations = sort(costs, decreasing = TRUE, index=TRUE);
      sortedPopulation  = c(pop[sortedEvaluations$ix])
      newConf            <-c(confPop[sortedEvaluations$ix])
      
      costs      <-sortedEvaluations$x
      pop       <-sortedPopulation
      confPop    <- newConf
      
      
      
      
      ##get the best member
      #best_member<-pop[[1]]
      best_cost<-mean(costs)
      
    }###end for maxgen_val
    
   
    ###apply best member to validation dataset
  #  best_cost<-calcFitTst(best_member, attType1, testdf)
    
     ##return to irace
    return(list(cost= best_cost*(-1)))  
    
  }###end my_target_runner
  
  
  
  ###call GA, fitness function target_runner
  x   <-my_GeneticAlg.int(genomeLen = 22, codonMin = 0, codonMax = 20,
                          evalFunc = my_target.runner, iterations = 200, popSize = 100,  mutationChance = 0.1)
  
  best_ch <-x$best$genome
  return(best_ch)
  
}


###############


opt_partition<-function(norm_rank, x_seq, pop, conf_pop){
  
  
  ###get global df
  instanceData<<-migrate_data
  ###init tune fold
  
  nfold_opt = 3
  foldtemp <-generateCVRuns(labels =  instanceData[, classIndex],
                          ntimes = 1,
                          nfold = nfold_opt,
                         stratified=TRUE)
  
  
  sum_curves<-0
  for(l in 1:nfold_opt){
    
    
  
    instance.test<-foldtemp $`Run  1`[[l]]
    instance.train<-instanceData[-instance.test,]
    instance.test<-instanceData[instance.test,]
    
    
    #opt_ch<-Curveopt(norm_rank, x_seq, pop, conf_pop, instance.train, instance.test)
    
    ###summ curves
    #sum_curves<-sum_curves + opt_ch
    
  }###end for l
 
  #sum_curves_avg<-sum_curves/nfold_opt
  
  opt_ch<-Curveopt(norm_rank, x_seq, pop, conf_pop, instanceData)
  sum_curves_avg<-opt_ch
  
  return(sum_curves_avg)
 
  
   
}###end opt_partition


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
