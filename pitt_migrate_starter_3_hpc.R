
source("Pitt_migrate_norm_3_1.R")


#a<-read.table("C:\\Drive E\\imbalance\\Data\\aus.csv",sep=",",header=F)
#a<-read.table("C:\\Drive E\\imbalance\\Data\\thyroidGland.csv",sep=",",header=F)
#a<-read.table("C:\\Drive E\\imbalance\\Data\\sattelite.csv",sep=",",header=F)
a<-read.table("page.csv",sep=",",header=F)
#a<-read.table("sonar.csv",sep=",",header=F)


#a<-read.table("C:\\Drive E\\imbalance\\Data\\wave.csv",sep=",",header=F)

#a<-read.arff("C:\\Drive E\\Fall16\\BBO-RMCode\\adult1.arff")
#a<-read.arff("C:/Drive E/Fall16/BBO-RMCode/cfs/optical_IG.arff")
#a<-read.table("C:\\Drive E\\imbalance\\Data\\german.csv",sep=",",header=F)
#a<-read.table("C:\\Drive E\\Fall 2017\\Data\\cmc.csv",sep=",",header=T)
#a<-read.table("C:\\Drive E\\imbalance\\Data\\sonar.csv",sep=",",header=F)
#a<-read.arff("C:\\Drive E\\Fall 2017\\Data\\thoracic.arff")
#a<-read.table("C:\\Drive E\\Fall 2017\\Data\\spect.csv",sep=",",header=T)
#a<-read.table("C:\\Drive E\\Fall 2017\\Data\\park.csv",sep=",",header=T)

rows<-nrow(a)

cols<-ncol(a)

##############
classIndex=11#23#--park #23--spect#17--thoracic#61--sonar #15#--adult##41#--wave#11 #--page#37--satt
################


### determine attType
f<-sapply(a, is.double)
f1<-sapply(a, is.factor)
f2<- sapply(f, function(x){if(x==TRUE) x=0 else x=1})
f3<-sapply(f1,function(x){if(x==TRUE) x=1 else x=0})
attType1<-rep(0, cols)

for(i in 1:cols){
  if(f3[[i]]==1)#factor
    attType1[i]<-1
  
  if(length(unique(a[,i])) <=20)#discrete var if <=20 unique values
    attType1[i]<-1
  
}


###for adult dataset only ..7/18
#attType1[1]<-0

nominal<-which(attType1==1)
##discard classIndex
nominal<-setdiff(nominal,classIndex)

y<-1
while(y<=length(nominal)){
  
  i<-nominal[y]
  if(i==classIndex)
    next
  #### here nominal values are encoded 1, 2,..R will consider them as numeric
  ### make 1, 2 as "1", "2"
  a[,i]<-as.factor(a[,i])
  
  colVal<-levels(a[,i])
  
  
  #d<-c(1:length(sets1[[3]]))
  d<-c(1:length(colVal))
  ###added new for binary feature 9/20
  
  if(length(colVal)==2){
    replaceNomVals<-as.numeric(colVal)
    # any_na<-any(is.na(replaceNomVals))
    
    is_binary<-sapply(colVal, function(x) x %in% c(0,1))
    res_logical_vector<-as.logical(is_binary)
    
    if(all(res_logical_vector)==TRUE)
      d<-replaceNomVals
    
    #binary_levels<-c(0,1)
    ###if binary levels, d is replaceNomVals..it is required for ordering of levels 0, 1
    # if(all(replaceNomVals) %in% binary_levels )
    #   d<-replaceNomVals
    
  }
  
  levels(a[,i]) <- c(d)
  
  y<-y+1
}



lowerRange<-c()
upperRange<-c()
for(y in 1:cols){
  
  i<-y
  colVal<-levels(a[,i])
  
  lowerRange[i]<-min(as.numeric((a[,i])))
  upperRange[i]<-max(as.numeric((a[,i])))
  
  
}

###apply cfs

#############################################################################################
### all migration curves seed 111
my_seed = 150

set.seed(my_seed)
a<-a[sample(nrow(a)),]


#attType1= c( 1, 0, 1, 1, 0, 1, 1,0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1);#ger

#sourceCpp("C:/Drive E/Fall16/BBO-RMCode/cppFiles/aus_num_2.cpp")
#sourceCpp("C:/Drive E/Fall16/BBO-RMCode/cppFiles/top_down.cpp")
#sourceCpp("C:/Drive E/Fall16/BBO-RMCode/cppFiles/thy.cpp")


sourceCpp("Pitt_mem_1.cpp")
#sourceCpp("C:\\Drive E\\Fall 2017\\migration_optimize\\HPC\\Pitt_mem_1.cpp")

outer_fold <-3

folds <-generateCVRuns(labels = a[, classIndex],
                       ntimes = 1,
                       nfold = outer_fold,
                       stratified=TRUE)


tuneFolds<-vector("list",3)
testfold  <-vector("list",2)

migrate_data<-NULL

sumAcc<-0



#ptm_allcv<-proc.time() - ptm

for(l in 1:1){
  
  
  ptm_cv_st<- proc.time() #ifelse(l==1, proc.time() - ptm, ptm_cv_end)
  #proc.time() - ptm
  
  
  cv.test<-folds$`Run  1`[[l]]
  cv.train<-a[-cv.test,]
  cv.test<-a[cv.test,]
  
  currDataSet<-cv.train
  
 
  
  migrate_data<-cv.train
  migrate_data <- data.frame(migrate_data)
  migrate_data <-data.matrix(migrate_data)
  
  
  
  ####change 10/1
  #set.seed(117)
  #cv.train<- cv.train[sample(nrow( cv.train)),]
  
  print(paste0("yey ", l))
  # attType1=c(1,1,1,1,1,1,1,1,1,1,1);
  
  
  cDist<-table(cv.train[,classIndex])
  c<-NULL
  for(j in 1:length(cDist))
    c<-c(c,cDist[[j]])
  
  
  df<-data.frame(currDataSet)
  df = data.matrix(df)
  
  classLabels<-allClass#c(1,2,3)
  
  defaultRule<-0
  defaultVal<-0
  for(i in 1:length(cDist)){
    if(cDist[[i]]>defaultVal){
      defaultVal<-cDist[[i]]
      defaultRule<-classLabels[i]
    }
  }
  print(paste0(classIndex))
  dfr<-setDefaultClass(defaultRule,classIndex,length(classLabels));
  
  x<-NULL
  
  
  # BBO_GA<-bbo2(fn=calcFit,genomeLen =6, rMin = 5, rMax = 20,  lower=lowerRange,#c(1,64,0.4,0.1,0.05,-.5), #lower index less than actual
  #              upper=upperRange,#c(3,145,26.3,12,58.4,58.3), 
  #              attType=attType1)
  genomeLen<-ncol(a)
  # satellite --prev 16, 25 rMin, rMax
  BBO_GA<-bbo2(fn=calcFit,genomeLen =37, rMin = 10, rMax = 20,  lower=lowerRange,#c(1,64,0.4,0.1,0.05,-.5), #lower index less than actual
               upper=upperRange, RandSeed = my_seed + 100,#c(3,145,26.3,12,58.4,58.3), 
               attType=attType1)
  
  #   BBO_GA<-bbo2(fn=calcFit,genomeLen =59, rMin = 25, rMax = 35,  lower=lowerRange,#c(1,64,0.4,0.1,0.05,-.5), #lower index less than actual
  #               upper=upperRange,#c(3,145,26.3,12,58.4,58.3), 
  #              attType=attType1)
  
  x<-BBO_GA$minCost$bestMember
  
  currDataSet<-cv.test
  ##changed 11/9
  df<-data.frame(currDataSet)
  df = data.matrix(df)
  evalRes_1 <-calcFitTst(x,attType1,df)
  testAcc<-evalRes_1[1]#calcFit(x,attType1)
  
  ptm_cv_end<-proc.time() - ptm_cv_st
  
  
  
  print(paste0("testAcc: ", l, " ", testAcc))
  sumAcc<-sumAcc+testAcc
  
  
}
#aus : 144--86.08, 177--86.08 seed 110--86.96, 120--85.94, 130--86.52

ptmElapse<-proc.time() - ptm
