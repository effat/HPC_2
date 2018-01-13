library(lattice)
library(robustbase)

source("Pitt_migrate_norm_3_1.R")


#a<-read.table("C:\\Drive E\\imbalance\\Data\\aus.csv",sep=",",header=F)

a<-read.table("page_transformed.csv",sep=",",header=F)# for recording lower and upper values
a_tr<-read.table("page_test.csv",sep=",",header=F) ###smaller file for training
a_tst<-read.table("page_tr.csv",sep=",",header=F)

#a<-read.table("C:\\Drive E\\imbalance\\Data\\thyroidGland.csv",sep=",",header=F)
#a<-read.table("C:\\Drive E\\imbalance\\Data\\sattelite.csv",sep=",",header=F)
#a<-read.table("C:\\Drive E\\imbalance\\Data\\page.csv",sep=",",header=F)
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
classIndex=11#23#--park #23--spect#17--thoracic#61--sonar #15#--adult##41#--wave#11 #--page#37--satt#--15: aus
################

attType1<-read.table("page_attTypes.csv",sep=",",header=F)
attType1 <- as.numeric(attType1)




lowerRange<-c()
upperRange<-c()
for(y in 1:cols){
  
  i<-y
  
  
  lowerRange[i]<-min(as.numeric((a[,i])))
  upperRange[i]<-max(as.numeric((a[,i])))
  
  
}

###apply cfs

#############################################################################################
### all migration curves seed 111
my_seed <-111
set.seed(155)
a<-a[sample(nrow(a)),]



#sourceCpp("C:/Drive E/Fall16/BBO-RMCode/cppFiles/Pitt_mem_1.cpp")
sourceCpp("Pitt_mem_1.cpp")

outer_fold <-3

folds <-generateCVRuns(labels = a[, classIndex],
                       ntimes = 1,
                       nfold = outer_fold,
                       stratified=TRUE)


tuneFolds<-vector("list",3)
testfold  <-vector("list",2)

migrate_data<-NULL

sumAcc<-0


allClass<-unique(a[, classIndex])


for(l in 1:1){
  
  
  ptm_cv_st<- proc.time() #ifelse(l==1, proc.time() - ptm, ptm_cv_end)
  #proc.time() - ptm
  
  
  #cv.test<-folds$`Run  1`[[l]]
  cv.train<-a_tr
  cv.test<-a_tst
  
  migrate_data<-cv.train
  migrate_data <- data.frame(migrate_data)
  migrate_data <-data.matrix(migrate_data)
  
  
  
  set.seed(my_seed)
  cv.train<- cv.train[sample(nrow( cv.train)),]
  
  currDataSet<-cv.train
  
  
  print(paste0("yey ", l))
  
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
  
  
  genomeLen<-ncol(a)
  # satellite --prev 16, 25 rMin, rMax
  BBO_GA<-bbo2(fn=calcFit,genomeLen =37, rMin = 10, rMax = 20,  lower=lowerRange,#c(1,64,0.4,0.1,0.05,-.5), #lower index less than actual
               upper=upperRange, RandSeed = 100 + my_seed,#c(3,145,26.3,12,58.4,58.3), 
               attType=attType1)
  
  
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
