
library('MASS')
library('dplyr')

a<-read.table("data1_clean_2.csv",sep=",",header=T)




##arguments to stepwise regression
classIndex  <- 20
remove_atts<-c(20, 1, 4) # classindex, stdID, MCAS
y           <-a[,classIndex]
df_subset <-as.data.frame(cbind(a[-remove_atts]))

full<-glm(y ~ ., family = binomial,data = df_subset)
full_summary<-summary(full)

##apply stepwise regression, default direction is both
step<-stepAIC(full, trace = FALSE)
##apply stepwise backward elimination
backward<-step<-stepAIC(full,direction="backward", trace = FALSE)


### final model chosen by stepwise logistic regression
final_default<-names(step$coefficients)

final_back<-names(backward$coefficients)

###print indices of final_back
final_back_ind<-c()
feature_name<-names(a)

for (i in 1:ncol(a)){
  
  if(feature_name[i]%in% final_back)
    final_back_ind<-c(final_back_ind, i-1)###python index 0 based
  
}

cat("att selected ", final_back,"\n")
cat("att selected indices ", final_back_ind,"\n")


###call from another function


do_step<-function(a){
  ##arguments to stepwise regression
  classIndex  <- 20
  remove_atts<-c(20, 1, 4) # classindex, stdID, MCAS
  y           <-a[,classIndex]
  df_subset <-as.data.frame(cbind(a[-remove_atts]))
  
  full<-glm(y ~ ., family = binomial,data = df_subset)
  full_summary<-summary(full)
  
  ##apply stepwise regression, default direction is both
  step<-stepAIC(full, trace = FALSE)
  ##apply stepwise backward elimination
  backward<-step<-stepAIC(full,direction="backward", trace = FALSE)
  
  
  ### final model chosen by stepwise logistic regression
  final_default<-names(step$coefficients)
  
  final_back<-names(backward$coefficients)
  
  ###print indices of final_back
  final_back_ind<-c()
  feature_name<-names(a)
  
  for (i in 1:ncol(a)){
    
    if(feature_name[i]%in% final_back)
      final_back_ind<-c(final_back_ind, i-1)###python index 0 based
    
  }
  
  cat("att selected ", final_back,"\n")
  cat("att selected indices ", final_back_ind,"\n")
  
}
