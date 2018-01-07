//#include <Rcpp.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//



int defaultClass=4;
int classIndex1=0;
int totalClass=0; 

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  defaultClass=6;
  Rcout<<defaultClass<<std::endl;
  return x * 2;
}


// [[Rcpp::export]]
NumericVector setDefaultClass(NumericVector x, NumericVector y, NumericVector z) {
  defaultClass=x[0];
  classIndex1=y[0]-1;
  totalClass=z[0]+1;
  Rcout<<defaultClass<<std::endl;
  return x * 1;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
#timesTwo(42)
*/





// [[Rcpp::export]]
NumericVector coverCpp(NumericVector x, NumericVector data, NumericVector attType) {
  
  int  i=0, rLen=x.size()-1, attId;
  Rcpp::NumericVector res(1);
  res[0]=1;
  
  while(i<rLen){
    
    attId=x[i];
    attId=attId-1;
    if(attType[attId]==0)
    {
      if((data[attId]< x[i+1]) || (x[i+2]<data[attId])){
        res[0]=0;  
        break;
      }
    }
    
    else
    {
      if(data[attId]!=(x[i+2])){
        res[0]=0;
        break;
      }
    }
    
    i=i+3;
    
    
  }
  
  return res;
}



NumericVector coverCpp1(NumericVector x, NumericMatrix data, int index, NumericVector attType) {
  
  int  i=0, rLen=x.size()-1, attId;
  Rcpp::NumericVector res(1);
  res[0]=1;
  //index=index-1;// C++ index 0 based
  
  while(i<rLen){
    
    attId=x[i];
    attId=attId-1;//C++ zero based indexing
    
    if(attType[attId]==0)
    {
      if((data(index, attId)< x[i+1]) || (x[i+2]<data(index, attId))){
        res[0]=0;  
        break;
      }
    }
    
    else
    {
      
      int  dataVal=data(index, attId), flag=1;
      
      
      //  Rcout<<"nom "<<x[i+dataVal]<<"data "<<dataVal<<std::endl;
      
      if(x[i+1]==1)//equal sign, break if val mismatch, continue if matches
      {
        if(x[i+2]!=dataVal){
          res[0]=0;
          flag=0;
          // break;
        }
      }
      else// not equal sign, break if matches, continue if mismatches
      {
        if(x[i+2]==dataVal){
          res[0]=0;
          flag=0;
          //break;
        }
      }  
      
      
      if(flag==0)
        break;
      
      
      
    }// else ends for checking nom att
    
    i=i+3;
    
    
  }
  
  return res;
}


// [[Rcpp::export]]
NumericVector calc_num_att(NumericVector x,  NumericVector attType) {
  
  int  i=0, rLen=x.size()-1, attId, attCnt;
  NumericVector res(rLen);

  //index=index-1;// C++ index 0 based
  attCnt=0;
  
  while(i<rLen){
    
    attId=x[i];
    attId=attId-1;//C++ zero based indexing
    
    if(attType[attId]==0)
    {
      // insert index in rule of the numerical att
      res[attCnt]=i+1;//to match with R index count
      attCnt+=1;
    }
    
   
    i=i+3;
    
    
  }//end while
  
  return res;
}

// [[Rcpp::export]]
NumericVector calcFitTst(List x, NumericVector attType, NumericMatrix currData){
  
  int  totalRow=currData.nrow(), stIndex, endIndex, classIndex, n, count=0 ;
  n=x.size();
  
  NumericVector covRes(1), fitVal(1);
  
  for(int i=0;i<totalRow;i++){
    stIndex=0;
    endIndex=1;
    classIndex=-1;
    covRes[0]=0;
    
    while(stIndex<n){
      
      SEXP ll = x[stIndex];
      Rcpp::NumericVector y(ll);
      // int rLen=y.size()-1, attId, j=0;
      
      //Rcout<<xa(i,0)<<std::endl;
      
      covRes=coverCpp1(y, currData, i, attType);
      
      if(covRes[0]==1) {
        
        endIndex=y[y.size()-1];
        classIndex=endIndex;
        //  Rcout<<"yey ="<<y<<" classin= "<<classIndex<<std::endl;
        //  Rcout<<"curr"<<currData(i,0)<<std::endl;
        break;
        
      }
      
      stIndex=stIndex+1;
    }
    
    if(stIndex>=n){
      classIndex=defaultClass;
      
    }
    
    
    
    if(classIndex==currData(i,classIndex1)){
      count=count+1;
      //  Rcout<<"curr "<<currData(i,0)<<" class "<<classIndex<<" count "<<count<<std::endl;
    }
    
  }
  
  fitVal[0]=(1)*((double)count/(double)totalRow);
  // Rcout<<"fit  "<<fitVal[0]<<std::endl;
  return fitVal;
  
  
}

// [[Rcpp::export]]
NumericVector rule_generalize(NumericVector y, NumericVector attType, NumericMatrix currData){
  
  int  totalRow=currData.nrow(),endIndex, classIndex,count=0, covered=0 ;
  
  int attributes=(y.size()-1)/3, pos_missed_index, neg_missed_index;
  
  NumericVector covRes(1), fitVal(1), output(2*attributes+4);
  
  endIndex=y[y.size()-1];
  classIndex=endIndex;
  
  
  for(int i=0;i<totalRow;i++){
    
    covRes[0]=0;
    //   while(stIndex<n)
    
    
    // int rLen=y.size()-1, attId, j=0;
    
    //Rcout<<xa(i,0)<<std::endl;
    
    // covRes=coverCpp1(y, currData, i, attType);
    int  i1=0, rLen=y.size()-1, attId, mismatch, attChecked, mismatch_index;
    Rcpp::NumericVector res(1);
    res[0]=1;
    //index=index-1;// C++ index 0 based
    
    mismatch=0; //break if mismatch >0
    attChecked=0;
    mismatch_index=0;
    
    while(i1<rLen){
      
      attId=y[i1];
      attId=attId-1;//C++ zero based indexing
      attChecked=attChecked+1;//1 more att being checked
      
      if(attType[attId]==0)
      {
        if((currData(i, attId)< y[i1+1]) || (y[i1+2]<currData(i, attId))){
          mismatch=mismatch+1;
          
          if(mismatch==1){
            res[0]=0;
            mismatch_index=attChecked;
          }
          else
            break;
          //end else 
          
        }//end if for numeric condition
      }
      
      else
      {
        if(currData(i, attId)!=(y[i1+2])){
          
          mismatch=mismatch+1;
          if(mismatch==1){
            res[0]=0;
            mismatch_index=attChecked;
          }
          else
            break;
          //end else 
          
        }//end nominal cond check
      }
      
      i=i+3;
      
      
    }
    
    if(res[0]==1) {
      
      covered+=1;
      //  Rcout<<"yey ="<<y<<" classin= "<<classIndex<<std::endl;
      //  Rcout<<"curr"<<currData(i,0)<<std::endl;
      // break;
      if(classIndex==currData(i,classIndex1)){
        count=count+1;
        //  Rcout<<"curr "<<currData(i,0)<<" class "<<classIndex<<" count "<<count<<std::endl;
      }
    }
    
    if((res[0]==0) &&(mismatch==1)){
      if(classIndex==currData(i,classIndex1)){
        pos_missed_index=(mismatch_index-1)*2; //get positioned after mismatch_index att*2 cells, don't need to inc for next cell, C++ 0 based
        output[pos_missed_index]+=1;
      }
      else{
        neg_missed_index=(mismatch_index-1)*2+1;
        output[neg_missed_index]+=1;
      }
      
    }
    
    
    
    
    
  }
  
  
  if(covered==0)
    fitVal[0]=0;
  else
    fitVal[0]=((double)count/(double)covered);
  // Rcout<<"fit  "<<fitVal[0]<<std::endl;
  
  output[2*attributes]=count;
  output[2*attributes+1]=covered;
  output[2*attributes+2]=fitVal[0];
  
  return output;
  
  
}

// [[Rcpp::export]]
NumericVector confCount(NumericVector y, NumericVector attType, NumericMatrix currData){
  
  int  totalRow=currData.nrow(),endIndex, classIndex,count=0, covered=0 ;
  
  
  NumericVector covRes(1), fitVal(1);
  
  for(int i=0;i<totalRow;i++){
    
    endIndex=1;
    classIndex=-1;
    covRes[0]=0;
    
    //   while(stIndex<n)
    
    
    // int rLen=y.size()-1, attId, j=0;
    
    //Rcout<<xa(i,0)<<std::endl;
    
    covRes=coverCpp1(y, currData, i, attType);
    
    if(covRes[0]==1) {
      
      endIndex=y[y.size()-1];
      classIndex=endIndex;
      covered+=1;
      //  Rcout<<"yey ="<<y<<" classin= "<<classIndex<<std::endl;
      //  Rcout<<"curr"<<currData(i,0)<<std::endl;
      // break;
      if(classIndex==currData(i,classIndex1)){
        count=count+1;
        //  Rcout<<"curr "<<currData(i,0)<<" class "<<classIndex<<" count "<<count<<std::endl;
      }
    }
    
    
    
    
    
    
    
  }
  
  if(covered==0)
    fitVal[0]=0;
  else
    fitVal[0]=((double)count/(double)covered);
  // Rcout<<"fit  "<<fitVal[0]<<std::endl;
  return fitVal;
  
  
}

// [[Rcpp::export]]
NumericVector confCountRand(NumericVector y, NumericVector attType, NumericMatrix currData){
  
  int  totalRow=currData.nrow(),covered=0 ;
  
  
  NumericVector covRes(1), fitVal(3), classMayBe(14);
  
  for(int i=0;i<totalRow;i++){
    
    covRes[0]=0;
    
    
    covRes=coverCpp1(y, currData, i, attType);
    
    if(covRes[0]==1) {
      
      covered+=1;
      
      
      classMayBe[currData(i,classIndex1)]+=1;
      
      
    }
    
    
  }
  
  
  int max=0,classIndex=-1,i; 
  for(i=0;i<=totalClass;i++){
    if(classMayBe[i]>max){
      max=classMayBe[i];
      classIndex=i;
      
    }
    
  }
  
  fitVal[0]=((double)max/(double)covered);
  fitVal[1]=(double)classIndex;
  // Rcout<<"fit  "<<fitVal[0]<<std::endl;
  return fitVal;
  
  
}


// [[Rcpp::export]]
NumericVector covByRuleset(List x, NumericVector attType, NumericMatrix currData){
  
  int  totalRow=currData.nrow(), stIndex, endIndex, classIndex, n, count=0 ;
  n=x.size();
  
  NumericVector covRes(1), covCountVector(totalRow+1);
  
  for(int i=0;i<totalRow;i++){
    stIndex=0;
    endIndex=1;
    classIndex=-1;
    covRes[0]=0;
    
    while(stIndex<n){
      
      SEXP ll = x[stIndex];
      Rcpp::NumericVector y(ll);
      // int rLen=y.size()-1, attId, j=0;
      
      //Rcout<<xa(i,0)<<std::endl;
      
      covRes=coverCpp1(y, currData, i, attType);
      
      if(covRes[0]==1) {
        
        endIndex=y[y.size()-1];
        classIndex=endIndex;
        covCountVector[i]=1;
        //  Rcout<<"yey ="<<y<<" classin= "<<classIndex<<std::endl;
        //  Rcout<<"curr"<<currData(i,0)<<std::endl;
        break;
        
      }
      
      stIndex=stIndex+1;
    }
    
    
    // do not need this block, we are intersted how many data items truely covered by this ruleset, not by default rule
    //   if(stIndex>=n){
    //    classIndex=defaultClass;
    
    // }
    
    
    
    if(classIndex==currData(i,classIndex1)){
      count=count+1;
      //  Rcout<<"curr "<<currData(i,0)<<" class "<<classIndex<<" count "<<count<<std::endl;
    }
    
  }
  
  // Rcout<<"fit  "<<fitVal[0]<<std::endl;
  return covCountVector;
  
  
}