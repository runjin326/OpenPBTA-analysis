# UQ transformation
UQ<-function(X){

    uq<-function(y){
      quantile(y, 0.75)
    }
    X<-X+0.1
    upperQ<-apply(X,2,uq)
    f<-upperQ/mean(upperQ) # calculate normalization factor
    res<-scale(X,center=FALSE,scale=f)
    return(res)
  }


# per gene normalization by Median: pgQ2
# X: a matrix of data with the multiplication of factor (f) as:
# f=100 as a default
uq.pgQ2<-function(X, f=100){

      uq.res<-UQ(X) #  perform UQ normalization per sample

    ## perform per gene nomralization
    m<-apply(uq.res,1, median)

    idx<-which(m==0) # avoid 0 median
    m[idx]=0.1

    si<-m/f # calculate normalization factor
    X1<-scale(t(uq.res),center=FALSE,scale=si)
    res<-t(X1)
    return(res)
}


