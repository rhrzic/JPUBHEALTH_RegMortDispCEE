segmented_wrapper <- function(y,x=NULL,se=NULL,npsi=0:5){
  
  require(segmented)
  
  if (is.null(x)) {
    n<-length(y)
    x<-1:n
  } else {
    x = x
  }
  
  initial_model<-glm(y~x, weights = 1/(se^2), family=gaussian)

  ris.ok<-NULL
  ris<-vector("list", length(npsi))  
  id.ris<-rep(FALSE, length(npsi))
  
  cat("Running  ... ")
  
  
  for(i in npsi){
    
    if(npsi[1] == 0){
      j = i+1
    } else {j = i}
    
    ris[[j]]<-suppressWarnings(try(segmented(initial_model, npsi=i, control=seg.control(n.boot=100, seed = 1234, it.max = 1000)), silent=TRUE))
    if(inherits(ris[[j]], "segmented")) {
      id.ris[j]<-TRUE
      ris.ok[[length(ris.ok)+1]]<-ris[[j]]
    }
  }
  
  n.ok<-sum(id.ris)
  r<-sapply(ris.ok, BIC)
  id.ok<-which.min(r)
  
  #browser()
  npsi<-npsi[id.ris] #omit the errors..
  n.psi.ok<-npsi[id.ok]
  o<-ris.ok[[id.ok]]
  cat(paste(n.psi.ok, "breakpoints selected by BIC\n"))
  names(r)<-npsi
  o$bic<-r
  class(o)<- c("segmented", "glm", "lm")
  return(o)
}