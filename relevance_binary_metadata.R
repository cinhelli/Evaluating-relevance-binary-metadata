library(igraph)
#Angles are in radians, not degrees, 
#for the standard versions (i.e., a right angle is ??/2), and in 'half-rotations' for cospi etc.
Relevance <- function(g, binary_vector){
  
  ############################computation of n1, m11 (1-1), m10 (1-0)
  n1 <- sum(binary_vector)
  binary_meta <- binary_vector
  c_0 <- which(binary_meta==0)
  m00_subgraph <- induced.subgraph(g, c_0)
  m00 <- ecount(m00_subgraph)
  
  c_1 <- which(binary_meta==1)
  m11_subgraph <- induced.subgraph(g, c_1)
  m11 <- ecount(m11_subgraph)
  m10 <- ecount(g) - m11 - m00
  ############################ expected values and bounds
  m11.med<-function(n1,d){
    choose(n1,2)*d
  }
  m10.med<-function(n1,N,d){
    n1*(N-n1)*d
  }
  m11.max.new <- function(n1,M,deg_seq){ #deg_seq deve essere ordinata dal pi? grande al pi? piccolo
    ds1 <- deg_seq[1:n1]
    ds1[ds1 > (n1-1)] <- n1-1 
    min(M,choose(n1,2),ceiling(sum(ds1)/2))
  }
  m10.max.new <- function(n1,N,M,deg_seq){
    ds1 <- deg_seq[1:n1]
    ds1[ds1 > (N - n1)] <- N - n1
    ds0 <- deg_seq[1:(N - n1)]
    ds0[ds0 > n1] <- n1
    min(M,n1*(N-n1),min(sum(ds1),sum(ds0)))
  }
  m11.min.new.new <- function(n1,N,deg_seq){
    n0 <- N - n1
    if(n1!=N){
      firstpart <- deg_seq[N:(N-n1+1)]
      #firstpart[firstpart>n0] <- n0
      secondpart <- deg_seq[1:(N-n1)]
      secondpart[secondpart>n1] <- n1
      max(0,floor((sum(firstpart)-sum(secondpart))/2))#il +1 ? per la vettorizzazione
    }else(sum(deg_seq)/2)
  }
  m10.min.new.new <- function(n1,N,deg_seq){
    t_f <- 0
    n0 <- N - n1
    if(n1 != 0 & n1 != N ){
      t1 <- 1
      t2 <- sum(deg_seq[N:(N-n1+1)])-n1*(n1-1)
      t3 <- sum(deg_seq[N:(N-n0+1)])-n0*(n0-1)
      if(n1 > n0){
        tmp_t4 <- deg_seq[1:n1]-(n1-1)
        tmp_t4[tmp_t4<0] <- 0
        t4 <- sum(tmp_t4)/2
        t_f <- max(t1,t2,t3,t4)
      }
      if(n1 <= n0){
        tmp_t4 <- deg_seq[1:n0]-(n0-1)
        tmp_t4[tmp_t4<0] <- 0
        t4 <- sum(tmp_t4)/2
        t_f <- max(t1,t2,t3,t4)
      }
    }
    return(t_f)
  }
  ############################
  Em11 <- m11.med(n1, graph.density(g))
  Em10 <- m10.med(n1,vcount(g), graph.density(g))
  UBm11 <- m11.max.new(n1,ecount(g),sort(degree(g),decreasing = T))
  UBm10 <- m10.max.new(n1,vcount(g),ecount(g),sort(degree(g),decreasing = T))
  LBm11 <- m11.min.new.new(n1,vcount(g),sort(degree(g),decreasing = T))
  LBm10 <- m10.min.new.new(n1,vcount(g),sort(degree(g),decreasing = T))
  Hmin <- LBm10/Em10
  Dmin <- LBm11/Em11
  Hmax <- UBm10/Em10
  Dmax <- UBm11/Em11
  
  H <- m10/Em10
  D <- m11/Em11
  
  vA <- sqrt((1-Hmin)^2+(Dmax-1)^2)
  OA <- acos((1-Hmin)/vA)#in radianti
  vB <- sqrt((Hmax-1)^2+(Dmax-1)^2)
  OB <- acos((Hmax-1)/vB)
  vC <- sqrt((Hmax-1)^2+(1-Dmin)^2)
  OC <- acos((Hmax-1)/vC)
  vD <- sqrt((1-Hmin)^2+(1-Dmin)^2)
  OD <- acos((1-Hmin)/vD)
  
  if(H==1 && D ==1){s <- 0}
  
  if(H < 1 && D > 1){
    print("a")
    viA <- sqrt((1-H)^2+(D-1)^2)
    OiA <- acos((1-H)/viA)
    
    if(OiA<=OA){
      pviA <- viA*cos(OA-OiA)
    }else{
      pviA <- viA*cos(OiA-OA)
    }
    
    s <- pviA/vA
  }
  if(H > 1 && D > 1){
    print("B")
    viB <- sqrt((H-1)^2+(D-1)^2)
    OiB <- acos((H-1)/viB)
    
    if(OiB<=OB){
      pviB <- viB*cos(OB-OiB)
    }else{
      pviB <- viB*cos(OiB-OB)
    }
    
    s <- pviB/vB
  }
  if(H > 1 && D < 1){
    print("C")
    viC <- sqrt((H-1)^2+(1-D)^2)
    OiC <- acos((H-1)/viC)
    
    if(OiC<=OC){
      pviC <- viC*cos(OC-OiC)
    }else{
      pviC <- viC*cos(OiC-OC)
    }
    
    s <- pviC/vC
  }
  if(H < 1 && D < 1){
    print("D")
    viD <- sqrt((1-H)^2+(1-D)^2)
    OiD <- acos((1-H)/viD)
    
    if(OiD<=OD){
      pviD <- viD*cos(OD-OiD)
    }else{
      pviD <- viD*cos(OiD-OD)
    }
    
    s <- pviD/vD
  }
  return(list("n1"=n1,"m11"=m11,"m10"=m10,"D"=D,"H"=H,"Relevance"=s))
}
#####################
#working example
#####################
library(igraph)
library(igraphdata)
#
data("yeast")
#sum(table(V(yeast)$Class, useNA="ifany"))
yeastBridge <- clusters(yeast)
yeastBridge <- which(yeastBridge$membership == which.max(yeastBridge$csize))
#
yeast1 <- induced.subgraph(yeast, yeastBridge)
#sum(table(V(yeast1)$Class, useNA="ifany"))
rm(yeastBridge)
rm(yeast)
#
# N <- vcount(yeast1)
# M <- ecount(yeast1)
# d <- graph.density(yeast1)
# DS <- sort(degree(yeast1), decreasing = T)
#
V(yeast1)$Class
V(yeast1)$Class[is.na(V(yeast1)$Class)] <- "MISS"
classi <- as.vector(table(V(yeast1)$Class))
names(classi) <- names(table(V(yeast1)$Class))

BinaryVectors <- matrix(NA, nrow = vcount(yeast1), ncol = length(classi), dimnames = list(c(1:vcount(yeast1)), names(classi)))
for(i in 1:length(classi)){
  nome <- names(classi[i])
  
  vettore <- V(yeast1)$Class
  vettore[which(vettore==nome)] <- 1
  vettore[which(vettore != 1) ] <- 0
  vettore[which(is.na(vettore))]<- 0
  vettore <- as.integer(vettore)
  
  BinaryVectors[,nome] <- vettore
}
BinaryVectors <- BinaryVectors[,-which(colnames(BinaryVectors)=="MISS")]
#
####Functional class P
#
Relevance(yeast1,BinaryVectors[,10])
