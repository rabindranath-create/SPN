#!/usr/bin/env Rscript


cat("Working directory:", getwd(), "\n")

# Set up and confirm output folder
output_dir <- file.path(getwd(), "outputs/script2")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cat("Created directory:", output_dir, "\n")

# Confirm contents before saving
print("Files in 'outputs' before saving:")
print(list.files("outputs", recursive = TRUE))




library(spatial)
library(spatstat)
library(igraph)

# 1. discretize graph
Graph_Discretized<-function(x,y){
  n=(x+1)*(y+1)
  temp<-graph.empty(n,directed=FALSE)
  G<-set.vertex.attribute(temp,"name",value=c(1:n))
  A<-get.adjacency(G,sparse=FALSE)
  #dim(A)=n by n
  #assume that x,y greater than 1.
  for(j in 1:(y-1)){
    for(i in 1:(x-1)){
      
      A[1+i+j*(x+1),i+j*(x+1)]<-1 #left
      A[1+i+j*(x+1),2+i+j*(x+1)]<-1 #right
      A[1+i+j*(x+1),1+i+(j-1)*(x+1)]<-1 #down
      A[1+i+j*(x+1),1+i+(j+1)*(x+1)]<-1 #up
      
      A[1+i+j*(x+1),i+(j+1)*(x+1)]<-sqrt(2) #left up corner
      A[1+i+j*(x+1),2+i+(j+1)*(x+1)]<-sqrt(2) #right up corner
      A[1+i+j*(x+1),i+(j-1)*(x+1)]<-sqrt(2) #left down corner
      A[1+i+j*(x+1),2+i+(j-1)*(x+1)]<-sqrt(2) #right down corner
      
    }#end of inner loop
  }#end of outer loop
  
  #j=0
  for(i in 1:(x-1)){
    A[1+i,i]<-1; A[1+i,2+i]<-1; A[1+i,1+i+(x+1)]<-1;
    A[1+i,i+(x+1)]<-sqrt(2); A[1+i,2+i+(x+1)]<-sqrt(2); 
  }
  #j=y
  for(i in 1:(x-1)){
    A[1+i+y*(x+1),i+y*(x+1)]<-1; A[1+i+y*(x+1),2+i+y*(x+1)]<-1; A[1+i+y*(x+1),1+i+(y-1)*(x+1)]<-1;
    A[1+i+y*(x+1),i+(y-1)*(x+1)]<-sqrt(2); A[1+i+y*(x+1),2+i+(y-1)*(x+1)]<-sqrt(2); 
  }
  #i=0
  for(j in 1:(y-1)){
    A[1+j*(x+1),1+(j-1)*(x+1)]<-1; A[1+j*(x+1),2+j*(x+1)]<-1; A[1+j*(x+1),1+(j+1)*(x+1)]<-1;
    A[1+j*(x+1),2+(j-1)*(x+1)]<-sqrt(2); A[1+j*(x+1),2+(j+1)*(x+1)]<-sqrt(2); 
  }
  #i=x
  for(j in 1:(y-1)){
    A[1+x+j*(x+1),1+x+(j-1)*(x+1)]<-1; A[1+x+j*(x+1),1+x+(j+1)*(x+1)]<-1; A[1+x+j*(x+1),x+j*(x+1)]<-1;
    A[1+x+j*(x+1),x+(j-1)*(x+1)]<-sqrt(2); A[1+x+j*(x+1),x+(j+1)*(x+1)]<-sqrt(2); 
  }
  #4 Corners
  A[1,2]<-1; A[1,2+x]<-1; A[1,3+x]<-sqrt(2);
  A[1+y*(x+1),1+(y-1)*(x+1)]<-1; A[1+y*(x+1),2+y*(x+1)]<-1; A[1+y*(x+1),2+(y-1)*(x+1)]<-sqrt(2);
  A[1+x,x]<-1; A[1+x,2*(1+x)]<-1; A[1+x,2*(1+x)-1]<-sqrt(2);
  A[(1+x)*(1+y),(1+x)*y]<-1; A[(1+x)*(1+y),(1+x)*(1+y)-1]<-1; A[(1+x)*(1+y),(1+x)*y-1]<-sqrt(2);
  
  G<-graph.adjacency(A,mode=c("undirected"),weighted=TRUE)
  return(G)
}
# 2. intersected edges
Intersect_Obs <- function(c,r,x,y){
  # input is c=(cx,cy) x,y coordinate of a circle. r-radius
  cx<-c[1]; cy<-c[2]
  cx1<-floor(cx); cy1<-floor(cy)
  r1<-ceiling(r)+2
  coor_info <- Lattice_Vertices(x,y)
  
  # points outside the boundary of circle
  X_temp <- matrix(c(rep(-r1:r1,times=length(seq(-r1,r1,by=1))),
                     rep(-r1:r1,each=length(seq(-r1,r1,by=1)))),ncol=2)
  Intersect_temp1 <- function(vector_ij){
    x_i<-vector_ij[1];x_j<-vector_ij[2]
    if(r^2< (cx-(cx1+x_i))^2+(cy-(cy1+x_j))^2 & (cx-(cx1+x_i))^2+(cy-(cy1+x_j))^2 <=(r+sqrt(2))^2){
      return(c(cx1+x_i,cy1+x_j))	
    }
  }
  temp<-as.numeric(unlist(apply(X_temp,1,Intersect_temp1)))  
  case2<-matrix(temp,ncol=2,byrow=TRUE)
  
  # points inside and on the boundary of circle
  Intersect_temp2 <- function(vector_ij){
    x_i<-vector_ij[1];x_j<-vector_ij[2]
    if((r-sqrt(2))^2<(cx-(cx1+x_i))^2+(cy-(cy1+x_j))^2 &  (cx-(cx1+x_i))^2+(cy-(cy1+x_j))^2 <= r^2){
      return(c(cx1+x_i,cy1+x_j))	
    }
  }
  temp<-as.numeric(unlist(apply(X_temp,1,Intersect_temp2)))
  case3<-matrix(temp,ncol=2,byrow=TRUE)
  
  #I am going to form edgelist matrix intersecting the circle
  el<-matrix(0,ncol=2)
  n2=nrow(case2); n3=nrow(case3)
  
  #temp1<-matrix(0,ncol=4)
  
  if(n2!=0 & n3!=0){
    X_temp2 <- matrix(c(rep(1:n3,times=length(seq(1,n2,by=1))),
                        rep(1:n2,each=length(seq(1,n3,by=1)))),ncol=2)
    Intersect_temp3 <- function(vector_ij){
      x_i<-vector_ij[1];x_j<-vector_ij[2]
      if(Dist_Euclidean(case2[x_j,],case3[x_i,])==1 || Dist_Euclidean(case2[x_j,],case3[x_i,])==sqrt(2)){
        e1<-which(coor_info[,1]==case2[x_j,1]&coor_info[,2]==case2[x_j,2])
        e2<-which(coor_info[,1]==case3[x_i,1]&coor_info[,2]==case3[x_i,2])
        return(sort(c(e1,e2),decreasing=F))
      }
    }
    el_vector <- as.numeric(unlist(apply(X_temp2,1,Intersect_temp3)))
    el <- matrix(el_vector,ncol=2,byrow=T)
  }#if statement
  
  #for(i in 2:89){points(x=c(temp1[i,1],temp1[i,3]),y=c(temp1[i,2],temp1[i,4]),type="l",col="red")}
  return(el)
}
# 3. get index for certain coordinates
Index_Coordinates<-function(m,x,y){
  # m is the coordinates of points of interests
  temp<-1+m[1]+m[2]*(x+1)
  return(temp)
}
# 4. distance between two points
Dist_Euclidean <- function(point1,point2){
  return(sqrt((point1[1]-point2[1])^2+(point1[2]-point2[2])^2))
}
# 5. generate coordinates for lattice
Lattice_Vertices <- function(x,y){
  LatticeCoordinates<-matrix(nrow=(x+1)*(y+1),ncol=2)
  temp<-rep(0:x,y+1)	
  LatticeCoordinates[,1]<-temp
  temp<-rep(0,0)
  for(i in 0:y){
    temp<-c(temp,rep(i,x+1))
  }
  LatticeCoordinates[,2]<-temp
  temp<-rep(0,0)
  return(LatticeCoordinates)
}









#20-4

for(kk in c(20, 40, 80)){


  
  obs_info_all1 <- read.csv(paste0('obs_info_all_', kk, '.csv'))


  obs_info_all1[, "cost"] <- 1
  for(i in 1:99){
    obs_info_all1[, paste0("cost.", i)] <- 1
  }




obs_info_all <- list()

for(i in 1:100){
  obs_info_all[[i]] <- obs_info_all1[(5*(i-1)+1):(5*(i-1)+5)]
  colnames(obs_info_all[[i]]) <- c('x','y','cost','prob','status')
}

# RD algorithm
Update_graph_intersect_DT<-function(g,x,y,circle_info,r){
  #read circle center x,y coordinate c-cost,p-prabability,True or False Obstacles
  #circles=read.csv("example1.csv",header=FALSE)
  n <- nrow(circle_info)
  elg <- get.data.frame(g,what="edges")
  colnames(elg) <- c("From","To","Cost")
  int_info <- matrix(0,ncol=nrow(circle_info),nrow=nrow(elg))
  for(i in 1:n){
    el<-Intersect_Obs(t(circle_info[i,1:2]),r,x,y)
    n1=nrow(el)
    dt <- Dist_Euclidean(as.numeric(circle_info[i,1:2]),c(50,1))
    for(k in 1:n1){el[k,]<-sort(el[k,],decreasing=FALSE)} #sort the elements
    
    for(j in 1:n1){
      index=which((elg[,1]==el[j,1] & elg[,2]==el[j,2]))
      elg[index,3] <- elg[index,3]+0.5*(circle_info[i,3]+(dt/(1-circle_info[i,4]))^(-log(1-circle_info[i,4])))
      int_info[index, i] <- 1
    }#inner loop
  }#outer loop
  
  updateg=graph.data.frame(elg,directed=0)
  output <- list(G_info=updateg, Int_info=int_info)
  return(output)
}



for(jj in c(1, 2, 3)){
  



DT_Alg <- function(obs_info){
  W <- jj
  x <- 100; y <- 50; r <- 5
  # begin the loop to travel from s to t
  s <- 5101
  t <- 152
  # create graph - based on W
  vertice_list <- Lattice_Vertices(x,y)
  G_original <- Graph_Discretized(x,y)
  if (W==0){
    obs_info$prob <- Inf
  }
  output_Ginfo <- Update_graph_intersect_DT(G_original, x, y, obs_info, r)
  G_ed <- output_Ginfo$G_info
  Int_info <- output_Ginfo$Int_info
  df_edge_ed <- get.data.frame(G_ed, what="edges")
  # some record vectors
  length_total <- 0 # record Euclidean length
  cost_total <- 0 # record disambiguation cost 
  reach_t <- F
  path_record <- s
  D_record <- c()
  while(reach_t!=T){
    # implement shortest path algorithm
    output <- get.shortest.paths(G_ed,
                                 which(vertex.attributes(G_ed)$name==as.character(s)),
                                 which(vertex.attributes(G_ed)$name==as.character(t)),
                                 weights=df_edge_ed$Cost, output="both",algorithm="dijkstra")
    V_list <- c(as.numeric(attributes(output$vpath[[1]])$names))
    # follow the path until the disambiguation state
    for (i in 2:length(V_list)){
      edge_ind_temp <- which(df_edge_ed$from==min(V_list[(i-1):i])&df_edge_ed$to==max(V_list[(i-1):i]))
      edge_length <- Dist_Euclidean(as.numeric(vertice_list[V_list[i-1],1:2]),as.numeric(vertice_list[V_list[i],1:2]))
      if (sum(Int_info[edge_ind_temp,])!=0){
        D_state <- V_list[i-1]
        D_record <- c(D_record, D_state)
        break
      } else{
        length_total <- length_total+edge_length
        path_record <- c(path_record,V_list[i])
        D_state <- NULL
      }
    }
    if(is.null(D_state)){
      # if didn't run into any obstacle and reach target
      reach_t=T
      output_final <- list(Optimal_found=T,Length_total=length_total,Cost_total=cost_total,
                           Optimal_path=path_record, Disambiguate_state=D_record)
    } else{
      # run into obstacle
      # update start to current disambiguation state
      # subtract one disambiguation 
      reach_t=F
      s <- D_state
      # adjust graph based on W
      if(W==0){
        # identify edges still intersect unclear obstacle
        edge_ind_temp <- which(rowSums(Int_info)>0)
        df_edge_ed[edge_ind_temp,3] <- Inf
      } else{
        # determine which obstacle
        obs_ind_temp <- which(Int_info[edge_ind_temp,]==1)
        if (length(obs_ind_temp)==1){
          if(W<obs_info$cost[obs_ind_temp]){
            df_edge_ed[which(Int_info[,obs_ind_temp]==1),3] <- Inf
          } else{
            W <- W - obs_info$cost[obs_ind_temp]
            # add cost of disambiguation
            cost_total <- cost_total+obs_info[obs_ind_temp,3]
            if(obs_info$status[obs_ind_temp]==1){
              # adjust based on true obstalce
              df_edge_ed[which(Int_info[,obs_ind_temp]==1),3] <- Inf
            } else{
              dt <- Dist_Euclidean(as.numeric(obs_info[obs_ind_temp,1:2]),c(50,1))
              # adjust based on false obstacle
              df_edge_ed[which(Int_info[,obs_ind_temp]==1),3] <- df_edge_ed[which(Int_info[,obs_ind_temp]==1),3]-
                0.5*(obs_info[obs_ind_temp,3]+(dt/(1-obs_info[obs_ind_temp,4]))^(-log(1-obs_info[obs_ind_temp,4])))
              Int_info[which(Int_info[,obs_ind_temp]==1),obs_ind_temp] <- 0 
            }
          }
        } else{
          dist_temp <- rep(0,length(obs_ind_temp))
          for(i in 1:length(obs_ind_temp)){
            dist_temp[i] <- Dist_Euclidean(as.numeric(vertice_list[D_state,1:2]),obs_info[obs_ind_temp[i],1:2])
          }
          obs_ind_temp2 <- obs_ind_temp[which.min(dist_temp)]
          # add cost of disambiguation
          if(W<obs_info$cost[obs_ind_temp2]){
            df_edge_ed[which(Int_info[,obs_ind_temp2]==1),3] <- Inf
          } else{
            W <- W - obs_info$cost[obs_ind_temp2]
            cost_total <- obs_info[obs_ind_temp2,3]
            if (obs_info$status[obs_ind_temp2]==1){
              # true obstacle
              df_edge_ed[which(Int_info[,obs_ind_temp2]==1),3] <- Inf
            } else{
              dt <- Dist_Euclidean(as.numeric(obs_info[obs_ind_temp2,1:2]),c(50,1))
              # false obstacle
              df_edge_ed[which(Int_info[,obs_ind_temp2]==1),3] <- df_edge_ed[which(Int_info[,obs_ind_temp2]==1),3]-
                0.5*(obs_info[obs_ind_temp2,3]+(dt/(1-obs_info[obs_ind_temp2,4]))^(-log(1-obs_info[obs_ind_temp2,4])))
              Int_info[which(Int_info[,obs_ind_temp2]==1),obs_ind_temp2] <- 0
            }
          }
        }
      }
      G_ed <- graph.data.frame(df_edge_ed,directed = F)
    }
  }
  return(output_final)
}
result_DT <- matrix(NA,ncol=3,nrow=100)
for (i in 1:100){
  result <- DT_Alg(obs_info_all[[i]])
  result_DT[i,1] <- result$Length_total
  result_DT[i,2] <- result$Cost_total
  result_DT[i,3] <- length(result$Disambiguate_state)
  #write.csv(result_DT,"result_DT.csv")
  write.csv(result_DT, file = file.path(output_dir, paste0("result_DT_", kk, "_", jj, ".csv") ))
}






}
}










