#!/usr/bin/env Rscript
#!/usr/bin/env Rscript
cat("Working directory:", getwd(), "\n")
# Set up and confirm output folder
output_dir <- file.path(getwd(), "outputs/script21")
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
# 7. simple node elimination algorithm_modified
WCSPP_Initial_modified <- function(G, s, t, W){
  #--------------------------------------------#
  # INPUT: graph, start node and target node   #
  #        weight constraint                   #
  # OUTPUT: Intial upper bound or trivial case #
  #--------------------------------------------#
  
  # edge list
  df_edge <- get.data.frame(G,what="edges")
  G_new <- graph.data.frame(df_edge,directed = F)
  #--------------------------------------------#
  # minimum weight path
  Q_inf <- get.shortest.paths(G,
                              which(vertex.attributes(G)$name==as.character(s)),
                              which(vertex.attributes(G)$name==as.character(t)),
                              output="both",
                              weights=df_edge$Weight,algorithm="dijkstra")
  W_Qinf <- sum(df_edge[Q_inf[[2]][[1]],4]) # W(min weight path)
  C_Qinf <- sum(df_edge[Q_inf[[2]][[1]],3]) # C(min weight path)
  # check if feasible
  feasibility <- W_Qinf<=W
  if(feasibility==T){
    #--------------------------------------------#
    # minimum cost path
    Q_0 <- get.shortest.paths(G,
                              which(vertex.attributes(G)$name==as.character(s)),
                              which(vertex.attributes(G)$name==as.character(t)),
                              output="both",
                              weights=df_edge$Cost,algorithm="dijkstra")
    W_Q0 <- sum(df_edge[Q_0[[2]][[1]],4]) # W(min cost path)
    C_Q0 <- sum(df_edge[Q_0[[2]][[1]],3]) # C(min cost path)
    # check if optimal found
    optimality <- W_Q0<=W
    if (optimality==FALSE){
      Upper <- C_Qinf
      P_upper <- Q_inf
      for(rep in 1:2){
        G_new <- graph.data.frame(df_edge,directed = F)
        # remove not feasible vertices
        # if minimum weight path through a vertex is over the max weight
        C_inf_fb <- rep(NA,length(V(G_new)))
        W_inf_fb <- rep(NA,length(V(G_new)))
        V_delete <- c()
        E_delete <- c()
        for(k in 1:length(V(G_new))){
          # forward and back trees (in terms of Gamma value)
          Q_inf_f <- get.shortest.paths(G_new,
                                        which(vertex.attributes(G_new)$name==as.character(s)),
                                        V(G_new)[k],output="both",
                                        weights=df_edge$Weight,algorithm="dijkstra")
          Q_inf_b <- get.shortest.paths(G_new,V(G_new)[k],
                                        which(vertex.attributes(G_new)$name==as.character(t)),
                                        output="both",
                                        weights=df_edge$Weight,algorithm="dijkstra")
          if(length(Q_inf_f$vpath[[1]])==0|length(Q_inf_b$vpath[[1]])==0){
            V_delete <- c(V_delete,as.numeric(vertex.attributes(G_new)$name[k]))
            E_delete <- rbind(E_delete,
                              df_edge[which(df_edge$from==as.numeric(vertex.attributes(G_new)$name[k])|
                                              df_edge$to==as.numeric(vertex.attributes(G_new)$name[k])),1:2])
          } else{
            C_inf_fb[k] <- sum(df_edge[Q_inf_f[[2]][[1]],3])+sum(df_edge[Q_inf_b[[2]][[1]],3])
            W_inf_fb[k] <- sum(df_edge[Q_inf_f[[2]][[1]],4])+sum(df_edge[Q_inf_b[[2]][[1]],4])
            # Update Upper bound
            if(C_inf_fb[k]<Upper&W_inf_fb[k]<=W){
              Upper <- C_inf_fb[k]
              P_upper <- list(forward=Q_inf_f, backward=Q_inf_b)
            }
            # remove vertices
            if(W_inf_fb[k]>W){
              V_delete <- c(V_delete,as.numeric(vertex.attributes(G_new)$name[k]))
              E_delete <- rbind(E_delete,
                                df_edge[which(df_edge$from==as.numeric(vertex.attributes(G_new)$name[k])|
                                                df_edge$to==as.numeric(vertex.attributes(G_new)$name[k])),1:2])
            }
          }
        }
        G_update <- delete_vertices(G_new,as.character(V_delete))
        df_edge <- get.data.frame(G_update,what="edges")
        G_new <- graph.data.frame(df_edge,directed = F)
        # remove minimum cost path over the upper bound
        C_0_fb <- c()
        W_0_fb <- c()
        V_delete <- c()
        E_delete <- c()
        for(k in 1:length(V(G_new))){
          # forward and back trees (in terms of Gamma value)
          Q_0_f <- get.shortest.paths(G_new,
                                      which(vertex.attributes(G_new)$name==as.character(s)),
                                      V(G_new)[k],output="both",
                                      weights=df_edge$Cost,algorithm="dijkstra")
          Q_0_b <- get.shortest.paths(G_new,V(G_new)[k],
                                      which(vertex.attributes(G_new)$name==as.character(t)),
                                      output="both",
                                      weights=df_edge$Cost,algorithm="dijkstra")
          C_0_fb[k] <- sum(df_edge[Q_0_f[[2]][[1]],3])+sum(df_edge[Q_0_b[[2]][[1]],3])
          W_0_fb[k] <- sum(df_edge[Q_0_f[[2]][[1]],4])+sum(df_edge[Q_0_b[[2]][[1]],4])
          # Update Upper bound
          if(C_0_fb[k]<Upper&W_0_fb[k]<=W){
            Upper <- C_0_fb[k]
            P_upper <- list(forward=Q_0_f, backward=Q_0_b)
          }
          # remove vertices
          if(round(C_0_fb[k],5)>round(Upper,5)){
            print(k)
            V_delete <- c(V_delete,as.numeric(vertex.attributes(G_new)$name[k]))
            E_delete <- rbind(E_delete,
                              df_edge[which(df_edge$from==as.numeric(vertex.attributes(G_new)$name[k])|
                                              df_edge$to==as.numeric(vertex.attributes(G_new)$name[k])),1:2])
          }
        }
        G_update <- delete_vertices(G_new,as.character(V_delete))
        df_edge <- get.data.frame(G_update,what="edges")
      }
      G_new <- graph.data.frame(df_edge,directed = F)
      # update minimum cost path
      Q_0 <- get.shortest.paths(G_new,
                                which(vertex.attributes(G_new)$name==as.character(s)),
                                which(vertex.attributes(G_new)$name==as.character(t)),
                                output="both",
                                weights=df_edge$Cost,algorithm="dijkstra")
      W_Q0 <- sum(df_edge[Q_0[[2]][[1]],4]) # W(min cost path)
      C_Q0 <- sum(df_edge[Q_0[[2]][[1]],3]) # C(min cost path)
      # update minimum weight path
      Q_inf <- get.shortest.paths(G_new,
                                  which(vertex.attributes(G_new)$name==as.character(s)),
                                  which(vertex.attributes(G_new)$name==as.character(t)),
                                  output="both",
                                  weights=df_edge$Weight,algorithm="dijkstra")
      W_Qinf <- sum(df_edge[Q_inf[[2]][[1]],4]) # W(min weight path)
      C_Qinf <- sum(df_edge[Q_inf[[2]][[1]],3]) # C(min weight path)
      # check if optimal found
      optimality <- W_Q0<=W
      if(optimality==T){
        output <- list(optimal=optimality, Q_inf=Q_inf, value_Qinf=c(C_Qinf,W_Qinf),
                       Q_0=Q_0, value_Q0=c(C_Q0,W_Q0), Value_upper=C_Q0, 
                       Info_path=Q_0, df_edge_update=df_edge, Graph=G_new)
      } else{
        output <- list(optimal=F, Q_inf=Q_inf, value_Qinf=c(C_Qinf,W_Qinf),
                       Q_0=Q_0, value_Q0=c(C_Q0,W_Q0), Value_upper=Upper, 
                       Info_path=P_upper, df_edge_update=df_edge, Graph=G_new) 
      }
    } else{
      output <- list(optimal=optimality, Q_inf=Q_inf, value_Qinf=c(C_Qinf,W_Qinf),
                     Q_0=Q_0, value_Q0=c(C_Q0,W_Q0), Value_upper=C_Q0, 
                     Info_path=Q_0, df_edge_update=df_edge, Graph=G_new)
    }
    return(output)
  } else{
    return(NULL)
  }
} 
Simple_Node_Eliminate_modified2 <- function(G, s, t, W){
  # do WCSPP initialization
  output_initial <- WCSPP_Initial_modified(G, s, t, W)
  if(is.null(output_initial)){
    # no feasible solution
    output <- NA
    return(output)
  } else if(output_initial$optimal==T){
    # minimum cost solution is optimal
    output <- list(Optimal_upper=T, Optimal_multiplier=F, Value_upper=output_initial$Value_upper, 
                   Info_path=output_initial$Info_path, Lambda_pos=0, Lambda_neg=Inf,
                   Lambda_new=0, Graph=output_initial$Graph, Step=1,Value_lower=output_initial$Value_upper) 
    return(output)
  } else{
    # -------------- STEP 1 --------------- #
    # edge list
    df_edge <- output_initial$df_edge_update
    # -------------- STEP 2 --------------- #
    # no trivial case
    lambda_pos <- 0
    lambda_neg <- Inf
    L <- sort(c(lambda_pos,lambda_neg))
    Upper <- output_initial$Value_upper
    P_upper <- output_initial$Info_path
    C_lambda_pos <- output_initial$value_Q0[1]
    C_lambda_neg <- output_initial$value_Qinf[1]
    W_lambda_pos <- output_initial$value_Q0[2]
    W_lambda_neg <- output_initial$value_Qinf[2]
    phi_prime_lambda_pos <- W_lambda_pos-W
    phi_prime_lambda_neg <- W_lambda_neg-W
    # create loop from step3 - step8
    output <- NULL
    while(is.null(output)){
      # -------------- STEP 3 --------------- #
      # new lambda and L_new
      lambda_new <- (C_lambda_pos-C_lambda_neg)/(W_lambda_neg-W_lambda_pos)
      L_new <- C_lambda_pos+lambda_new*(W_lambda_pos-W)
      # shortest path trees for new lambda
      df_edge$Gamma_lambda <- (df_edge$Cost)+lambda_new*(df_edge$Weight)
      G_new <- graph.data.frame(df_edge,directed = F)
      Q_lambda <- get.shortest.paths(G_new,
                                     which(vertex.attributes(G_new)$name==as.character(s)),
                                     which(vertex.attributes(G_new)$name==as.character(t)),
                                     output="both",
                                     weights=df_edge$Gamma_lambda,algorithm="dijkstra")
      C_lambda <- sum(df_edge[Q_lambda[[2]][[1]],3])
      W_lambda <- sum(df_edge[Q_lambda[[2]][[1]],4])
      phi_prime_lambda <- W_lambda-W
      phi_lambda <- C_lambda+lambda_new*phi_prime_lambda
      if(phi_prime_lambda==0){
        output <- list(Optimal_upper=T, Optimal_multiplier=T, Value_upper=phi_lambda, 
                       Info_path=Q_lambda, Lambda_pos=lambda_pos, Lambda_neg=lambda_neg,
                       Lambda_new=lambda_new, Graph=G_new, Step=3,Value_lower=phi_lambda)
        break
      }
      if(!lambda_new%in%L) L <- sort(c(L,lambda_new))
      # -------------- STEP 4 --------------- #
      # loop for each node to update upper bound
      C_lambda_fb <- c()
      W_lambda_fb <- c()
      for(k in 1:length(V(G_new))){
        # forward and back trees (in terms of Gamma value)
        Q_lambda_f <- get.shortest.paths(G_new,
                                         which(vertex.attributes(G_new)$name==as.character(s)),
                                         V(G_new)[k],output="both",
                                         weights=df_edge$Gamma_lambda,algorithm="dijkstra")
        Q_lambda_b <- get.shortest.paths(G_new,V(G_new)[k],
                                         which(vertex.attributes(G_new)$name==as.character(t)),
                                         output="both",
                                         weights=df_edge$Gamma_lambda,algorithm="dijkstra")
        C_lambda_fb[k] <- sum(df_edge[Q_lambda_f[[2]][[1]],3])+sum(df_edge[Q_lambda_b[[2]][[1]],3])
        W_lambda_fb[k] <- sum(df_edge[Q_lambda_f[[2]][[1]],4])+sum(df_edge[Q_lambda_b[[2]][[1]],4])
        # check if new upper bound available
        if(C_lambda_fb[k]<Upper&W_lambda_fb[k]<=W){
          Upper <- C_lambda_fb[k]
          P_upper <- list(forward=Q_lambda_f, backward=Q_lambda_b)
        }
      }
      # -------------- STEP 5 --------------- #
      # define vector to store node and edge to be deleted
      V_delete <- c()
      E_delete <- c()
      phi_lambda_fb <- c()
      index_delete <- c()
      # loop over node to complete V_delete and A_delete
      for (k in 1:length(V(G_new))){
        phi_lambda_fb[k] <- C_lambda_fb[k]+lambda_new*(W_lambda_fb[k]-W)
        if (round(phi_lambda_fb[k],5)>round(Upper,5)){
          index_delete <- c(index_delete,k)
          V_delete <- c(V_delete,as.numeric(vertex.attributes(G_new)$name[k]))
          E_delete <- rbind(E_delete,
                            df_edge[which(df_edge$from==as.numeric(vertex.attributes(G_new)$name[k])|
                                            df_edge$to==as.numeric(vertex.attributes(G_new)$name[k])),1:2])
        }
      }
      # update graph vertice set and edge set
      G_update <- delete_vertices(G_new,as.character(V_delete))
      df_edge <- get.data.frame(G_update,what="edges")
      # -------------- STEP 6 --------------- #
      # new path based on lambda_new and new graph
      if(is.null(index_delete)){
        output_step6 <- round(min(phi_lambda_fb),5)==round(Upper,5)
      } else{
        output_step6 <- round(min(phi_lambda_fb[setdiff(seq_along(phi_lambda_fb), index_delete)]),5) == round(Upper,5)
      }
      if(output_step6==T){
        # check if stop with empty graph
        # current upper bound is optimal
        output <- list(Optimal_upper=T, Optimal_multiplier=F, Value_upper=Upper, 
                       Info_path=P_upper, Lambda_pos=lambda_pos, Lambda_neg=lambda_neg,
                       Lambda_new=lambda_new, Graph=G_update, Step=6,Value_lower=Upper)
      } else{
        Q_lambda <- get.shortest.paths(G_update,
                                       which(vertex.attributes(G_update)$name==as.character(s)),
                                       which(vertex.attributes(G_update)$name==as.character(t)),
                                       output="both",
                                       weights=df_edge$Gamma_lambda,algorithm="dijkstra")
        C_lambda <- sum(df_edge[Q_lambda[[2]][[1]],3])
        W_lambda <- sum(df_edge[Q_lambda[[2]][[1]],4])
        phi_prime_lambda <- W_lambda-W
        phi_lambda <- C_lambda+lambda_new*phi_prime_lambda
        if(round(L_new,5)==round(phi_lambda,5)){
          # optimal multiplier found
          output <- list(Optimal_upper=T, Optimal_multiplier=T, Value_upper=Upper, 
                         Info_path=P_upper, Lambda_pos=lambda_pos, Lambda_neg=lambda_neg,
                         Lambda_new=lambda_new, Graph=G_update,Step=6, Value_lower=min(phi_lambda_fb[setdiff(seq_along(phi_lambda_fb), index_delete)]))
        } else if(phi_prime_lambda<0){
          # -------------- STEP 7 --------------- #
          # continue to search for new lambda  
          # lambda_new as new lambda_neg
          lambda_neg <- lambda_new
          C_lambda_neg <- C_lambda
          W_lambda_neg <- W_lambda
          phi_prime_lambda_neg <- phi_prime_lambda
          phi_lambda_neg <- phi_lambda
          # if tree out-of-date, recaculate for lambda_pos
          if (!is.empty(V_delete)){
            # shortest path tree update for lambda_pos
            df_edge$Gamma_lambda <- (df_edge$Cost)+lambda_pos*(df_edge$Weight)
            G_update <- graph.data.frame(df_edge,directed = F)
            Q_lambda_pos <- get.shortest.paths(G_update,
                                               which(vertex.attributes(G_update)$name==as.character(s)),
                                               which(vertex.attributes(G_update)$name==as.character(t)),
                                               output="both",
                                               weights=df_edge$Gamma_lambda,algorithm="dijkstra")
            C_lambda_pos <- sum(df_edge[Q_lambda_pos[[2]][[1]],3])
            W_lambda_pos <- sum(df_edge[Q_lambda_pos[[2]][[1]],4])
            phi_prime_lambda_pos <- W_lambda_pos-W
            phi_lambda_pos <- C_lambda_pos+lambda_pos*phi_prime_lambda_pos
            # check if still positive - stop once reach 0 or becomes positive
            while (phi_prime_lambda_pos<=0 & lambda_pos!=0){
              ind_pos <- which(L==lambda_pos)
              lambda_pos <- L[ind_pos-1]
              df_edge$Gamma_lambda <- (df_edge$Cost)+lambda_pos*(df_edge$Weight)
              G_update <- graph.data.frame(df_edge,directed = F)
              Q_lambda_pos <- get.shortest.paths(G_update,
                                                 which(vertex.attributes(G_update)$name==as.character(s)),
                                                 which(vertex.attributes(G_update)$name==as.character(t)),
                                                 output="both",
                                                 weights=df_edge$Gamma_lambda,algorithm="dijkstra")
              C_lambda_pos <- sum(df_edge[Q_lambda_pos[[2]][[1]],3])
              W_lambda_pos <- sum(df_edge[Q_lambda_pos[[2]][[1]],4])
              phi_prime_lambda_pos <- W_lambda_pos-W
              phi_lambda_pos <- C_lambda_pos+lambda_pos*phi_prime_lambda_pos
            }
          }
          # check if need to stop
          if (lambda_pos==0&phi_prime_lambda_pos<=0){
            # check if current upper bound is optimal
            output <- list(Optimal_upper=T, Optimal_multiplier=F, Value_upper=Upper, 
                           Info_path=P_upper, Lambda_pos=lambda_pos, Lambda_neg=lambda_neg,
                           Lambda_new=lambda_new, Graph=G_update, Step=7,Value_lower=Upper)
          }
        } else{
          # lambda_new as new lambda_pos
          lambda_pos <- lambda_new
          C_lambda_pos <- C_lambda
          W_lambda_pos <- W_lambda
          phi_prime_lambda_pos <- phi_prime_lambda
          phi_lambda_pos <- phi_lambda
          # if tree out-of-date, recaculate for lambda_neg
          if (!is.empty(V_delete)){
            # shortest path tree update for lambda_neg
            if(lambda_neg==Inf){
              df_edge$Gamma_lambda <- (df_edge$Weight)
              G_update <- graph.data.frame(df_edge,directed = F)
              Q_lambda_neg <- get.shortest.paths(G_update,
                                                 which(vertex.attributes(G_update)$name==as.character(s)),
                                                 which(vertex.attributes(G_update)$name==as.character(t)),
                                                 output="both",
                                                 weights=df_edge$Gamma_lambda,algorithm="dijkstra")
              C_lambda_neg <- sum(df_edge[Q_lambda_neg[[2]][[1]],3])
              W_lambda_neg <- sum(df_edge[Q_lambda_neg[[2]][[1]],4])
              phi_prime_lambda_neg <- W_lambda_neg-W
            } else{
              df_edge$Gamma_lambda <- (df_edge$Cost)+lambda_neg*(df_edge$Weight)
              G_update <- graph.data.frame(df_edge,directed = F)
              Q_lambda_neg <- get.shortest.paths(G_update,
                                                 which(vertex.attributes(G_update)$name==as.character(s)),
                                                 which(vertex.attributes(G_update)$name==as.character(t)),
                                                 output="both",
                                                 weights=df_edge$Gamma_lambda,algorithm="dijkstra")
              C_lambda_neg <- sum(df_edge[Q_lambda_neg[[2]][[1]],3])
              W_lambda_neg <- sum(df_edge[Q_lambda_neg[[2]][[1]],4])
              phi_prime_lambda_neg <- W_lambda_neg-W
              phi_lambda_neg <- C_lambda_neg+lambda_neg*phi_prime_lambda_neg
            }
            # check if still negative - stop once reach inf or becomes negative
            while (phi_prime_lambda_neg>0 & lambda_neg!=Inf){
              ind_neg <- which(L==lambda_neg)
              lambda_neg <- L[ind_neg+1]
              if (lambda_neg==Inf){
                Q_lambda_neg <- get.shortest.paths(G_update,
                                                   which(vertex.attributes(G_update)$name==as.character(s)),
                                                   which(vertex.attributes(G_update)$name==as.character(t)),
                                                   output="both",
                                                   weights=df_edge$Weight,algorithm="dijkstra")
                C_lambda_neg <- sum(df_edge[Q_lambda_neg[[2]][[1]],3])
                W_lambda_neg <- sum(df_edge[Q_lambda_neg[[2]][[1]],4])
                phi_prime_lambda_neg <- W_lambda_neg-W
                phi_lambda_neg <- C_lambda_neg+lambda_neg*phi_prime_lambda_neg
              } else{
                df_edge$Gamma_lambda <- (df_edge$Cost)+lambda_neg*(df_edge$Weight)
                G_update <- graph.data.frame(df_edge,directed = F)
                Q_lambda_neg <- get.shortest.paths(G_update,
                                                   which(vertex.attributes(G_update)$name==as.character(s)),
                                                   which(vertex.attributes(G_update)$name==as.character(t)),
                                                   output="both",
                                                   weights=df_edge$Gamma_lambda,algorithm="dijkstra")
                C_lambda_neg <- sum(df_edge[Q_lambda_neg[[2]][[1]],3])
                W_lambda_neg <- sum(df_edge[Q_lambda_neg[[2]][[1]],4])
                phi_prime_lambda_neg <- W_lambda_neg-W
                phi_lambda_neg <- C_lambda_neg+lambda_neg*phi_prime_lambda_neg
              }
            }
          }
          # check if need to stop
          if (lambda_neg==Inf&W_lambda_neg>W){
            # check if feasible solution exist in current
            output <- list(Optimal_upper=T, Optimal_multiplier=F, Value_upper=Upper, 
                           Info_path=P_upper, Lambda_pos=lambda_pos, Lambda_neg=lambda_neg,
                           Lambda_new=lambda_new, Graph=G_update, Step=7,Value_lower=Upper)
          }
        } 
      }
    }
    return(output)
  }
}

# Generate Obstacle information




alpha <- 21


obs_info_all1 <- read.csv('obs_info_all_80.csv')

obs_info_all1[, "cost"] <- 1
for(i in 1:99){
  obs_info_all1[, paste0("cost.", i)] <- 1
}


obs_info_all <- list()
for(i in 1:100){
  obs_info_all[[i]] <- obs_info_all1[(5*(i-1)+1):(5*(i-1)+5)]
  colnames(obs_info_all[[i]]) <- c('x','y','cost','prob','status')
}

# SR risk function - alpha=30
Update_graph_intersect<-function(g,x,y,circle_info,r){
  #read circle center x,y coordinate c-cost,p-prabability,True or False Obstacles
  #circles=read.csv("example1.csv",header=FALSE)
  n <- nrow(circle_info)
  elg <- get.data.frame(g,what="edges")
  colnames(elg) <- c("From","To","Cost")
  elg$Weight <- rep(0,nrow(elg))
  int_info <- matrix(0,ncol=nrow(circle_info),nrow=nrow(elg))
  for(i in 1:n){
    el<-Intersect_Obs(t(circle_info[i,1:2]),r,x,y)
    n1=nrow(el)
    for(k in 1:n1){el[k,]<-sort(el[k,],decreasing=FALSE)} #sort the elements
    
    for(j in 1:n1){
      index=which((elg[,1]==el[j,1] & elg[,2]==el[j,2]))
      elg[index,3] <- elg[index,3]-0.5*alpha*log(1-circle_info[i,4])
      elg[index,4] <- elg[index,4]+0.5*circle_info$cost[i]
      int_info[index, i] <- 1
    }#inner loop
  }#outer loop
  
  updateg=graph.data.frame(elg,directed=0)
  output <- list(G_info=updateg, Int_info=int_info)
  return(output)
}



WCSPP_Node_risk_30 <- function(obs_info){
  W <- 2
  x <- 100; y <- 50; r <- 5
  # create graph
  vertice_list <- Lattice_Vertices(x,y)
  G_original <- Graph_Discretized(x,y)
  output_Ginfo <- Update_graph_intersect(G_original, x, y, obs_info, r)
  G_ed <- output_Ginfo$G_info
  Int_info <- output_Ginfo$Int_info
  df_edge_ed <- get.data.frame(G_ed, what="edges")
  # begin the loop to travel from s to t
  s <- 5101
  t <- 152
  length_total <- 0 # record Euclidean length
  cost_total <- 0 # record disambiguation cost 
  reach_t <- F
  path_record <- s
  D_record <- c()
  LU_diff <- c()
  while(reach_t!=T){
    # implement simple node elimination algorithm
    output <- Simple_Node_Eliminate_modified2(G_ed, s, t, W)
    if (length(output)==1){
      output_final <- list(Optimal_found=F)
      break
    }
    if (output$Optimal_upper==T){
      # if optimal solution found
      if (output$Step==1){
        # if minimum cost path is optimal
        P_optimal <- output$Info_path
        V_list <- c(as.numeric(attributes(P_optimal$vpath[[1]])$names))
        LU_diff <- c(LU_diff,(output$Value_upper-output$Value_lower)/output$Value_lower)
      } else if(output$Step==3){
        P_optimal <- output$Info_path
        V_list <- c(as.numeric(attributes(P_optimal$vpath[[1]])$names))
        LU_diff <- c(LU_diff,(output$Value_upper-output$Value_lower)/output$Value_lower)
      } else{
        # optimal found in the SNE process
        P_optimal <- output$Info_path
        V_list <- c(as.numeric(attributes(P_optimal$forward$vpath[[1]])$names),as.numeric(attributes(P_optimal$backward$vpath[[1]])$names)[-1])
        LU_diff <- c(LU_diff,0)
      }
      #follow the path until the first disambiguation state
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
                             Optimal_path=path_record, Disambiguate_state=D_record, LU_diff=LU_diff)
      } else{
        # run into obstacle
        # update start to current disambiguation state
        # subtract one disambiguation 
        reach_t=F
        s <- D_state
        # determine which obstacle
        obs_ind_temp <- which(Int_info[edge_ind_temp,]==1)
        if (length(obs_ind_temp)==1){
          W <- W - obs_info$cost[obs_ind_temp]
          # add cost of disambiguation
          cost_total <- cost_total+obs_info[obs_ind_temp,3]
          if (obs_info$status[obs_ind_temp]==1){
            # asjust based on true obstacle
            df_edge_ed[which(Int_info[,obs_ind_temp]==1),3] <- Inf
          } else{
            # adjust based on false obstacle
            df_edge_ed[which(Int_info[,obs_ind_temp]==1),3] <- df_edge_ed[which(Int_info[,obs_ind_temp]==1),3]+
              0.5*alpha*log(1-obs_info[obs_ind_temp,4])
            df_edge_ed[which(Int_info[,obs_ind_temp]==1),4] <- df_edge_ed[which(Int_info[,obs_ind_temp]==1),4]-
              0.5*obs_info[obs_ind_temp,3]
            Int_info[which(Int_info[,obs_ind_temp]==1),obs_ind_temp] <- 0
          }
        } else{
          dist_temp <- rep(0,length(obs_ind_temp))
          for(i in 1:length(obs_ind_temp)){
            dist_temp[i] <- Dist_Euclidean(as.numeric(vertice_list[D_state,1:2]),obs_info[obs_ind_temp[i],1:2])
          }
          obs_ind_temp2 <- obs_ind_temp[which.min(dist_temp)]
          W <- W-obs_info$cost[obs_ind_temp2]
          # add cost of disambiguation
          cost_total <- cost_total+obs_info[obs_ind_temp2,3]
          if(obs_info$status[obs_ind_temp2]==1){
            # true obstacle
            df_edge_ed[which(Int_info[,obs_ind_temp2]==1),3] <- Inf
          } else{
            # false obstacle
            df_edge_ed[which(Int_info[,obs_ind_temp2]==1),3] <- df_edge_ed[which(Int_info[,obs_ind_temp2]==1),3]+
              0.5*alpha*log(1-obs_info[obs_ind_temp2,4])
            df_edge_ed[which(Int_info[,obs_ind_temp2]==1),4] <- df_edge_ed[which(Int_info[,obs_ind_temp2]==1),4]-
              0.5*obs_info[obs_ind_temp2,3]
            Int_info[which(Int_info[,obs_ind_temp2]==1),obs_ind_temp2] <- 0
          }
        }
        G_ed <- graph.data.frame(df_edge_ed,directed = F)
      }
    } else{
      output_final <- list(Optimal_found=F)
    }
  }
  return(output_final)
}





result_WCSPP_risk_30 <- matrix(NA,ncol=7,nrow=100)
write.csv(result_WCSPP_risk_30, file = file.path(output_dir, paste0("result_WCSPP_risk_", alpha, "_80_1.csv")))
for (i in 1:10){
    obs_info_all_use <- obs_info_all[[i]]
    result <- WCSPP_Node_risk_30(obs_info_all_use)
    result_WCSPP_risk_30[i, 1] <- result$Length_total
    result_WCSPP_risk_30[i,2] <- result$Cost_total
    result_WCSPP_risk_30[i,3] <- length(result$Disambiguate_state)
    result_WCSPP_risk_30[i,4] <- result$LU_diff[1]
    result_WCSPP_risk_30[i,5] <- result$LU_diff[2]
    result_WCSPP_risk_30[i,6] <- result$LU_diff[3]
    result_WCSPP_risk_30[i,7] <- result$LU_diff[4]
    write.csv(result_WCSPP_risk_30, file = file.path(output_dir, paste0("result_WCSPP_risk_", alpha, "_80_1.csv")))
  }




