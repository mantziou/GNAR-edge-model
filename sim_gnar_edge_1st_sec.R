library(igraph)
library(vars)
library(forecast)
library(reshape2)
library(ggplot2)
library(dplyr)

########################################################
####### SIMULATION EXPERIMENTS FOR SECTION 5.1 #######
########################################################

####### SECTION 5.1.1 #######

# Global alpha
# Four cases for each network structure: GNAR(1,[1]), GNAR(1,[2]), GNAR(3,[1,1,1]), GNAR(3,[2,2,2])

# Regimes : Parameter specification

# Regime 1: GNAR(1,[1])
alpha_par_1 <- .2
beta_par_1 <- list(c(0.3))

# Regime 2: GNAR(1,[2])
alpha_par_2 <- .2
beta_par_2 <- list(c(0.3,0.4))

# Regime 3: GNAR(3,[1,1,1])
alpha_par_3 <- c(.2,.4,-0.6)
beta_par_3 <- list(c(0.2),c(0.1),c(-.2))

# Regime 4: GNAR(3,[2,2,2])
alpha_par_4 <- c(.2,.4,-0.6)
beta_par_4 <- list(c(0.3,.1),c(0.1,.1),c(-0.2,.3))

# Regime 5: GNAR(3,[2,0,0])
alpha_par_5 <- c(.2,.4,-0.6)
beta_par_5 <- list(c(0.3,.4),c(0),c(0))


set.seed(44)
seeds <- sample(1:500,50,replace=FALSE)
#SBM parameters
probmat <- matrix(c(.7,.2,.1,.7),nrow = 2)


alphaOrder_reg <- c(1,1,3,3,3)
betaOrder_reg <- list(c(1),c(2),rep(1,3),rep(2,3),c(2,0,0))
num_param <- c(2,3,6,9,5)

for (reg in 1:5){
  for (mod in c("grg","er","sbm")){
    print(c("regime ",reg, " model ",mod))
    cov_mat <- matrix(nrow=num_param[reg],ncol = 50)
    rmse_mat <- matrix(nrow=num_param[reg],ncol = 50)
    for (i in 1:50){
      print(c("seed ",i))
      set.seed(seeds[i])
      
      if (mod=="sbm"){# SBM model
        net_sbm <- sample_sbm(20,probmat,c(10,10),directed = TRUE)
        V(net_sbm)$name <- as.character(seq(1,20,1))# name nodes (characters)
        edgelist_sbm <- get.edgelist(net_sbm)
        nedges_sbm <- ecount(net_sbm)
        if (reg %in% c(1,2)){
          assign(paste("alpha_par_",reg,"_sbm",sep=""),list(rep(get(paste("alpha_par_",reg,sep=""))[1],nedges_sbm)))
        }else{
          assign(paste("alpha_par_",reg,"_sbm",sep=""),list(rep(get(paste("alpha_par_",reg,sep=""))[1],nedges_sbm),
                                                           rep(get(paste("alpha_par_",reg,sep=""))[2],nedges_sbm),
                                                           rep(get(paste("alpha_par_",reg,sep=""))[3],nedges_sbm)))
        }
      }
      if (mod=="er"){# ER model
        net_er <- erdos.renyi.game(20,p.or.m = 168,type = "gnm",directed = TRUE) # p=.4 fix the number of edges in graph ("gnm") rather than give prob of edge ("gnp")
        V(net_er)$name <- as.character(seq(1,20,1))# name nodes (characters)
        edgelist_er <- get.edgelist(net_er)
        nedges_er <- ecount(net_er)
        if (reg %in% c(1,2)){
          assign(paste("alpha_par_",reg,"_er",sep=""),list(rep(get(paste("alpha_par_",reg,sep=""))[1],nedges_er)))
        }else{
          assign(paste("alpha_par_",reg,"_er",sep=""),list(rep(get(paste("alpha_par_",reg,sep=""))[1],nedges_er),
                                                            rep(get(paste("alpha_par_",reg,sep=""))[2],nedges_er),
                                                            rep(get(paste("alpha_par_",reg,sep=""))[3],nedges_er)))
        }
      }
      if (mod=="grg"){# Random Geometric Graph
        lpvs2 <- sample_sphere_surface(dim = 2, n = 20,radius = .7)#.35 for 0.1 density, .7 for density 0.4
        net_grg <- sample_dot_product(lpvs2,directed = TRUE)
        #net_grg <- sample_grg(20, 0.45)
        V(net_grg)$name <- as.character(seq(1,20,1))# name nodes (characters)
        edgelist_grg <- get.edgelist(net_grg)
        nedges_grg <- ecount(net_grg)
        if (reg %in% c(1,2)){
          assign(paste("alpha_par_",reg,"_grg",sep=""),list(rep(get(paste("alpha_par_",reg,sep=""))[1],nedges_grg)))
        }else{
          assign(paste("alpha_par_",reg,"_grg",sep=""),list(rep(get(paste("alpha_par_",reg,sep=""))[1],nedges_grg),
                                                            rep(get(paste("alpha_par_",reg,sep=""))[2],nedges_grg),
                                                            rep(get(paste("alpha_par_",reg,sep=""))[3],nedges_grg)))
        }
      }
      

      
      simdata <- gnar_edge_sim(n=200,net=get(paste("net_",mod,sep = "")),alphaParams=get(paste("alpha_par_",reg,"_",mod,sep = "")),
                                                                         betaParams=get(paste("beta_par_",reg,sep = "")),
                                                                         sigma = 1, meann=0,nedges=get(paste("nedges_",mod,sep = "")),
                                                                        data_edges = get(paste("edgelist_",mod,sep = "")))
      
      fit <- gnar_edge_fit(simdata,get(paste("edgelist_",mod,sep = "")),alphaOrder_reg[reg],
                           betaOrder_reg[[reg]],get(paste("net_",mod,sep = "")),lead_lag_mat = NULL,globalalpha = TRUE,lead_lag_weights = FALSE)

      comp_truevec <- c()
      for (ord in 1:alphaOrder_reg[reg]){
        if (betaOrder_reg[[reg]][ord]>0){
          comp_truevec <- c(comp_truevec,get(paste("alpha_par_",reg,sep = ""))[ord],get(paste("beta_par_",reg,sep = ""))[[ord]])
        }else{
          comp_truevec <- c(comp_truevec,get(paste("alpha_par_",reg,sep = ""))[ord])
        }
      }
      
      for(tr in 1:length(comp_truevec)){
        cov_mat[tr,i] <- between(comp_truevec[tr],confint(fit$mod)[tr,1],confint(fit$mod)[tr,2])
      }
      rmse_mat[,i] <- (fit$mod$coefficients-comp_truevec)^2

    }
    rownames(cov_mat) <- rownames(confint(fit$mod,level=.95))
    rownames(rmse_mat) <- names(fit$mod$coefficients)
    assign(paste("coverage_reg_",reg,"_mod_",mod,sep = ""),apply(cov_mat,1,sum)/50)
    assign(paste("rmse_reg_",reg,"_mod_",mod,sep = ""),apply(rmse_mat,1,function(x) sqrt(mean(x))))
  }
}

###### simulation scenario relevant to real data ######

# ER model
set.seed(12)
net_er <- erdos.renyi.game(86,p.or.m = 6858,type = "gnm",directed = TRUE) 
# name nodes (characters)
V(net_er)$name <- as.character(seq(1,86,1))
edgelist_er <- get.edgelist(net_er)
plot(net_er,vertex.size=15,edge.color="slategrey", vertex.label.color="black",vertex.color="lightblue",
     vertex.label.font=2, edge.arrow.size=0.18,vertex.label.degree=3.7,cex.sub=15,vertex.label.cex=.8) 

nedges <- ecount(net_er)

# Global alpha
# Two cases for ER network: GNAR(4,[1,1,1,1]), GNAR(4,[2,2,2,2]) 

# Regimes : Parameter specification

# Regime 1: GNAR(4,[1,1,1,1])
alpha_par_1 <- c(-0.6,-0.4,-0.2,-0.1)
beta_par_1 <- list(c(0.2),c(0.1),c(0.3),c(0.05))

# Regime 2: GNAR(4,[2,2,2,2]) 
alpha_par_2 <- c(-0.6,-0.4,-0.2,-0.1)
beta_par_2 <- list(c(0.4,-0.4),c(0.3,-0.4),c(0.5,-0.3),c(0.05,-0.1))

# Simulate data for each regime
seed_data <- 100
for (i in 1){
  set.seed(seed_data)
  assign(paste("gnar_edge_sim_data_",i,"_er_real",sep = ""),gnar_edge_sim(n=90,net = net_er,alphaParams = get(paste("alpha_par_",i,sep = "")),
                                                                          betaParams = get(paste("beta_par_",i,sep = "")),sigma = 1, meann=0,
                                                                          nedges=nedges,data_edges=edgelist_er))
  seed_data <- seed_data+10
}


# Fit the model on Simulated data
alphaOrder_reg <- c(4,4)
betaOrder_reg <- list(rep(1,4),rep(2,4))
for (i in 1){
  assign(paste("fit_er_",i,sep = ""),gnar_edge_fit(get(paste("gnar_edge_sim_data_",i,"_er_real",sep = "")),edgelist_er,alphaOrder_reg[i],
                                                   betaOrder_reg[[i]],net_er,lead_lag_mat = NULL,globalalpha = TRUE,lead_lag_weights = FALSE))
}

confint(fit_er_1$mod,  level=0.95)
confint(fit_er_2$mod,  level=0.95)


