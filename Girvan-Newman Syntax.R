library(igraph)
library(sna)
library(network)
library(statnet)
setwd("C:/Users/Allison Barnard/Dropbox/SNH (1)/AlopeciaProject2016/Community Structure/")
getwd()

##########           ##########
########## NC9 Groom ##########
##########           ##########

#Read your edgelist into R; create a graph object
NC9_groom_initialEdgelist=read.csv("9D_GroomExHPEdgelist.csv", header = TRUE)
#Change the items in groom to characters because it automatically reads in as factor regardless of type
NC9_groom_stringEdgelist<-data.frame(lapply(NC9_groom_initialEdgelist, as.character), stringsAsFactors = FALSE)
#Change the data frame to a graph object
NC9_groom_graphDF<-graph.data.frame(NC9_groom_stringEdgelist)
#Since we're using this for network measures, it should be undirected no matter what
NC9_groom_graph_undirected <- as.undirected(NC9_groom_graphDF, mode='collapse')


#Run the Girvan-Newman and Walk-Trap algorithms;
##Try Walk-Trap with 3 and 5 steps choose the larger modularity score 
###If mods are equal, choose the one that returns the smaller number of groups
NC9_groom_community_eb <- edge.betweenness.community(NC9_groom_graph_undirected, weights = E(NC9_groom_graphDF)$Weight)
NC9_groom_community_wk_5 <- walktrap.community(NC9_groom_graph_undirected, weights = E(NC9_groom_graphDF)$Weight, steps = 5, merges=TRUE, modularity=TRUE, membership=TRUE)
NC9_groom_community_wk_3 <- walktrap.community(NC9_groom_graph_undirected, weights = E(NC9_groom_graphDF)$Weight, steps = 3, merges=TRUE, modularity=TRUE, membership=TRUE)

#Return the modularity for the best-fit cluster
NC9_groom_community_eb_mod<-(NC9_groom_community_eb$modularity)
NC9_groom_community_eb_mod_max<-max(NC9_groom_community_eb_mod)

NC9_groom_community_wk_5_mod<-(NC9_groom_community_wk_5$modularity)
NC9_groom_community_wk_5_mod_max<-max(NC9_groom_community_wk_5_mod)

NC9_groom_community_wk_3_mod<-(NC9_groom_community_wk_3$modularity)
NC9_groom_community_wk_3_mod_max<-max(NC9_groom_community_wk_3_mod)

wk5<-NC9_groom_community_wk_5_mod_max
wk3<-NC9_groom_community_wk_3_mod_max

## Choose the better Walk-Trap model (3 or 5) and output that one as "...wk_bestFit" to use for output
if (wk5 > wk3) {
  if((wk5-wk3)>0.01){
    NC9_groom_community_wk_bestFit <- NC9_groom_community_wk_5
    NC9_groom_community_wk_mod_max_bestFit <-NC9_groom_community_wk_5_mod_max
    NC9_groom_communities_wk_cBind=cbind(NC9_groom_community_wk_bestFit$names, NC9_groom_community_wk_bestFit$membership, NC9_groom_community_wk_mod_max_bestFit,"5")
  } else
    if (max(NC9_groom_community_wk_5$membership) < max(NC9_groom_community_wk_3$membership)) {
      NC9_groom_community_wk_bestFit <- NC9_groom_community_wk_5
      NC9_groom_community_wk_mod_max_bestFit <-NC9_groom_community_wk_5_mod_max
      NC9_groom_communities_wk_cBind=cbind(NC9_groom_community_wk_bestFit$names, NC9_groom_community_wk_bestFit$membership, NC9_groom_community_wk_mod_max_bestFit,"5")
    } else {
      NC9_groom_community_wk_bestFit <- NC9_groom_community_wk_3
      NC9_groom_community_wk_mod_max_bestFit <-NC9_groom_community_wk_3_mod_max
      NC9_groom_communities_wk_cBind=cbind(NC9_groom_community_wk_bestFit$names, NC9_groom_community_wk_bestFit$membership, NC9_groom_community_wk_mod_max_bestFit,"3")
    }
} else {
  if((wk3-wk5)>0.01) {
    NC9_groom_community_wk_bestFit <- NC9_groom_community_wk_3
    NC9_groom_community_wk_mod_max_bestFit <-NC9_groom_community_wk_3_mod_max
    NC9_groom_communities_wk_cBind=cbind(NC9_groom_community_wk_bestFit$names, NC9_groom_community_wk_bestFit$membership, NC9_groom_community_wk_mod_max_bestFit,"3")
  } else
    if (max(NC9_groom_community_wk_3$membership) < max(NC9_groom_community_wk_5$membership)) {
      NC9_groom_community_wk_bestFit <- NC9_groom_community_wk_3
      NC9_groom_community_wk_mod_max_bestFit <-NC9_groom_community_wk_3_mod_max
      NC9_groom_communities_wk_cBind=cbind(NC9_groom_community_wk_bestFit$names, NC9_groom_community_wk_bestFit$membership, NC9_groom_community_wk_mod_max_bestFit,"3")
    } else {
      NC9_groom_community_wk_bestFit <- NC9_groom_community_wk_5
      NC9_groom_community_wk_mod_max_bestFit <-NC9_groom_community_wk_5_mod_max
      NC9_groom_communities_wk_cBind=cbind(NC9_groom_community_wk_bestFit$names, NC9_groom_community_wk_bestFit$membership, NC9_groom_community_wk_mod_max_bestFit,"5")
    }
 }

#Write out animal ID's, their community memberships, the modularity score for the best-fit cluster, and the number of steps used for the best Walk-Trap cluster (3 or 5)
NC9_groom_communities_eb_cBind=cbind(NC9_groom_community_eb$names, NC9_groom_community_eb$membership, NC9_groom_community_eb_mod_max)
write.csv(NC9_groom_communities_eb_cBind,"9D_GroomExHP_Communities_EdgeBetweenness.csv")
write.csv(NC9_groom_communities_wk_cBind,"9D_GroomExHP_Communities_WalkTrap.csv")
