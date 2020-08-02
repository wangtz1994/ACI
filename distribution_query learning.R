rm(list = ls())
library("pcalg")
poss = 0.6
for (n in c(2500)) {
  for (total_cost in c(1:4)) {
    res = c()
    rres = c()
    cost_a1 = 0  #cost of once intervention
    cost_a2 = 1  #cost of once observation (observing the response variable) under intervention
    node_num = 30 # variable number
    for (asd in c(1:100)) {
      set.seed(asd)
      
      
      setwd("C:/code/active learning/R code/ACI_release/")
      source("./distribution_function.R")
      source("./generate_graph.R") #generate random causal graph and simulate the response variable
      source("./distribution_max_main.R") #main program
      
      
      rres = rbind(rres,
                   criterion(adj_mat_truth, adj_mat_pre, adj_mat_pri, goal))
    }
    name = paste(
      "C:/code/active learning/R code/ACI_release/data/v2_ACI_max_thre_1-100_",
      poss,
      "_node_num_",
      node_num,
      "_sample_",
      n,
      "_total_cost_",
      total_cost,
      ".Rdata",
      sep = ""
    )
    save(rres, file = name)
    
  }
}




result = NULL
for (i in c(1:4)) {
  name = paste(
    "C:/code/active learning/R code/ACI_release/data/v2_ACI_max_thre_1-100_0.6_node_num_30_sample_2500_total_cost_",
    i,
    ".Rdata",
    sep = ""
  )
  load(name)
  result = rbind(result, colSums(rres)[1])
}
