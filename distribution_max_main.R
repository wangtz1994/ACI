#adj_mat: adjacent matrix describing the graph which is updated continuously, [i,j]=1 means i->j
#adj_mat_pri: adjacent matrix describing the Essential graph, [i,j]=1 means i->j
#adj_mat_truth: adjacent matrix describing the causal graph, [i,j]=1 means i->j

To_be_explored = as(Essgraph, "graphNEL") #To_be_explored: the graph which is updated continuously
To_be_explored2 = To_be_explored #To_be_explored2: Essential graph

adj_mat = (as(To_be_explored, "matrix") != 0) * 1
adj_mat_pri = (as(To_be_explored, "matrix") != 0) * 1
adj_mat_truth = ((as(myDAG, "matrix") != 0) * 1)


all_G = cclist(adj_mat)  #divide the chain graph into several chain_component
all_G = filtrate(adj_mat, all_G) # delete the chain components without directed paths to Y
len = length(all_G)
max_length = 10
stopcount = 0
query_number = 0 #record the total number of observed variables (in our paper, query_number equels to the intervention number because only one variable Y is observed under each interveniton)
Inter_number = 0 #record the intervention times
while (len > 0) {
  chain_component = all_G[[1]] # select one chain component to identify
  adj_mat = (as(To_be_explored, "matrix") != 0) * 1
  if (goal %in% chain_component) {
    # if Y is in this chain component
    Chosen_set <-
      as.numeric(which((adj_mat[, goal] == 1 &
                          adj_mat[goal, ] == 1))) # which variables are siblings of Y
    query_node <-
      query_index_in_subgraph(adj_mat, Chosen_set) #select the variable to be intervened
    Direct_set = goal
    all_G = chain_component_update(query_node, chain_component, Direct_set) # learn graph
  } else{
    Directed_graph = adj_mat - ((adj_mat + t(adj_mat)) == 2) * 1 #delete undirected paths
    Chosen_set = possAn(t(Directed_graph), goal, type = "dag") #find which variables are ancestors of Y
    Direct_set = intersect(Chosen_set, chain_component) #find which variables are ancestors of Y in this chain component
    if (length(Direct_set) == 0) {
      #No directed paths is impossible because we have deleted the chain components without directed paths to Y
      print('error')
    }
    adj_subgraph = adj_mat[Direct_set, Direct_set] #get the subgraph comprised of Direct_set
    
    tem_adj_sub = ((adj_subgraph + t(adj_subgraph)) == 2) * 1
    if (is.null(dim(tem_adj_sub))) {
      # there is no directed edges between Direct_set
      tem_adj_sub = t(matrix(tem_adj_sub))
    }
    
    if (length(tem_adj_sub) == 1) {
      #only one variable in Direct_set
      undirect_degree = tem_adj_sub
    } else{
      undirect_degree = colSums(tem_adj_sub)
    }
    degree_query <- max(undirect_degree)
    
    
    # chosen_set is the intervention variable selection set obtained by Fig.2 in main paper
    if (degree_query >= 1) {
      Chosen_set <- Direct_set
      query_node <-
        query_index_in_subgraph(adj_mat, Chosen_set) # select the variable with the maximum siblings
      all_G = chain_component_update(query_node, chain_component, Direct_set)
    } else{
      Chosen_set = c()
      for (index in 1:length(Direct_set)) {
        tem_set = which(tem_adj_mat[Direct_set[index], ] == 1) # which variable is the sibling of the variable in Direct_set
        Chosen_set = union(Chosen_set, tem_set)
      }
      Chosen_set = setdiff(Chosen_set, Direct_set)
      query_node <-
        query_index_in_subgraph(adj_mat, Chosen_set) # select the variable with the maximum siblings
      all_G = chain_component_update(query_node, chain_component, Direct_set)
    }
  }
  len = length(all_G)
  if (len == 0) {
    break
  }
  if (identical(all_G, 'break')) {
    break
  }
  if (identical(chain_component, all_G[[1]])) {
    stopcount = stopcount + 1
  } else{
    stopcount = 0
  }
  if (stopcount > 9) {
    browser()
  }
  # if(len>0){
  #   if(chain_component==all_G[[1]])
  #     break
  # }
}

if (require(Rgraphviz)) {
  plot(To_be_explored)
}

adj_mat_pre = (as(To_be_explored, "matrix") != 0) * 1
pa_goal_truth = possAn(t(adj_mat_truth), goal, type = "dag")
pa_goal_pre = possAn(t(adj_mat_pre), goal, type = "dag")
