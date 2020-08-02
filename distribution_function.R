# The following code is finished by Tian-Zuo Wang
query_distribution <-
  function(query_node, chain_component, Direct_set) {
    #this function learn the causal structure by the interventional data of the response vairable
    #query_node is the intervention variable
    
    
    adj_mat = (as(To_be_explored, "matrix") != 0) * 1
    query_index = which(chain_component == query_node)
    tem_adj = adj_mat[chain_component, chain_component]
    tem_adj_sub = ((tem_adj + t(tem_adj)) == 2) * 1
    index = which(tem_adj_sub[query_index,] == 1)  #find siblings
    index = as.numeric(index)
    
    index = chain_component[index]  #siblings
    
    ind1 = NULL
    n.ind <- length(index)
    for (ww in 1:n.ind) {
      ind1 <- rbind(ind1, c(query_node, index[ww]))
    }
    
    Struc = generate_candidate(n.ind) #1 in, 0 out.
    res_chosen = chosen_potential(query_node, index, Struc, goal)
    potential = res_chosen[[1]]
    group = res_chosen[[2]]
    Backvariable = res_chosen[[3]] #the last element of backvariable is query_node
    
    Est_effect = c()
    # Est_effect=NULL
    intervention_value = 2
    for (i in c(1:length(potential))) {
      parent = Backvariable[[i]]
      Est_effect = c(Est_effect,
                     effect_estimation(parent, query_node, intervention_value))
      # Est_effect=rbind(Est_effect,effect_estimation(parent,query_node,intervention_value))
    }
    print(Est_effect)
    query_result = mean(real_distribution(query_node, intervention_value, goal))
    # query_result=real_distribution(query_node,intervention_value,goal)
    if (Check_cost(total_cost, cost_a1, cost_a2)) {
      query_number2 <<- query_number - 1
      return('break')
    }
    
    
    chosen_index = which.min(abs(Est_effect - query_result))
    
    result = potential[[chosen_index]]
    if (is.vector(group[[chosen_index]])) {
      result = potential[[chosen_index]]
    } else{
      result = rep(-1, length(potential[[chosen_index]]))
      tem = colSums(group[[chosen_index]]) / (dim(group[[chosen_index]])[1])
      result[which(tem == 1)] = 1
      result[which(tem == 0)] = 0
    }
    
    out_index = which(result == 1)
    in_index = which(result == 0)
    
    
    for (i in in_index) {
      # debug_me(ind1[i,2],ind1[i,1],query_node,index)
      To_be_explored <<-
        addBgKnowledge(To_be_explored,
                       ind1[i, 2],
                       ind1[i, 1],
                       verbose = TRUE,
                       checkInput = TRUE)
    }
    for (i in out_index) {
      # debug_me(ind1[i,1],ind1[i,2],query_node,index)
      To_be_explored <<-
        addBgKnowledge(To_be_explored,
                       ind1[i, 1],
                       ind1[i, 2],
                       verbose = TRUE,
                       checkInput = TRUE)
    }
    adj_mat = (as(To_be_explored, "matrix") != 0) * 1
    return(adj_mat)
  }

Check_cost <- function(total_cost, a1, a2) {
  # check whether intervention cost limit has ran out
  if (a1 * Inter_number + a2 * query_number > total_cost) {
    return(1)
  } else{
    return(0)
  }
}



chain_component_update <-
  function(query_node, chain_component, Direct_set) {
    #this function update the graph divide the new chain components
    
    
    new_mat = query_distribution(query_node, chain_component, Direct_set) #update the mat and grpah
    
    #if the cost limit run out
    if (identical(new_mat, 'break')) {
      return('break')
    }
    
    #Next, we will divide tht new chain graph into chain components
    tem_component = new_mat[chain_component, chain_component]
    all_new_G = list()
    adj_mat = (as(To_be_explored, "matrix") != 0) * 1
    if (isValidGraph(tem_component, "dag") == FALSE) {
      all_new_G = cclist(tem_component)
    }
    tem_len = length(all_new_G)
    if (tem_len > 0) {
      for (tt in c(1:tem_len)) {
        all_new_G[[tt]] = as.numeric(all_new_G[[tt]])
        all_new_G[[tt]] = chain_component[all_new_G[[tt]]]
      }
    }
    if (tem_len > 1) {
      all_G[[1]] = NULL
      filtrated_component = filtrate(adj_mat, all_new_G)
      len_of_filt = length(filtrated_component)
      if (len_of_filt > 0) {
        all_G[(length(all_G) + 1):(length(all_G) + len_of_filt)] = filtrated_component
      }
    } else if (tem_len == 0) {
      all_G[[1]] = NULL
    } else if (tem_len == 1) {
      fil = filtrate(adj_mat, all_new_G)
      if (length(fil) == 0) {
        all_G[[1]] = NULL
      } else{
        all_G[1] = fil
      }
    }
    return(all_G)
  }


real_distribution <- function(query_goal, query_value, goal) {
  # this function generate the interventional data
  
  inter_data = modified_rmvnorm.ivent(n = 1000,
                                      myDAG,
                                      target = query_goal,
                                      target.value = query_value)
  
  Inter_number <<- Inter_number + 1
  query_number <<- query_number + 1
  return(inter_data[, goal])
}


effect_estimation <- function(parent, query_goal, query_value) {
  # this function estimates the causal effect given each minimal parental back-door admissible set
  
  if (identical(parent, -1)) {
    return(mean(Data[, goal]))
  }
  if (length(parent) == 1) {
    model_1 = NULL
    Data_1 = cbind(Data[, goal], Data[, query_goal])
  } else{
    Data_1 = cbind(Data[, goal], Data[, parent])
    model_1 = Linear_model(Data_1)
  }
  dataset = Data_1
  if (is.null(model_1)) {
    for (ts in seq(0.1, 5, 0.2)) {
      index = which(dataset[, dim(dataset)[2]] < (query_value + ts) &
                      dataset[, dim(dataset)[2]] > (query_value - ts))
      n = dim(dataset)[1]
      thre = n / 100
      if (length(index) > thre) {
        predictions = dataset[index, 1]
        return(mean(predictions))
        # return(predictions)
      }
    }
  } else{
    n = dim(dataset)[1]
    # est_sample_index=sample(n,size=n/10,replace = TRUE)
    est_sample_index = sample(n, size = n / 2, replace = TRUE)
    est_sample = dataset[est_sample_index, ]
    est_sample[, dim(est_sample)[2]] = query_value
    predictions = predict(model_1, data.frame(est_sample[, -1]))
    return(mean(predictions))
    # return(predictions)
  }
  
}




add_edge <- function(tem_graph,
                     query_goal,
                     index,
                     edge_information) {
  len = length(edge_information)
  for (j in 1:len) {
    if (candidate[j] == 1) {
      tem_graph = addBgKnowledge(tem_graph,
                                 query_goal,
                                 index[j],
                                 verbose = FALSE,
                                 checkInput = TRUE)
    } else{
      tem_graph = addBgKnowledge(tem_graph,
                                 index[j],
                                 query_goal,
                                 verbose = FALSE,
                                 checkInput = TRUE)
    }
  }
  mat = (as(tem_graph, "matrix") != 0) * 1
  return()
}


chosen_potential <- function(query_goal, index, potential, goal) {
  # this function records the minimal parental back-door admissible set
  
  #if there exists a backvariable==-1, it is located at the last one
  #Backvariable records the minimal parental back-door admissible set
  #tem records the orientation of each possible causal structure
  
  
  candidate <- NULL
  group <- NULL
  Backvariable <- NULL
  for (i in 1:length(potential)) {
    tem_graph = To_be_explored
    tem = potential[[i]]
    len = length(tem)
    activate = 0
    for (j in 1:len) {
      # at first check whether such orientation is legal (no cycles, no new v-structures) or not
      if (tem[j] == 1) {
        tem_graph = addBgKnowledge(tem_graph,
                                   query_goal,
                                   index[j],
                                   verbose = FALSE,
                                   checkInput = TRUE)
        if (is.null(tem_graph)) {
          activate = 1
          break
        }
      } else{
        tem_graph = addBgKnowledge(tem_graph,
                                   index[j],
                                   query_goal,
                                   verbose = FALSE,
                                   checkInput = TRUE)
        if (is.null(tem_graph)) {
          activate = 1
          break
        }
      }
    }
    if (activate == 1) {
      # such orientation is not legal
      next
    }
    
    mat = (as(tem_graph, "matrix") != 0) * 1
    mat = Construct_DAG(mat)
    De = possDe(t(mat), query_node, type = "cpdag")
    if (!goal %in% De) {
      # the response variable is not the ancestor of intervention variable
      Back = -1
    } else{
      # find the minimal parental back-door admissible set
      
      new_DAG = as(mat, "graphNEL")
      Back = backdoor(mat, query_goal, goal, type = "cpdag")
      
      qq = 1
      end_node = as.character(goal)
      
      while (length(Back) > 0 & qq <= length(Back)) {
        start_node = as.character(Back[qq])
        S = Back[-qq]
        S = as.character(c(S, query_node))
        if (dsep(start_node, end_node, S, new_DAG) == TRUE) {
          Back = Back[-qq]
          next
        } else{
          qq = qq + 1
        }
      }
      if (length(Back) == 0) {
        Back = NULL
      }
    }
    if (!identical(Back, -1)) {
      # add the intervention variable index behind the minimal parental back-door admissible set for the convenience to calculate the expectation of causal effect
      Back = c(Back, query_node)
    }
    
    
    tlen = length(Backvariable)
    if (tlen == 0) {
      Backvariable = c(Backvariable, list(Back))
      candidate = c(candidate, list(tem))
      group = c(group, list(tem))
    } else{
      activate2 = 0
      for (ss in 1:tlen) {
        if (identical(Back, Backvariable[[ss]])) {
          # if this MPS has been added to the MPS set, then add this orientation (with this MPS) to candidate
          if (sum(tem) < sum(candidate[[ss]])) {
            candidate[[ss]] = tem
          }
          group[[ss]] = rbind(group[[ss]], tem)
          activate2 = 1
          break
        }
      }
      
      #if there is no back-door variables same as previous
      if (activate2 == 0) {
        Backvariable = c(Backvariable, list(Back))
        candidate = c(candidate, list(tem))
        group = c(group, list(tem))
      }
    }
  }
  tlen = length(Backvariable)
  for (ss in 1:tlen) {
    if (identical(-1, Backvariable[[ss]])) {
      tem_candidate = candidate[[ss]]
      tem_group = group[[ss]]
      tem_Backvariable = Backvariable[[ss]]
      candidate[[ss]] = NULL
      group[[ss]] = NULL
      Backvariable[[ss]] = NULL
      candidate[[tlen]] = tem_candidate
      group[[tlen]] = tem_group
      Backvariable[[tlen]] = tem_Backvariable
      break
    }
  }
  return(list(candidate, group, Backvariable))
}


Construct_DAG <- function(mat) {
  # randomly generate a new DAG based on a PDAG on the premise that generates no new cycles or v-structures
  DAG_adj = mat
  tem_graph = as(DAG_adj, "graphNEL")
  while (isValidGraph(DAG_adj, "dag") == FALSE) {
    aa = which(DAG_adj + t(DAG_adj) == 2, arr.ind = TRUE)[1, ]
    # DAG_adj[aa]=0
    # tem_graph=as(DAG_adj,"graphNEL")
    tem_graph = addBgKnowledge(tem_graph, aa[1], aa[2])
    DAG_adj = as(tem_graph, "matrix")
  }
  return(DAG_adj)
}



# debug_me<-function(a,b,query_node,index){
#   if(adj_mat_truth[a,b]==0){
#     lindex=sort(c(query_node,index))
#     # print(which(lindex==query_node))
#     plot(as(adj_mat_truth, "graphNEL"))
#     plot(as(adj_mat, "graphNEL"))
#     print(c(a,b,goal))
#     tem_mat_truth=adj_mat_truth[lindex,lindex]
#   }
# }



Linear_model <- function(dataset) {
  #this function trains a linear model to estimate the expectation E[Y|X,Z]
  
  Values = dataset[, -1]
  # a1=dataset[,2]
  # a2=dataset[,3]
  Targets = dataset[, 1]
  total_data = data.frame(dataset)
  len = dim(total_data)[2]
  fm = as.formula(paste(
    colnames(total_data)[1],
    "~",
    paste(colnames(total_data)[-1], collapse = "+"),
    sep = ""
  ))
  model_linear = lm(fm, total_data)
  # test_data=data.frame(Values[1:100,])
  # predict(model_linear,test_data)
  return(model_linear)
}


estimation_est <- function(model, dataset, query_value) {
  #this function estimates the causal effect given a specific intervention value
  if (is.null(model)) {
    for (ts in c(seq(0.1, 5, 0.2))) {
      index = which(dataset[, dim(dataset)[2]] < (query_value + ts) &
                      dataset[, dim(dataset)[2]] > (query_value - ts))
      if (length(index) > 10) {
        predictions = dataset[index, 1]
        return(predictions)
      }
    }
  } else{
    n = dim(dataset)[1]
    est_sample_index = sample(n, size = n / 10, replace = TRUE)
    est_sample = dataset[est_sample_index, ]
    est_sample[, dim(est_sample)[2]] = query_value
    predictions = predict(model, data.frame(est_sample[, -1]))
    return(predictions)
  }
}

criterion <- function(adj_mat_truth,
                      adj_mat_pre,
                      adj_mat_pri,
                      goal) {
  # this funciton return the accuracy of estimated ancestor causal graph and real ancestor causal graph
  pa_goal_truth = possAn(t(adj_mat_truth), goal, type = "dag")
  pa_goal_pre = possAn(t(adj_mat_pre), goal, type = "dag")
  
  # #past
  # tem=adj_mat_pri[pa_goal_truth,pa_goal_truth]
  # ind<-which((tem==t(tem)& tem==1),arr.ind=TRUE)
  
  ind <-
    which((adj_mat_pri == t(adj_mat_pri) &
             adj_mat_pri == 1), arr.ind = TRUE)
  is.legal = apply(ind, 1, function(edge) {
    if (edge[1] %in% pa_goal_truth || edge[2] %in% pa_goal_truth)
      1
    else
      0
  })
  ind = ind[which(is.legal == 1), ]
  
  
  edge_to_be_pre <- ind[(ind[, 1] < ind[, 2]), ]
  edge_to_be_pre <- matrix(edge_to_be_pre, ncol = 2)
  edge_num_to_be_pre <-
    length(edge_to_be_pre) / 2 #number of edge to be discovered
  true_edge_found = 0
  unfound = 0
  wrong_found = 0
  # i=i+1
  # true_edge_found_old=0
  if (edge_num_to_be_pre > 0) {
    for (i in c(1:edge_num_to_be_pre)) {
      tt = edge_to_be_pre[i, ]
      # if(adj_mat_truth[tt[1],tt[2]]==0)
      
      if ((adj_mat_pre[tt[1], tt[2]] == adj_mat_truth[tt[1], tt[2]]) &&
          (adj_mat_pre[tt[2], tt[1]] == adj_mat_truth[tt[2], tt[1]])) {
        true_edge_found = true_edge_found + 1
      }
      if ((adj_mat_pre[tt[1], tt[2]] == 1) &&
          (adj_mat_pre[tt[2], tt[1]] == 1)) {
        unfound = unfound + 1
      }
      # if(((adj_mat_pre[tt[1],tt[2]]==0))&&((adj_mat_pre[tt[2],tt[1]]==1)||(adj_mat_pre[tt[2],tt[1]]==0))){
      if ((adj_mat_pre[tt[2], tt[1]] == adj_mat_truth[tt[1], tt[2]]) &&
          (adj_mat_pre[tt[1], tt[2]] == adj_mat_truth[tt[2], tt[1]])) {
        wrong_found = wrong_found + 1
      }
    }
  }
  return(c(true_edge_found, unfound, wrong_found, edge_num_to_be_pre))
}




filtrate <- function(mat, component_index) {
  # this function delete the chain components without directed paths to Y
  return_index = list()
  tem_len = length(component_index)
  if (tem_len == 0) {
    return(NULL)
  }
  for (i in 1:tem_len) {
    tem = as.numeric(component_index[[i]])
    if (judge_path(tem, goal) == TRUE) {
      #judge whether there is a directed path from tem to goal
      return_index = c(return_index, list(tem))
    }
  }
  return(return_index)
}

judge_path = function(X, Y) {
  # this function judges whether there is a directed path from set X to Y
  if (length(intersect(X, Y)) > 0) {
    # Y\in X
    return(TRUE)
  }
  tem_num = length(X)
  adj_mat = (as(To_be_explored, "matrix") != 0) * 1
  tem_adj_mat = adj_mat - ((adj_mat + t(adj_mat)) == 2) * 1
  
  
  chosen = possAn(t(tem_adj_mat), Y, type = "dag")
  chosen = chosen[-which(chosen == Y)]
  
  if (length(intersect(chosen, X)) > 0) {
    return(TRUE)
  } else{
    return(FALSE)
  }
}


query_index_in_subgraph <- function(mat, chosen_index) {
  # return the variable to be intervened in the intervention variable selection set
  adj_subgraph = mat[chosen_index, chosen_index]
  tem_adj_sub = ((adj_subgraph + t(adj_subgraph)) == 2) * 1
  if (is.null(dim(tem_adj_sub))) {
    tem_adj_sub = matrix(tem_adj_sub)
  }
  
  
  undirect_degree = colSums(tem_adj_sub)
  
  #----- if select the intervention variable with the maximum siblings-----
  return(chosen_index[which.max(undirect_degree)])
  
  #----- if randomly select the intervention variable-----
  # if(length(chosen_index)==1){
  #   return(chosen_index)
  # }else{
  #   return(sample(chosen_index,1))
  # }
}


generate_candidate <- function(n) {
  # this function generate all possible orientations of the edges between the intervention variable and its siblings
  # 0 and 1 represent two possible orientations of one edge
  
  if (n == 1) {
    return(list(0, 1))
  }
  candidate = NULL
  for (i in 1:(2 ^ n)) {
    tem_array = c()
    count = 1
    s = i
    while (count <= n) {
      if (s > 2 ^ (n - count)) {
        tem_array[count] = 1
        s = s - 2 ^ (n - count)
      } else{
        tem_array[count] = 0
      }
      count = count + 1
    }
    candidate[[i]] = tem_array
  }
  return(candidate)
}

#-----------------------------------------------------------------#
# this part aims to generate the data with the graph from randomDAG


modified_get_adj <-
  function (object, target = integer(0))
    #returned result is like the new relation matrix
  {
    p <-  length(object@nodes)
    target <- as.integer(sort(target))
    result <- matrix(0, p, p)
    cumul = 0
    for (i in 1:p) {
      if (length(object@edgeL[[i]]$edges) > 0) {
        for (k in 1:length(object@edgeL[[i]]$edges)) {
          ss = object@edgeL[[i]]$edges[k]
          result[i, ss] <- object@edgeData@data[[cumul + k]]$weight
        }
        cumul = cumul + length(object@edgeL[[i]]$edges)
      }
    }
    for (i in 1:p) {
      if (as.integer(i) %in% target)
        result[i, ]
      
    }
    rownames(result) <- object@nodes
    colnames(result) <- object@nodes
    return(result)
  }
modified_rmvnorm.ivent <-
  function (n,
            object = myDAG2,
            target = integer(0),
            target.value = numeric(0))
  {
    p <- length(object@nodes)
    stopifnot(length(target) == 0 ||
                (1 <= min(target) && max(target) <=
                   p))
    stopifnot((
      is.vector(target.value) && length(target.value) ==
        length(target)
    ) || (is.matrix(target.value) && dim(target.value) ==
            c(n, length(target))))
    sigma <- 1
    mu <- 0
    Y <- matrix(rnorm(n * p, 0, 1), nrow = p, ncol = n)
    Y[target,] <- target.value
    A <- -t(modified_get_adj(object, target))
    diag(A) <- 1
    t(solve(A, Y))
  }




cclist <- function(gm) {
  # Purpose: Find all chain components from a chain graph
  # This function is finished by He & Geng 2008
  
  ind1 <- which((gm == t(gm) & gm == 1), arr.ind = TRUE)
  ind1 <- ind1[(ind1[, 1] < ind1[, 2]), ]
  n.ind <- length(ind1) / 2
  cclist <- vector("list", n.ind)
  if (n.ind == 1)
    cclist[[1]] <- ind1
  if (n.ind > 1) {
    for (i in 1:n.ind) {
      cclist[[i]] <- ind1[i, ]
    }
    
    i = 1
    j = 2
    ccfind <- FALSE
    repeat {
      if (length(intersect(cclist[[i]], cclist[[j]]))) {
        cclist[[i]] <- union(cclist[[i]], cclist[[j]])
        cclist = cclist[-j]
        ccfind <- TRUE
      }
      else {
        j = j + 1
      }#if (length
      if (j > length(cclist)) {
        if (!ccfind) {
          i = i + 1
          j = i + 1
        } else {
          j = i + 1
          ccfind <- FALSE
        }
      }
      if (i >= length(cclist))
        break
    }#repeat
    
    for (i in 1:length(cclist)) {
      cclist[[i]] <- sort(cclist[[i]])
    }
  }#if (n.
  cclist
}
