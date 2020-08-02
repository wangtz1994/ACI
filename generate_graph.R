myDAG <- randomDAG(node_num, poss)
Data <- rmvDAG(n = n, myDAG)
Essgraph = dag2essgraph(myDAG, targets = list(integer(0)))
if (require(Rgraphviz)) {
  plot(Essgraph)
}
adj_mat_Ess = (as(Essgraph, "matrix") != 0) * 1
goal = which.max(colSums(adj_mat_Ess))
