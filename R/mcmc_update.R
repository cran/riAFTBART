blocked.mcmc.update <- function(status, y.train, x.train, x.test, cluster.id, cur,aft_model_scale,SA = SA,trt.train,prior_c_function_used,gps) {
  cur <- sample.tree(status, y.train, x.train, x.test, cluster.id, cur,aft_model_scale,SA = SA,trt.train,prior_c_function_used,gps)
  cur <- sample.b(status, y.train, x.train, x.test, cluster.id, cur)
  cur <- sample.tau(status, y.train, x.train, x.test, cluster.id, cur)
  cur <- sample.alpha(status, y.train, x.train, x.test, cluster.id, cur)
  cur <- sample.time.censored(status, y.train, x.train, cluster.id, cur)
  return(cur)
}
