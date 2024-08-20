#' RRC method to obtain Cox model coefficients
#'
#'@import spatstat
#'@import survival
#'@import nleqslv
#'@import sandwich
#'@import magic
#'@import boot
#' @param data dataset.
#' @param Iname name of the indicator on the main/validation study (I=1 for validation study sample, I=0 for main study samples).
#' @param Tname name of the observed failure time.
#' @param Dname name of the failure indicator.
#' @returns RRC estimates of the Cox regression coefficients.
cox_rrc.f=function(data,Time_sep,Iname,Tname,Dname,Aname,Asname,Wname) {
  Wf = paste(Wname,collapse="+")
  Outcome.formula.rrc = as.formula(paste("survival::Surv(tstart,",Tname,",",Dname,") ~", Aname, "+", Wf))
  Outcome.formula = as.formula(paste("survival::Surv(",Tname,",",Dname,") ~", Aname, "+", Wf))
  ME.formula = as.formula(paste(Aname, "~", Asname, "+", Wf))

  data.main = data[which(data[,Iname]==0),]
  data.main[,Aname] = data.main[,Asname]
  data.main = data.main[order(data.main[,Tname]),]
  data.va   = data[which(data[,Iname]==1),]
  data.va = data.va[order(data.va[,Tname]),]
  data=rbind(data.main,data.va)
  n=dim(data)[1]
  n1 = dim(data.main)[1]
  n2=n-n1
  #######################################
  # split the datasets based on Time_sep
  #######################################
  index_t_main = index_t_va = c()
  n_rrc = length(Time_sep)
  for (i in (1:n_rrc)) {
    index_t_main[i] = which(data.main[,Tname]>=Time_sep[i])[1]
    index_t_va[i] = which(data.va[,Tname]>=Time_sep[i])[1]
  }
  #index_t_main[i+1] = n1;index_t_va[i+1] = n2;
  ##########################
  # measurement error models
  ##########################
  gamma2 = c()
  for (j in (1:n_rrc)) {
    gamma2.now = lm(ME.formula,data=data.va[index_t_va[j]:n2,])$coefficients
    gamma2 = c(gamma2,gamma2.now)
  }
  # design matrices
  Xe2 = model.matrix(ME.formula,data.va)
  Ae = data.va[,all.vars(ME.formula)[1]]
  ### all parameters in ME models
  #gamma2 = gamma[(gamma_sep[2]+1):gamma_sep[3]]
  n_gamma2 = length(gamma2)/n_rrc
  # estimating equation for gamma
  # type 1, bread; type 2, meat
  ee_gamma = function(gamma,type=1) {
    gamma2 = gamma
    # ee for gamma2
    for (j in (1:n_rrc)) {
      index.now = index_t_va[j]:n2
      gamma2.now = gamma2[((j-1)*n_gamma2+1):(n_gamma2*j)]
      resi=rep(0,n2)
      resi[index.now] = c(Ae - c(Xe2 %*% gamma2.now))[index.now]
      if (j == 1) {
        out5 = Xe2 * resi
      } else {
        out5 = cbind(out5,Xe2 * resi)
      }
    }
    out6=matrix(0,ncol=dim(out5)[2],nrow=n)
    out6[(n1+1):n,] = out5
    out=cbind(out6)
    if (type==1) return(apply(out,2,sum))
    if (type==2) return(out)
  }
  Ug = ee_gamma(gamma=gamma2,type=2)
  Ig = numDeriv::jacobian(ee_gamma, x=gamma2, method="simple", type=1)
  # Outcome model
  split.formula = as.formula(paste("Surv(",Tname,",",Dname,")~."))
  Xt0 = survSplit(split.formula, data=as.data.frame(data.main), cut=Time_sep[-1],episode="timegroup")
  index_t_main=c(index_t_main,n1)
  for (j in (1:n_rrc)) {
    gamma2.now = gamma2[((j-1)*n_gamma2+1):(n_gamma2*j)]
    index.now = which(Xt0$timegroup==j)
    Xt0[index.now,Aname] = c(model.matrix(ME.formula,Xt0[index.now,]) %*% gamma2.now)
  }
  beta = survival::coxph(Outcome.formula.rrc,data=Xt0,ties="breslow")$coefficients
  # now we calculate the sandwich variance
  param = c(beta,gamma2)
  param_sep = c(0,cumsum(c(length(beta),length(gamma2))))
  # joint estimating equation for alpha and beta
  # type 1, bread; type 2, meat
  Time = data.main[,Tname]
  D = data.main[,Dname]
  unique_t<-unique(Time)
  index_t1 = num_t <- c()
  for (i in (1:length(unique_t))) {
    index=which(Time==unique_t[i])
    num_t[i] = length(index)
    index_t1[i] = index[1]
  }
  ee_beta = function(param, type=1) {
    beta = param[(param_sep[1]+1):param_sep[2]]
    gamma2 = param[(param_sep[2]+1):param_sep[3]]
    # ee for outcome model
    Xt0 = list()
    out_cox = matrix(0,ncol=length(beta),nrow = n1)
    for (j in (1:n_rrc)) {
      Xt0[[j]] = data.main
      gamma2.now = gamma2[((j-1)*n_gamma2+1):(n_gamma2*j)]
      Xt0[[j]][,Aname] = c(model.matrix(ME.formula,data.main) %*% gamma2.now)
      Xt0[[j]] = model.matrix(Outcome.formula,Xt0[[j]])[,-1]
      exp_LP = exp(c(Xt0[[j]] %*% beta))
      E_no = rev(cumsum(rev(exp_LP)))
      E_deno = apply(apply(apply(Xt0[[j]] * exp_LP,2,rev),2,cumsum),2,rev)
      E = E_deno/E_no
      # This is for the ties based on Breslow's approximation
      E = matrix(unlist(rep(as.data.frame(t(E[index_t1,])),num_t)),ncol=length(beta),byrow=T)
      out_cox_now = (Xt0[[j]] - E) * (D %*% t(rep(1,dim(Xt0[[j]])[2])))
      if (j<n_rrc) index_now = which(Time>Time_sep[j] & Time<Time_sep[j+1])
      if (j==n_rrc) index_now = which(Time>Time_sep[j])
      out_cox[index_now,] = out_cox_now[index_now,]
    }
    out3=rbind(out_cox,matrix(0,ncol=length(beta),nrow=n2))
    out = out3
    if (type==1) return(apply(out,2,sum))
    if (type==2) return(out)
  }
  Uab = ee_beta(param=param,type=2)
  Iabg = numDeriv::jacobian(ee_beta, x=param, method="simple", type=1)/n
  Iabg_ab = Iabg[,(param_sep[1]+1):param_sep[2]]
  Iabg_g  = Iabg[,(param_sep[2]+1):param_sep[3]]
  Meat = crossprod(Uab + t(Iabg_g %*% solve(Ig) %*% t(Ug)))
  cov_ab = (solve(Iabg_ab) %*% (Meat/n) %*% t(solve(Iabg_ab)))/n

  Point=c(beta)
  SE = sqrt(abs(diag(cov_ab)))
  CI_low = Point - 1.96*SE
  CI_up = Point + 1.96*SE

  res = rbind(Point,SE,CI_low,CI_up)
  rownames(res) = c("point","SE","CI.lower","CI.upper")
  return(res)
}



#' ORC method to obtain the Cox model coefficients
#'@import spatstat
#'@import survival
#'@import nleqslv
#'@import sandwich
#'@import magic
#'@import boot
#' @param data dataset.
#' @param Iname name of the indicator on the main/validation study (I=1 for validation study sample, I=0 for main study samples).
#' @param Tname name of the observed failure time.
#' @param Dname name of the failure indicator.
#' @returns ORC estimates of the Cox regression coefficients.
cox_orc.f=function(data,Iname,Tname,Dname,Aname,Asname,Wname) {
  Wf = paste(Wname,collapse="+")
  Outcome.formula = as.formula(paste("survival::Surv(",Tname,",",Dname,") ~", Aname, "+", Wf))
  ME.formula = as.formula(paste(Aname, "~", Asname, "+", Wf))

  data.main = data[which(data[,Iname]==0),]
  data.main[,Aname] = data.main[,Asname]
  data.main = data.main[order(data.main[,Tname]),]
  data.va   = data[which(data[,Iname]==1),]
  data=rbind(data.main,data.va)
  n=dim(data)[1]
  n1 = dim(data.main)[1]
  n2=n-n1
  ##########################
  # measurement error models
  ##########################
  ## measurement error model 2
  gamma2.lm = lm(ME.formula,data=data.va)
  gamma2.coef = gamma2.lm$coefficient
  #s2gamma2 = mean((gamma2.lm$residuals)^2)
  # design matrices
  #Xe1 = model.matrix(gamma1.lm,data.va)
  Xe2 = model.matrix(gamma2.lm,data.va)
  Ae = data.va[,all.vars(ME.formula)[1]]
  ### all parameters in ME models
  #gamma = c(gamma2.coef)
  gamma2 = gamma2.coef
  # estimating equation for gamma
  # type 1, bread; type 2, meat
  ee_gamma = function(gamma,type=1) {
    gamma2 = gamma
    # ee for gamma2
    resi=Ae - c(Xe2 %*% gamma2)
    out5=Xe2 * resi
    out7 = out5
    out8=matrix(0,ncol=dim(out7)[2],nrow=n)
    out8[(n1+1):n,] = out7
    out=out8
    if (type==1) return(apply(out,2,sum))
    if (type==2) return(out)
  }
  Ug = ee_gamma(gamma=gamma2,type=2)
  Ig = numDeriv::jacobian(ee_gamma, x=gamma2, method="simple", type=1)
  # Outcome model
  Xt0 = as.matrix(data.main)
  Xt0[,Aname] = c(model.matrix(ME.formula,data.main) %*% gamma2)
  Xt0=as.data.frame(Xt0)
  beta = survival::coxph(Outcome.formula,data=Xt0,ties="breslow")$coefficients
  # now we calculate the sandwich variance
  param = c(beta,gamma2)
  param_sep = c(0,cumsum(c(length(beta),length(gamma2))))
  # joint estimating equation for alpha and beta
  # type 1, bread; type 2, meat
  ee_beta = function(param, type=1) {
    beta = param[(param_sep[1]+1):param_sep[2]]
    gamma2 = param[(param_sep[2]+1):param_sep[3]]
    # ee for outcome model
    Xt0 = as.matrix(data.main)
    Xt0[,Aname] = c(model.matrix(ME.formula,data.main) %*% gamma2)
    Xt0=as.data.frame(Xt0)
    Xt = model.matrix(Outcome.formula,Xt0)[,-1]
    Time = data.main[,Tname]
    D = data.main[,Dname]
    unique_t<-unique(Time)
    index_t1 = num_t <- c()
    for (i in (1:length(unique_t))) {
      index=which(Time==unique_t[i])
      num_t[i] = length(index)
      index_t1[i] = index[1]
    }
    exp_LP = exp(c(Xt %*% beta))
    E_no = rev(cumsum(rev(exp_LP)))
    E_deno = apply(apply(apply(Xt * exp_LP,2,rev),2,cumsum),2,rev)
    E = E_deno/E_no
    # This is for the ties based on Breslow's approximation
    E = matrix(unlist(rep(as.data.frame(t(E[index_t1,])),num_t)),ncol=length(beta),byrow=T)
    out3=rbind((Xt - E) * (D %*% t(rep(1,dim(Xt)[2]))),matrix(0,ncol=length(beta),nrow=n2))
    out = out3
    if (type==1) return(apply(out,2,sum))
    if (type==2) return(out)
  }
  Uab = ee_beta(param=param,type=2)
  Iabg = numDeriv::jacobian(ee_beta, x=param, method="simple", type=1)/n
  Iabg_ab = Iabg[,(param_sep[1]+1):param_sep[2]]
  Iabg_g  = Iabg[,(param_sep[2]+1):param_sep[3]]
  Meat = crossprod(Uab + t(Iabg_g %*% solve(Ig) %*% t(Ug)))
  cov_ab = (solve(Iabg_ab) %*% (Meat/n) %*% t(solve(Iabg_ab)))/n


  Point=c(beta)
  SE = sqrt(abs(diag(cov_ab)))
  CI_low = Point - 1.96*SE
  CI_up = Point + 1.96*SE

  res = rbind(Point,SE,CI_low,CI_up)
  rownames(res) = c("point","SE","CI.lower","CI.upper")
  return(res)
}
