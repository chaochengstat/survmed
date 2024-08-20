#' Unadjusted method that mistakenly treats the mismeasured exposure as the true exposure
#'
#' @param Outcome.formula formula for the Cox outcome model.
#' @param Mediator.formula formula for the linear mediator model.
#' @param data dataset.
#' @param Iname name of the indicator on the main/validation study (I=1 for validation study sample, I=0 for main study samples).
#' @param Tname name of the observed failure time.
#' @param Dname name of the failure indicator.
#' @returns Unadjusted estimates of the regression coefficients.
unadj.f=function(Outcome.formula,Mediator.formula,data,Iname,Tname,Dname,is.boot=0) {
  data[,Aname] = data[,Asname]
  data.main = data[data[,Iname]==0,]
  data.main = data.main[order(data.main[,Tname]),]
  data=rbind(data.main,data[data[,Iname]==1,])
  n=dim(data)[1]
  n1 = dim(data.main)[1]
  n2=n-n1
  # mediator model
  alpha.lm = lm(Mediator.formula,data)
  alpha = alpha.lm$coefficient
  s2alpha = mean((alpha.lm$residuals)^2)
  Xm = model.matrix(Mediator.formula,data)
  M = data[,all.vars(Mediator.formula)[1]]
  # estimating equation for alpha
  # type 1, bread; type 2, meat
  ee_alpha = function(alpha,type=1) {
    alpha.coef = alpha[-length(alpha)]
    s2alpha= alpha[length(alpha)]
    resi=M - c(Xm %*% alpha.coef)
    out1=Xm * resi
    out2= s2alpha - resi^2
    out = cbind(out1,out2)
    if (type==1) return(apply(out,2,sum))
    if (type==2) return(out)
  }
  alpha = c(alpha,s2alpha)
  Ua = ee_alpha(alpha=alpha,type=2)
  Ia = numDeriv::jacobian(ee_alpha, x=alpha, method="simple", type=1)/n
  cov_a = (solve(Ia) %*% (crossprod(Ua)/n) %*% t(solve(Ia)))/n
  #outcome model
  beta = survival::coxph(Outcome.formula,data=data.main,ties="breslow")$coefficients
  Xt = model.matrix(Outcome.formula,data.main)[,-1]
  Time = data.main[,Tname]
  D = data.main[,Dname]
  # create some quantities for all unique time (for ties)
  unique_t<-unique(Time)
  index_t1 = num_t <- c()
  for (i in (1:length(unique_t))) {
    index=which(Time==unique_t[i])
    num_t[i] = length(index)
    index_t1[i] = index[1]
  }
  # estimating equation for beta
  ee_beta = function(beta,type=1) {
    exp_LP = exp(c(Xt %*% beta))
    E_no = rev(cumsum(rev(exp_LP)))
    E_deno = apply(apply(apply(Xt * exp_LP,2,rev),2,cumsum),2,rev)
    E = E_deno/E_no
    # This is for the ties based on Breslow's approximation
    E = matrix(unlist(rep(as.data.frame(t(E[index_t1,])),num_t)),ncol=length(beta),byrow=T)
    out = (Xt - E) * (D %*% t(rep(1,dim(Xt)[2])))
    if (type==1) return(apply(out,2,sum))
    if (type==2) return(out)
  }
  beta=nleqslv::nleqslv(x=beta,ee_beta)$x
  if (is.boot==0) {
    Ub = ee_beta(beta=beta,type=2)
    Ib = numDeriv::jacobian(ee_beta, x=beta, method="simple", type=1)/n1
    cov_b = (solve(Ib) %*% (crossprod(Ub)/n1) %*% t(solve(Ib)))/n1
  } else {
    cov_b=diag(1,length(beta))
  }
  return(list(Point=c(alpha,beta),COV=magic::adiag(cov_a,cov_b)))
}

#' ORC1 method to obtain the mediator and outcome coefficients
#' @import spatstat
#' @import survival
#' @import nleqslv
#' @import sandwich
#' @import magic
#' @import boot
#' @param Outcome.formula formula for the Cox outcome model.
#' @param Mediator.formula formula for the linear mediator model.
#' @param data dataset.
#' @param Iname name of the indicator on the main/validation study (I=1 for validation study sample, I=0 for main study samples).
#' @param Tname name of the observed failure time.
#' @param Dname name of the failure indicator.
#' @returns ORC1 estimates of the regression coefficients.
orc1.f=function(Outcome.formula,Mediator.formula,ME1.formula,
                data,Iname,Tname,Dname,Aname,Asname,Mname,Wname,is.boot=0) {
  data.main = data[data[,Iname]==0,]
  data.main[,Aname] = data.main[,Asname]
  data.main = data.main[order(data.main[,Tname]),]
  data.va   = data[data[,Iname]==1,]
  data=rbind(data.main,data.va)
  n=dim(data)[1]
  n1 = dim(data.main)[1]
  n2=n-n1
  # measurement error model
  gamma.lm = lm(ME1.formula,data=data.va)
  gamma.coef = gamma.lm$coefficient
  s2gamma = mean((gamma.lm$residuals)^2)
  Xe = model.matrix(gamma.lm,data.va)
  Ae = data.va[,all.vars(ME1.formula)[1]]
  # estimating equation for alpha
  # type 1, bread; type 2, meat
  ee_gamma = function(gamma,type=1) {
    gamma.coef = gamma[-length(gamma)]
    s2gamma= gamma[length(gamma)]
    resi=Ae - c(Xe %*% gamma.coef)
    out1=Xe * resi
    out2= s2gamma - resi^2
    out3 = cbind(out1,out2)
    out=matrix(0,ncol=dim(out3)[2],nrow=n)
    out[(n1+1):n,] = out3
    if (type==1) return(apply(out,2,sum))
    if (type==2) return(out)
  }
  gamma = c(gamma.coef,s2gamma)
  Ug = ee_gamma(gamma=gamma,type=2)
  Ig = numDeriv::jacobian(ee_gamma, x=gamma, method="simple", type=1)
  #cov_g = (solve(Ig) %*% (crossprod(Ug)/n2) %*% t(solve(Ig)))/n2
  # mediator model
  Xm = model.matrix(Mediator.formula,data)
  Xm[1:n1,Aname] = c(model.matrix(ME1.formula,data.main) %*% gamma[-length(gamma)])
  M = data[,all.vars(Mediator.formula)[1]]
  alpha.lm = lm(M~Xm-1)
  alpha.coef = alpha.lm$coefficients
  s2alpha = (sum((alpha.lm$residuals)^2) - n1*(alpha.coef[2]^2)*gamma[length(gamma)])/n
  alpha = c(alpha.coef,s2alpha)
  # Outcome model
  Xt0 = as.matrix(data.main)
  Xt0[,Aname] = (alpha[2]*(Xt0[,Mname]-alpha[1]- c(Xt0[,Wname,drop=F] %*% alpha[-c(1:2,length(alpha))]))*gamma[length(gamma)]
                 +c(model.matrix(ME1.formula,as.data.frame(Xt0)) %*% gamma[-length(gamma)])*alpha[length(alpha)])/(
                   alpha[2]^2*gamma[length(gamma)]+alpha[length(alpha)])
  Xt0=as.data.frame(Xt0)
  beta = survival::coxph(Outcome.formula,data=Xt0,ties="breslow")$coefficients
  # now we calculate the sandwich variance
  param = c(alpha,beta,gamma)
  param_sep = c(0,cumsum(c(length(alpha),length(beta),length(gamma))))
  # joint estimating equation for alpha and beta
  # type 1, bread; type 2, meat
  ee_alpha_beta = function(param, type=1) {
    alpha = param[(param_sep[1]+1):param_sep[2]]
    beta = param[(param_sep[2]+1):param_sep[3]]
    gamma = param[(param_sep[3]+1):param_sep[4]]
    # ee for mediator model
    Xm[1:n1,Aname] = c(model.matrix(ME1.formula,data.main) %*% gamma[-length(gamma)])
    M = data[,all.vars(Mediator.formula)[1]]
    alpha.coef = alpha[-length(alpha)]
    s2alpha= alpha[length(alpha)]
    resi=M - c(Xm %*% alpha.coef)
    out1=Xm * resi
    offset.term = c(rep((alpha.coef[2]^2)*gamma[length(gamma)],n1),rep(0,n2))
    out2= s2alpha+offset.term - resi^2
    out = cbind(out1,out2)
    # ee for outcome model
    Xt0 = as.matrix(data.main)
    Xt0[,Aname] = (alpha[2]*(Xt0[,Mname]-alpha[1]- c(Xt0[,Wname,drop=F] %*% alpha[-c(1:2,length(alpha))]))*gamma[length(gamma)]
                   +c(model.matrix(ME1.formula,as.data.frame(Xt0)) %*% gamma[-length(gamma)])*alpha[length(alpha)])/(
                     alpha[2]^2*gamma[length(gamma)]+alpha[length(alpha)])
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
    out = cbind(out,out3)
    if (type==1) return(apply(out,2,sum))
    if (type==2) return(out)
  }
  if (is.boot==0) {
    Uab = ee_alpha_beta(param=param,type=2)
    Iabg = numDeriv::jacobian(ee_alpha_beta, x=param, method="simple", type=1)/n
    Iabg_ab = Iabg[,(param_sep[1]+1):param_sep[3]]
    Iabg_g  = Iabg[,(param_sep[3]+1):param_sep[4]]
    Meat = crossprod(Uab + t(Iabg_g %*% solve(Ig) %*% t(Ug)))
    cov_ab = (solve(Iabg_ab) %*% (Meat/n) %*% t(solve(Iabg_ab)))/n
  } else {
    cov_ab = diag(1,length(c(alpha,beta)))
  }
  return(list(Point=c(alpha,beta),COV=cov_ab))
}

#' ORC2 method to obtain the mediator and outcome coefficients
#'@import spatstat
#'@import survival
#'@import nleqslv
#'@import sandwich
#'@import magic
#'@import boot
#' @param Outcome.formula formula for the Cox outcome model.
#' @param Mediator.formula formula for the linear mediator model.
#' @param data dataset.
#' @param Iname name of the indicator on the main/validation study (I=1 for validation study sample, I=0 for main study samples).
#' @param Tname name of the observed failure time.
#' @param Dname name of the failure indicator.
#' @returns ORC2 estimates of the regression coefficients.
orc2.f=function(Outcome.formula,Mediator.formula,ME1.formula,ME2.formula,data,
                Iname,Tname,Dname,Aname,Asname,Mname,Wname,is.boot=0) {
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
  ## measurement error model 1
  gamma1.lm = lm(ME1.formula,data=data.va)
  gamma1.coef = gamma1.lm$coefficient
  s2gamma1 = mean((gamma1.lm$residuals)^2)
  ## measurement error model 2
  gamma2.lm = lm(ME2.formula,data=data.va)
  gamma2.coef = gamma2.lm$coefficient
  s2gamma2 = mean((gamma2.lm$residuals)^2)
  # design matrices
  Xe1 = model.matrix(gamma1.lm,data.va)
  Xe2 = model.matrix(gamma2.lm,data.va)
  Ae = data.va[,all.vars(ME1.formula)[1]]
  ### all parameters in ME models
  gamma = c(gamma1.coef,s2gamma1,gamma2.coef,s2gamma2)
  gamma_sep = c(0,cumsum(c(length(gamma1.coef)+1,length(gamma2.coef)+1)))
  gamma1 = gamma[(gamma_sep[1]+1):gamma_sep[2]]
  gamma2 = gamma[(gamma_sep[2]+1):gamma_sep[3]]
  # estimating equation for gamma
  # type 1, bread; type 2, meat
  ee_gamma = function(gamma,type=1) {
    gamma1 = gamma[(gamma_sep[1]+1):gamma_sep[2]]
    gamma2 = gamma[(gamma_sep[2]+1):gamma_sep[3]]
    gamma1.coef = gamma1[-length(gamma1)]
    s2gamma1= gamma1[length(gamma1)]
    gamma2.coef = gamma2[-length(gamma2)]
    s2gamma2= gamma2[length(gamma2)]
    # ee for gamma1
    resi=Ae - c(Xe1 %*% gamma1.coef)
    out1=Xe1 * resi
    out2= s2gamma1 - resi^2
    out3 = cbind(out1,out2)
    out4=matrix(0,ncol=dim(out3)[2],nrow=n)
    out4[(n1+1):n,] = out3
    # ee for gamma2
    resi=Ae - c(Xe2 %*% gamma2.coef)
    out5=Xe2 * resi
    out6= s2gamma2 - resi^2
    out7 = cbind(out5,out6)
    out8=matrix(0,ncol=dim(out7)[2],nrow=n)
    out8[(n1+1):n,] = out7
    out=cbind(out4,out8)
    if (type==1) return(apply(out,2,sum))
    if (type==2) return(out)
  }
  Ug = ee_gamma(gamma=gamma,type=2)
  Ig = numDeriv::jacobian(ee_gamma, x=gamma, method="simple", type=1)
  # mediator model
  Xm = model.matrix(Mediator.formula,data)
  Xm[1:n1,Aname] = c(model.matrix(ME1.formula,data.main) %*% gamma1[-length(gamma1)])
  M = data[,all.vars(Mediator.formula)[1]]
  alpha.lm = lm(M~Xm-1)
  alpha.coef = alpha.lm$coefficients
  s2alpha = (sum((alpha.lm$residuals)^2) - n1*(alpha.coef[2]^2)*gamma1[length(gamma1)])/n
  alpha = c(alpha.coef,s2alpha)
  # Outcome model
  Xt0 = as.matrix(data.main)
  Xt0[,Aname] = c(model.matrix(ME2.formula,data.main) %*% gamma2[-length(gamma2)])
  Xt0=as.data.frame(Xt0)
  beta = survival::coxph(Outcome.formula,data=Xt0,ties="breslow")$coefficients
  # now we calculate the sandwich variance
  param = c(alpha,beta,gamma)
  param_sep = c(0,cumsum(c(length(alpha),length(beta),length(gamma1),length(gamma2))))
  # joint estimating equation for alpha and beta
  # type 1, bread; type 2, meat
  ee_alpha_beta = function(param, type=1) {
    alpha = param[(param_sep[1]+1):param_sep[2]]
    beta = param[(param_sep[2]+1):param_sep[3]]
    gamma1 = param[(param_sep[3]+1):param_sep[4]]
    gamma2 = param[(param_sep[4]+1):param_sep[5]]
    # ee for mediator model
    Xm[1:n1,Aname] = c(model.matrix(ME1.formula,data.main) %*% gamma1[-length(gamma1)])
    M = data[,all.vars(Mediator.formula)[1]]
    alpha.coef = alpha[-length(alpha)]
    s2alpha= alpha[length(alpha)]
    resi=M - c(Xm %*% alpha.coef)
    out1=Xm * resi
    offset.term = c(rep((alpha.coef[2]^2)*gamma1[length(gamma1)],n1),rep(0,n2))
    out2= s2alpha+offset.term - resi^2
    out = cbind(out1,out2)
    # ee for outcome model
    Xt0 = as.matrix(data.main)
    Xt0[,Aname] = c(model.matrix(ME2.formula,data.main) %*% gamma2[-length(gamma2)])
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
    out = cbind(out,out3)
    if (type==1) return(apply(out,2,sum))
    if (type==2) return(out)
  }
  if (is.boot==0) {
    Uab = ee_alpha_beta(param=param,type=2)
    Iabg = numDeriv::jacobian(ee_alpha_beta, x=param, method="simple", type=1)/n
    Iabg_ab = Iabg[,(param_sep[1]+1):param_sep[3]]
    Iabg_g  = Iabg[,(param_sep[3]+1):param_sep[5]]
    Meat = crossprod(Uab + t(Iabg_g %*% solve(Ig) %*% t(Ug)))
    cov_ab = (solve(Iabg_ab) %*% (Meat/n) %*% t(solve(Iabg_ab)))/n
  } else {
    cov_ab = diag(1,length(c(alpha,beta)))
  }
  return(list(Point=c(alpha,beta),COV=cov_ab))
}


#' RRC method to obtain the mediator and outcome coefficients
#'
#'@import spatstat
#'@import survival
#'@import nleqslv
#'@import sandwich
#'@import magic
#'@import boot
#' @param Outcome.formula formula for the Cox outcome model.
#' @param Mediator.formula formula for the linear mediator model.
#' @param data dataset.
#' @param Iname name of the indicator on the main/validation study (I=1 for validation study sample, I=0 for main study samples).
#' @param Tname name of the observed failure time.
#' @param Dname name of the failure indicator.
#' @returns RRC estimates of the regression coefficients.
rrc.f=function(Outcome.formula,Outcome.formula.rrc,Mediator.formula,ME1.formula,ME2.formula,data,Time_sep,
               Iname,Tname,Dname,Aname,Asname,Mname,Wname,is.boot=0) {
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
  ## measurement error model 1
  gamma1.lm = lm(ME1.formula,data=data.va)
  gamma1.coef = gamma1.lm$coefficient
  s2gamma1 = mean((gamma1.lm$residuals)^2)
  ## measurement error model 2
  gamma2 = c()
  for (j in (1:n_rrc)) {
    gamma2.now = lm(ME2.formula,data=data.va[index_t_va[j]:n2,])$coefficients
    gamma2 = c(gamma2,gamma2.now)
  }
  # design matrices
  Xe1 = model.matrix(ME1.formula,data.va)
  Xe2 = model.matrix(ME2.formula,data.va)
  Ae = data.va[,all.vars(ME1.formula)[1]]
  ### all parameters in ME models
  gamma = c(gamma1.coef,s2gamma1,gamma2)
  gamma_sep = c(0,cumsum(c(length(gamma1.coef)+1,length(gamma2))))
  gamma1 = gamma[(gamma_sep[1]+1):gamma_sep[2]]
  gamma2 = gamma[(gamma_sep[2]+1):gamma_sep[3]]
  n_gamma2 = length(gamma2)/n_rrc
  # estimating equation for gamma
  # type 1, bread; type 2, meat
  ee_gamma = function(gamma,type=1) {
    gamma1 = gamma[(gamma_sep[1]+1):gamma_sep[2]]
    gamma2 = gamma[(gamma_sep[2]+1):gamma_sep[3]]
    gamma1.coef = gamma1[-length(gamma1)]
    s2gamma1= gamma1[length(gamma1)]
    # ee for gamma1
    resi=Ae - c(Xe1 %*% gamma1.coef)
    out1=Xe1 * resi
    out2= s2gamma1 - resi^2
    out3 = cbind(out1,out2)
    out4=matrix(0,ncol=dim(out3)[2],nrow=n)
    out4[(n1+1):n,] = out3
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
    out=cbind(out4,out6)
    if (type==1) return(apply(out,2,sum))
    if (type==2) return(out)
  }
  Ug = ee_gamma(gamma=gamma,type=2)
  Ig = numDeriv::jacobian(ee_gamma, x=gamma, method="simple", type=1)
  # mediator model
  Xm = model.matrix(Mediator.formula,data)
  Xm[1:n1,Aname] = c(model.matrix(ME1.formula,data.main) %*% gamma1[-length(gamma1)])
  M = data[,all.vars(Mediator.formula)[1]]
  alpha.lm = lm(M~Xm-1)
  alpha.coef = alpha.lm$coefficients
  s2alpha = (sum((alpha.lm$residuals)^2) - n1*(alpha.coef[2]^2)*gamma1[length(gamma1)])/n
  alpha = c(alpha.coef,s2alpha)
  # Outcome model
  split.formula = as.formula(paste("Surv(",Tname,",",Dname,")~."))
  Xt0 = survSplit(split.formula, data=as.data.frame(data.main), cut=Time_sep[-1],episode="timegroup")
  index_t_main=c(index_t_main,n1)
  for (j in (1:n_rrc)) {
    gamma2.now = gamma2[((j-1)*n_gamma2+1):(n_gamma2*j)]
    index.now = which(Xt0$timegroup==j)
    Xt0[index.now,Aname] = c(model.matrix(ME2.formula,Xt0[index.now,]) %*% gamma2.now)
  }
  beta = survival::coxph(Outcome.formula.rrc,data=Xt0,ties="breslow")$coefficients
  # now we calculate the sandwich variance
  param = c(alpha,beta,gamma)
  param_sep = c(0,cumsum(c(length(alpha),length(beta),length(gamma1),length(gamma2))))
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
  ee_alpha_beta = function(param, type=1) {
    alpha = param[(param_sep[1]+1):param_sep[2]]
    beta = param[(param_sep[2]+1):param_sep[3]]
    gamma1 = param[(param_sep[3]+1):param_sep[4]]
    gamma2 = param[(param_sep[4]+1):param_sep[5]]
    # ee for mediator model
    Xm[1:n1,Aname] = c(model.matrix(ME1.formula,data.main) %*% gamma1[-length(gamma1)])
    M = data[,all.vars(Mediator.formula)[1]]
    alpha.coef = alpha[-length(alpha)]
    s2alpha= alpha[length(alpha)]
    resi=M - c(Xm %*% alpha.coef)
    out1=Xm * resi
    offset.term = c(rep((alpha.coef[2]^2)*gamma1[length(gamma1)],n1),rep(0,n2))
    out2= s2alpha+offset.term - resi^2
    out = cbind(out1,out2)
    # ee for outcome model
    Xt0 = list()
    out_cox = matrix(0,ncol=length(beta),nrow = n1)
    for (j in (1:n_rrc)) {
      Xt0[[j]] = data.main
      gamma2.now = gamma2[((j-1)*n_gamma2+1):(n_gamma2*j)]
      Xt0[[j]][,Aname] = c(model.matrix(ME2.formula,data.main) %*% gamma2.now)
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
    out = cbind(out,out3)
    if (type==1) return(apply(out,2,sum))
    if (type==2) return(out)
  }
  if (is.boot==0) {
    Uab = ee_alpha_beta(param=param,type=2)
    Iabg = numDeriv::jacobian(ee_alpha_beta, x=param, method="simple", type=1)/n
    Iabg_ab = Iabg[,(param_sep[1]+1):param_sep[3]]
    Iabg_g  = Iabg[,(param_sep[3]+1):param_sep[5]]
    Meat = crossprod(Uab + t(Iabg_g %*% solve(Ig) %*% t(Ug)))
    cov_ab = (solve(Iabg_ab) %*% (Meat/n) %*% t(solve(Iabg_ab)))/n
  } else {
    cov_ab = diag(1,length(c(alpha,beta)))
  }
  return(list(Point=c(alpha,beta),COV=cov_ab))
}

#' Obtain the estimated mediation effect measures based on the approximate expressions under a rare
#' outcome assumption
#'
#'@import spatstat
#'@import survival
#'@import nleqslv
#'@import sandwich
#'@import magic
#'@import boot
#' @param res results from the unadj.f, orc1.f, orc2.f, or rrc.f
#' @param Mname name of the mediator
#' @param Aname name of the true exposure
#' @param Wname names of the covariates
#' @param a0 baseline exposure level in the mediation effect measures
#' @param a1 active exposure level in the mediation effect measures
#' @param cw levels of the covariates
#' @returns estimation of mediation effect measures
get_estimated_effect = function(res,Mname,Aname,Wname,a0=0,a1=1,cw=0) {
  n_alpha = length(Wname)+3
  n_beta  = length(Wname)+3
  param = res$Point
  param_cov = res$COV
  names(param)=c(paste("alpha",0:(n_alpha-2),sep=""),"s2alpha",paste("beta",1:n_beta,sep=""))
  colnames(param_cov)=rownames(param_cov) = names(param)
  name_alpha_w = names(param)[3:(n_alpha-1)]
  name_beta_w   = names(param)[(n_alpha+4):(n_alpha+n_beta)]

  alpha0=param[1];alpha1=param[2];s2alpha = param[n_alpha]
  beta1 = param[n_alpha+1];beta2 = param[n_alpha+2];beta3 = param[n_alpha+3]
  for (i in (1:length(Wname))) {
    assign(name_alpha_w[i],param[2+i])
    assign(name_beta_w[i],param[n_alpha+3+i])
  }

  # calculate the mediation measures
  .A = "(beta2*alpha1+beta3*alpha1*a1)*(a1-a0)"
  .B = paste("beta1*(a1-a0)+beta3*(alpha0+alpha1*a0+beta2*s2alpha)*(a1-a0)+0.5*(beta3^2)*s2alpha*(a1^2-a0^2)",
             paste("+beta3*(",paste(paste(name_alpha_w,"*",cw),collapse = "+"),")*(a1-a0)"))
  # NIE
  NIE_fun = function() {output = .A;return(output)}
  variable=names(param)
  NIE_D=deriv(parse(text=NIE_fun()),variable)
  NIE_D = eval(NIE_D)
  NIE_p = NIE_D[1]
  lambda= t(attr(NIE_D,"gradient"))
  V_NIE = as.vector(t(lambda) %*% param_cov %*% lambda)

  # NDE
  NDE_fun = function() {output = .B;return(output)}
  variable=names(param)
  NDE_D=deriv(parse(text=NDE_fun()),variable)
  NDE_D = eval(NDE_D)
  NDE_p = NDE_D[1]
  lambda= t(attr(NDE_D,"gradient"))
  V_NDE = as.vector(t(lambda) %*% param_cov %*% lambda)

  # TE
  TE_fun = function() {output = paste(.A,"+",.B,sep="");return(output)}
  variable=names(param)
  TE_D=deriv(parse(text=TE_fun()),variable)
  TE_D = eval(TE_D)
  TE_p = TE_D[1]
  lambda= t(attr(TE_D,"gradient"))
  V_TE = as.vector(t(lambda) %*% param_cov %*% lambda)

  # MP
  MP_fun = function() {output = paste("(",.A,")/(",.A,"+",.B,")",sep="");return(output)}
  variable=names(param)
  MP_D=deriv(parse(text=MP_fun()),variable)
  MP_D = eval(MP_D)
  MP_p = MP_D[1]
  lambda= t(attr(MP_D,"gradient"))
  V_MP = as.vector(t(lambda) %*% param_cov %*% lambda)
  list(Point=c(NIE=NIE_p,NDE=NDE_p,TE=TE_p,MP=MP_p),
       SE=c(NIE=sqrt(V_NIE),NDE=sqrt(V_NDE),TE=sqrt(V_TE),MP=sqrt(V_MP)))
}

#' Mediation analysis based on the approximate expressions under a rare outcome assumption.
#'
#' @import spatstat
#' @import survival
#' @import nleqslv
#' @import sandwich
#' @import magic
#' @import boot
#' @param data dataset.
#' @param Iname name of the indicator on the main/validation study (I=1 for validation study sample, I=0 for main study samples).
#' @param Tname name of the observed failure time.
#' @param Dname name of the failure indicator.
#' @param Mname name of the mediator
#' @param Aname name of the true exposure
#' @param Asname name of the mismeasured exposure
#' @param Wname names of the covariates
#' @param Time_sep a vector of time split points for RRC method.
#' @param a0 baseline exposure level in the mediation effect measures
#' @param a1 active exposure level in the mediation effect measures
#' @param cw levels of the covariates
#' @returns estimation of mediation effect measures based on the unadjusted method, ORC1, ORC2, and RRC
#' @examples
#' ## attach example dataset
#' attach(example_data)
#'
#' ## specify the exposure, mediator, outcome, and covariates
#' Iname="I" # name of the indicator of main/validation study samples
#' Tname="Time" # name of observed failure time
#' Dname="Event" # name of failure indicator
#' Aname="A" # name of true exposure
#' Asname = "As" # name of mismeasured exposure
#' Mname = "M" # name of mediator
#' Wname="W" # name of covariates
#'
#' ## specify split points for RRC method
#' Time_sep=c(0,15,25,35)
#'
#' ## specify the conditional values in the mediation effect measures
#' a0=0; a1=1; cw=0
#'
#' ## run the following code to perform the mediation analysis with approximate expressions
#' ## mediate_approx(example_data,Iname,Tname,Dname,Aname,
#' ##                Asname,Mname,Wname,Time_sep,a0,a1,cw)
mediate_approx = function(data,
                          Iname,
                          Tname,
                          Dname,
                          Aname,
                          Asname,
                          Mname,
                          Wname,
                          Time_sep,
                          a0,
                          a1,
                          cw) {
  Wf = paste(Wname,collapse="+")
  Outcome.formula = as.formula(paste("survival::Surv(",Tname,",",Dname,") ~", Aname, "+", Mname, "+ I(", Aname,"*", Mname, ") +", Wf))
  Outcome.formula.rrc = as.formula(paste("survival::Surv(tstart,",Tname,",",Dname,") ~", Aname, "+", Mname, "+ I(", Aname,"*", Mname, ") +", Wf))
  Mediator.formula = as.formula(paste(Mname, "~", Aname, "+", Wf))
  ME1.formula = as.formula(paste(Aname, "~", Asname, "+", Wf))
  ME2.formula = as.formula(paste(Aname, "~", Asname, "+", Mname, "+", Wf))
  unadj.e = unadj.f(Outcome.formula,Mediator.formula,data,Iname,Tname,Dname,is.boot=0)
  orc1.e=orc1.f(Outcome.formula,Mediator.formula,ME1.formula,data,Iname,Tname,Dname,Aname,Asname,Mname,Wname,is.boot=0)
  orc2.e=orc2.f(Outcome.formula,Mediator.formula,ME1.formula,ME2.formula,data,Iname,Tname,Dname,Aname,Asname,Mname,Wname,is.boot=0)
  rrc.e=rrc.f(Outcome.formula,Outcome.formula.rrc,Mediator.formula,ME1.formula,ME2.formula,
              data,Time_sep,Iname,Tname,Dname,Aname,Asname,Mname,Wname,is.boot=0)

  unadj.res = get_estimated_effect(res=unadj.e,Mname,Aname,Wname,a0,a1,cw)
  orc1.res = get_estimated_effect(res=orc1.e,Mname,Aname,Wname,a0,a1,cw)
  orc2.res = get_estimated_effect(res=orc2.e,Mname,Aname,Wname,a0,a1,cw)
  rrc.res = get_estimated_effect(res=rrc.e,Mname,Aname,Wname,a0,a1,cw)

  Unadj.res = rbind(unadj.res$Point,unadj.res$Point-1.96*unadj.res$SE,unadj.res$Point+1.96*unadj.res$SE)
  ORC1.res = rbind(orc1.res$Point,orc1.res$Point-1.96*orc1.res$SE,orc1.res$Point+1.96*orc1.res$SE)
  ORC2.res = rbind(orc2.res$Point,orc2.res$Point-1.96*orc2.res$SE,orc2.res$Point+1.96*orc2.res$SE)
  RRC.res = rbind(rrc.res$Point,rrc.res$Point-1.96*rrc.res$SE,rrc.res$Point+1.96*rrc.res$SE)

  rownames(Unadj.res)=rownames(ORC1.res)=rownames(ORC2.res)=rownames(RRC.res)=c("point","CI.lower","CI.upper")
  colnames(Unadj.res)=colnames(ORC1.res)=colnames(ORC2.res)=colnames(RRC.res)=c("NIE","NDE","TE","MP")
  list(Unadj=Unadj.res,ORC1=ORC1.res,ORC2=ORC2.res,RRC=RRC.res)
}

