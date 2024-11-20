#' Unadjusted method for regression coefficients (no measurement error)
#'
#' @param Outcome.formula formula for the Cox outcome model.
#' @param Mediator.formula formula for the linear mediator model.
#' @param data dataset
#' @param Tname name of the observed failure time.
#' @param Dname name of the failure indicator.
#' @returns Unadjusted estimates of the regression coefficients.
param_unadj_approx.f=function(Outcome.formula,Mediator.formula,data,Tname,Dname,is.boot=0) {
  # mediator model
  alpha.lm = lm(Mediator.formula,data)
  alpha = alpha.lm$coefficient
  s2alpha = mean((alpha.lm$residuals)^2)
  Xm = model.matrix(Mediator.formula,data)
  M = data[,all.vars(Mediator.formula)[1]]
  n=dim(data)[1]
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
  beta = survival::coxph(Outcome.formula,data=data,ties="breslow")$coefficients
  Xt = model.matrix(Outcome.formula,data)[,-1]
  Time = data[,Tname]
  D = data[,Dname]
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
    Ib = numDeriv::jacobian(ee_beta, x=beta, method="simple", type=1)/n
    cov_b = (solve(Ib) %*% (crossprod(Ub)/n) %*% t(solve(Ib)))/n
  } else {
    cov_b=diag(1,length(beta))
  }
  return(list(Point=c(alpha,beta),COV=magic::adiag(cov_a,cov_b)))
}

#' Mediation analysis based on the approximate expressions under a rare outcome assumption (no measurement error).
#'
#' @import spatstat
#' @import survival
#' @import nleqslv
#' @import sandwich
#' @import magic
#' @import boot
#' @param data dataset.
#' @param Tname name of the observed failure time.
#' @param Dname name of the failure indicator.
#' @param Mname name of the mediator
#' @param Aname name of the true exposure
#' @param Wname names of the covariates
#' @param a0 baseline exposure level in the mediation effect measures
#' @param a1 active exposure level in the mediation effect measures
#' @param cw levels of the covariates
#' @returns estimation of mediation effect measures based on the unadjusted method
mediate_unadj_approx = function(data,
                                Tname,
                                Dname,
                                Aname,
                                Mname,
                                Wname,
                                a0,
                                a1,
                                cw) {
  Wf = paste(Wname,collapse="+")
  Outcome.formula = as.formula(paste("survival::Surv(",Tname,",",Dname,") ~", Aname, "+", Mname, "+ I(", Aname,"*", Mname, ") +", Wf))
  Mediator.formula = as.formula(paste(Mname, "~", Aname, "+", Wf))
  unadj.e = param_unadj_approx.f(Outcome.formula,Mediator.formula,data,Tname,Dname,is.boot=0)
  unadj.res = get_estimated_effect(res=unadj.e,Mname,Aname,Wname,a0,a1,cw)
  Unadj.res = rbind(unadj.res$Point,unadj.res$Point-1.96*unadj.res$SE,unadj.res$Point+1.96*unadj.res$SE)
  rownames(Unadj.res)=c("point","CI.lower","CI.upper")
  colnames(Unadj.res)=c("NIE","NDE","TE","MP")
  Unadj.res
}



#' Unadjusted method for regression coefficients and baseline hazards (no measurement error)
#'
#' @param Outcome.formula formula for the Cox outcome model.
#' @param Mediator.formula formula for the linear mediator model.
#' @param data dataset.
#' @param Iname name of the indicator on the main/validation study (I=1 for validation study sample, I=0 for main study samples).
#' @param Tname name of the observed failure time.
#' @param Dname name of the failure indicator.
#' @returns Unadjusted estimates of the regression coefficients.
param_unadj_exact.f=function(Outcome.formula,Mediator.formula,data,Tname,Dname,is.boot=0) {
  n=dim(data)[1]
  # mediator model
  alpha.lm = lm(Mediator.formula,data)
  alpha = alpha.lm$coefficient
  s2alpha = mean((alpha.lm$residuals)^2)
  Xm = model.matrix(Mediator.formula,data)
  M = data[,all.vars(Mediator.formula)[1]]
  # estimating equation for alpha
  # type 1, bread; type 2, meat
  alpha = c(alpha,s2alpha)
  #outcome model
  outcome.m = survival::coxph(Outcome.formula,data=data,ties="breslow")
  beta = outcome.m$coefficients
  # cumulative baseline hazard
  Lambda0 = basehaz(outcome.m,centered=FALSE)
  return(list(Point=c(alpha,beta),Lambda0=Lambda0))
}


#' Point estimation of the mediation effect based on exact expressions (no measurement error)
mediate_unadj_point_exact.f=function(Outcome.formula,Mediator.formula,data,Tname,Dname,Aname,Mname,Wname,a0,a1,cw,t.list) {
  unadj.e = param_unadj_exact.f(Outcome.formula,Mediator.formula,data,Tname,Dname,is.boot=0)
  for (i in (1:length(t.list))) {
    t=t.list[i]
    unadj.res = get_estimated_exact_effect(res=unadj.e,Mname,Aname,Wname,a0=a0,a1=a1,cw=cw,t=t)
    if (i==1) {out=c(unadj.res)}
    if (i>1) {out=rbind(out,unadj.res)}
  }
  colnames(out) = c("NIE_unadj","NDE_unadj","TE_unadj","MP_unadj")
  rownames(out) = t.list
  return(out)
}

#' Bootstrap confidence interval of the mediation effect based on exact expressions (no measurement error)
boot_unadj_point_exact.f=function(Outcome.formula,Mediator.formula,data,Tname,Dname,Aname,Mname,Wname,a0,a1,cw,t.list=c(10,20,30,40),R=50) {
  #a0=0;a1=1;cw=0;t.list=c(10,20,30,40)
  parameter.f=function(data,indices) {
    data=data[indices,]
    unadj.e = param_unadj_exact.f(Outcome.formula,Mediator.formula,data,Tname,Dname,is.boot=0)
    for (i in (1:length(t.list))) {
      t=t.list[i]
      unadj.res = get_estimated_exact_effect(res=unadj.e,Mname,Aname,Wname,a0=a0,a1=a1,cw=cw,t=t)
      if (i==1) {out=c(unadj.res)}
      if (i>1) {out=rbind(out,unadj.res)}
    }
    colnames(out) = c("NIE_unadj","NDE_unadj","TE_unadj","MP_unadj")
    rownames(out) = t.list
    out=c(out)
  }
  boot.res <- boot::boot(data=data, statistic=parameter.f, R=R)
  bootci.f=function(res=boot.res) {
    mat=matrix(NA,ncol=4*length(t.list),nrow=2)
    for (j in (1:(4*length(t.list)))) {
      mat[,j]=boot::boot.ci(res, type="perc",index=j)$percent[4:5]
    }
    mat
  }
  ci=bootci.f(res=boot.res)

  return(ci)
}




#' Mediation analysis based on the exact expressions (no measurement error)
#'
#' @import spatstat
#' @import survival
#' @import nleqslv
#' @import sandwich
#' @import magic
#' @import boot
#' @param data dataset.
#' @param Tname name of the observed failure time.
#' @param Dname name of the failure indicator.
#' @param Mname name of the mediator
#' @param Aname name of the true exposure
#' @param Wname names of the covariates
#' @param t.list a list of time points for the mediation effect measures
#' @param a0 baseline exposure level in the mediation effect measures
#' @param a1 active exposure level in the mediation effect measures
#' @param cw levels of the covariates
#' @param R number of bootstrap numbers (default 50)
#' @returns estimation of mediation effect measures based on the unadjusted method
mediate_unadj_exact = function(data,
                         Tname,
                         Dname,
                         Aname,
                         Mname,
                         Wname,
                         t.list,
                         a0,
                         a1,
                         cw,
                         R=50) {
  Wf = paste(Wname,collapse="+")
  Outcome.formula = as.formula(paste("survival::Surv(",Tname,",",Dname,") ~", Aname, "+", Mname, "+ I(", Aname,"*", Mname, ") +", Wf))
  Mediator.formula = as.formula(paste(Mname, "~", Aname, "+", Wf))

  allres =  mediate_unadj_point_exact.f(Outcome.formula,Mediator.formula,data,Tname,Dname,Aname,Mname,Wname,a0=a0,a1=a1,cw=cw,t.list=t.list)
  bootres = boot_unadj_point_exact.f(Outcome.formula,Mediator.formula,data,Tname,Dname,Aname,Mname,Wname,a0,a1,cw,t.list,R=R)

  out = matrix(paste(sprintf("%0.3f",allres), " (",sprintf("%0.3f",bootres[1,]),",",sprintf("%0.3f",bootres[2,]),")",sep=""),ncol=4)
  out=as.data.frame(out)
  colnames(out) = colnames(allres)
  rownames(out) = rownames(allres)

  out
}

