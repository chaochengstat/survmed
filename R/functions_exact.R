#' Gauss-Hermite Quadrature Approximation (copied from spatstat package)
my.ghq = function(f, mu = 0, sd = 1, ..., order = 5) {
  Hn <- spatstat.random::HermiteCoefs(order)
  Hn1 <- spatstat.random::HermiteCoefs(order - 1)
  x <- sort(Re(polyroot(Hn)))
  Hn1x <- matrix(Hn1, nrow = 1) %*% t(outer(x, 0:(order - 1),"^"))
  w <- 2^(order - 1) * factorial(order) * sqrt(pi)/(order * Hn1x)^2
  ww <- w/sqrt(pi)
  xx <- mu + sqrt(2)*sd %*% t(x)
  fall=f(xx[,1],...)
  for (i in (2:order)) {
    fall=cbind(fall,f(xx[,i],...))
  }
  c(fall %*% t(ww))
}

#' Unadjusted method (for exact expressions)
#'
#' @param Outcome.formula formula for the Cox outcome model.
#' @param Mediator.formula formula for the linear mediator model.
#' @param data dataset.
#' @param Iname name of the indicator on the main/validation study (I=1 for validation study sample, I=0 for main study samples).
#' @param Tname name of the observed failure time.
#' @param Dname name of the failure indicator.
#' @returns Unadjusted estimates of the regression coefficients.
unadj.ef=function(Outcome.formula,Mediator.formula,data,Iname,Tname,Dname,is.boot=0) {
  data[,Aname] = data[,Asname]
  data.main <<- data[data[,Iname]==0,]
  data.main <<- data.main[order(data.main[,Tname]),]
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
  alpha = c(alpha,s2alpha)
  #outcome model
  outcome.m = survival::coxph(Outcome.formula,data=data.main,ties="breslow")
  beta = outcome.m$coefficients
  # cumulative baseline hazard
  Lambda0 = basehaz(outcome.m,centered=FALSE)
  return(list(Point=c(alpha,beta),Lambda0=Lambda0))
}

#' ORC1 method (for exact expressions)
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
orc1.ef=function(Outcome.formula,Mediator.formula,ME1.formula,data,Iname,Tname,Dname,Aname,Asname,Mname,Wname,is.boot=0) {
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
  gamma = c(gamma.coef,s2gamma)
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
  outcome.m = survival::coxph(Outcome.formula,data=Xt0,ties="breslow")
  beta = outcome.m$coefficients
  Xt0 <<- Xt0
  # cumulative baseline hazard
  Lambda0 = basehaz(outcome.m,centered=FALSE)

  return(list(Point=c(alpha,beta),Lambda0=Lambda0))
}

#' ORC2 method (for exact expressions)
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
orc2.ef=function(Outcome.formula,Mediator.formula,ME1.formula,ME2.formula,data,
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
  outcome.m = survival::coxph(Outcome.formula,data=Xt0,ties="breslow")
  beta = outcome.m$coefficients
  # cumulative baseline hazard
  Xt0 <<- Xt0
  Lambda0 = basehaz(outcome.m,centered=FALSE)
  return(list(Point=c(alpha,beta),Lambda0=Lambda0))
}

#' RRC method (for exact expressions)
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
#' @returns ORC2 estimates of the regression coefficients.
rrc.ef=function(Outcome.formula,Outcome.formula.rrc,Mediator.formula,ME1.formula,ME2.formula,data,Time_sep,
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
  # mediator model
  Xm = model.matrix(Mediator.formula,data)
  Xm[1:n1,Aname] = c(model.matrix(ME1.formula,data.main) %*% gamma1[-length(gamma1)])
  M = data[,all.vars(Mediator.formula)[1]]
  alpha.lm = lm(M~Xm-1)
  alpha.coef = alpha.lm$coefficients
  s2alpha = (sum((alpha.lm$residuals)^2) - n1*(alpha.coef[2]^2)*gamma1[length(gamma1)])/n
  alpha = c(alpha.coef,s2alpha)
  # Outcome model
  Xt0 = survSplit(as.formula(paste("Surv(",Tname,",",Dname,") ~.")), data=data.main, cut=Time_sep[-1],episode="timegroup")
  index_t_main=c(index_t_main,n1)
  for (j in (1:n_rrc)) {
    #index.now = index_t_main[j]:index_t_main[j+1]
    gamma2.now = gamma2[((j-1)*n_gamma2+1):(n_gamma2*j)]
    index.now = which(Xt0$timegroup==j)
    Xt0[index.now,Aname] = c(model.matrix(ME2.formula,Xt0[index.now,]) %*% gamma2.now)
  }
  #Xt0[,Aname]= c(model.matrix(ME2.formula,data.main) %*% gamma2[1:4]
  Xt0 <<- Xt0
  outcome.m = survival::coxph(Outcome.formula.rrc,data=Xt0,ties="breslow")
  beta = outcome.m$coefficients
  # cumulative baseline hazard
  Lambda0 = basehaz(outcome.m,centered=FALSE)

  return(list(Point=c(alpha,beta),Lambda0=Lambda0))
}


#' Obtain the estimated mediation effect measures based on exact expressions
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
get_estimated_exact_effect = function(res,Mname,Aname,Wname,a0=0,a1=1,cw=0,t=0.5) {
  n_alpha = length(Wname)+3
  n_beta  = length(Wname)+3
  param = res$Point
  names(param)=c(paste("alpha",0:(n_alpha-2),sep=""),"s2alpha",paste("beta",1:n_beta,sep=""))
  #colnames(param_cov)=rownames(param_cov) = names(param)
  name_alpha_w = names(param)[3:(n_alpha-1)]
  name_beta_w   = names(param)[(n_alpha+4):(n_alpha+n_beta)]
  alpha0=param[1];alpha1=param[2];s2alpha = param[n_alpha]
  beta1 = param[n_alpha+1];beta2 = param[n_alpha+2];beta3 = param[n_alpha+3]
  for (i in (1:length(Wname))) {
    assign(name_alpha_w[i],param[2+i])
    assign(name_beta_w[i],param[n_alpha+3+i])
  }
  Lambda.t=res$Lambda0[which.min(abs(res$Lambda0$time-t))[1],"hazard"]
  alpha = param[1:(n_alpha-1)];beta = param[-c(1:n_alpha)]
  # mediation effect on log hazard ratio scale
  lambda.s=function(a0,a1) {
    fa =function(m,a) exp((beta[2]+beta[3]*a)*m)*exp(-Lambda.t*exp(beta[1]*a + beta[2]*m + beta[3]*a*m + sum(beta[-c(1:3)]*cw)))
    fb =function(m,a) exp(-Lambda.t*exp(beta[1]*a + beta[2]*m + beta[3]*a*m + sum(beta[-c(1:3)]*cw)))
    a=my.ghq(fa, mu = alpha[1]+alpha[2]*a0 + sum(alpha[-c(1:2)]*cw), sd = sqrt(s2alpha),order =25,a=a1)
    b=my.ghq(fb, mu = alpha[1]+alpha[2]*a0 + sum(alpha[-c(1:2)]*cw), sd = sqrt(s2alpha),order =25,a=a1)
    exp(beta[1]*a1+beta[4]*cw)*a/b
  }
  lambda11 = lambda.s(a0=a1,a1=a1)
  lambda10 = lambda.s(a0=a0,a1=a1)
  lambda00 = lambda.s(a0=a0,a1=a0)
  NIE.l = log(lambda11) -log(lambda10); TE.l = log(lambda11) -log(lambda00); MP.l = NIE.l/TE.l

  return(c(NIE=NIE.l,NDE=TE.l-NIE.l,TE=TE.l,MP=MP.l))
}




#' Point estimation of the mediation effect based on exact expressions by all proposed estimators
combined.f=function(Outcome.formula,Outcome.formula.rrc,Mediator.formula,ME1.formula,ME2.formula,data,Time_sep,
                    Iname,Tname,Dname,Aname,Asname,Mname,Wname,a0,a1,cw,t.list) {

  unadj.e = unadj.ef(Outcome.formula,Mediator.formula,data,Iname,Tname,Dname,is.boot=0)
  orc1.e=orc1.ef(Outcome.formula,Mediator.formula,ME1.formula,data,Iname,Tname,Dname,Aname,Asname,Mname,Wname,is.boot=0)
  orc2.e=orc2.ef(Outcome.formula,Mediator.formula,ME1.formula,ME2.formula,data,Iname,Tname,Dname,Aname,Asname,Mname,Wname,is.boot=0)
  rrc.e=rrc.ef(Outcome.formula,Outcome.formula.rrc,Mediator.formula,ME1.formula,ME2.formula,
              data,Time_sep,Iname,Tname,Dname,Aname,Asname,Mname,Wname,is.boot=0)

  for (i in (1:length(t.list))) {
    t=t.list[i]
    unadj.res = get_estimated_exact_effect(res=unadj.e,Mname,Aname,Wname,a0=a0,a1=a1,cw=cw,t=t)
    orc1.res = get_estimated_exact_effect(res=orc1.e,Mname,Aname,Wname,a0=a0,a1=a1,cw=cw,t=t)
    orc2.res = get_estimated_exact_effect(res=orc2.e,Mname,Aname,Wname,a0=a0,a1=a1,cw=cw,t=t)
    rrc.res = get_estimated_exact_effect(res=rrc.e,Mname,Aname,Wname,a0=a0,a1=a1,cw=cw,t=t)
    if (i==1) {out=c(unadj.res,orc1.res,orc2.res,rrc.res)}
    if (i>1) {out=rbind(out,c(unadj.res,orc1.res,orc2.res,rrc.res))}
  }
  colnames(out) = paste(c("NIE","NDE","TE","MP","NIE","NDE","TE","MP","NIE","NDE","TE","MP","NIE","NDE","TE","MP"
                          ),"_",rep(c("unadj","orc1","orc2","rrc"),each=4),sep="")
  rownames(out) = t.list
  return(out)
}

#' Bootstrap confidence interval of the mediation effect based on exact expressions by all proposed estimators
bootci.f=function(Outcome.formula,Outcome.formula.rrc,Mediator.formula,ME1.formula,ME2.formula,data,Time_sep,
                  Iname,Tname,Dname,Aname,Asname,Mname,Wname,a0,a1,cw,t.list=c(10,20,30,40),R=50) {
  #a0=0;a1=1;cw=0;t.list=c(10,20,30,40)
  parameter.f=function(data,indices) {
    data=data[indices,]
    unadj.e = unadj.ef(Outcome.formula,Mediator.formula,data,Iname,Tname,Dname,is.boot=0)
    orc1.e=orc1.ef(Outcome.formula,Mediator.formula,ME1.formula,data,Iname,Tname,Dname,Aname,Asname,Mname,Wname,is.boot=0)
    orc2.e=orc2.ef(Outcome.formula,Mediator.formula,ME1.formula,ME2.formula,data,Iname,Tname,Dname,Aname,Asname,Mname,Wname,is.boot=0)
    rrc.e=rrc.ef(Outcome.formula,Outcome.formula.rrc,Mediator.formula,ME1.formula,ME2.formula,
                data,Time_sep,Iname,Tname,Dname,Aname,Asname,Mname,Wname,is.boot=0)
    for (i in (1:length(t.list))) {
      t=t.list[i]
      unadj.res = get_estimated_exact_effect(res=unadj.e,Mname,Aname,Wname,a0=a0,a1=a1,cw=cw,t=t)
      orc1.res = get_estimated_exact_effect(res=orc1.e,Mname,Aname,Wname,a0=a0,a1=a1,cw=cw,t=t)
      orc2.res = get_estimated_exact_effect(res=orc2.e,Mname,Aname,Wname,a0=a0,a1=a1,cw=cw,t=t)
      rrc.res = get_estimated_exact_effect(res=rrc.e,Mname,Aname,Wname,a0=a0,a1=a1,cw=cw,t=t)
      if (i==1) {out=c(unadj.res,orc1.res,orc2.res,rrc.res)}
      if (i>1) {out=rbind(out,c(unadj.res,orc1.res,orc2.res,rrc.res))}
    }
    colnames(out) = paste(c("NIE","NDE","TE","MP","NIE","NDE","TE","MP","NIE","NDE","TE","MP","NIE","NDE","TE","MP"
                            ),"_",rep(c("unadj.res","orc1","orc2","rrc"),each=4),sep="")
    rownames(out) = t.list
    out=c(out)
  }
  boot.res <- boot::boot(data=data, statistic=parameter.f, R=R, strata = data[,Iname])
  bootci.f=function(res=boot.res) {
    mat=matrix(NA,ncol=16*length(t.list),nrow=2)
    for (j in (1:(16*length(t.list)))) {
      mat[,j]=boot::boot.ci(res, type="perc",index=j)$percent[4:5]
    }
    mat
  }
  ci=bootci.f(res=boot.res)

  return(ci)
}



#' Mediation analysis based on the exact expressions
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
#' @param t.list a list of time points for the mediation effect measures
#' @param Time_sep a vector of time split points for RRC method.
#' @param a0 baseline exposure level in the mediation effect measures
#' @param a1 active exposure level in the mediation effect measures
#' @param cw levels of the covariates
#' @param R number of bootstrap numbers (default 50)
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
#' t.list=c(10,15,20,25,30,35,40,45); a0=0; a1=1; cw=0
#'
#' ## run the following code to perform the mediation analysis with approximate expressions
#' ## Here, R=20 bootstrap replications were used for illustrative purposes
#' ## In practice, one should try a large number (e.g., R=500)
#'
#' ## mediate_exact(example_data,Iname,Tname,Dname,Aname,
#' ##                Asname,Mname,Wname,t.list,Time_sep,a0,a1,cw,R=20)
mediate_exact = function(data,
                          Iname,
                          Tname,
                          Dname,
                          Aname,
                          Asname,
                          Mname,
                          Wname,
                          t.list,
                          Time_sep,
                          a0,
                          a1,
                          cw,
                          R=50) {
  Wf = paste(Wname,collapse="+")
  Outcome.formula = as.formula(paste("survival::Surv(",Tname,",",Dname,") ~", Aname, "+", Mname, "+ I(", Aname,"*", Mname, ") +", Wf))
  Outcome.formula.rrc = as.formula(paste("survival::Surv(tstart,",Tname,",",Dname,") ~", Aname, "+", Mname, "+ I(", Aname,"*", Mname, ") +", Wf))
  Mediator.formula = as.formula(paste(Mname, "~", Aname, "+", Wf))
  ME1.formula = as.formula(paste(Aname, "~", Asname, "+", Wf))
  ME2.formula = as.formula(paste(Aname, "~", Asname, "+", Mname, "+", Wf))

  allres =  combined.f(Outcome.formula,Outcome.formula.rrc,Mediator.formula,ME1.formula,ME2.formula,data,Time_sep,
                       Iname,Tname,Dname,Aname,Asname,Mname,Wname,a0=a0,a1=a1,cw=cw,t.list=t.list)
  bootres = bootci.f(Outcome.formula,Outcome.formula.rrc,Mediator.formula,ME1.formula,ME2.formula,data,Time_sep,
                     Iname,Tname,Dname,Aname,Asname,Mname,Wname,a0=a0,a1=a1,cw=cw,t.list=t.list,R=R)

  out = matrix(paste(sprintf("%0.3f",allres), " (",sprintf("%0.3f",bootres[1,]),",",sprintf("%0.3f",bootres[2,]),")",sep=""),ncol=16)
  out=as.data.frame(out)
  colnames(out) = colnames(allres)
  rownames(out) = rownames(allres)

  list(Unadjusted=out[,1:4],ORC1=out[,5:8],ORC2=out[,9:12],RRC=out[,13:16])
}
