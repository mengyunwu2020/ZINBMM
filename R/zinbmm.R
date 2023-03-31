zinbmm=function(mat,Bmat,lib=NULL,K,phi_init=NULL,phi_global=T,pi_ind=F,
                tune_set=NULL,ntune=10,c='BIC',ncores=1,maxit=500,tol=1e-05,vrbs=F){

  n=dim(mat)[1]
  p=dim(mat)[2]

  if(!length(lib)) lib=rep(1,n)
  if(!length(tune_set)) tune_set=seq(0.01,20,length.out=ntune)
  phi_old=phi_init
  if (!length(phi_old)) phi_old=rep(1,p)
  pi_old=colSums(mat==0)/n
  gamma_old=rep(0, dim(Bmat)[2])
  beta_old=log(colMeans(mat)+1e-03)

  iter_init=1
  repeat{

    if(vrbs) print(paste0('The ',iter_init,'-th initialization iteration'))
    set_mat=matrix(rep(exp(Bmat%*%gamma_old+log(lib)),p),nrow=n,byrow=F)
    pi_mat=matrix(rep(pi_old,n),nrow=n,byrow=T)
    phi_mat=matrix(rep(phi_old,n),nrow=n,byrow=T)
    beta_mat=matrix(rep(exp(beta_old),n),nrow=n,byrow=T)
    mu_mat=set_mat*beta_mat

    ## unobserved m
    m=pi_mat/(pi_mat+(1-pi_mat)*(phi_mat/(mu_mat+phi_mat))^phi_mat)*(mat==0)

    ## update pi
    pi_new=colSums(m)/n

    ## update gamma
    if(length(gamma_old)==1){gamma_new=gamma_old}else{
      gamma_new=gamma_update(mat,lib,Bmat,gamma_old,beta_old,phi_old,k=1,z=NULL,m,maxit=maxit)}

    ## update beta
    beta_new=beta_vec_update(mat,lib,Bmat,gamma_new,beta_old,phi_old,beta_global=NULL,
                             tune=NULL,z_vec=rep(1,n),m,penalize=F,maxit=300)
    ## update phi
    if (length(phi_init)){phi_new=phi_old}else{
      psi_new=psi_update(mat,lib,Bmat,psi_old=log(phi_old),beta=beta_new,
                         gamma=gamma_new,m,psi_global=phi_global,maxit=500)
    }
    phi_new=exp(psi_new)

    iter_init=iter_init+1

    if(iter_init>300){break}

    diff=max(mean(abs(beta_new-beta_old),trim=5e-03),mean(abs(pi_new-pi_old),trim=5e-03),
             max(abs(gamma_new-gamma_old)), max(abs(phi_new- phi_old)))

    if(vrbs) print(paste0('The difference is ', round(diff,4)))
    if(diff<tol){break}
    beta_old=beta_new
    pi_old=pi_new
    gamma_old=gamma_new
    phi_old=phi_new

  }
  if(vrbs) print('initalization end')

  init_set=list('beta_global'=beta_new,'pi'=pi_new,'phi'=phi_new,'gamma'=gamma_new)

  mc.cores= ncores
  models=mclapply(1:length(tune_set),function(i){
    res=zinbmm_onetune(mat,Bmat,lib,K,init_set,tune=tune_set[i],pi_ind,maxit,tol,vrbs)
    return(res)
  },mc.cores=mc.cores)

  if (c=='BIC'){
    cset=sapply(models,function(x)x$bic) }else if(c=='AIC'){
      cset=sapply(models,function(x)x$aic)}else{
        cset=sapply(models,function(x)x$bic_modified)
      }

  model0=models[[which.min(cset)]]
  clust=apply(model0$result$zmatrix, 1, which.max)
  beta=which(apply(model0$result$beta,1,function(x)length(unique(x)))!=1)
  gamma=model0$result$gamma
  pclust=model0$result$pclust
  pi=model0$result$pi

  output=list(clust,pclust,beta,gamma,pi)
  names(output)=c('clust','pclust','beta_info','batch_effect','zero_infl')

  return(output)
}

zinbmm_onetune=function(mat,Bmat,lib,K,init_set,tune,pi_ind,
                        maxit,tol,vrbs=F){

  n=dim(mat)[1]
  p=dim(mat)[2]

  if(K>1){
    iter_em=1

    beta_global=init_set$beta_global
    phi=init_set$phi
    beta_init=NULL
    for(k in 1:K){
      beta_init=cbind(beta_init,beta_global+0.01*(k-1))
    }
    pi_old=matrix(rep(init_set$pi,K),nrow=p,byrow=F)
    beta_old=beta_init
    gamma_old=init_set$gamma
    pc_old=rep(1/K,K)
    repeat{
      if(vrbs) print(paste0('The ',iter_em,'-th iteration'))

      mult_pdf=matrix(NA,nrow=n,ncol=K)
      set_vec=exp(Bmat%*%gamma_old)*lib

      for(i in 1:n){
        for(l in 1:K){
          mult_pdf[i,l]=as.numeric(mult_density(x=mat[i,],mu=set_vec[i]*exp(beta_old[,l]),
                                                phi, pi=pi_old[,l]))
        }
      }

      ###Updating z matrix
      z=apply(mult_pdf,1,function(x){
        d=brob(x)*pc_old
        return(as.numeric(d/sum(d)))
      })
      z=t(z)

      m=list()
      set_mat=matrix(rep(exp(Bmat%*% gamma_old+log(lib)),p),nrow=n,byrow=F)
      phi_mat=matrix(rep(phi,n),nrow=n,byrow=T)

      for (l in 1:K){
        pi_mat=matrix(rep(pi_old[,l],n),nrow=n,byrow=T)
        beta_mat<-matrix(rep(exp(beta_old[,l]),n),nrow=n,byrow=T)
        mu_mat=set_mat*beta_mat
        m[[l]]=pi_mat/(pi_mat+(1-pi_mat)*(phi_mat/(mu_mat+phi_mat))^phi_mat)*(mat==0)
      }
      ## update pc
      pc_new=colSums(z)/n

      ## update pi
      pi_new=matrix(NA,nrow=p,ncol=K)
      if(pi_ind==F){
        for ( l in 1:K){
          s=z[,l]*m[[l]]
          s=colSums(s)/sum(z[,l])
          pi_new[,l]=s
        }}else{
          s=0
          for( l in 1:K){
            s=s+z[,l]*m[[l]]
          }
          s=colSums(s)/n
          pi_new=matrix(data=s,nrow=p,ncol=K,byrow=F)
        }
      pi_new[which(pi_new>=0.95)]=0.95

      ## update gamma IRLS
      if(length(gamma_old)==1){gamma_new=gamma_old}else{
        gamma_new=gamma_update(mat,lib,Bmat,gamma_old,beta_old,phi,
                               k=K,z=z,m,maxit=maxit)
      }

      ## update beta IRLS
      beta_new=beta_mat_update(mat,lib,Bmat,k=K,gamma_new,beta_old,phi,
                               beta_global,tune,z,m,penalize=T,maxit=maxit)

      iter_em=iter_em+1
      if(iter_em>maxit){break}

      diff=max(mean(abs(beta_new-beta_old),trim=5e-03),
               mean(abs(pi_new-pi_old),trim=5e-03),
               max(abs(gamma_new-gamma_old)),
               max(abs(pc_new- pc_old)))

      if(vrbs) print(paste0('The difference is ', round(diff,4)))
      if(diff<tol){break}

      beta_old=beta_new
      gamma_old=gamma_new
      pi_old=pi_new
      pc_old=pc_new
    }

    result_list=list('beta'=beta_new,'pi'=pi_new,'gamma'=gamma_new,'zmatrix'=z,'pclust'=pc_new)
    log_lik_vector=LogLikFunc(mat,lib,Bmat,k=K,gamma_new,beta_new,phi,pi_new,pc_new,z)
    num_beta=apply(beta_new,1,function(x)length(unique(x)))
    s=sum(num_beta)
  }else{
    result_list=list('beta'=init_set$beta_global,'pi'=init_set$pi,'gamma'=init_set$gamma)
    log_lik_vector=LogLikFunc(mat,lib,Bmat,k=1,gamma=init_set$gamma,beta=init_set$beta_global,phi=init_set$phi,pi=init_set$pi)
    s=length(init_set$beta_global)
  }

  BIC=-2*log_lik_vector+log(n)*(K-1+s+dim(Bmat)[2]+p*K)
  AIC=-2*log_lik_vector+2*(K-1+s+dim(Bmat)[2]+p*K)
  BIC_modified=-2*log_lik_vector+log(n)*log(log(p))*(K-1+s+dim(Bmat)[2]+p*K)

  return(list('result'=result_list,'bic'=BIC,'aic'=AIC,'bic_modified'=BIC_modified))
}

gamma_update=function(mat,lib,Bmat,gamma_old,beta,phi,
                      k,z,m,maxit=500){
  n=dim(mat)[1]
  p=dim(mat)[2]

  phi_mat=matrix(rep(phi,n),nrow=n,byrow=T)
  iter_gamma=0
  repeat{
    set_mat=matrix(rep(exp(Bmat%*%gamma_old)*lib,p),nrow=n,byrow=F)

    A_mat=B_mat=matrix(0,nrow=n,ncol=p)

    if(k == 1){
      beta_mat=matrix(rep(exp(beta),n),nrow=n,byrow=T)
      mu_mat=set_mat*beta_mat

      A_mat=(1-m)*phi_mat*(mat-mu_mat)/(mu_mat*(mu_mat+phi_mat))* mu_mat
      B_mat=A_mat+(1-m)*(-2*phi_mat*mat*mu_mat+phi_mat*mu_mat^2-phi_mat^2*mat)/
        (mu_mat^2*(mu_mat+phi_mat)^2)*mu_mat^2
    }else{
      for (l in 1:k){

        beta_mat=matrix( rep(exp(beta[,l]),n),nrow=n,byrow=T)
        mu_mat=set_mat*beta_mat

        dmumat=z[,l]*(1-m[[l]])*phi_mat*(mat-mu_mat)/(mu_mat*(mu_mat+phi_mat))*mu_mat
        dmu2mat=dmumat+z[,l]*(1-m[[l]])*(-2*phi_mat*mat*mu_mat+phi_mat*mu_mat^2-phi_mat^2*mat)/
          (mu_mat^2*(mu_mat+phi_mat)^2)*mu_mat^2
        A_mat=A_mat+dmumat
        B_mat=B_mat+dmu2mat
      }
    }
    dmu_sumj=rowSums(A_mat)
    dmu2_sumj=rowSums(B_mat)

    dgamma=colSums(dmu_sumj*Bmat)
    dgamma2=colSums(dmu2_sumj*Bmat)

    gamma_new=gamma_old-dgamma/dgamma2
    iter_gamma=iter_gamma+1

    if(iter_gamma>maxit){break}
    if(max(abs(gamma_new-gamma_old))<10^(-6)){break}
    gamma_old=gamma_new
  }

  return(gamma_new)
}

beta_vec_update=function(mat,lib,Bmat,gamma,beta_old,phi,beta_global,tune,z_vec, m_mat,penalize,maxit=500){
  n=dim(mat)[1]
  p=dim(mat)[2]

  set_mat=matrix(rep(exp(Bmat%*%gamma)*lib,p),nrow=n,byrow=F)
  phi_mat=matrix(rep(phi,n),nrow=n,byrow=T)

  iter_beta=0
  repeat{
    beta_mat=matrix(rep(beta_old,n),nrow=n,byrow=T)

    mu_mat=set_mat*exp(beta_mat)

    A_mat=phi_mat*(mat - mu_mat)/(mu_mat*(mu_mat+phi_mat))*mu_mat
    B_mat=A_mat+(-2*phi_mat*mat*mu_mat+phi_mat*mu_mat^2-phi_mat^2*mat)/
      (mu_mat^2*(mu_mat+phi_mat)^2)*mu_mat^2

    omega_mat=B_mat
    tau_mat=beta_mat-A_mat/B_mat
    a_mat=matrix(0,nrow=n,ncol=p)

    loc=which(B_mat==0,arr.ind=T)

    if(length(loc)){
      tau_mat[loc]=0
      a_mat[loc]=A_mat[loc]
      index=unique(loc[,2])
    }

    beta_tilde=colSums(z_vec*(1-m_mat)*(omega_mat*tau_mat-a_mat))
    beta_tilde=beta_tilde/colSums((z_vec*(1-m_mat)*omega_mat))
    beta_control=log(colSums(z_vec*mat)/sum(z_vec)+1e-02)
    beta_tilde[which(is.nan(beta_tilde))]=beta_control[which(is.nan(beta_tilde))]
    beta_tilde[which(beta_tilde<=log(1e-02)|beta_tilde>=log(1e03))]=
      beta_control[which(beta_tilde<=log(1e-02)|beta_tilde>=log(1e03))]


    if(penalize){
      betanew=tune*sign(beta_tilde-beta_global)+colSums((z_vec*(1-m_mat)*omega_mat*tau_mat))


      betanew=betanew/colSums((z_vec*(1-m_mat)*omega_mat))

      index_select=which(sign(beta_tilde-beta_global)*(betanew-beta_global)<=0)
      betanew[index_select]=beta_global[index_select]
    }else{
      betanew=beta_tilde
    }
    betanew[which(is.nan(betanew))]=beta_control[which(is.nan(betanew))]
    betanew[which(betanew<=log(1e-02)|betanew>=log(1e03))]=beta_control[which(betanew<=log(1e-02)|betanew>=log(1e03))]
    betanew[which(is.na(betanew))]=beta_control[which(is.na(betanew))]

    iter_beta=iter_beta +1
    if(iter_beta>maxit){break}
    if(max(abs(betanew-beta_old))<10^(-6)){break}
    beta_old=betanew
  }
  return(betanew)
}

beta_mat_update=function(mat,lib,Bmat,k,gamma,beta_old,phi,
                         beta_global,tune,z,m,penalize=T,maxit=500){
  n=dim(mat)[1]
  p=dim(mat)[2]
  beta_new_mat=matrix(NA,nrow=p,ncol=k)
  for ( l in 1:k){
    beta_new_mat[,l]=beta_vec_update(mat,lib,Bmat,gamma,beta_old[,l],phi,beta_global,tune,
                                     z_vec=z[,l], m[[l]],penalize,maxit)
  }
  return(beta_new_mat)
}

psi_update=function(mat,lib,Bmat,psi_old,beta,gamma,m,psi_global=T,maxit=500){

  n=dim(mat)[1]
  p=dim(mat)[2]

  set_mat=matrix(rep(exp(Bmat%*%gamma)*lib,p),nrow=n,byrow=F)
  beta_mat=matrix(rep(exp(beta),n),nrow=n,byrow=T)
  mu_mat=set_mat*beta_mat
  iter_psi=0
  repeat{
    A_mat=B_mat=matrix(0,nrow=n,ncol=p)
    psi_mat=matrix(rep(psi_old,n),nrow=n,byrow=T)

    phi_mat=exp(psi_mat)

    A_mat=(1-m)*(digamma(mat+phi_mat)-digamma(phi_mat)-
                   (mat+phi_mat)/(mu_mat+phi_mat)+log(phi_mat) -
                   log(mu_mat+phi_mat)+1)*phi_mat
    B_mat=A_mat+(1-m)*(trigamma(mat+phi_mat)-trigamma(phi_mat)-
                         (mu_mat-mat)/(mu_mat+phi_mat)^2+
                         1.0/phi_mat-1.0/(mu_mat+phi_mat))*phi_mat^2
    if(psi_global){
      dpsi=sum(A_mat)
      dpsi2=sum(B_mat)
    }else{
      dpsi=colSums(A_mat)
      dpsi2=colSums(B_mat)
    }

    psi_new=psi_old-dpsi/dpsi2

    iter_psi=iter_psi+1
    if(iter_psi>maxit){break}
    if(max(abs(psi_new-psi_old))<10^(-7)){break}
    psi_old=psi_new
  }
  return(psi_new)
}


mult_density=function(x,mu,phi,pi){
  density_log=log(pi*(x==0)+ (1-pi)*dnbinom(x,mu=mu,size=phi,log=F))
  density_log[which(density_log==-Inf)]=log(1e-20)
  sum=sum(density_log)
  return(sum)
}

LogLikFunc=function(mat,lib,Bmat,k,gamma,beta,phi,pi,pc,z){
  n=dim(mat)[1]
  p=dim(mat)[2]

  set_vec=exp(Bmat%*%gamma)*lib

  if( k > 1){
    mult_pdf=matrix(NA,nrow=n,ncol=k)
    for(i in 1:n){
      for(l in 1:k){
        mult_pdf[i,l]=as.numeric(mult_density(x=mat[i,],mu=set_vec[i]*exp(beta[,l]),
                                              phi,pi=pi[,l]))
      }
    }
    pc_mat=matrix(rep(pc,n),nrow=n,byrow=T)
    loglik=sum(z*(log(pc_mat)+mult_pdf))
  }else{

    mult_pdf=c()
    for(i in 1:n){
      mult_pdf[i]=as.numeric(mult_density(x=mat[i,],mu=set_vec[i]*exp(beta),phi,pi=pi))
    }
    loglik=sum(mult_pdf)
  }

  return(loglik)
}

