########################################
##### FULL PARAMETRIC SIMULATIONs ###### 
########################################

# Thesis seed: 11,2,3
#              16032021, 17032021, 18032021
flush.console()
data.seed = 03042021

## Library                          -----
library(readr)
library(ggplot2)
library(ggpubr)
library(latex2exp)
library(gridExtra)
library(extrafont)
library(ggdendro)
library(reshape2)
library(grid)
library(mvtnorm)

### Objects                         -----
J = 25  # Framework components
L = 3   # Tessiue type 
S = 2   # States 

n.chain = 2
seeds = 160
count = 0
Names = rep(c("High_Low",
              "High_Mid",
              "High_High",
              "Mid_Low",
              "Mid_Mid",
              "Mid_High",
              "Low_Low",
              "Low_Mid",
              "Low_High"),2)


for(sigma.test in c(1,.1,.01)){
  for(ES in c(0.5,1,3)){
    count = count +1 
    for(N in c(5,10)){
      for(chain in 1: n.chain){
        Effect.Size = ES
        
        FileName = paste("05042021_Unscaled_Normal_chain",chain,"N",N,"Sigma2_ES",Names[count], sep = "_")
        
        
        sub = paste0("sub.",1:N)
        loc = paste0("loc.",1:L)
        measure = paste0("comp.",1:J)
        states = c("diseased","healthy")
        
        ### Effects & Dataset               -----
        set.seed(data.seed)
        
        # Xi
        Xi.j    = c(2,2,2,2,5,5,5,15,15,15,
                    2,2,2,2,5,5,5,15,15,15,
                    2,2,2,2,10)
        names(Xi.j) = measure
        
        Omega.j = rep(1,J)**2
        names(Omega.j) = measure
        
        Mu.jl   = matrix(data = rnorm(J*L, mean = Xi.j , sd = sqrt(Omega.j)), ncol = L, nrow = J)
        colnames(Mu.jl) = loc
        rownames(Mu.jl) = measure
        
        lambda = 0.3^2
        delta.i = rnorm(N, mean = 0, sd = sqrt(lambda))
        names(delta.i) = sub
        
        rho.j2 = paste0("grp.", c(rep(1,5),rep(2,12),rep(3,8))) 
        names(rho.j2) = measure
        
        Theta_Star = list()
        m = c(0,0,0) * Effect.Size
        names(m) = loc
        Theta_Star[["grp.1"]] = m
        m = c(1,1,1) * Effect.Size
        names(m) = loc 
        Theta_Star[["grp.2"]] = m
        m =-c(1,1,1) * Effect.Size
        names(m) = loc
        Theta_Star[["grp.3"]] = m
        
        sigma2 = sigma.test**2
        
        ID = c(sapply(sub, function(x) rep(x,L*S)))
        State = rep(c(sapply(states, function(x) rep(x,3))),N)
        Location = rep(loc,S*N)
        Matrix = NULL
        names = NULL
        
        for(j in measure){
          grp = rho.j2[j]
          thet = Theta_Star[[grp]][Location]
          thet[State == "healthy"] = 0
          effects = Mu.jl[j,Location] + thet  + delta.i + rnorm(N*L*S, mean = 0, sd = sqrt(sigma2))
          Matrix = cbind(Matrix, effects)
          names = c(names,j)
          colnames(Matrix) = names
        }
        
        Data.Simulated = data.frame(ID,Matrix, State, Location)
        remove(effects, grp, ID, j, Matrix, names, rho.j2, State, thet , Location)     
        ### Prior Setting                   ------
        
        # Xij (effetto fisso citochina)
        m0j   =  Xi.j
        s20j  =  rep(5,25)/10
        
        # Omegaj (varianza effetto random)
        a0wj = rep(10, J)
        b0wj = rep(2 , J)
        # plot(seq(0,5,by=0.01),dgamma(1/seq(0,5,by=0.01), 10/2,2/2), type ="l")
        
        # eta_sigma
        eta0  = c(0,0,0)
        m0 = 3
        SIGMA0 = matrix(diag(rep(5,3)), ncol = 3, nrow = 3)
        d0 = 3
        
        # lambda
        l0       = 10
        lambda0  = 2
        
        # plot(seq(0,5,by=0.01), dgamma(x=1/seq(0,5,by=0.01), l0/2, lambda0/2), type="l")
        
        lambda_par = c(l0, lambda0)
        names(lambda_par) = c("shape","rate")
        
        # sigma
        n0    = 5
        sig20 = n0*2
        sigma_par = c(n0/2, sig20*n0/2)
        names(sigma_par) = c("shape","rate")
        
        
        PRIORS = list(sigma  = sigma_par,
                      lambda = lambda_par,
                      SIGMA_eta = list(eta0 = eta0, m0 = m0, Sigma = SIGMA0, d0 = d0 ),
                      xi = matrix(c(m0j,  s20j), ncol=2 , dimnames = list(measure,c("mj0","s2j0"))),
                      omega = matrix(c(a0wj, b0wj), ncol=2 , dimnames = list(measure,c("a0j","b0j"))))
        
        remove(m0j ,s20j,a0wj ,b0wj ,eta0,m0,SIGMA0 ,d0 ,l0,lambda0,n0, sig20,  sigma_par,lambda_par)
        
        ### Starting Values                 -----
        seeds = seeds + 1
        
        set.seed(seeds)
        
        mu_start     = t(apply(na.omit(Data.Simulated[Data.Simulated$State=="healthy",2:26]), 2, function(x)tapply(x, na.omit(Data.Simulated[Data.Simulated$State=="healthy",])$Location, mean)))
        
        
        if(chain %%2 == 1){
          Thetas = matrix(0, nrow = J, ncol=L)}else{
            Thetas = - mu_start + t(apply(na.omit(Data.Simulated[Data.Simulated$State=="diseased",2:26]), 2,
                                          function(x)tapply(x, na.omit(Data.Simulated[Data.Simulated$State=="diseased",])$Location, mean)))}
        
        rownames(Thetas) = measure
        colnames(Thetas) = loc
        
        omega_start  = apply(na.omit(Data.Simulated[2:26]), 2, var)
        xi_start     = apply(na.omit(Data.Simulated[2:26]), 2, mean)
        
        delta_start  = vector(mode = "numeric", length = N) 
        names(delta_start) = unique(Data.Simulated$ID)
        
        sigma_start  =  .1
        lambda_start =  .1
        
        SIGMA_start  = matrix(diag(c(1/5,1/5,1/5)), ncol = 3, nrow = 3)
        eta_start    = c(0,0,0)
        
        STATE = list(
          betas  = Thetas,
          omega  = omega_start,
          xi     = xi_start,
          mu     = mu_start,
          delta  = delta_start,
          sigma  = sigma_start, 
          lambda = lambda_start,
          SIGMA  = SIGMA_start,
          eta    = eta_start,
          Iter = 0,
          sample = 0)
        
        remove(omega_start, xi_start ,mu_start ,delta_start,sigma_start,lambda_start,SIGMA_start,eta_start, Thetas)
        
        ### SAMPLER -----
        data = Data.Simulated
        prior = PRIORS
        starting_values = STATE
        sample =  1000
        burn.in = 1000
        thining = 50
        verbouse = TRUE
        diagnostic = TRUE
        
        sample.MCMC = list()
        next.sample = 1
        
        ## Prepare data 
        Y1 = reshape2::melt(Data.Simulated,  id.vars = c("ID", "State", "Location"))
        Y3 = reshape2::dcast(Y1, formula = ID+State+variable~Location )
        
        ## Prepare Gibbs
        Max.iter = burn.in + sample*thining
        STATE = starting_values
        
        if(verbouse){
          width = options()$width  
          if(width>=114){
            barwidth = 100
          }else{
            barwidth = width-14
          }
          step = round(seq(from = 1, to= Max.iter, length.out = barwidth))
          perc = round(step/Max.iter*100)
          len = as.numeric(perc>=10)
          n.cat = 1
        }
        
        SIGMA1 = solve(STATE$SIGMA)
        t0 = Sys.time()
        for(iter in 1:Max.iter){ 
          #### STEP I     - Betas                  -----
          
          for(j in rownames(STATE$betas)){
            
            data.j =  Y3[Y3$variable == j & Y3$State == 'diseased' ,]
            
            M.j  = data.j[, loc]
            ID.j = data.j$ID
            Nj = nrow(M.j)
            
            M.mus.j    = matrix(STATE$mu[j,] , ncol= L , nrow= J, byrow=TRUE)
            M.deltas.j = matrix(rep(delta.i[ID.j],3), ncol=3, nrow=Nj, byrow = FALSE)
            
            
            ga.j = M.j - M.mus.j - M.deltas.j 
            
            DELTAj1 = diag(c(1,1,1))/STATE$sigma
            
            ThetaStarInvCov =  SIGMA1  + Nj*DELTAj1
            ThetaStarMean   =  c(SIGMA1%*%STATE$eta) + c(DELTAj1%*%colSums(ga.j))
            
            ThetaStarVar  = solve(ThetaStarInvCov) 
            ThetaStarMean = c(ThetaStarVar%*%ThetaStarMean)
            
            STATE$betas[j,] = c(rmvnorm(1, mean= ThetaStarMean, sigma = ThetaStarVar))
          }
          
          #### STEP II    -  Xi                    ----
          for(j in names(STATE$xi)){
            xi_var  = 1/(1/PRIORS$xi[j,"s2j0"] + 3/STATE$omega[j])
            xi_mean = xi_var*(PRIORS$xi[j,"mj0"]/PRIORS$xi[j,"s2j0"]
                              + sum(STATE$mu[j,])/STATE$omega[j])
            
            STATE$xi[j] = rnorm(1, mean = xi_mean, sd = sqrt(xi_var))
          }
          
          #### STEP III   - Omega                  -----
          for(j in names(STATE$omega)){
            omega_shape = PRIORS$omega[j,"a0j"]/2 + L/2
            omega_rate  = PRIORS$omega[j,"b0j"]/2 + sum((STATE$mu[j,]-STATE$xi[j])^2)/2
            STATE$omega[j] = 1/rgamma(1, shape = omega_shape, rate = omega_rate)
          }
          
          
          #### STEP VI    - Mu                     ----
          for(j in rownames(STATE$mu)){
            for(l in colnames(STATE$mu)){
              yjl = Y1[Y1$Location == l & Y1$variable == j ,]
              group_effect = rep(STATE$betas[j,l], nrow(yjl))
              group_effect[yjl$State=="healthy"] = 0
              
              obs_mu = c(yjl$value - STATE$delta[yjl$ID] - group_effect)
              
              mu_var = 1/(1/STATE$omega[j] + N*S / STATE$sigma)
              mu_mean = (STATE$xi[j]/STATE$omega[j] + sum(obs_mu)/STATE$sigma)*mu_var
              
              STATE$mu[j,l] = rnorm(1, mean = mu_mean, sd = sqrt(mu_var))
            }
          }
          
          #### STEP VII  - delta                   ----
          for(i in names(STATE$delta)){
            
            Y_delta.i = Y1[Y1$ID==i,]
            
            Y_delta_theta = apply(Y_delta.i, 1, function(x) STATE$betas[x[4], x[3]])
            Y_delta_theta[Y_delta.i$State=="healthy"] = 0
            
            Y_delta_mujl =  apply(Y_delta.i, 1, function(x) STATE$mu[x[4],x[3]])
            
            obs_delta = Y_delta.i$value - Y_delta_mujl  - Y_delta_theta 
            
            delta_var = 1/(1/STATE$lambda + S*L*J/STATE$sigma)
            delta_mean = sum(obs_delta)/STATE$sigma*delta_var
            
            STATE$delta[i] = rnorm(1, mean = delta_mean, sd = sqrt(delta_var))
          }
          
          #### STEP VIII - sigma2                  -----
          
          Y1_theta =  apply(Y1, 1, function(x) STATE$betas[x[4], x[3]])
          Y1_theta[Y1$State=="healthy"] = 0
          
          Y1_mujl =  apply(Y1, 1, function(x) STATE$mu[x[4],x[3]])
          Y1_delta = apply(Y1, 1, function(x) STATE$delta[x[1]])
          obs_sigma = Y1$value - ( Y1_mujl + Y1_theta + Y1_delta)
          
          sigma_shape = PRIORS$sigma["shape"]/2 + N*J*L*S/2
          sigma_rate  = PRIORS$sigma["rate"]/2 + sum(obs_sigma**2)/2
          
          STATE$sigma = 1/rgamma(1, shape = sigma_shape, rate = sigma_rate)
          
          
          #### STEP IX   - lambda                  -----
          
          lambda_shape = PRIORS$lambda["shape"]/2 + N/2
          lambda_rate =  PRIORS$lambda["rate"]/2 + sum(STATE$delta^2)/2 
          
          STATE$lambda = 1/rgamma(1, shape = lambda_shape, rate = lambda_rate)
          
          #### STEP IX - DP prior SIGMA.eta      -----   
          
          bar_theta = colSums(STATE$betas)/J
          
          S_theta = matrix( rowSums(apply(STATE$betas, 1, function(x) (x-bar_theta)%*%t(x-bar_theta) ))/J, ncol=3)
          
          SIGMA_df = PRIORS$SIGMA_eta$d0 + J
          n_h = PRIORS$SIGMA_eta$m0 + J
          eta_mean = (PRIORS$SIGMA_eta$m0*PRIORS$SIGMA_eta$eta0 + J*bar_theta) / n_h
          SIGMA_var = PRIORS$SIGMA_eta$Sigma + S_theta + J*PRIORS$SIGMA_eta$m0 / n_h * 
            (bar_theta - PRIORS$SIGMA_eta$eta0)%*%t(bar_theta - PRIORS$SIGMA_eta$eta0)
          
          SIGMA1= rWishart(1, df = SIGMA_df, Sigma = solve(SIGMA_var))[,,1]
          
          STATE$SIGMA = solve(SIGMA1)
          
          STATE$eta = c(rmvnorm(1, mean = eta_mean, checkSymmetry = FALSE , sigma = STATE$SIGMA/SIGMA_df, method = "chol"))
          
          ### STORE RESULT ###
          
          if(iter>burn.in & iter%%thining==0){
            sample.MCMC[[paste0("sample.",next.sample)]] = STATE
            next.sample = next.sample + 1
            
          }
          
          if(verbouse)
            if(step[n.cat]==iter){
              t1 = Sys.time()
              iter.time = difftime(t1,t0,units='secs')/iter
              time.to.end = as.numeric(iter.time)*(Max.iter-iter)
              if(time.to.end>=3600*24){
                time.print = "+ 1d"
              }else if(time.to.end>3600){
                time.print =  paste0(round((time.to.end)/3600,1)," h")
              }else if(time.to.end>60){
                time.print =  paste0(round((time.to.end)/60,1)," m")
              }else{        
                time.print =  paste0(round((time.to.end),1)," s")
              }
              
              
              cat('\f')
              cat(paste0("|",paste0(rep('=', n.cat), collapse = ''), paste0(rep(' ', barwidth-n.cat), collapse = ''),
                         "|",paste0(rep(' ', 2-len[n.cat]), collapse = ''),perc[n.cat],"%")," (",time.print,")")
              n.cat= n.cat + 1
            }
          
          ### STTUFS
          
        }
        
        saveRDS(sample.MCMC, file = paste0(FileName,".rds"))
      }
    }
  }
}
