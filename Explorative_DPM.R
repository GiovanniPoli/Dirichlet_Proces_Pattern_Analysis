########################################
######### UNSCALED SIMULATION ########## 
########################################

# Thesis seed: 11,2,3
#              16032021, 17032021, 18032021

data.seed = 18032021


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
J = 10  # Framework components
L = 3   # Tessiue type 
S = 2   # States 

seeds = 14
count = 0
Names = c("High_Low",
          "High_Mid",
          "High_High",
          "Mid_Low",
          "Mid_Mid",
          "Mid_High",
          "Low_Low",
          "Low_Mid",
          "Low_High")
n.chain = 2

for(sigma.test in c(1,.1,.01)){
  for(Effect.Size in c(0.5,1,3)){
    count = count +1 
    for(N in c(5,10)){
      for(chain in 1:n.chain){
        
        FileName = paste("26032021_Unscaled_DPMM_chain",chain,"N",N,"Sigma2_ES",Names[count], sep = "_")
        
        sub = paste0("sub.",1:N)
        loc = paste0("loc.",1:L)
        measure = paste0("comp.",1:J)
        states = c("diseased","healthy")
        
        ### Effects & Dataset               -----
        set.seed(data.seed)
        
        # Xi
        Xi.j    = c(2,2,2,2,5,5,5,15,15,15)
        names(Xi.j) = measure
        
        Omega.j = rep(1,J)**2
        names(Omega.j) = measure
        
        Mu.jl   = matrix(data = rnorm(30, mean = Xi.j , sd = sqrt(Omega.j)), ncol = L, nrow = J)
        colnames(Mu.jl) = loc
        rownames(Mu.jl) = measure
        
        lambda = 0.3^2
        delta.i = rnorm(N, mean = 0, sd = sqrt(lambda))
        names(delta.i) = sub
        
        rho.j2 = paste0("grp.", c(1,2,2,3,2,2,3,2,3,1)) 
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
        m0j   =  c(2,2,2,2,5,5,5,15,15,15)
        s20j  =  c(5,5,5,5,5,5,5,5,5,5)/10
        
        # Omegaj (varianza effetto random)
        a0wj = rep(10, J)
        b0wj = rep(2 , J)

        # eta_sigma
        eta0  = c(0,0,0)
        m0 = 3
        SIGMA0 = matrix(diag(rep(5,3)), ncol = 3, nrow = 3)
        d0 = 3
        
        ThetaStart = list()
        
        m =  c(0,0,0) 
        names(m) = loc
        ThetaStart[["grp.1"]] = m
        
        
        # lambda
        l0       = 10
        lambda0  = 2
        
        lambda_par = c(l0, lambda0)
        names(lambda_par) = c("shape","rate")
        
        # sigma
        n0    = 5
        sig20 = n0*2
        sigma_par = c(n0/2, sig20*n0/2)
        names(sigma_par) = c("shape","rate")
        sigma2= 0.8
        
        # M 
        a0M = 5
        b0M = 7.5
        
        M = c(a0M, b0M)
        names(M) = c("shape","rate")
        
        PRIORS = list(sigma  = sigma_par,
                      lambda = lambda_par,
                      SIGMA_eta = list(eta0 = eta0, m0 = m0, Sigma = SIGMA0, d0 = d0 ),
                      xi = matrix(c(m0j,  s20j), ncol=2 , dimnames = list(measure,c("mj0","s2j0"))),
                      omega = matrix(c(a0wj, b0wj), ncol=2 , dimnames = list(measure,c("a0j","b0j"))),
                      M= M)
        
        remove(m0j ,s20j,a0wj ,b0wj ,a0M ,b0M ,eta0,m0,SIGMA0 ,d0 ,l0,lambda0,n0, sig20,  sigma_par,lambda_par,M)
        
        ### Starting Values                 -----
        theta_start = list()
        
        if(chain%%2==1){
          m=c(0,0,0)
          names(m) = loc
          theta_start[["grp.1"]] = m
          
          rho_start    =c(rep("grp.1",10))
          
        }else{
          for(ifg in 1:J){
            grp = paste0("grp."  ,ifg)
            cp  = paste0("comp." ,ifg)
            
            theta_start[[grp]] =  tapply(Data.Simulated[Data.Simulated$State=="diseased",cp], 
                                         Data.Simulated$Location[Data.Simulated$State=="diseased"], mean) -
              tapply(Data.Simulated[Data.Simulated$State=="healthy",cp],
                     Data.Simulated$Location[Data.Simulated$State=="healthy"], mean)
            
            
          }
          
          rho_start = paste0("grp.",1:10)
        }
        
        names(rho_start) =  colnames(Data.Simulated[2:11])
        
        omega_start  = apply(na.omit(Data.Simulated[2:11]), 2, var)
        xi_start     = apply(na.omit(Data.Simulated[2:11]), 2, mean)
        
        mu_start     = t(apply(na.omit(Data.Simulated[Data.Simulated$State=="healthy",2:11]), 2, function(x)tapply(x, na.omit(Data.Simulated[Data.Simulated$State=="healthy",])$Location, mean)))
        delta_start  = vector(mode = "numeric", length = N) 
        names(delta_start) = unique(Data.Simulated$ID)
        
        sigma_start  =  .1
        lambda_start =  .1
        
        SIGMA_start  = matrix(diag(c(1/5,1/5,1/5)), ncol = 3, nrow = 3)
        eta_start    = c(0,0,0)
        M_start      = 0.5
        
        STATE = list(
          theta.star = theta_start,
          rho    = rho_start,
          omega  = omega_start,
          xi     = xi_start,
          mu     = mu_start,
          delta  = delta_start,
          sigma  = sigma_start, 
          lambda = lambda_start,
          SIGMA  = SIGMA_start,
          eta    = eta_start,
          M = M_start,
          Iter = 0,
          sample = 0)
        
        remove(m,rho_start,omega_start,xi_start ,mu_start ,delta_start,sigma_start,lambda_start,SIGMA_start,eta_start,M_start, theta_start)
        
        
        ### SAMPLER -----
        data = Data.Simulated
        prior = PRIORS
        starting_values = STATE
        sample =  1000
        burn.in = 1000
        thining = 50
        verbouse = TRUE
        
        sample.MCMC = list()
        next.sample = 1
        
        Theta_Integrated_Log_Likelihood = function(data , STATE, Cytokine){
          
          Nj = nrow(data)
          M. = as.matrix(data[,loc])
          
          M.mus    = matrix(rep(STATE$mu[Cytokine,],       Nj), ncol=3, nrow=Nj, byrow = TRUE)
          M.deltas = matrix(rep(STATE$delta[data$ID],3)       , ncol=3, nrow=Nj, byrow = FALSE)
          
          M.gammas = M. - M.mus - M.deltas
          
          DELTAj1 = diag(c(1,1,1))/STATE$sigma 
          SIGMA1 = solve(STATE$SIGMA)
          
          sumgj = colSums(M.gammas)
          A1. = DELTAj1*Nj+SIGMA1
          A. = solve(A1.)                                             
          a. = c(A. %*% (DELTAj1%*%sumgj + SIGMA1%*%STATE$eta))
          
          logk1 = -0.5*Nj*3 * log(2*pi*STATE$sigma)
          logk2 = c(-0.5* ( sum(apply(M.gammas, 1 , function(x) t(x)%*%DELTAj1%*%x)) + t(STATE$eta) %*% SIGMA1 %*% STATE$eta ))
          logk3 = c(0.5 * t(a.)%*% A1. %*% a.)
          logk4 = 1.5*log(2*pi) + 0.5*determinant(A., logarithm = TRUE)$modulus[1]
          
          sum(c(logk1,logk2,logk3,logk4))
        }
        seeds = seeds + 1
        
        set.seed(seeds)
        
        ## Prepare data 
        # (6*8)*27 + (6*6)*6 
        Y1 = reshape2::melt(Data.Simulated,  id.vars = c("ID", "State", "Location"))
        Y3 = reshape2::dcast(Y1, formula = ID+State+variable~Location )
        ## Usefull Constant & Functions
        labels = paste0("grp.",1:28)
        unseen_groups =  labels[! labels %in% names(STATE$theta.star)]
        
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
          
          #### STEP I    - rho                     ----
          for(j in names(STATE$rho)){
            n.obs = table(STATE$rho)
            old = STATE$rho[j]
            Yj.2. = Y1[Y1$State=="diseased" & Y1$variable == j,]
            Y3j.2. = Y3[Y3$State=="diseased" & Y3$variable == j,]
            
            # Remove state
            if(n.obs[old]==1){
              STATE$theta.star = STATE$theta.star[names(STATE$theta.star) != old]
              unseen_groups = c(c(old, use.names = FALSE),unseen_groups) 
            }else{
              n.obs[old] = n.obs[old] - 1
            }
            
            probs = NULL
            
            for(h in names(STATE$theta.star)){
              means =  STATE$mu[j, Yj.2.$Location] + STATE$delta[Yj.2.$ID] + STATE$theta.star[[h]][Yj.2.$Location]
              
              probs = c(probs, log(n.obs[h]) + sum(dnorm(Yj.2.$value, mean = means, sd = sqrt(STATE$sigma), log = TRUE)))
            }
            
            probs = c(probs,
                      log(STATE$M) + Theta_Integrated_Log_Likelihood(data = Y3j.2., STATE = STATE, Cytokine =  j))
            names(probs)[length(probs)] = unseen_groups[1]
            
            probs = probs - max(probs)
            
            new = sample(names(probs),1, prob = exp(probs))
            
            STATE$rho[j] = new
            
            if(new==unseen_groups[1]){
              unseen_groups = unseen_groups[2:length(unseen_groups)]
              
              Nj = nrow(Y3j.2.)
              M. = as.matrix(Y3j.2.[,loc])
              
              M.mus    = matrix(rep(STATE$mu[j,],       Nj),   ncol=3, nrow=Nj, byrow = TRUE)
              M.deltas = matrix(rep(STATE$delta[Y3j.2.$ID],3), ncol=3, nrow=Nj, byrow = FALSE)
              
              M.gammas = M. - M.mus - M.deltas
              
              DELTAj1 = diag(c(1,1,1))/STATE$sigma 
              
              sumgj = colSums(M.gammas)
              A1. = DELTAj1*Nj+SIGMA1
              A. = solve(A1.)                                             
              a. = c(A. %*% (DELTAj1%*%sumgj + SIGMA1%*%STATE$eta))
              
              STATE$theta.star[[new]] = c(rmvnorm(1, mean= a., sigma = A.))
              names(STATE$theta.star[[new]]) = loc
            }
            
          }
          
          #### STEP II   - theta.star              -----
          
          for(h in names(STATE$theta.star)){
            group = names(STATE$rho[STATE$rho==h])
            indexs = Y3[Y3$variable%in%group & Y3$State == 'diseased',"variable"]
            
            M.j =  Y3[Y3$variable%in%group & Y3$State == 'diseased' ,loc]
            Nj = nrow(M.j)
            
            M.mus.j    = STATE$mu[indexs,]
            M.deltas.j = matrix(rep(delta.i[Y3j.2.$ID],3), ncol=3, nrow=Nj, byrow = FALSE)
            
            
            ga.j = M.j - M.mus.j - M.deltas.j 
            
            DELTAj1 = diag(c(1,1,1))/STATE$sigma
            
            ThetaStarInvCov =  SIGMA1  + Nj*DELTAj1
            ThetaStarMean   =  c(SIGMA1%*%STATE$eta) + c(DELTAj1%*%colSums(ga.j))
            
            ThetaStarVar  = solve(ThetaStarInvCov) 
            ThetaStarMean = c(ThetaStarVar%*%ThetaStarMean)
            
            STATE$theta.star[[h]] = c(rmvnorm(1, mean= ThetaStarMean, sigma = ThetaStarVar))
            names(STATE$theta.star[[h]]) = loc
          }
          
          #### STEP III  - Xi                      ----
          for(j in names(STATE$xi)){
            xi_var  = 1/(1/PRIORS$xi[j,"s2j0"] + 3/STATE$omega[j])
            xi_mean = xi_var*(PRIORS$xi[j,"mj0"]/PRIORS$xi[j,"s2j0"]
                              + sum(STATE$mu[j,])/STATE$omega[j])
            
            STATE$xi[j] = rnorm(1, mean = xi_mean, sd = sqrt(xi_var))
          }
          
          #### STEP IV   - Omega                   -----
          for(j in names(STATE$omega)){
            omega_shape = PRIORS$omega[j,"a0j"]/2 + L/2
            omega_rate  = PRIORS$omega[j,"b0j"]/2 + sum((STATE$mu[j,]-STATE$xi[j])^2)/2
            STATE$omega[j] = 1/rgamma(1, shape = omega_shape, rate = omega_rate)
          }
          
          
          #### STEP VI   - Mu                     ----
          for(j in rownames(STATE$mu)){
            for(l in colnames(STATE$mu)){
              yjl = Y1[Y1$Location == l & Y1$variable == j ,]
              group_effect = rep(STATE$theta.star[[STATE$rho[j]]][l], nrow(yjl))
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
            
            Y_delta_theta = apply(Y_delta.i, 1, function(x) STATE$theta.star[[STATE$rho[x[4]]]][x[3]])
            Y_delta_theta[Y_delta.i$State=="healthy"] = 0 
            
            Y_delta_mujl =  apply(Y_delta.i, 1, function(x) STATE$mu[x[4],x[3]])
            
            
            obs_delta = Y_delta.i$value - Y_delta_mujl  - Y_delta_theta 
            
            delta_var = 1/(1/STATE$lambda + S*L*J/STATE$sigma)
            delta_mean = sum(obs_delta)/STATE$sigma*delta_var
            
            STATE$delta[i] = rnorm(1, mean = delta_mean, sd = sqrt(delta_var))
          }
          
          #### STEP VIII - sigma2                  -----
          
          Y1_theta =  apply(Y1, 1, function(x) STATE$theta.star[[STATE$rho[x[4]]]][x[3]])
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
          
          #### STEP x XI - DP prior SIGMA.eta      -----   
          Hi = length(names(STATE$theta.star)) 
          
          bar_theta = rowSums(sapply(STATE$theta.star, function(x) x))/Hi
          S_theta = Reduce("+", x = lapply(STATE$theta.star, function(x) (x-bar_theta)%*%t(x-bar_theta))) / Hi
          
          SIGMA_df = PRIORS$SIGMA_eta$d0 + Hi
          n_h = PRIORS$SIGMA_eta$m0 + Hi
          eta_mean = (PRIORS$SIGMA_eta$m0*PRIORS$SIGMA_eta$eta0 + Hi*bar_theta) / n_h
          SIGMA_var = PRIORS$SIGMA_eta$Sigma + S_theta + Hi*PRIORS$SIGMA_eta$m0 / n_h * 
            (bar_theta - PRIORS$SIGMA_eta$eta0)%*%t(bar_theta - PRIORS$SIGMA_eta$eta0)
          
          SIGMA1= rWishart(1, df = SIGMA_df, Sigma = solve(SIGMA_var))[,,1]
          
          STATE$SIGMA = solve(SIGMA1)
          
          STATE$eta = c(rmvnorm(1, mean = eta_mean, checkSymmetry = FALSE , sigma = STATE$SIGMA/SIGMA_df, method = "chol"))
          
          #### STEP XIi - M                        ----- 
          phi = rbeta(1, shape1 = STATE$M + 1, shape2 = 27) # 27 n citochine
          odds = (PRIORS$M["shape"] + Hi -1)/(27*(PRIORS$M["rate"]-log(phi)))  
          pi_G1 = odds/(1+odds)
          
          if(sample(c(TRUE,FALSE), size = 1, prob = c(pi_G1,1-pi_G1))){
            STATE$M = rgamma(1, shape = PRIORS$M["shape"] + Hi,     rate =  PRIORS$M["rate"]- log(phi))
          }else{
            STATE$M = rgamma(1, shape = PRIORS$M["shape"] + Hi - 1, rate =  PRIORS$M["rate"]- log(phi))
          }
          
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
