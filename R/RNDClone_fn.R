
###################################################################################
# 1. RNDClone: Function for fixed-dimensional MCMC
###################################################################################
RNDClone_PT = function(n, N, m, M, C, 
  g_fun = NULL, K_min = 1, K_max = 3,
  niter = 5000, burnin = 20000, thin = 2, Delta = 1.15^(9:0),
  a_w = 1, b_w = 1, d = 1, d0 = 0.03, a_lambda = 1, b_lambda = 1, 
  a_pai = NULL, b_pai = NULL, a_zeta = NULL, b_zeta = NULL,
  a_gamma_D = 1, b_gamma_D = NULL, a_nu_D = 1, b_nu_D = NULL,
  a_gamma_R = 1, b_gamma_R = NULL, a_nu_R = 1, b_nu_R = NULL,
  a_phi = NULL, b_phi = NULL, a_psi = NULL, b_psi = NULL,
  verbose = FALSE){
  
  # Input:
  # niter: number of iterations
  # burnin: number of burn-in iterations to initialize chains for each temprature
  # thin: 
  # n, N, m, M: data, each one is a S*T matrix
  # g_fun: from 1 to G. Need to use (g_fun - 1) when called by C functions, 
  #        because vectors in C start from 0
  # C: number of subclones
  # Delta: PT temperatures


  cat("RNDClone MCMC (with parallel tempering, fixed C) started.\n")
  cat(sprintf("Date: %s.\n\n", date()))
  
  S = dim(n)[1]
  T = dim(n)[2]
  
  if(is.null(g_fun)) g_fun = 1:S

  G = max(g_fun)
  
  
  ############################################################################
  # Setting Hyperparameters
  ############################################################################
  if(is.null(b_phi)) b_phi = rep(10, T)
  if(is.null(b_psi)) b_psi = rep(10, T)
  if(is.null(a_phi)) a_phi = colMeans(N) * b_phi
  if(is.null(a_psi)) a_psi = colMeans(M) * b_psi
  if(is.null(b_gamma_D)) b_gamma_D = a_gamma_D * 10 * mean(N)
  if(is.null(b_gamma_R)) b_gamma_R = a_gamma_R * 10 * mean(M)
  if(is.null(b_nu_D)) b_nu_D = a_nu_D * 10 * mean(N)
  if(is.null(b_nu_R)) b_nu_R = a_nu_R * 10 * mean(M)
  
  if(is.null(a_pai) | is.null(b_pai)){
    a_pai = C - 1
    b_pai = 1
  }
  
  if(is.null(a_zeta) | is.null(b_zeta)){
    a_zeta = C - 1
    b_zeta = 1
  }
  
  if (verbose) {
    print("a_phi = ")
    print(a_phi)
    print("b_phi = ")
    print(b_phi)
    print("a_psi = ")
    print(a_psi)
    print("b_psi = ")
    print(b_psi)
  }
  
  
  
  # Delta = c(xx, xx, ..., 1), last one is 1
  nThread = length(Delta)
  
  
  ############################################################################
  # Setting storage space for two result lists
  # PT_state_list and sample_list
  ############################################################################
  # the samples we are interested, only the last one with temp = 1
  L_spls = array(0, c(S, C, niter))
  Z_spls = array(0, c(S, C, niter))
  W_spls = array(0, c(T, C+1, niter))
  Lambda_spls = array(0, c(G, C, niter))
  
  pai_spls = array(0, c(C, niter))
  zeta_spls = array(0, c(C, niter))
  
  phi_spls = array(0, c(T, niter))
  psi_spls = array(0, c(T, niter))
  gamma_D_spls = array(0, c(T, niter))
  nu_D_spls = array(0, c(T, niter))
  gamma_R_spls = array(0, c(T, niter))
  nu_R_spls = array(0, c(T, niter))
  
  l_D0_spls = rep(0, niter)
  z_D0_spls = rep(0, niter)
  l_R0_spls = rep(0, niter)
  z_R0_spls = rep(0, niter)
  
  loglik_spls = rep(0, niter)
  logpost_spls = rep(0, niter)
  
  
  # the samples we are not interested in, other tempratures. 
  # So only save their current state
  L_PTR = array(2, c(S, C, nThread))
  Z_PTR = array(0, c(S, C, nThread))
  W_PTR = array(0, c(T, C+1, nThread))
  W_PTR[ , 1:C, ] = 0.99 / C
  W_PTR[ , C+1, ] = 0.01
  Lambda_PTR = array(a_lambda / b_lambda, c(G, C, nThread))
  
  pai_PTR = array(a_pai / (a_pai + b_pai), c(C, nThread))
  zeta_PTR = array(a_zeta / (a_zeta + b_zeta), c(C, nThread))
  
  phi_PTR = matrix(rep(a_phi / b_phi, nThread), T, nThread)
  psi_PTR = matrix(rep(a_psi / b_psi, nThread), T, nThread)
  gamma_D_PTR = matrix(a_gamma_D / b_gamma_D, T, nThread)
  nu_D_PTR = matrix(a_nu_D / b_nu_D, T, nThread)
  gamma_R_PTR = matrix(a_gamma_R / b_gamma_R, T, nThread)
  nu_R_PTR = matrix(a_nu_R / b_nu_R, T, nThread)
  
  l_D0_PTR = rep((K_min + K_max) / 2, nThread)
  z_D0_PTR = rep((K_min + K_max) / 4, nThread)
  l_R0_PTR = rep((K_min + K_max) / 2, nThread)
  z_R0_PTR = rep((K_min + K_max) / 4, nThread)
  
  loglik_PTR = rep(0, nThread)
  logpost_PTR = rep(0, nThread)

  J1 = as.integer(burnin / 4)
  ############################################################################
  # 1. burn-in for DNA
  ############################################################################    
  output = .C("DClone_MCMC", L_PTR = as.integer(L_PTR),
              Z_PTR = as.integer(Z_PTR),
              W_PTR = as.double(W_PTR),
              pai_PTR = as.double(pai_PTR),
              zeta_PTR = as.double(zeta_PTR),
              phi_PTR = as.double(phi_PTR),
              gamma_D_PTR = as.double(gamma_D_PTR),
              nu_D_PTR = as.double(nu_D_PTR),
              l_D0_PTR = as.double(l_D0_PTR),
              z_D0_PTR = as.double(z_D0_PTR),
              loglik_PTR = as.double(loglik_PTR),
              logpost_PTR = as.double(logpost_PTR),
              N = as.double(N), n = as.double(n),
              S = as.integer(S), T = as.integer(T), 
              C = as.integer(C), K_min = as.integer(K_min), K_max = as.integer(K_max),
              a_w = as.double(a_w), b_w = as.double(b_w), 
              d = as.double(d), d0 = as.double(d0),
              a_pai = as.double(a_pai), b_pai = as.double(b_pai), 
              a_zeta = as.double(a_zeta), b_zeta = as.double(b_zeta),
              a_phi = as.double(a_phi), b_phi = as.double(b_phi),
              a_gamma_D = as.double(a_gamma_D), b_gamma_D = as.double(b_gamma_D),
              a_nu_D = as.double(a_nu_D), b_nu_D = as.double(b_nu_D),
              niter = as.integer(J1), Delta = as.double(Delta), 
              nThread = as.integer(nThread))
  
  L_PTR = array(output$L_PTR, c(S, C, nThread))
  Z_PTR = array(output$Z_PTR, c(S, C, nThread))
  W_PTR = array(output$W_PTR, c(T, C+1, nThread))
  
  pai_PTR = matrix(output$pai_PTR, C, nThread)
  zeta_PTR = matrix(output$zeta_PTR, C, nThread)
  
  phi_PTR = matrix(output$phi_PTR, T, nThread)
  gamma_D_PTR = matrix(output$gamma_D_PTR, T, nThread)
  nu_D_PTR = matrix(output$nu_D_PTR, T, nThread)
  l_D0_PTR = output$l_D0_PTR
  z_D0_PTR = output$z_D0_PTR
  loglik_PTR = output$loglik_PTR
  logpost_PTR = output$logpost_PTR
  
  L_PTR[] = L_PTR[ , , nThread]
  Z_PTR[] = Z_PTR[ , , nThread]
  W_PTR[] = W_PTR[ , , nThread]
  
  pai_PTR[] = pai_PTR[ , nThread]
  zeta_PTR[] = zeta_PTR[ , nThread]
  
  phi_PTR[] = phi_PTR[ , nThread]
  gamma_D_PTR[] = gamma_D_PTR[ , nThread]
  nu_D_PTR[] = nu_D_PTR[ , nThread]
  l_D0_PTR[] = l_D0_PTR[nThread]
  z_D0_PTR[] = z_D0_PTR[nThread]
  loglik_PTR[] = loglik_PTR[nThread]
  logpost_PTR[] = logpost_PTR[nThread]
  
  J2 = as.integer(burnin / 4)
  ############################################################################
  # 2. burn-in for RNA
  ############################################################################    
  output = .C("RClone_MCMC_fixed_W", L = as.integer(L_PTR[ , , nThread]),
              Z = as.integer(Z_PTR[ , , nThread]),
              W = as.double(W_PTR[ , , nThread]),
              Lambda_PTR = as.double(Lambda_PTR),
              psi_PTR = as.double(psi_PTR),
              gamma_R_PTR = as.double(gamma_R_PTR),
              nu_R_PTR = as.double(nu_R_PTR),
              l_R0_PTR = as.double(l_R0_PTR),
              z_R0_PTR = as.double(z_R0_PTR),
              loglik_PTR = as.double(loglik_PTR),
              logpost_PTR = as.double(logpost_PTR),
              M = as.double(M), m = as.double(m),
              S = as.integer(S), T = as.integer(T), 
              C = as.integer(C), G = as.integer(G),
              K_min = as.integer(K_min), K_max = as.integer(K_max),
              a_lambda = as.double(a_lambda), b_lambda = as.double(b_lambda), 
              a_psi = as.double(a_psi), b_psi = as.double(b_psi),
              a_gamma_R = as.double(a_gamma_R), b_gamma_R = as.double(b_gamma_R),
              a_nu_R = as.double(a_nu_R), b_nu_R = as.double(b_nu_R),
              niter = as.integer(J2), g_fun = as.integer(g_fun - 1),
              Delta = as.double(Delta), nThread = as.integer(nThread))
  
  Lambda_PTR = array(output$Lambda_PTR, c(G, C, nThread))
  psi_PTR = matrix(output$psi_PTR, T, nThread)
  gamma_R_PTR = matrix(output$gamma_R_PTR, T, nThread)
  nu_R_PTR = matrix(output$nu_R_PTR, T, nThread)
  l_R0_PTR = output$l_R0_PTR
  z_R0_PTR = output$z_R0_PTR
  loglik_PTR = output$loglik_PTR
  logpost_PTR = output$logpost_PTR
  
  Lambda_PTR[] = Lambda_PTR[ , , nThread]
  psi_PTR[] = psi_PTR[ , nThread]
  gamma_R_PTR[] = gamma_R_PTR[ , nThread]
  nu_R_PTR[] = nu_R_PTR[ , nThread]
  l_R0_PTR[] = l_R0_PTR[nThread]
  z_R0_PTR[] = z_R0_PTR[nThread]
  loglik_PTR[] = loglik_PTR[nThread]
  logpost_PTR[] = logpost_PTR[nThread]
  
  J0 = as.integer(burnin / 2)
  ############################################################################
  # 3. burn-in for DNA+RNA
  ############################################################################
  output = .C("RNDClone_MCMC", L_PTR = as.integer(L_PTR),
              Z_PTR = as.integer(Z_PTR),
              W_PTR = as.double(W_PTR),
              Lambda_PTR = as.double(Lambda_PTR),
              pai_PTR = as.double(pai_PTR),
              zeta_PTR = as.double(zeta_PTR),
              phi_PTR = as.double(phi_PTR),
              psi_PTR = as.double(psi_PTR),
              gamma_D_PTR = as.double(gamma_D_PTR),
              nu_D_PTR = as.double(nu_D_PTR),
              gamma_R_PTR = as.double(gamma_R_PTR),
              nu_R_PTR = as.double(nu_R_PTR),
              l_D0_PTR = as.double(l_D0_PTR),
              z_D0_PTR = as.double(z_D0_PTR),
              l_R0_PTR = as.double(l_R0_PTR),
              z_R0_PTR = as.double(z_R0_PTR),
              loglik_PTR = as.double(loglik_PTR),
              logpost_PTR = as.double(logpost_PTR),
              N = as.double(N), n = as.double(n),
              M = as.double(M), m = as.double(m),
              S = as.integer(S), T = as.integer(T), 
              C = as.integer(C), G = as.integer(G),
              K_min = as.integer(K_min), K_max = as.integer(K_max),
              a_w = as.double(a_w), b_w = as.double(b_w), 
              d = as.double(d), d0 = as.double(d0),
              a_lambda = as.double(a_lambda), b_lambda = as.double(b_lambda), 
              a_pai = as.double(a_pai), b_pai = as.double(b_pai), 
              a_zeta = as.double(a_zeta), b_zeta = as.double(b_zeta),
              a_phi = as.double(a_phi), b_phi = as.double(b_phi),
              a_psi = as.double(a_psi), b_psi = as.double(b_psi),
              a_gamma_D = as.double(a_gamma_D), b_gamma_D = as.double(b_gamma_D),
              a_gamma_R = as.double(a_gamma_R), b_gamma_R = as.double(b_gamma_R),
              a_nu_D = as.double(a_nu_D), b_nu_D = as.double(b_nu_D),
              a_nu_R = as.double(a_nu_R), b_nu_R = as.double(b_nu_R),
              niter = as.integer(J0), g_fun = as.integer(g_fun - 1),
              Delta = as.double(Delta), nThread = as.integer(nThread))
  
  L_PTR = array(output$L_PTR, c(S, C, nThread))
  Z_PTR = array(output$Z_PTR, c(S, C, nThread))
  W_PTR = array(output$W_PTR, c(T, C+1, nThread))
  Lambda_PTR = array(output$Lambda_PTR, c(G, C, nThread))
  
  pai_PTR = matrix(output$pai_PTR, C, nThread)
  zeta_PTR = matrix(output$zeta_PTR, C, nThread)
  
  phi_PTR = matrix(output$phi_PTR, T, nThread)
  psi_PTR = matrix(output$psi_PTR, T, nThread)
  gamma_D_PTR = matrix(output$gamma_D_PTR, T, nThread)
  nu_D_PTR = matrix(output$nu_D_PTR, T, nThread)
  gamma_R_PTR = matrix(output$gamma_R_PTR, T, nThread)
  nu_R_PTR = matrix(output$nu_R_PTR, T, nThread)
  l_D0_PTR = output$l_D0_PTR
  z_D0_PTR = output$z_D0_PTR
  l_R0_PTR = output$l_R0_PTR
  z_R0_PTR = output$z_R0_PTR
  loglik_PTR = output$loglik_PTR
  logpost_PTR = output$logpost_PTR
  
  L_spls[ , , 1] = L_PTR[ , , nThread]
  Z_spls[ , , 1] = Z_PTR[ , , nThread]
  W_spls[ , , 1] = W_PTR[ , , nThread]
  Lambda_spls[ , , 1] = Lambda_PTR[ , , nThread]
  
  pai_spls[ , 1] = pai_PTR[ , nThread]
  zeta_spls[ , 1] = zeta_PTR[ , nThread]
  
  phi_spls[ , 1] = phi_PTR[ , nThread]
  psi_spls[ , 1] = psi_PTR[ , nThread]
  gamma_D_spls[ , 1] = gamma_D_PTR[ , nThread]
  nu_D_spls[ , 1] = nu_D_PTR[ , nThread]
  gamma_R_spls[ , 1] = gamma_R_PTR[ , nThread]
  nu_R_spls[ , 1] = nu_R_PTR[ , nThread]
  l_D0_spls[1] = l_D0_PTR[nThread]
  z_D0_spls[1] = z_D0_PTR[nThread]
  l_R0_spls[1] = l_R0_PTR[nThread]
  z_R0_spls[1] = z_R0_PTR[nThread]

  loglik_spls[1] = loglik_PTR[nThread]
  logpost_spls[1] = logpost_PTR[nThread]
  
  cat(sprintf("RNDClone MCMC burn-in finished. Date: %s.\n", date()))
  
  
  ############################################################################
  # Start MCMC iteration
  ############################################################################
  for(iter in 2:niter){
    
    output = .C("RNDClone_MCMC", L_PTR = as.integer(L_PTR),
                Z_PTR = as.integer(Z_PTR),
                W_PTR = as.double(W_PTR),
                Lambda_PTR = as.double(Lambda_PTR),
                pai_PTR = as.double(pai_PTR),
                zeta_PTR = as.double(zeta_PTR),
                phi_PTR = as.double(phi_PTR),
                psi_PTR = as.double(psi_PTR),
                gamma_D_PTR = as.double(gamma_D_PTR),
                nu_D_PTR = as.double(nu_D_PTR),
                gamma_R_PTR = as.double(gamma_R_PTR),
                nu_R_PTR = as.double(nu_R_PTR),
                l_D0_PTR = as.double(l_D0_PTR),
                z_D0_PTR = as.double(z_D0_PTR),
                l_R0_PTR = as.double(l_R0_PTR),
                z_R0_PTR = as.double(z_R0_PTR),
                loglik_PTR = as.double(loglik_PTR),
                logpost_PTR = as.double(logpost_PTR),
                N = as.double(N), n = as.double(n),
                M = as.double(M), m = as.double(m),
                S = as.integer(S), T = as.integer(T), 
                C = as.integer(C), G = as.integer(G),
                K_min = as.integer(K_min), K_max = as.integer(K_max),
                a_w = as.double(a_w), b_w = as.double(b_w), 
                d = as.double(d), d0 = as.double(d0),
                a_lambda = as.double(a_lambda), b_lambda = as.double(b_lambda), 
                a_pai = as.double(a_pai), b_pai = as.double(b_pai), 
                a_zeta = as.double(a_zeta), b_zeta = as.double(b_zeta),
                a_phi = as.double(a_phi), b_phi = as.double(b_phi),
                a_psi = as.double(a_psi), b_psi = as.double(b_psi),
                a_gamma_D = as.double(a_gamma_D), b_gamma_D = as.double(b_gamma_D),
                a_gamma_R = as.double(a_gamma_R), b_gamma_R = as.double(b_gamma_R),
                a_nu_D = as.double(a_nu_D), b_nu_D = as.double(b_nu_D),
                a_nu_R = as.double(a_nu_R), b_nu_R = as.double(b_nu_R),
                niter = as.integer(thin), g_fun = as.integer(g_fun - 1),
                Delta = as.double(Delta), nThread = as.integer(nThread))
    
    L_PTR = array(output$L_PTR, c(S, C, nThread))
    Z_PTR = array(output$Z_PTR, c(S, C, nThread))
    W_PTR = array(output$W_PTR, c(T, C+1, nThread))
    Lambda_PTR = array(output$Lambda_PTR, c(G, C, nThread))
    
    pai_PTR = matrix(output$pai_PTR, C, nThread)
    zeta_PTR = matrix(output$zeta_PTR, C, nThread)
    
    phi_PTR = matrix(output$phi_PTR, T, nThread)
    psi_PTR = matrix(output$psi_PTR, T, nThread)
    gamma_D_PTR = matrix(output$gamma_D_PTR, T, nThread)
    nu_D_PTR = matrix(output$nu_D_PTR, T, nThread)
    gamma_R_PTR = matrix(output$gamma_R_PTR, T, nThread)
    nu_R_PTR = matrix(output$nu_R_PTR, T, nThread)
    l_D0_PTR = output$l_D0_PTR
    z_D0_PTR = output$z_D0_PTR
    l_R0_PTR = output$l_R0_PTR
    z_R0_PTR = output$z_R0_PTR
    loglik_PTR = output$loglik_PTR
    logpost_PTR = output$logpost_PTR
    
    L_spls[ , , iter] = L_PTR[ , , nThread]
    Z_spls[ , , iter] = Z_PTR[ , , nThread]
    W_spls[ , , iter] = W_PTR[ , , nThread]
    Lambda_spls[ , , iter] = Lambda_PTR[ , , nThread]
    
    pai_spls[ , iter] = pai_PTR[ , nThread]
    zeta_spls[ , iter] = zeta_PTR[ , nThread]
    
    phi_spls[ , iter] = phi_PTR[ , nThread]
    psi_spls[ , iter] = psi_PTR[ , nThread]
    gamma_D_spls[ , iter] = gamma_D_PTR[ , nThread]
    nu_D_spls[ , iter] = nu_D_PTR[ , nThread]
    gamma_R_spls[ , iter] = gamma_R_PTR[ , nThread]
    nu_R_spls[ , iter] = nu_R_PTR[ , nThread]
    l_D0_spls[iter] = l_D0_PTR[nThread]
    z_D0_spls[iter] = z_D0_PTR[nThread]
    l_R0_spls[iter] = l_R0_PTR[nThread]
    z_R0_spls[iter] = z_R0_PTR[nThread]

    loglik_spls[iter] = loglik_PTR[nThread]
    logpost_spls[iter] = logpost_PTR[nThread]
    
    if((iter %% 500) == 0){
      cat(sprintf("RNDClone MCMC %.2f %% finished. ", round(i/niter*100, 2)))
      cat(sprintf("Date: %s.\n", date()))
    }
  }
  
  
  
  sample_list = list()
  sample_list$L_spls = L_spls
  sample_list$Z_spls = Z_spls
  sample_list$W_spls = W_spls
  sample_list$Lambda_spls = Lambda_spls
  sample_list$pai_spls = pai_spls
  sample_list$zeta_spls = zeta_spls
  sample_list$phi_spls = phi_spls
  sample_list$psi_spls = psi_spls
  sample_list$gamma_D_spls = gamma_D_spls
  sample_list$nu_D_spls = nu_D_spls
  sample_list$gamma_R_spls = gamma_R_spls
  sample_list$nu_R_spls = nu_R_spls
  sample_list$l_D0_spls = l_D0_spls
  sample_list$z_D0_spls = z_D0_spls
  sample_list$l_R0_spls = l_R0_spls
  sample_list$z_R0_spls = z_R0_spls
  sample_list$loglik_spls = loglik_spls
  sample_list$logpost_spls = logpost_spls
  
  PT_state_list = list()
  PT_state_list$L_PTR = L_PTR
  PT_state_list$Z_PTR = Z_PTR
  PT_state_list$W_PTR = W_PTR
  PT_state_list$Lambda_PTR = Lambda_PTR
  PT_state_list$pai_PTR = pai_PTR
  PT_state_list$zeta_PTR = zeta_PTR
  PT_state_list$phi_PTR = phi_PTR
  PT_state_list$psi_PTR = psi_PTR
  PT_state_list$gamma_D_PTR = gamma_D_PTR
  PT_state_list$nu_D_PTR = nu_D_PTR
  PT_state_list$gamma_R_PTR = gamma_R_PTR
  PT_state_list$nu_R_PTR = nu_R_PTR
  PT_state_list$l_D0_PTR = l_D0_PTR
  PT_state_list$z_D0_PTR = z_D0_PTR
  PT_state_list$l_R0_PTR = l_R0_PTR
  PT_state_list$z_R0_PTR = z_R0_PTR
  PT_state_list$loglik_PTR = loglik_PTR
  PT_state_list$logpost_PTR = logpost_PTR
  
  result_list = list()
  result_list$sample_list = sample_list
  result_list$PT_state_list = PT_state_list
  
  cat("\nRNDClone MCMC (with parallel tempering, fixed C) finished.\n")
  cat(sprintf("Date: %s.\n\n", date()))
  
  return(result_list)
  
}
















###################################################################################
# 2. RNDClone: Function for Trans-dimensional (Reversible Jump) MCMC
###################################################################################
RNDClone_RJMCMC = function(n, N, m, M, 
  C_min = 2, C_max = 7, 
  g_fun = NULL, K_min = 1, K_max = 3,
  niter = 5000, burnin = 20000, thin = 2, Delta = 1.15^(9:0), tau = 0.99,
  alpha = 0.8, a_w = 1, b_w = 1, d = 1, d0 = 0.03, a_lambda = 1, b_lambda = 1, 
  a_pai = NULL, b_pai = NULL, a_zeta = NULL, b_zeta = NULL,
  a_gamma_D = 1, b_gamma_D = NULL, a_nu_D = 1, b_nu_D = NULL,
  a_gamma_R = 1, b_gamma_R = NULL, a_nu_R = 1, b_nu_R = NULL,
  a_phi = NULL, b_phi = NULL, a_psi = NULL, b_psi = NULL,
  verbose = FALSE){
  
  # tau: power prior

  cat("RNDClone RJ-MCMC (with parallel tempering, jump across different C) started.\n")
  cat(sprintf("C_min = %d, C_max = %d. Date: %s.\n\n", C_min, C_max, date()))
  
  S = dim(n)[1]
  T = dim(n)[2]

  if(is.null(g_fun)) g_fun = 1:S
  
  G = max(g_fun)
  
  ############################################################################
  # Setting Hyperparameters
  ############################################################################
  if(is.null(b_phi)) b_phi = rep(10, T)
  if(is.null(b_psi)) b_psi = rep(10, T)
  if(is.null(a_phi)) a_phi = colMeans(N) * b_phi
  if(is.null(a_psi)) a_psi = colMeans(M) * b_psi
  if(is.null(b_gamma_D)) b_gamma_D = a_gamma_D * 10 * mean(N)
  if(is.null(b_gamma_R)) b_gamma_R = a_gamma_R * 10 * mean(M)
  if(is.null(b_nu_D)) b_nu_D = a_nu_D * 10 * mean(N)
  if(is.null(b_nu_R)) b_nu_R = a_nu_R * 10 * mean(M)
  
  if (verbose) {
    print("a_phi = ")
    print(a_phi)
    print("b_phi = ")
    print(b_phi)
    print("a_psi = ")
    print(a_psi)
    print("b_psi = ")
    print(b_psi)
  }
  
  # Delta = c(xx, xx, ..., 1), last one is 1
  nThread = length(Delta)
  
  # thus, the last temperature corresponds to power prior with power tau
  Delta_tau = Delta / tau
  
  
  ############################################################################
  # Setting storage space for two result lists
  # PT_AC_state_list and sample_list
  ############################################################################
  # save the current state for all PT temperatures (PT) and all C (AC) in C_min -- C_max
  PT_AC_state_list = list()
  PT_AC_state_list$L_PTR_AC = list()
  PT_AC_state_list$Z_PTR_AC = list()
  PT_AC_state_list$W_PTR_AC = list()
  PT_AC_state_list$Lambda_PTR_AC = list()
  PT_AC_state_list$pai_PTR_AC = list()
  PT_AC_state_list$zeta_PTR_AC = list()
  PT_AC_state_list$phi_PTR_AC = array(0, c(T, nThread, C_max))
  PT_AC_state_list$psi_PTR_AC = array(0, c(T, nThread, C_max))
  PT_AC_state_list$gamma_D_PTR_AC = array(0, c(T, nThread, C_max))
  PT_AC_state_list$nu_D_PTR_AC = array(0, c(T, nThread, C_max))
  PT_AC_state_list$gamma_R_PTR_AC = array(0, c(T, nThread, C_max))
  PT_AC_state_list$nu_R_PTR_AC = array(0, c(T, nThread, C_max))
  PT_AC_state_list$l_D0_PTR_AC = array(0, c(nThread, C_max))
  PT_AC_state_list$z_D0_PTR_AC = array(0, c(nThread, C_max))
  PT_AC_state_list$l_R0_PTR_AC = array(0, c(nThread, C_max))
  PT_AC_state_list$z_R0_PTR_AC = array(0, c(nThread, C_max))
  PT_AC_state_list$loglik_PTR_AC = array(0, c(nThread, C_max))
  PT_AC_state_list$logpost_PTR_AC = array(0, c(nThread, C_max))
  
  # save all RJ samples
  # each sample_list$x_spls is a list with niter elements
  # x_spls[[iter1]] probably has different dimension with x_spls[[iter2]]
  sample_list = list()
  sample_list$L_spls = list()
  sample_list$Z_spls = list()
  sample_list$W_spls = list()
  sample_list$Lambda_spls = list()
  sample_list$pai_spls = list()
  sample_list$zeta_spls = list()
  
  sample_list$phi_spls = array(0, c(T, niter))
  sample_list$psi_spls = array(0, c(T, niter))
  sample_list$gamma_D_spls = array(0, c(T, niter))
  sample_list$nu_D_spls = array(0, c(T, niter))
  sample_list$gamma_R_spls = array(0, c(T, niter))
  sample_list$nu_R_spls = array(0, c(T, niter))
  
  sample_list$l_D0_spls = rep(0, niter)
  sample_list$z_D0_spls = rep(0, niter)
  sample_list$l_R0_spls = rep(0, niter)
  sample_list$z_R0_spls = rep(0, niter)
  
  sample_list$loglik_spls = rep(0, niter)
  sample_list$logpost_spls = rep(0, niter)
  
  sample_list$C_spls = rep(0, niter)
  
  J1 = as.integer(burnin / 3)
  J2 = as.integer(burnin / 3)
  J0 = as.integer(burnin / 3)
  
  ############################################################################
  # Burnin for each C in C_min:C_max
  ############################################################################
  for(C_iter in C_min:C_max){
    
    cat(sprintf("RNDClone RJ-MCMC burn-in started for C = %d. Date: %s.\n", C_iter, date()))
    
    if(is.null(a_pai) | is.null(b_pai)){
      a_pai1 = C_iter - 1
      b_pai1 = 1
    } else{
      a_pai1 = a_pai
      b_pai1 = b_pai
    }
    
    if(is.null(a_zeta) | is.null(b_zeta)){
      a_zeta1 = C_iter - 1
      b_zeta1 = 1
    } else{
      a_zeta1 = a_zeta
      b_zeta1 = b_zeta
    }
    
    ## initialize
    PT_AC_state_list$L_PTR_AC[[C_iter]] = array(2, c(S, C_iter, nThread))
    PT_AC_state_list$Z_PTR_AC[[C_iter]] = array(0, c(S, C_iter, nThread))
    PT_AC_state_list$W_PTR_AC[[C_iter]] = array(0, c(T, C_iter + 1, nThread))
    PT_AC_state_list$W_PTR_AC[[C_iter]][ , 1:C_iter, ] = 0.99 / C_iter
    PT_AC_state_list$W_PTR_AC[[C_iter]][ , C_iter+1, ] = 0.01
    PT_AC_state_list$Lambda_PTR_AC[[C_iter]] = array(a_lambda / b_lambda, c(G, C_iter, nThread))
    
    PT_AC_state_list$pai_PTR_AC[[C_iter]] = array(a_pai1 / (a_pai1 + b_pai1), c(C_iter, nThread))
    PT_AC_state_list$zeta_PTR_AC[[C_iter]] = array(a_zeta1 / (a_zeta1 + b_zeta1), c(C_iter, nThread))
    
    PT_AC_state_list$phi_PTR_AC[ , , C_iter] = matrix(rep(a_phi / b_phi, nThread), T, nThread)
    PT_AC_state_list$psi_PTR_AC[ , , C_iter] = matrix(rep(a_psi / b_psi, nThread), T, nThread)
    PT_AC_state_list$gamma_D_PTR_AC[ , , C_iter] = a_gamma_D / b_gamma_D
    PT_AC_state_list$nu_D_PTR_AC[ , , C_iter] = a_nu_D / b_nu_D
    PT_AC_state_list$gamma_R_PTR_AC[ , , C_iter] = a_gamma_R / b_gamma_R
    PT_AC_state_list$nu_R_PTR_AC[ , , C_iter] = a_nu_R / b_nu_R
    PT_AC_state_list$l_D0_PTR_AC[ , C_iter] = (K_min + K_max) / 2
    PT_AC_state_list$z_D0_PTR_AC[ , C_iter] = (K_min + K_max) / 4
    PT_AC_state_list$l_R0_PTR_AC[ , C_iter] = (K_min + K_max) / 2
    PT_AC_state_list$z_R0_PTR_AC[ , C_iter] = (K_min + K_max) / 4
    # PT_AC_state_list$loglik_PTR_AC[ , C_iter] = rep(0, nThread)
    # PT_AC_state_list$logpost_PTR_AC[ , C_iter] = rep(0, nThread)
    
    
    # Burnin using DNA
    output = .C("DClone_MCMC", L_PTR = as.integer(PT_AC_state_list$L_PTR_AC[[C_iter]]),
              Z_PTR = as.integer(PT_AC_state_list$Z_PTR_AC[[C_iter]]),
              W_PTR = as.double(PT_AC_state_list$W_PTR_AC[[C_iter]]),
              pai_PTR = as.double(PT_AC_state_list$pai_PTR_AC[[C_iter]]),
              zeta_PTR = as.double(PT_AC_state_list$zeta_PTR_AC[[C_iter]]),
              phi_PTR = as.double(PT_AC_state_list$phi_PTR_AC[ , , C_iter]),
              gamma_D_PTR = as.double(PT_AC_state_list$gamma_D_PTR_AC[ , , C_iter]),
              nu_D_PTR = as.double(PT_AC_state_list$nu_D_PTR_AC[ , , C_iter]),
              l_D0_PTR = as.double(PT_AC_state_list$l_D0_PTR_AC[ , C_iter]),
              z_D0_PTR = as.double(PT_AC_state_list$z_D0_PTR_AC[ , C_iter]),
              loglik_PTR = as.double(rep(0, nThread)),
              logpost_PTR = as.double(rep(0, nThread)),
              N = as.double(N), n = as.double(n),
              S = as.integer(S), T = as.integer(T), 
              C = as.integer(C_iter), K_min = as.integer(K_min), K_max = as.integer(K_max),
              a_w = as.double(a_w), b_w = as.double(b_w), 
              d = as.double(d), d0 = as.double(d0),
              a_pai = as.double(a_pai1), b_pai = as.double(b_pai1), 
              a_zeta = as.double(a_zeta1), b_zeta = as.double(b_zeta1),
              a_phi = as.double(a_phi), b_phi = as.double(b_phi),
              a_gamma_D = as.double(a_gamma_D), b_gamma_D = as.double(b_gamma_D),
              a_nu_D = as.double(a_nu_D), b_nu_D = as.double(b_nu_D),
              niter = as.integer(J1), Delta = as.double(Delta_tau), 
              nThread = as.integer(nThread))
    
    PT_AC_state_list$L_PTR_AC[[C_iter]] = array(output$L_PTR, c(S, C_iter, nThread))
    PT_AC_state_list$Z_PTR_AC[[C_iter]] = array(output$Z_PTR, c(S, C_iter, nThread))
    PT_AC_state_list$W_PTR_AC[[C_iter]] = array(output$W_PTR, c(T, C_iter + 1, nThread))
    
    PT_AC_state_list$pai_PTR_AC[[C_iter]] = matrix(output$pai_PTR, C_iter, nThread)
    PT_AC_state_list$zeta_PTR_AC[[C_iter]] = matrix(output$zeta_PTR, C_iter, nThread)
    
    PT_AC_state_list$phi_PTR_AC[ , , C_iter] = matrix(output$phi_PTR, T, nThread)
    PT_AC_state_list$gamma_D_PTR_AC[ , , C_iter] = matrix(output$gamma_D_PTR, T, nThread)
    PT_AC_state_list$nu_D_PTR_AC[ , , C_iter] = matrix(output$nu_D_PTR, T, nThread)
    PT_AC_state_list$l_D0_PTR_AC[ , C_iter] = output$l_D0_PTR
    PT_AC_state_list$z_D0_PTR_AC[ , C_iter] = output$z_D0_PTR
    # PT_AC_state_list$loglik_PTR_AC[ , C_iter] = output$loglik_PTR
    # PT_AC_state_list$logpost_PTR_AC[ , C_iter] = output$logpost_PTR
  
    PT_AC_state_list$L_PTR_AC[[C_iter]][] = PT_AC_state_list$L_PTR_AC[[C_iter]][ , , nThread]
    PT_AC_state_list$Z_PTR_AC[[C_iter]][] = PT_AC_state_list$Z_PTR_AC[[C_iter]][ , , nThread]
    PT_AC_state_list$W_PTR_AC[[C_iter]][] = PT_AC_state_list$W_PTR_AC[[C_iter]][ , , nThread]
  
    PT_AC_state_list$pai_PTR_AC[[C_iter]][] = PT_AC_state_list$pai_PTR_AC[[C_iter]][ , nThread]
    PT_AC_state_list$zeta_PTR_AC[[C_iter]][] = PT_AC_state_list$zeta_PTR_AC[[C_iter]][ , nThread]
  
    PT_AC_state_list$phi_PTR_AC[ , , C_iter] = PT_AC_state_list$phi_PTR_AC[ , nThread, C_iter]
    PT_AC_state_list$gamma_D_PTR_AC[ , , C_iter] = PT_AC_state_list$gamma_D_PTR_AC[ , nThread, C_iter]
    PT_AC_state_list$nu_D_PTR_AC[ , , C_iter] = PT_AC_state_list$nu_D_PTR_AC[ , nThread, C_iter]
    PT_AC_state_list$l_D0_PTR_AC[ , C_iter] = PT_AC_state_list$l_D0_PTR_AC[nThread, C_iter]
    PT_AC_state_list$z_D0_PTR_AC[ , C_iter] = PT_AC_state_list$z_D0_PTR_AC[nThread, C_iter]
    # PT_AC_state_list$loglik_PTR_AC[ , C_iter] = PT_AC_state_list$loglik_PTR_AC[nThread, C_iter]
    # PT_AC_state_list$logpost_PTR_AC[ , C_iter] = PT_AC_state_list$logpost_PTR_AC[nThread, C_iter]
    
    if (verbose) {
      cat(sprintf("RNDClone DNA burn-in finished for C = %d. Date: %s.\n", C_iter, date()))
    }
    
    # Burnin using RNA
    output = .C("RClone_MCMC_fixed_W", L = as.integer(PT_AC_state_list$L_PTR_AC[[C_iter]][ , , nThread]),
              Z = as.integer(PT_AC_state_list$Z_PTR_AC[[C_iter]][ , , nThread]),
              W = as.double(PT_AC_state_list$W_PTR_AC[[C_iter]][ , , nThread]),
              Lambda_PTR = as.double(PT_AC_state_list$Lambda_PTR_AC[[C_iter]]),
              psi_PTR = as.double(PT_AC_state_list$psi_PTR_AC[ , , C_iter]),
              gamma_R_PTR = as.double(PT_AC_state_list$gamma_R_PTR_AC[ , , C_iter]),
              nu_R_PTR = as.double(PT_AC_state_list$nu_R_PTR_AC[ , , C_iter]),
              l_R0_PTR = as.double(PT_AC_state_list$l_R0_PTR_AC[ , C_iter]),
              z_R0_PTR = as.double(PT_AC_state_list$z_R0_PTR_AC[ , C_iter]),
              loglik_PTR = as.double(rep(0, nThread)),
              logpost_PTR = as.double(rep(0, nThread)),
              M = as.double(M), m = as.double(m),
              S = as.integer(S), T = as.integer(T), 
              C = as.integer(C_iter), G = as.integer(G),
              K_min = as.integer(K_min), K_max = as.integer(K_max),
              a_lambda = as.double(a_lambda), b_lambda = as.double(b_lambda), 
              a_psi = as.double(a_psi), b_psi = as.double(b_psi),
              a_gamma_R = as.double(a_gamma_R), b_gamma_R = as.double(b_gamma_R),
              a_nu_R = as.double(a_nu_R), b_nu_R = as.double(b_nu_R),
              niter = as.integer(J2), g_fun = as.integer(g_fun - 1),
              Delta = as.double(Delta_tau), nThread = as.integer(nThread))
  
    PT_AC_state_list$Lambda_PTR_AC[[C_iter]] = array(output$Lambda_PTR, c(G, C_iter, nThread))
    PT_AC_state_list$psi_PTR_AC[ , , C_iter] = matrix(output$psi_PTR, T, nThread)
    PT_AC_state_list$gamma_R_PTR_AC[ , , C_iter] = matrix(output$gamma_R_PTR, T, nThread)
    PT_AC_state_list$nu_R_PTR_AC[ , , C_iter] = matrix(output$nu_R_PTR, T, nThread)
    PT_AC_state_list$l_R0_PTR_AC[ , C_iter] = output$l_R0_PTR
    PT_AC_state_list$z_R0_PTR_AC[ , C_iter] = output$z_R0_PTR
    # PT_AC_state_list$loglik_PTR_AC[ , C_iter] = output$loglik_PTR
    # PT_AC_state_list$logpost_PTR_AC[ , C_iter] = output$logpost_PTR
  
    PT_AC_state_list$Lambda_PTR_AC[[C_iter]][] = PT_AC_state_list$Lambda_PTR_AC[[C_iter]][ , , nThread]
    PT_AC_state_list$psi_PTR_AC[ , , C_iter] = PT_AC_state_list$psi_PTR_AC[ , nThread, C_iter]
    PT_AC_state_list$gamma_R_PTR_AC[ , , C_iter] = PT_AC_state_list$gamma_R_PTR_AC[ , nThread, C_iter]
    PT_AC_state_list$nu_R_PTR_AC[ , , C_iter] = PT_AC_state_list$nu_R_PTR_AC[ , nThread, C_iter]
    PT_AC_state_list$l_R0_PTR_AC[ , C_iter] = PT_AC_state_list$l_R0_PTR_AC[nThread, C_iter]
    PT_AC_state_list$z_R0_PTR_AC[ , C_iter] = PT_AC_state_list$z_R0_PTR_AC[nThread, C_iter]
    # PT_AC_state_list$loglik_PTR_AC[ , C_iter] = PT_AC_state_list$loglik_PTR_AC[nThread, C_iter]
    # PT_AC_state_list$logpost_PTR_AC[ , C_iter] = PT_AC_state_list$logpost_PTR_AC[nThread, C_iter]
  
    if (verbose) {
      cat(sprintf("RNDClone RNA burn-in finished for C = %d. Date: %s.\n", C_iter, date()))
    }    
    # Burnin using DNA + RNA
    output = .C("RNDClone_MCMC", L_PTR = as.integer(PT_AC_state_list$L_PTR_AC[[C_iter]]),
              Z_PTR = as.integer(PT_AC_state_list$Z_PTR_AC[[C_iter]]),
              W_PTR = as.double(PT_AC_state_list$W_PTR_AC[[C_iter]]),
              Lambda_PTR = as.double(PT_AC_state_list$Lambda_PTR_AC[[C_iter]]),
              pai_PTR = as.double(PT_AC_state_list$pai_PTR_AC[[C_iter]]),
              zeta_PTR = as.double(PT_AC_state_list$zeta_PTR_AC[[C_iter]]),
              phi_PTR = as.double(PT_AC_state_list$phi_PTR_AC[ , , C_iter]),
              psi_PTR = as.double(PT_AC_state_list$psi_PTR_AC[ , , C_iter]),
              gamma_D_PTR = as.double(PT_AC_state_list$gamma_D_PTR_AC[ , , C_iter]),
              nu_D_PTR = as.double(PT_AC_state_list$nu_D_PTR_AC[ , , C_iter]),
              gamma_R_PTR = as.double(PT_AC_state_list$gamma_R_PTR_AC[ , , C_iter]),
              nu_R_PTR = as.double(PT_AC_state_list$nu_R_PTR_AC[ , , C_iter]),
              l_D0_PTR = as.double(PT_AC_state_list$l_D0_PTR_AC[ , C_iter]),
              z_D0_PTR = as.double(PT_AC_state_list$z_D0_PTR_AC[ , C_iter]),
              l_R0_PTR = as.double(PT_AC_state_list$l_R0_PTR_AC[ , C_iter]),
              z_R0_PTR = as.double(PT_AC_state_list$z_R0_PTR_AC[ , C_iter]),
              loglik_PTR = as.double(rep(0, nThread)),
              logpost_PTR = as.double(rep(0, nThread)),
              N = as.double(N), n = as.double(n),
              M = as.double(M), m = as.double(m),
              S = as.integer(S), T = as.integer(T), 
              C = as.integer(C_iter), G = as.integer(G),
              K_min = as.integer(K_min), K_max = as.integer(K_max),
              a_w = as.double(a_w), b_w = as.double(b_w), 
              d = as.double(d), d0 = as.double(d0),
              a_lambda = as.double(a_lambda), b_lambda = as.double(b_lambda), 
              a_pai = as.double(a_pai1), b_pai = as.double(b_pai1), 
              a_zeta = as.double(a_zeta1), b_zeta = as.double(b_zeta1),
              a_phi = as.double(a_phi), b_phi = as.double(b_phi),
              a_psi = as.double(a_psi), b_psi = as.double(b_psi),
              a_gamma_D = as.double(a_gamma_D), b_gamma_D = as.double(b_gamma_D),
              a_gamma_R = as.double(a_gamma_R), b_gamma_R = as.double(b_gamma_R),
              a_nu_D = as.double(a_nu_D), b_nu_D = as.double(b_nu_D),
              a_nu_R = as.double(a_nu_R), b_nu_R = as.double(b_nu_R),
              niter = as.integer(J0), g_fun = as.integer(g_fun - 1),
              Delta = as.double(Delta_tau), nThread = as.integer(nThread))
    
    PT_AC_state_list$L_PTR_AC[[C_iter]] = array(output$L_PTR, c(S, C_iter, nThread))
    PT_AC_state_list$Z_PTR_AC[[C_iter]] = array(output$Z_PTR, c(S, C_iter, nThread))
    PT_AC_state_list$W_PTR_AC[[C_iter]] = array(output$W_PTR, c(T, C_iter + 1, nThread))
    PT_AC_state_list$Lambda_PTR_AC[[C_iter]] = array(output$Lambda_PTR, c(G, C_iter, nThread))
    
    PT_AC_state_list$pai_PTR_AC[[C_iter]] = matrix(output$pai_PTR, C_iter, nThread)
    PT_AC_state_list$zeta_PTR_AC[[C_iter]] = matrix(output$zeta_PTR, C_iter, nThread)
    
    PT_AC_state_list$phi_PTR_AC[ , , C_iter] = matrix(output$phi_PTR, T, nThread)
    PT_AC_state_list$psi_PTR_AC[ , , C_iter] = matrix(output$psi_PTR, T, nThread)
    PT_AC_state_list$gamma_D_PTR_AC[ , , C_iter] = matrix(output$gamma_D_PTR, T, nThread)
    PT_AC_state_list$nu_D_PTR_AC[ , , C_iter] = matrix(output$nu_D_PTR, T, nThread)
    PT_AC_state_list$gamma_R_PTR_AC[ , , C_iter] = matrix(output$gamma_R_PTR, T, nThread)
    PT_AC_state_list$nu_R_PTR_AC[ , , C_iter] = matrix(output$nu_R_PTR, T, nThread)
    PT_AC_state_list$l_D0_PTR_AC[ , C_iter] = output$l_D0_PTR
    PT_AC_state_list$z_D0_PTR_AC[ , C_iter] = output$z_D0_PTR
    PT_AC_state_list$l_R0_PTR_AC[ , C_iter] = output$l_R0_PTR
    PT_AC_state_list$z_R0_PTR_AC[ , C_iter] = output$z_R0_PTR
    
    PT_AC_state_list$loglik_PTR_AC[ , C_iter] = output$loglik_PTR
    PT_AC_state_list$logpost_PTR_AC[ , C_iter] = output$logpost_PTR
    
    if (verbose) {
      print("Log-likelihood after burn-in:")
      print(PT_AC_state_list$loglik_PTR_AC)
      print("Log-posterior after burn-in:")
      print(PT_AC_state_list$logpost_PTR_AC)
    }
    
    cat(sprintf("RNDClone RJ-MCMC burn-in finished for C = %d. Date: %s.\n", C_iter, date()))
  }
  ## End burnin
  
  
  # initial value of C, and loglik on test set (for evaluating p_acc)
  C_cur = sample(C_min:C_max, 1)
  
            
  loglik_cur = (1 - tau) * PT_AC_state_list$loglik_PTR_AC[nThread, C_cur]
  
  logpost_cur = loglik_cur + (C_cur - 1) * log(1 - alpha)
  
  for(i in 1:niter){
    
    C_pro = sample(C_min:C_max, 1)
    
    if(is.null(a_pai) | is.null(b_pai)){
      a_pai1 = C_pro - 1
      b_pai1 = 1
    } else{
      a_pai1 = a_pai
      b_pai1 = b_pai
    }
    
    if(is.null(a_zeta) | is.null(b_zeta)){
      a_zeta1 = C_pro - 1
      b_zeta1 = 1
    } else{
      a_zeta1 = a_zeta
      b_zeta1 = b_zeta
    }
    
    output = .C("RNDClone_MCMC", L_PTR = as.integer(PT_AC_state_list$L_PTR_AC[[C_pro]]),
              Z_PTR = as.integer(PT_AC_state_list$Z_PTR_AC[[C_pro]]),
              W_PTR = as.double(PT_AC_state_list$W_PTR_AC[[C_pro]]),
              Lambda_PTR = as.double(PT_AC_state_list$Lambda_PTR_AC[[C_pro]]),
              pai_PTR = as.double(PT_AC_state_list$pai_PTR_AC[[C_pro]]),
              zeta_PTR = as.double(PT_AC_state_list$zeta_PTR_AC[[C_pro]]),
              phi_PTR = as.double(PT_AC_state_list$phi_PTR_AC[ , , C_pro]),
              psi_PTR = as.double(PT_AC_state_list$psi_PTR_AC[ , , C_pro]),
              gamma_D_PTR = as.double(PT_AC_state_list$gamma_D_PTR_AC[ , , C_pro]),
              nu_D_PTR = as.double(PT_AC_state_list$nu_D_PTR_AC[ , , C_pro]),
              gamma_R_PTR = as.double(PT_AC_state_list$gamma_R_PTR_AC[ , , C_pro]),
              nu_R_PTR = as.double(PT_AC_state_list$nu_R_PTR_AC[ , , C_pro]),
              l_D0_PTR = as.double(PT_AC_state_list$l_D0_PTR_AC[ , C_pro]),
              z_D0_PTR = as.double(PT_AC_state_list$z_D0_PTR_AC[ , C_pro]),
              l_R0_PTR = as.double(PT_AC_state_list$l_R0_PTR_AC[ , C_pro]),
              z_R0_PTR = as.double(PT_AC_state_list$z_R0_PTR_AC[ , C_pro]),
              loglik_PTR = as.double(rep(0, nThread)),
              logpost_PTR = as.double(rep(0, nThread)),
              N = as.double(N), n = as.double(n),
              M = as.double(M), m = as.double(m),
              S = as.integer(S), T = as.integer(T), 
              C = as.integer(C_pro), G = as.integer(G),
              K_min = as.integer(K_min), K_max = as.integer(K_max),
              a_w = as.double(a_w), b_w = as.double(b_w), 
              d = as.double(d), d0 = as.double(d0),
              a_lambda = as.double(a_lambda), b_lambda = as.double(b_lambda), 
              a_pai = as.double(a_pai1), b_pai = as.double(b_pai1), 
              a_zeta = as.double(a_zeta1), b_zeta = as.double(b_zeta1),
              a_phi = as.double(a_phi), b_phi = as.double(b_phi),
              a_psi = as.double(a_psi), b_psi = as.double(b_psi),
              a_gamma_D = as.double(a_gamma_D), b_gamma_D = as.double(b_gamma_D),
              a_gamma_R = as.double(a_gamma_R), b_gamma_R = as.double(b_gamma_R),
              a_nu_D = as.double(a_nu_D), b_nu_D = as.double(b_nu_D),
              a_nu_R = as.double(a_nu_R), b_nu_R = as.double(b_nu_R),
              niter = as.integer(thin), g_fun = as.integer(g_fun - 1),
              Delta = as.double(Delta_tau), nThread = as.integer(nThread))
    
    PT_AC_state_list$L_PTR_AC[[C_pro]] = array(output$L_PTR, c(S, C_pro, nThread))
    PT_AC_state_list$Z_PTR_AC[[C_pro]] = array(output$Z_PTR, c(S, C_pro, nThread))
    PT_AC_state_list$W_PTR_AC[[C_pro]] = array(output$W_PTR, c(T, C_pro + 1, nThread))
    PT_AC_state_list$Lambda_PTR_AC[[C_pro]] = array(output$Lambda_PTR, c(G, C_pro, nThread))
    
    PT_AC_state_list$pai_PTR_AC[[C_pro]] = matrix(output$pai_PTR, C_pro, nThread)
    PT_AC_state_list$zeta_PTR_AC[[C_pro]] = matrix(output$zeta_PTR, C_pro, nThread)
    
    PT_AC_state_list$phi_PTR_AC[ , , C_pro] = matrix(output$phi_PTR, T, nThread)
    PT_AC_state_list$psi_PTR_AC[ , , C_pro] = matrix(output$psi_PTR, T, nThread)
    PT_AC_state_list$gamma_D_PTR_AC[ , , C_pro] = matrix(output$gamma_D_PTR, T, nThread)
    PT_AC_state_list$nu_D_PTR_AC[ , , C_pro] = matrix(output$nu_D_PTR, T, nThread)
    PT_AC_state_list$gamma_R_PTR_AC[ , , C_pro] = matrix(output$gamma_R_PTR, T, nThread)
    PT_AC_state_list$nu_R_PTR_AC[ , , C_pro] = matrix(output$nu_R_PTR, T, nThread)
    PT_AC_state_list$l_D0_PTR_AC[ , C_pro] = output$l_D0_PTR
    PT_AC_state_list$z_D0_PTR_AC[ , C_pro] = output$z_D0_PTR
    PT_AC_state_list$l_R0_PTR_AC[ , C_pro] = output$l_R0_PTR
    PT_AC_state_list$z_R0_PTR_AC[ , C_pro] = output$z_R0_PTR
    
    PT_AC_state_list$loglik_PTR_AC[ , C_pro] = output$loglik_PTR
    PT_AC_state_list$logpost_PTR_AC[ , C_pro] = output$logpost_PTR
    
    # loglik on test set
    loglik_pro = (1 - tau) * PT_AC_state_list$loglik_PTR_AC[nThread, C_pro]
    
    logpost_pro = loglik_pro + (C_pro - 1) * log(1 - alpha)
    
    
    if(is.na(logpost_cur)  | is.na(logpost_pro)){
      # if current posterior is -INF, accept proposal. (Not likely to happen)
      if(is.na(logpost_cur)){
        C_cur = C_pro
        loglik_cur = loglik_pro
        logpost_cur = logpost_pro
      }
      # otherwise, if proposal posterior is -INF, reject
    } else{
      u = runif(1, 0, 1)
      # accept the proposal with min(1, p_acc)
      if(log(u) < logpost_pro - logpost_cur){
        C_cur = C_pro
        loglik_cur = loglik_pro
        logpost_cur = logpost_pro
      }
    } # end if(is.na(logpost_cur)  | is.na(logpost_pro))
    
    
    
    sample_list$C_spls[i] = C_cur
    
    sample_list$L_spls[[i]] = PT_AC_state_list$L_PTR_AC[[C_cur]][ , , nThread]
    sample_list$Z_spls[[i]] = PT_AC_state_list$Z_PTR_AC[[C_cur]][ , , nThread]
    sample_list$W_spls[[i]] = PT_AC_state_list$W_PTR_AC[[C_cur]][ , , nThread]
    sample_list$Lambda_spls[[i]] = PT_AC_state_list$Lambda_PTR_AC[[C_cur]][ , , nThread]
    
    sample_list$pai_spls[[i]] = PT_AC_state_list$pai_PTR_AC[[C_cur]][ , nThread]
    sample_list$zeta_spls[[i]] = PT_AC_state_list$zeta_PTR_AC[[C_cur]][ , nThread]
  
    sample_list$phi_spls[ , i] = PT_AC_state_list$phi_PTR_AC[ , nThread, C_cur]
    sample_list$psi_spls[ , i] = PT_AC_state_list$psi_PTR_AC[ , nThread, C_cur]
    sample_list$gamma_D_spls[ , i] = PT_AC_state_list$gamma_D_PTR_AC[ , nThread, C_cur]
    sample_list$nu_D_spls[ , i] = PT_AC_state_list$nu_D_PTR_AC[ , nThread, C_cur]
    sample_list$gamma_R_spls[ , i] = PT_AC_state_list$gamma_R_PTR_AC[ , nThread, C_cur]
    sample_list$nu_R_spls[ , i] = PT_AC_state_list$nu_R_PTR_AC[ , nThread, C_cur]
  
    sample_list$l_D0_spls[i] = PT_AC_state_list$l_D0_PTR_AC[nThread, C_cur]
    sample_list$z_D0_spls[i] = PT_AC_state_list$z_D0_PTR_AC[nThread, C_cur]
    sample_list$l_R0_spls[i] = PT_AC_state_list$l_R0_PTR_AC[nThread, C_cur]
    sample_list$z_R0_spls[i] = PT_AC_state_list$z_R0_PTR_AC[nThread, C_cur]
    
    sample_list$loglik_spls[i] = PT_AC_state_list$loglik_PTR_AC[nThread, C_cur]
    sample_list$logpost_spls[i] = PT_AC_state_list$logpost_PTR_AC[nThread, C_cur]
    
    
    # timer
    if((i%%200)==0){
      cat(sprintf("RNDClone RJ-MCMC %.2f %% finished. ", round(i/niter*100, 2)))
      cat(sprintf("Date: %s.\n", date()))
    }
    
  } # end for(i in 1:niter)
  
  
  result_list = list()
  result_list$sample_list = sample_list
  result_list$PT_AC_state_list = PT_AC_state_list
  
  cat("\nRNDClone RJ-MCMC (with parallel tempering, jump across different C) finished.\n")
  cat(sprintf("C_min = %d, C_max = %d. Date: %s.\n\n", C_min, C_max, date()))
  
  return(result_list)
  
}


























###################################################################################
# 3. DClone: Function for Trans-dimensional (Reversible Jump) MCMC,
#    if only DNA data are available
###################################################################################
DClone_RJMCMC = function(n, N,
  C_min = 2, C_max = 7, K_min = 1, K_max = 3,
  niter = 5000, burnin = 20000, thin = 2, Delta = 1.15^(9:0), tau = 0.99,
  alpha = 0.8, a_w = 1, b_w = 1, d = 1, d0 = 0.03,
  a_pai = NULL, b_pai = NULL, a_zeta = NULL, b_zeta = NULL,
  a_gamma_D = 1, b_gamma_D = NULL, a_nu_D = 1, b_nu_D = NULL,
  a_phi = NULL, b_phi = NULL,
  verbose = FALSE){
  
  # tau: power prior
  cat("DClone MCMC (with parallel tempering, jump across different C) started.\n")
  cat(sprintf("C_min = %d, C_max = %d. Date: %s.\n\n", C_min, C_max, date()))
  
  S = dim(n)[1]
  T = dim(n)[2]
  
  ############################################################################
  # Setting Hyperparameters
  ############################################################################
  if(is.null(b_phi)) b_phi = rep(10, T)
  if(is.null(a_phi)) a_phi = colMeans(N) * b_phi
  if(is.null(b_gamma_D)) b_gamma_D = a_gamma_D * 10 * mean(N)
  if(is.null(b_nu_D)) b_nu_D = a_nu_D * 10 * mean(N)
  
  if(verbose){
    print("a_phi = ")
    print(a_phi)
    print("b_phi = ")
    print(b_phi)
  }
  

  # Delta = c(xx, xx, ..., 1), last one is 1
  nThread = length(Delta)
  
  # thus, the last temperature corresponds to power prior with power tau
  Delta_tau = Delta / tau
  
  
  ############################################################################
  # Setting storage space for two result lists
  # PT_AC_state_list and sample_list
  ############################################################################
  # save the current state for all PT temperatures (PT) and all C (AC) in C_min -- C_max
  PT_AC_state_list = list()
  PT_AC_state_list$L_PTR_AC = list()
  PT_AC_state_list$Z_PTR_AC = list()
  PT_AC_state_list$W_PTR_AC = list()
  PT_AC_state_list$pai_PTR_AC = list()
  PT_AC_state_list$zeta_PTR_AC = list()
  PT_AC_state_list$phi_PTR_AC = array(0, c(T, nThread, C_max))
  PT_AC_state_list$gamma_D_PTR_AC = array(0, c(T, nThread, C_max))
  PT_AC_state_list$nu_D_PTR_AC = array(0, c(T, nThread, C_max))
  PT_AC_state_list$l_D0_PTR_AC = array(0, c(nThread, C_max))
  PT_AC_state_list$z_D0_PTR_AC = array(0, c(nThread, C_max))
  PT_AC_state_list$loglik_PTR_AC = array(0, c(nThread, C_max))
  PT_AC_state_list$logpost_PTR_AC = array(0, c(nThread, C_max))
  
  # save all RJ samples
  # each sample_list$x_spls is a list with niter elements
  # x_spls[[iter1]] probably has different dimension with x_spls[[iter2]]
  sample_list = list()
  sample_list$L_spls = list()
  sample_list$Z_spls = list()
  sample_list$W_spls = list()
  sample_list$pai_spls = list()
  sample_list$zeta_spls = list()
  
  sample_list$phi_spls = array(0, c(T, niter))
  sample_list$gamma_D_spls = array(0, c(T, niter))
  sample_list$nu_D_spls = array(0, c(T, niter))
  
  sample_list$l_D0_spls = rep(0, niter)
  sample_list$z_D0_spls = rep(0, niter)
  
  sample_list$loglik_spls = rep(0, niter)
  sample_list$logpost_spls = rep(0, niter)
  
  sample_list$C_spls = rep(0, niter)
  
  
  ############################################################################
  # Burnin for each C in C_min:C_max
  ############################################################################
  for(C_iter in C_min:C_max){
    
    cat("DClone RJ-MCMC burn-in started for ")
    cat(sprintf("C = %d. Date: %s.\n", C_iter, date()))
    
    if(is.null(a_pai) | is.null(b_pai)){
      a_pai1 = C_iter - 1
      b_pai1 = 1
    } else{
      a_pai1 = a_pai
      b_pai1 = b_pai
    }
    
    if(is.null(a_zeta) | is.null(b_zeta)){
      a_zeta1 = C_iter - 1
      b_zeta1 = 1
    } else{
      a_zeta1 = a_zeta
      b_zeta1 = b_zeta
    }
    
    ## initialize
    PT_AC_state_list$L_PTR_AC[[C_iter]] = array(2, c(S, C_iter, nThread))
    PT_AC_state_list$Z_PTR_AC[[C_iter]] = array(0, c(S, C_iter, nThread))
    PT_AC_state_list$W_PTR_AC[[C_iter]] = array(0, c(T, C_iter + 1, nThread))
    PT_AC_state_list$W_PTR_AC[[C_iter]][ , 1:C_iter, ] = 0.99 / C_iter
    PT_AC_state_list$W_PTR_AC[[C_iter]][ , C_iter+1, ] = 0.01

    PT_AC_state_list$pai_PTR_AC[[C_iter]] = array(a_pai1 / (a_pai1 + b_pai1), c(C_iter, nThread))
    PT_AC_state_list$zeta_PTR_AC[[C_iter]] = array(a_zeta1 / (a_zeta1 + b_zeta1), c(C_iter, nThread))
    
    PT_AC_state_list$phi_PTR_AC[ , , C_iter] = matrix(rep(a_phi / b_phi, nThread), T, nThread)
    PT_AC_state_list$gamma_D_PTR_AC[ , , C_iter] = a_gamma_D / b_gamma_D
    PT_AC_state_list$nu_D_PTR_AC[ , , C_iter] = a_nu_D / b_nu_D
    PT_AC_state_list$l_D0_PTR_AC[ , C_iter] = (K_min + K_max) / 2
    PT_AC_state_list$z_D0_PTR_AC[ , C_iter] = (K_min + K_max) / 4
    # PT_AC_state_list$loglik_PTR_AC[ , C_iter] = rep(0, nThread)
    # PT_AC_state_list$logpost_PTR_AC[ , C_iter] = rep(0, nThread)
    
    
    ## Burnin 
    output = .C("DClone_MCMC", L_PTR = as.integer(PT_AC_state_list$L_PTR_AC[[C_iter]]),
              Z_PTR = as.integer(PT_AC_state_list$Z_PTR_AC[[C_iter]]),
              W_PTR = as.double(PT_AC_state_list$W_PTR_AC[[C_iter]]),
              pai_PTR = as.double(PT_AC_state_list$pai_PTR_AC[[C_iter]]),
              zeta_PTR = as.double(PT_AC_state_list$zeta_PTR_AC[[C_iter]]),
              phi_PTR = as.double(PT_AC_state_list$phi_PTR_AC[ , , C_iter]),
              gamma_D_PTR = as.double(PT_AC_state_list$gamma_D_PTR_AC[ , , C_iter]),
              nu_D_PTR = as.double(PT_AC_state_list$nu_D_PTR_AC[ , , C_iter]),
              l_D0_PTR = as.double(PT_AC_state_list$l_D0_PTR_AC[ , C_iter]),
              z_D0_PTR = as.double(PT_AC_state_list$z_D0_PTR_AC[ , C_iter]),
              loglik_PTR = as.double(rep(0, nThread)),
              logpost_PTR = as.double(rep(0, nThread)),
              N = as.double(N), n = as.double(n),
              S = as.integer(S), T = as.integer(T), 
              C = as.integer(C_iter), K_min = as.integer(K_min), K_max = as.integer(K_max),
              a_w = as.double(a_w), b_w = as.double(b_w), 
              d = as.double(d), d0 = as.double(d0),
              a_pai = as.double(a_pai1), b_pai = as.double(b_pai1), 
              a_zeta = as.double(a_zeta1), b_zeta = as.double(b_zeta1),
              a_phi = as.double(a_phi), b_phi = as.double(b_phi),
              a_gamma_D = as.double(a_gamma_D), b_gamma_D = as.double(b_gamma_D),
              a_nu_D = as.double(a_nu_D), b_nu_D = as.double(b_nu_D),
              niter = as.integer(burnin), Delta = as.double(Delta_tau), 
              nThread = as.integer(nThread))
    
    PT_AC_state_list$L_PTR_AC[[C_iter]] = array(output$L_PTR, c(S, C_iter, nThread))
    PT_AC_state_list$Z_PTR_AC[[C_iter]] = array(output$Z_PTR, c(S, C_iter, nThread))
    PT_AC_state_list$W_PTR_AC[[C_iter]] = array(output$W_PTR, c(T, C_iter + 1, nThread))
    
    PT_AC_state_list$pai_PTR_AC[[C_iter]] = matrix(output$pai_PTR, C_iter, nThread)
    PT_AC_state_list$zeta_PTR_AC[[C_iter]] = matrix(output$zeta_PTR, C_iter, nThread)
    
    PT_AC_state_list$phi_PTR_AC[ , , C_iter] = matrix(output$phi_PTR, T, nThread)
    PT_AC_state_list$gamma_D_PTR_AC[ , , C_iter] = matrix(output$gamma_D_PTR, T, nThread)
    PT_AC_state_list$nu_D_PTR_AC[ , , C_iter] = matrix(output$nu_D_PTR, T, nThread)
    PT_AC_state_list$l_D0_PTR_AC[ , C_iter] = output$l_D0_PTR
    PT_AC_state_list$z_D0_PTR_AC[ , C_iter] = output$z_D0_PTR
    
    PT_AC_state_list$loglik_PTR_AC[ , C_iter] = output$loglik_PTR
    PT_AC_state_list$logpost_PTR_AC[ , C_iter] = output$logpost_PTR
    
    if (verbose) {
      print("Log-likelihood after burn-in:")
      print(PT_AC_state_list$loglik_PTR_AC)
      print("Log-posterior after burn-in:")
      print(PT_AC_state_list$logpost_PTR_AC)
    }
    

    cat("DClone RJ-MCMC burn-in finished for ")
    cat(sprintf("C = %d. Date: %s.\n", C_iter, date()))
    
  }
  ## End burnin
  
  
  # initial value of C, and loglik on test set (for evaluating p_acc)
  C_cur = sample(C_min:C_max, 1)
  
            
  loglik_cur = (1 - tau) * PT_AC_state_list$loglik_PTR_AC[nThread, C_cur]
  
  logpost_cur = loglik_cur + (C_cur - 1) * log(1 - alpha)
  
  for(i in 1:niter){
    
    C_pro = sample(C_min:C_max, 1)
    
    if(is.null(a_pai) | is.null(b_pai)){
      a_pai1 = C_pro - 1
      b_pai1 = 1
    } else{
      a_pai1 = a_pai
      b_pai1 = b_pai
    }
    
    if(is.null(a_zeta) | is.null(b_zeta)){
      a_zeta1 = C_pro - 1
      b_zeta1 = 1
    } else{
      a_zeta1 = a_zeta
      b_zeta1 = b_zeta
    }
    
    output = .C("DClone_MCMC", L_PTR = as.integer(PT_AC_state_list$L_PTR_AC[[C_pro]]),
              Z_PTR = as.integer(PT_AC_state_list$Z_PTR_AC[[C_pro]]),
              W_PTR = as.double(PT_AC_state_list$W_PTR_AC[[C_pro]]),
              pai_PTR = as.double(PT_AC_state_list$pai_PTR_AC[[C_pro]]),
              zeta_PTR = as.double(PT_AC_state_list$zeta_PTR_AC[[C_pro]]),
              phi_PTR = as.double(PT_AC_state_list$phi_PTR_AC[ , , C_pro]),
              gamma_D_PTR = as.double(PT_AC_state_list$gamma_D_PTR_AC[ , , C_pro]),
              nu_D_PTR = as.double(PT_AC_state_list$nu_D_PTR_AC[ , , C_pro]),
              l_D0_PTR = as.double(PT_AC_state_list$l_D0_PTR_AC[ , C_pro]),
              z_D0_PTR = as.double(PT_AC_state_list$z_D0_PTR_AC[ , C_pro]),
              loglik_PTR = as.double(rep(0, nThread)),
              logpost_PTR = as.double(rep(0, nThread)),
              N = as.double(N), n = as.double(n),
              S = as.integer(S), T = as.integer(T), 
              C = as.integer(C_pro),
              K_min = as.integer(K_min), K_max = as.integer(K_max),
              a_w = as.double(a_w), b_w = as.double(b_w), 
              d = as.double(d), d0 = as.double(d0),
              a_pai = as.double(a_pai1), b_pai = as.double(b_pai1), 
              a_zeta = as.double(a_zeta1), b_zeta = as.double(b_zeta1),
              a_phi = as.double(a_phi), b_phi = as.double(b_phi),
              a_gamma_D = as.double(a_gamma_D), b_gamma_D = as.double(b_gamma_D),
              a_nu_D = as.double(a_nu_D), b_nu_D = as.double(b_nu_D),
              niter = as.integer(thin),
              Delta = as.double(Delta_tau), nThread = as.integer(nThread))
    
    PT_AC_state_list$L_PTR_AC[[C_pro]] = array(output$L_PTR, c(S, C_pro, nThread))
    PT_AC_state_list$Z_PTR_AC[[C_pro]] = array(output$Z_PTR, c(S, C_pro, nThread))
    PT_AC_state_list$W_PTR_AC[[C_pro]] = array(output$W_PTR, c(T, C_pro + 1, nThread))

    PT_AC_state_list$pai_PTR_AC[[C_pro]] = matrix(output$pai_PTR, C_pro, nThread)
    PT_AC_state_list$zeta_PTR_AC[[C_pro]] = matrix(output$zeta_PTR, C_pro, nThread)
    
    PT_AC_state_list$phi_PTR_AC[ , , C_pro] = matrix(output$phi_PTR, T, nThread)
    PT_AC_state_list$gamma_D_PTR_AC[ , , C_pro] = matrix(output$gamma_D_PTR, T, nThread)
    PT_AC_state_list$nu_D_PTR_AC[ , , C_pro] = matrix(output$nu_D_PTR, T, nThread)
    PT_AC_state_list$l_D0_PTR_AC[ , C_pro] = output$l_D0_PTR
    PT_AC_state_list$z_D0_PTR_AC[ , C_pro] = output$z_D0_PTR

    PT_AC_state_list$loglik_PTR_AC[ , C_pro] = output$loglik_PTR
    PT_AC_state_list$logpost_PTR_AC[ , C_pro] = output$logpost_PTR
    
    # loglik on test set
    loglik_pro = (1 - tau) * PT_AC_state_list$loglik_PTR_AC[nThread, C_pro]
    
    logpost_pro = loglik_pro + (C_pro - 1) * log(1 - alpha)
    
    
    ###### Evaluate whether to accept or not
    if(is.na(logpost_cur)  | is.na(logpost_pro)){
      # if current posterior is -INF, accept proposal. (Not likely to happen)
      if(is.na(logpost_cur)){
        C_cur = C_pro
        loglik_cur = loglik_pro
        logpost_cur = logpost_pro
      }
      # otherwise, if proposal posterior is -INF, reject
    } else{
      u = runif(1, 0, 1)
      # accept the proposal with min(1, p_acc)
      if(log(u) < logpost_pro - logpost_cur){
        C_cur = C_pro
        loglik_cur = loglik_pro
        logpost_cur = logpost_pro
      }
    } # end if(is.na(logpost_cur)  | is.na(logpost_pro))
    
    
    ########## Update Sample list
    sample_list$C_spls[i] = C_cur
    
    sample_list$L_spls[[i]] = PT_AC_state_list$L_PTR_AC[[C_cur]][ , , nThread]
    sample_list$Z_spls[[i]] = PT_AC_state_list$Z_PTR_AC[[C_cur]][ , , nThread]
    sample_list$W_spls[[i]] = PT_AC_state_list$W_PTR_AC[[C_cur]][ , , nThread]
    
    sample_list$pai_spls[[i]] = PT_AC_state_list$pai_PTR_AC[[C_cur]][ , nThread]
    sample_list$zeta_spls[[i]] = PT_AC_state_list$zeta_PTR_AC[[C_cur]][ , nThread]
  
    sample_list$phi_spls[ , i] = PT_AC_state_list$phi_PTR_AC[ , nThread, C_cur]
    sample_list$gamma_D_spls[ , i] = PT_AC_state_list$gamma_D_PTR_AC[ , nThread, C_cur]
    sample_list$nu_D_spls[ , i] = PT_AC_state_list$nu_D_PTR_AC[ , nThread, C_cur]
  
    sample_list$l_D0_spls[i] = PT_AC_state_list$l_D0_PTR_AC[nThread, C_cur]
    sample_list$z_D0_spls[i] = PT_AC_state_list$z_D0_PTR_AC[nThread, C_cur]
    
    sample_list$loglik_spls[i] = PT_AC_state_list$loglik_PTR_AC[nThread, C_cur]
    sample_list$logpost_spls[i] = PT_AC_state_list$logpost_PTR_AC[nThread, C_cur]
    
    
    # timer
    if((i%%200)==0){
      cat(sprintf("DClone RJ-MCMC %.2f %% finished. ", round(i/niter*100, 2)))
      cat(sprintf("Date: %s.\n", date()))
    }
    
  } # end for(i in 1:niter)
  
  
  result_list = list()
  result_list$sample_list = sample_list
  result_list$PT_AC_state_list = PT_AC_state_list
  
  cat("\nDClone MCMC (with parallel tempering, jump across different C) finished.\n")
  cat(sprintf("C_min = %d, C_max = %d. Date: %s.\n\n", C_min, C_max, date()))
  
  return(result_list)
  
}

































###################################################################################
# 4. RClone: Function for Trans-dimensional (Reversible Jump) MCMC,
#    if only RNA data are available
###################################################################################
RClone_RJMCMC = function(M, 
  C_min = 2, C_max = 7, g_fun = NULL, K_min = 1, K_max = 3, 
  niter = 5000, burnin = 20000, thin = 2, Delta = 1.15^(9:0), tau = 0.99,
  alpha = 0.8, a_w = 1, b_w = 1, d = 1, d0 = 0.03, 
  a_lambda = 1, b_lambda = 1,
  a_gamma_R = 1, b_gamma_R = NULL, a_nu_R = 1, b_nu_R = NULL,
  a_psi = NULL, b_psi = NULL,
  verbose = FALSE){
  
  # tau: power prior
  cat("RClone MCMC (with parallel tempering, jump across different C) started.\n")
  cat(sprintf("C_min = %d, C_max = %d. Date: %s.\n\n", C_min, C_max, date()))  
  
  S = dim(M)[1]
  T = dim(M)[2]

  if(is.null(g_fun)) g_fun = 1:S

  m_dummy = matrix(0, S, T)
  
  G = max(g_fun)
  
  ############################################################################
  # Setting Hyperparameters
  ############################################################################
  if(is.null(b_psi)) b_psi = rep(10, T)
  if(is.null(a_psi)) a_psi = colMeans(M) * b_psi
  if(is.null(b_gamma_R)) b_gamma_R = a_gamma_R * 10 * mean(M)
  if(is.null(b_nu_R)) b_nu_R = a_nu_R * 10 * mean(M)
  
  if (verbose) {
    print("a_psi = ")
    print(a_psi)
    print("b_psi = ")
    print(b_psi)
  }
  
  
  # Delta = c(xx, xx, ..., 1), last one is 1
  nThread = length(Delta)
  
  # thus, the last temperature corresponds to power prior with power tau
  Delta_tau = Delta / tau
  
  
  ############################################################################
  # Setting storage space for two result lists
  # PT_AC_state_list and sample_list
  ############################################################################
  # save the current state for all PT temperatures (PT) and all C (AC) in C_min -- C_max
  PT_AC_state_list = list()
  PT_AC_state_list$W_PTR_AC = list()
  PT_AC_state_list$Lambda_PTR_AC = list()
  PT_AC_state_list$psi_PTR_AC = array(0, c(T, nThread, C_max))
  PT_AC_state_list$gamma_R_PTR_AC = array(0, c(T, nThread, C_max))
  PT_AC_state_list$nu_R_PTR_AC = array(0, c(T, nThread, C_max))
  PT_AC_state_list$l_R0_PTR_AC = array(0, c(nThread, C_max))
  PT_AC_state_list$z_R0_PTR_AC = array(0, c(nThread, C_max))
  PT_AC_state_list$loglik_PTR_AC = array(0, c(nThread, C_max))
  PT_AC_state_list$logpost_PTR_AC = array(0, c(nThread, C_max))
  
  # save all RJ samples
  # each sample_list$x_spls is a list with niter elements
  # x_spls[[iter1]] probably has different dimension with x_spls[[iter2]]
  sample_list = list()
  sample_list$W_spls = list()
  sample_list$Lambda_spls = list()
  
  sample_list$psi_spls = array(0, c(T, niter))
  sample_list$gamma_R_spls = array(0, c(T, niter))
  sample_list$nu_R_spls = array(0, c(T, niter))
  
  sample_list$l_R0_spls = rep(0, niter)
  sample_list$z_R0_spls = rep(0, niter)
  
  sample_list$loglik_spls = rep(0, niter)
  sample_list$logpost_spls = rep(0, niter)
  
  sample_list$C_spls = rep(0, niter)
  

  ############################################################################
  # Burnin for each C in C_min:C_max
  ############################################################################
  for(C_iter in C_min:C_max){
    
    cat("RClone RJ-MCMC burn-in started for ")
    cat(sprintf("C = %d. Date: %s.\n", C_iter, date())) 
    
    ## initialize
    L_dummy = matrix(2, S, C_iter)
    Z_dummy = matrix(0, S, C_iter)

    PT_AC_state_list$W_PTR_AC[[C_iter]] = array(0, c(T, C_iter + 1, nThread))
    PT_AC_state_list$W_PTR_AC[[C_iter]][ , 1:C_iter, ] = 0.99 / C_iter
    PT_AC_state_list$W_PTR_AC[[C_iter]][ , C_iter+1, ] = 0.01
    PT_AC_state_list$Lambda_PTR_AC[[C_iter]] = array(a_lambda / b_lambda, c(G, C_iter, nThread))

    PT_AC_state_list$psi_PTR_AC[ , , C_iter] = matrix(rep(a_psi / b_psi, nThread), T, nThread)
    PT_AC_state_list$gamma_R_PTR_AC[ , , C_iter] = a_gamma_R / b_gamma_R
    PT_AC_state_list$nu_R_PTR_AC[ , , C_iter] = a_nu_R / b_nu_R
    PT_AC_state_list$l_R0_PTR_AC[ , C_iter] = (K_min + K_max) / 2
    PT_AC_state_list$z_R0_PTR_AC[ , C_iter] = (K_min + K_max) / 4
    # PT_AC_state_list$loglik_PTR_AC[ , C_iter] = rep(0, nThread)
    # PT_AC_state_list$logpost_PTR_AC[ , C_iter] = rep(0, nThread)
    
    # Burnin using RNA
    output = .C("RClone_MCMC", L = as.integer(L_dummy),
              Z = as.integer(Z_dummy),
              W_PTR = as.double(PT_AC_state_list$W_PTR_AC[[C_iter]]),
              Lambda_PTR = as.double(PT_AC_state_list$Lambda_PTR_AC[[C_iter]]),
              psi_PTR = as.double(PT_AC_state_list$psi_PTR_AC[ , , C_iter]),
              gamma_R_PTR = as.double(PT_AC_state_list$gamma_R_PTR_AC[ , , C_iter]),
              nu_R_PTR = as.double(PT_AC_state_list$nu_R_PTR_AC[ , , C_iter]),
              l_R0_PTR = as.double(PT_AC_state_list$l_R0_PTR_AC[ , C_iter]),
              z_R0_PTR = as.double(PT_AC_state_list$z_R0_PTR_AC[ , C_iter]),
              loglik_PTR = as.double(rep(0, nThread)),
              logpost_PTR = as.double(rep(0, nThread)),
              M = as.double(M), m = as.double(m_dummy),
              S = as.integer(S), T = as.integer(T), 
              C = as.integer(C_iter), G = as.integer(G),
              K_min = as.integer(K_min), K_max = as.integer(K_max),
              a_w = as.double(a_w), b_w = as.double(b_w), 
              d = as.double(d), d0 = as.double(d0),
              a_lambda = as.double(a_lambda), b_lambda = as.double(b_lambda), 
              a_psi = as.double(a_psi), b_psi = as.double(b_psi),
              a_gamma_R = as.double(a_gamma_R), b_gamma_R = as.double(b_gamma_R),
              a_nu_R = as.double(a_nu_R), b_nu_R = as.double(b_nu_R),
              niter = as.integer(burnin), g_fun = as.integer(g_fun - 1),
              Delta = as.double(Delta_tau), nThread = as.integer(nThread))
  
    PT_AC_state_list$W_PTR_AC[[C_iter]] = array(output$W_PTR, c(T, C_iter + 1, nThread))
    PT_AC_state_list$Lambda_PTR_AC[[C_iter]] = array(output$Lambda_PTR, c(G, C_iter, nThread))
    PT_AC_state_list$psi_PTR_AC[ , , C_iter] = matrix(output$psi_PTR, T, nThread)
    PT_AC_state_list$gamma_R_PTR_AC[ , , C_iter] = matrix(output$gamma_R_PTR, T, nThread)
    PT_AC_state_list$nu_R_PTR_AC[ , , C_iter] = matrix(output$nu_R_PTR, T, nThread)
    PT_AC_state_list$l_R0_PTR_AC[ , C_iter] = output$l_R0_PTR
    PT_AC_state_list$z_R0_PTR_AC[ , C_iter] = output$z_R0_PTR
    
    PT_AC_state_list$loglik_PTR_AC[ , C_iter] = output$loglik_PTR
    PT_AC_state_list$logpost_PTR_AC[ , C_iter] = output$logpost_PTR
    
    if (verbose) {
      print("Log-likelihood after burn-in:")
      print(PT_AC_state_list$loglik_PTR_AC)
      print("Log-posterior after burn-in:")
      print(PT_AC_state_list$logpost_PTR_AC)
    }
    
    cat("RClone RJ-MCMC burn-in finished for ")
    cat(sprintf("C = %d. Date: %s.\n", C_iter, date()))

  }
  ## End burnin
  
  
  # initial value of C, and loglik on test set (for evaluating p_acc)
  C_cur = sample(C_min:C_max, 1)
  
            
  loglik_cur = (1 - tau) * PT_AC_state_list$loglik_PTR_AC[nThread, C_cur]
  
  logpost_cur = loglik_cur + (C_cur - 1) * log(1 - alpha)
  
  for(i in 1:niter){
    
    C_pro = sample(C_min:C_max, 1)

    L_dummy = matrix(2, S, C_pro)
    Z_dummy = matrix(0, S, C_pro)

    output = .C("RClone_MCMC", L = as.integer(L_dummy),
              Z = as.integer(Z_dummy),
              W_PTR = as.double(PT_AC_state_list$W_PTR_AC[[C_pro]]),
              Lambda_PTR = as.double(PT_AC_state_list$Lambda_PTR_AC[[C_pro]]),
              psi_PTR = as.double(PT_AC_state_list$psi_PTR_AC[ , , C_pro]),
              gamma_R_PTR = as.double(PT_AC_state_list$gamma_R_PTR_AC[ , , C_pro]),
              nu_R_PTR = as.double(PT_AC_state_list$nu_R_PTR_AC[ , , C_pro]),
              l_R0_PTR = as.double(PT_AC_state_list$l_R0_PTR_AC[ , C_pro]),
              z_R0_PTR = as.double(PT_AC_state_list$z_R0_PTR_AC[ , C_pro]),
              loglik_PTR = as.double(rep(0, nThread)),
              logpost_PTR = as.double(rep(0, nThread)),
              M = as.double(M), m = as.double(m_dummy),
              S = as.integer(S), T = as.integer(T), 
              C = as.integer(C_pro), G = as.integer(G),
              K_min = as.integer(K_min), K_max = as.integer(K_max),
              a_w = as.double(a_w), b_w = as.double(b_w), 
              d = as.double(d), d0 = as.double(d0),
              a_lambda = as.double(a_lambda), b_lambda = as.double(b_lambda), 
              a_psi = as.double(a_psi), b_psi = as.double(b_psi),
              a_gamma_R = as.double(a_gamma_R), b_gamma_R = as.double(b_gamma_R),
              a_nu_R = as.double(a_nu_R), b_nu_R = as.double(b_nu_R),
              niter = as.integer(thin), g_fun = as.integer(g_fun - 1),
              Delta = as.double(Delta_tau), nThread = as.integer(nThread))
    
    PT_AC_state_list$W_PTR_AC[[C_pro]] = array(output$W_PTR, c(T, C_pro + 1, nThread))
    PT_AC_state_list$Lambda_PTR_AC[[C_pro]] = array(output$Lambda_PTR, c(G, C_pro, nThread))
    
    PT_AC_state_list$psi_PTR_AC[ , , C_pro] = matrix(output$psi_PTR, T, nThread)
    PT_AC_state_list$gamma_R_PTR_AC[ , , C_pro] = matrix(output$gamma_R_PTR, T, nThread)
    PT_AC_state_list$nu_R_PTR_AC[ , , C_pro] = matrix(output$nu_R_PTR, T, nThread)
    PT_AC_state_list$l_R0_PTR_AC[ , C_pro] = output$l_R0_PTR
    PT_AC_state_list$z_R0_PTR_AC[ , C_pro] = output$z_R0_PTR
    
    PT_AC_state_list$loglik_PTR_AC[ , C_pro] = output$loglik_PTR
    PT_AC_state_list$logpost_PTR_AC[ , C_pro] = output$logpost_PTR
    
    # loglik on test set
    loglik_pro = (1 - tau) * PT_AC_state_list$loglik_PTR_AC[nThread, C_pro]
    
    logpost_pro = loglik_pro + (C_pro - 1) * log(1 - alpha)
    
    
    if(is.na(logpost_cur)  | is.na(logpost_pro)){
      # if current posterior is -INF, accept proposal. (Not likely to happen)
      if(is.na(logpost_cur)){
        C_cur = C_pro
        loglik_cur = loglik_pro
        logpost_cur = logpost_pro
      }
      # otherwise, if proposal posterior is -INF, reject
    } else{
      u = runif(1, 0, 1)
      # accept the proposal with min(1, p_acc)
      if(log(u) < logpost_pro - logpost_cur){
        C_cur = C_pro
        loglik_cur = loglik_pro
        logpost_cur = logpost_pro
      }
    } # end if(is.na(logpost_cur)  | is.na(logpost_pro))
    
    
    sample_list$C_spls[i] = C_cur
    
    sample_list$W_spls[[i]] = PT_AC_state_list$W_PTR_AC[[C_cur]][ , , nThread]
    sample_list$Lambda_spls[[i]] = PT_AC_state_list$Lambda_PTR_AC[[C_cur]][ , , nThread]
  
    sample_list$psi_spls[ , i] = PT_AC_state_list$psi_PTR_AC[ , nThread, C_cur]
    sample_list$gamma_R_spls[ , i] = PT_AC_state_list$gamma_R_PTR_AC[ , nThread, C_cur]
    sample_list$nu_R_spls[ , i] = PT_AC_state_list$nu_R_PTR_AC[ , nThread, C_cur]
  
    sample_list$l_R0_spls[i] = PT_AC_state_list$l_R0_PTR_AC[nThread, C_cur]
    sample_list$z_R0_spls[i] = PT_AC_state_list$z_R0_PTR_AC[nThread, C_cur]
    
    sample_list$loglik_spls[i] = PT_AC_state_list$loglik_PTR_AC[nThread, C_cur]
    sample_list$logpost_spls[i] = PT_AC_state_list$logpost_PTR_AC[nThread, C_cur]
    
    
    # timer
    if((i%%200)==0){
      cat(sprintf("RClone RJ-MCMC %.2f %% finished. ", round(i/niter*100, 2)))
      cat(sprintf("Date: %s.\n", date()))
    }
    
  } # end for(i in 1:niter)
  
  
  result_list = list()
  result_list$sample_list = sample_list
  result_list$PT_AC_state_list = PT_AC_state_list
  
  cat("\nRClone MCMC (with parallel tempering, across different C) finished.\n")
  cat(sprintf("C_min = %d, C_max = %d. Date: %s.\n\n", C_min, C_max, date()))  

  return(result_list)
  
}

