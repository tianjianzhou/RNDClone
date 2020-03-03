/*
Used for R
*/

#include "stdio.h"
#include "R.h"
#include "Rmath.h"
#include "assert.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"

#define MIN_COUNT 1.0E-10
#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))

extern "C"{
/**********************************************************************
 MCMC functions
**********************************************************************/

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// DNA + RNA integrated analysis
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////
// Update L
/////////////////////////////////////////
void update_L(int *L, int *Z, double *W, double *Lambda, double *pai, double *zeta, double *phi, double *psi,
    double *gamma_D, double *nu_D, double *gamma_R, double *nu_R,
    double *l_D0, double *z_D0, double *l_R0, double *z_R0,
    double *N, double *n, double *M, double *m, 
    int S, int T, int C, int K_min, int K_max, int G,
    int *g_fun, double Tmp){
    
    // element-wise update L, sequentially update l_sc
    
    int s, t, c, k, g, c_iter, z_sc, k1, max_z_sc_K_min;
    // int K = K_max - K_min + 1;
    
    
    double l_D0_num = *l_D0;
    double z_D0_num = *z_D0;
    double l_R0_num = *l_R0;
    double z_R0_num = *z_R0;
    
    double temp1;
    double A_st, B_st, A_tilde_st, B_tilde_st, p_st, q_st, lambda_gt0;
    double loglik, logpost;
    double max_log_prob_l_sc = 0.0, sum_prob_l_sc;
    double randomUnif;
    
    
    double *prob_l_sc; 
    prob_l_sc = (double *)calloc(K_max + 1, sizeof(double));
    
    for(s = 0; s < S; s++){
        g = g_fun[s];
        for(c = 1; c < C; c++){
            // update l_sc (c >= 1, c = 0 corresponds to the normal subclone)
            z_sc = Z[c * S + s];
            max_z_sc_K_min = max(z_sc, K_min);
            
            if(max_z_sc_K_min == K_max){
                L[c * S + s] = K_max;
            }
            else{
                for(k = max_z_sc_K_min; k <= K_max; k++){
                    // calculate p(l_sc = k | ...)
                    loglik = 0;
                    for(t = 0; t < T; t++){
                        A_st = 0;
                        B_st = 0;
                        A_tilde_st = 0;
                        B_tilde_st = 0;
                        lambda_gt0 = 0;
                        
                        for(c_iter = 0; c_iter < C; c_iter++){
                            if(c_iter != c){
                                A_st = A_st + W[c_iter * T + t] * L[c_iter * S + s];
                                B_st = B_st + W[c_iter * T + t] * L[c_iter * S + s] * Lambda[c_iter * G + g];
                            }
                            else{
                                A_st = A_st + W[c_iter * T + t] * k;
                                B_st = B_st + W[c_iter * T + t] * k * Lambda[c_iter * G + g];
                            }
                            
                            A_tilde_st = A_tilde_st + W[c_iter * T + t] * Z[c_iter * S + s];
                            B_tilde_st = B_tilde_st + W[c_iter * T + t] * Z[c_iter * S + s]* Lambda[c_iter * G + g];
                            lambda_gt0 = lambda_gt0 + W[c_iter * T + t] * Lambda[c_iter * G + g];
                        }
                    
                        A_st = A_st + W[C * T + t] * l_D0_num;
                        A_tilde_st = A_tilde_st + W[C * T + t] * z_D0_num;
                        p_st = A_tilde_st / A_st;
                        B_st = B_st + W[C * T + t] * l_R0_num * lambda_gt0;
                        B_tilde_st = B_tilde_st + W[C * T + t] * z_R0_num * lambda_gt0;
                        q_st = B_tilde_st / B_st;
                        
                        // if N_st = 0 and M_st = 0, done
                        loglik = loglik - (1 / gamma_D[t] + N[t * S + s]) * log(1 + gamma_D[t] * phi[t] * A_st / 2) - (1 / gamma_R[t] + M[t * S + s]) * log(1 + gamma_R[t] * psi[t] * B_st / 2);
                        
                        // if N_st > 0
                        if(N[t * S + s] > MIN_COUNT){
                            loglik = loglik + N[t * S + s] * log(gamma_D[t] * phi[t] * A_st / 2);
                            // if n_st = 0, done
                            
                            if((N[t * S + s] - n[t * S + s]) > MIN_COUNT){
                                loglik = loglik + lgamma((1 - p_st) / nu_D[t] + N[t * S + s] - n[t * S + s]) - lgamma((1 - p_st) / nu_D[t]);
                            }
                            
                            if(n[t * S + s] > MIN_COUNT){
                                loglik = loglik + lgamma(p_st / nu_D[t] + n[t * S + s]) - lgamma(p_st / nu_D[t]);
                            }
                        }
                        
                        // if M_st > 0
                        if(M[t * S + s] > MIN_COUNT){
                            loglik = loglik + M[t * S + s] * log(gamma_R[t] * psi[t] * B_st / 2);
                            // if m_st = 0, done
                            
                            if((M[t * S + s] - m[t * S + s]) > MIN_COUNT){
                                loglik = loglik + lgamma((1 - q_st) / nu_R[t] + M[t * S + s] - m[t * S + s]) - lgamma((1 - q_st) / nu_R[t]);
                            }
                            
                            if(m[t * S + s] > MIN_COUNT){
                                loglik = loglik + lgamma(q_st / nu_R[t] + m[t * S + s]) - lgamma(q_st / nu_R[t]);
                            }
                        }
                        
                    }
                    
                    temp1 = 0;
                    for(k1 = 0; k1 <= k; k1++){
                      temp1 = temp1 + pow(1 - zeta[c], k1);
                    }
                    
                    // setting pai[c, 0, .. K_min - 1] = 0
                    logpost = loglik / Tmp + (double) abs(k - 2) * log(1 - pai[c]) - log(temp1);
                    
                    // logpost = logpost / Tmp;
                
                    prob_l_sc[k] = logpost;
                    if(k == max_z_sc_K_min){
                        max_log_prob_l_sc = logpost;
                    }
                    else{
                        if(max_log_prob_l_sc < logpost) max_log_prob_l_sc = logpost;
                    }
                } // end for(k). Now we have p(l_sc = k | ...) for all k, stored in prob_l_sc
                
                sum_prob_l_sc = 0;
                
                for(k = max_z_sc_K_min; k <= K_max; k++){
                    prob_l_sc[k] = prob_l_sc[k] - max_log_prob_l_sc;
                    prob_l_sc[k] = exp(prob_l_sc[k]);
                    sum_prob_l_sc = sum_prob_l_sc + prob_l_sc[k];
                }
                
                for(k = 0; k < max_z_sc_K_min; k++){
                    prob_l_sc[k] = 0;
                }
                
                GetRNGstate();
                randomUnif = runif(0.0, 1.0);
                PutRNGstate();
  
                // l_sc ~ Discrete(0:K, prob_l_sc). Sample l_sc.
                k = max_z_sc_K_min;
                prob_l_sc[k] = prob_l_sc[k] / sum_prob_l_sc;
                while(randomUnif > prob_l_sc[k]){
  	                k++;
  	                prob_l_sc[k] = prob_l_sc[k] / sum_prob_l_sc + prob_l_sc[k - 1];
                }
  
                L[c * S + s] = k;
            
            }// end if(max_z_sc_K_min == K_max)
        }// end for(c)
    }// end for(s)
    
    free(prob_l_sc);
    return;
}



/////////////////////////////////////////
// Update Z
/////////////////////////////////////////
void update_Z(int *L, int *Z, double *W, double *Lambda, double *zeta,
    double *nu_D, double *nu_R,
    double *l_D0, double *z_D0, double *l_R0, double *z_R0,
    double *N, double *n, double *M, double *m, 
    int S, int T, int C, int K_min, int K_max, int G,
    int *g_fun, double Tmp){
    
    int s, t, c, k, g, c_iter, l_sc;
    
    double l_D0_num = *l_D0;
    double z_D0_num = *z_D0;
    double l_R0_num = *l_R0;
    double z_R0_num = *z_R0;
    
    double A_st, B_st, A_tilde_st, B_tilde_st, p_st, q_st, lambda_gt0;
    double loglik, logpost;
    double max_log_prob_z_sc, sum_prob_z_sc;
    double randomUnif;
    
    double *prob_z_sc; 
    prob_z_sc = (double *)calloc((K_max + 1), sizeof(double));
    
    for(s = 0; s < S; s++){
        g = g_fun[s];
        for(c = 1; c < C; c++){
            // update z_sc (c >= 1, c = 0 corresponds to the normal subclone)
            
            l_sc = L[c * S + s];
            // min(l_sc, K_max) is always l_sc so no need
            
            if(l_sc == 0){
                Z[c * S + s] = 0;
            }
            else{
                for(k = 0; k <= l_sc; k++){
                    // calculate p(z_sc = k | ...)
                    loglik = 0;
                    for(t = 0; t < T; t++){
                        A_st = 0;
                        B_st = 0;
                        A_tilde_st = 0;
                        B_tilde_st = 0;
                        lambda_gt0 = 0;
                        
                        for(c_iter = 0; c_iter < C; c_iter++){
                            if(c_iter != c){
                                A_tilde_st = A_tilde_st + W[c_iter * T + t] * Z[c_iter * S + s];
                                B_tilde_st = B_tilde_st + W[c_iter * T + t] * Z[c_iter * S + s]* Lambda[c_iter * G + g];
                            }
                            else{
                                A_tilde_st = A_tilde_st + W[c_iter * T + t] * k;
                                B_tilde_st = B_tilde_st + W[c_iter * T + t] * k * Lambda[c_iter * G + g];
                            }
                            A_st = A_st + W[c_iter * T + t] * L[c_iter * S + s];
                            B_st = B_st + W[c_iter * T + t] * L[c_iter * S + s] * Lambda[c_iter * G + g];
                            lambda_gt0 = lambda_gt0 + W[c_iter * T + t] * Lambda[c_iter * G + g];
                        }
                    
                        A_st = A_st + W[C * T + t] * l_D0_num;
                        A_tilde_st = A_tilde_st + W[C * T + t] * z_D0_num;
                        p_st = A_tilde_st / A_st;
                        B_st = B_st + W[C * T + t] * l_R0_num * lambda_gt0;
                        B_tilde_st = B_tilde_st + W[C * T + t] * z_R0_num * lambda_gt0;
                        q_st = B_tilde_st / B_st;
                        
                        // if N_st > 0
                        if(N[t * S + s] > MIN_COUNT){
                            if((N[t * S + s] - n[t * S + s]) > MIN_COUNT){
                                loglik = loglik + lgamma((1 - p_st) / nu_D[t] + N[t * S + s] - n[t * S + s]) - lgamma((1 - p_st) / nu_D[t]);
                            }
                            
                            if(n[t * S + s] > MIN_COUNT){
                                loglik = loglik + lgamma(p_st / nu_D[t] + n[t * S + s]) - lgamma(p_st / nu_D[t]);
                            }
                        }
                        
                        // if M_st > 0
                        if(M[t * S + s] > MIN_COUNT){
                            if((M[t * S + s] - m[t * S + s]) > MIN_COUNT){
                                loglik = loglik + lgamma((1 - q_st) / nu_R[t] + M[t * S + s] - m[t * S + s]) - lgamma((1 - q_st) / nu_R[t]);
                            }
                            
                            if(m[t * S + s] > MIN_COUNT){
                                loglik = loglik + lgamma(q_st / nu_R[t] + m[t * S + s]) - lgamma(q_st / nu_R[t]);
                            }
                        }
                        
                    }
                
                    logpost = loglik / Tmp + (double) k * log(1 - zeta[c]);
                    // logpost = logpost / Tmp;
                
                    prob_z_sc[k] = logpost;
                    if(k == 0){
                        max_log_prob_z_sc = logpost;
                    }
                    else{
                        if(max_log_prob_z_sc < logpost) max_log_prob_z_sc = logpost;
                    }
                } // end for(k). Now we have p(z_sc = k | ...) for all k, stored in prob_z_sc
                
                sum_prob_z_sc = 0;
                
                for(k = 0; k <= l_sc; k++){
                    prob_z_sc[k] = prob_z_sc[k] - max_log_prob_z_sc;
                    prob_z_sc[k] = exp(prob_z_sc[k]);
                    sum_prob_z_sc = sum_prob_z_sc + prob_z_sc[k];
                }
                
                for(k = l_sc + 1; k <= K_max; k++){
                    prob_z_sc[k] = 0;
                }
                
                /* for debug
                if(Tmp == 1.0){
                    fprintf(stdout, "Tmp = %.3lf; s = %d, c = %d, ", Tmp, s, c);
                    for(k = 0; k <= K_max; k++){
                        fprintf(stdout, "prob_z_sc[%d] = %.8lf. ", k, prob_z_sc[k]);
                    }
                }
                */
                
                GetRNGstate();
                randomUnif = runif(0.0, 1.0);
                PutRNGstate(); 
  
                // z_sc ~ Discrete(0:K_max, prob_z_sc). Sample z_sc.
                k = 0;
                prob_z_sc[k] = prob_z_sc[k] / sum_prob_z_sc;
                while(randomUnif > prob_z_sc[k]){
  	                k++;
  	                prob_z_sc[k] = prob_z_sc[k] / sum_prob_z_sc + prob_z_sc[k - 1];
                }
                
                /* for debug
                if(Tmp == 1.0){
                    fprintf(stdout, "sampled k = %d", k);
                    fprintf(stdout, "\n");
                }
                */
                
                Z[c * S + s] = k;
            }// end if(l_sc == 0)
        }// end for(c)
    }// end for(s)
    
    free(prob_z_sc);
    return;
}



/////////////////////////////////////////
// Update pai
/////////////////////////////////////////
void update_pai(int *L, double *pai, 
    int S, int C, int K_min, int K_max,
    double a_pai, double b_pai){
    
    // Update pai
    // pai: C dimensional vector
    
    int c, s, k;
    
    double pai_c_cur, pai_c_pro;  
    
    double norm_prop;
    double sd_prop = 0.1;
    double loglik_cur, loglik_pro, logpost_cur, logpost_pro, u;
    
    double temp1_cur, temp1_pro;
    
    for(c = 1; c < C; c++){
        pai_c_cur = pai[c];
        
        GetRNGstate();
        norm_prop = rnorm(0.0, sd_prop);
        PutRNGstate();
        
        pai_c_pro = pai_c_cur + norm_prop;
        
        if((pai_c_pro > 0) && (pai_c_pro < 1)){
            loglik_cur = 0.0;
            loglik_pro = 0.0;
            for(s = 0; s < S; s++){
                loglik_cur = loglik_cur + abs(L[c * S + s] - 2) * log(1 - pai_c_cur);
                loglik_pro = loglik_pro + abs(L[c * S + s] - 2) * log(1 - pai_c_pro);
            }
            
            temp1_cur = 0;
            temp1_pro = 0;
            
            for(k = K_min; k <= K_max; k++){
                temp1_cur = temp1_cur + pai_c_cur * pow(1 - pai_c_cur, abs(k - 2));
                temp1_pro = temp1_pro + pai_c_pro * pow(1 - pai_c_pro, abs(k - 2));
            }
            
            loglik_cur = loglik_cur + (double) S * log(pai_c_cur) - (double) S * log(temp1_cur);
            loglik_pro = loglik_pro + (double) S * log(pai_c_pro) - (double) S * log(temp1_pro);
            
            
            logpost_cur = loglik_cur + (a_pai - 1) * log(pai_c_cur) + (b_pai - 1) * log(1 - pai_c_cur);
            logpost_pro = loglik_pro + (a_pai - 1) * log(pai_c_pro) + (b_pai - 1) * log(1 - pai_c_pro);
            
            // logpost_cur = logpost_cur / Tmp;
            // logpost_pro = logpost_pro / Tmp;
            
            GetRNGstate();
            u = runif(0.0, 1.0);   
            PutRNGstate();
            
            if(log(u) < logpost_pro - logpost_cur){
                pai[c] = pai_c_pro;
            }
        }
    
    } //end for(c)
    
    return;
}




/////////////////////////////////////////
// Update zeta
/////////////////////////////////////////
void update_zeta(int *L, int *Z, double *zeta, 
    int S, int C,
    double a_zeta, double b_zeta){
    
    // Update zeta
    // zeta: C dimensional vector
    
    int c, s, k;
    
    double zeta_c_cur, zeta_c_pro;  
    
    double norm_prop;
    double sd_prop = 0.1;
    double loglik_cur, loglik_pro, logpost_cur, logpost_pro, u;
    
    double temp1_cur, temp1_pro;
    
    for(c = 1; c < C; c++){
        zeta_c_cur = zeta[c];
        
        GetRNGstate();
        norm_prop = rnorm(0.0, sd_prop);   
        PutRNGstate();
        
        zeta_c_pro = zeta_c_cur + norm_prop;
        
        if((zeta_c_pro > 0) && (zeta_c_pro < 1)){
            loglik_cur = 0.0;
            loglik_pro = 0.0;
            for(s = 0; s < S; s++){
                temp1_cur = 0;
                temp1_pro = 0;
                
                for(k = 0; k <= L[c * S + s]; k++){
                    temp1_cur = temp1_cur + zeta_c_cur * pow(1 - zeta_c_cur, k);
                    temp1_pro = temp1_pro + zeta_c_pro * pow(1 - zeta_c_pro, k);
                }
                
                loglik_cur = loglik_cur + Z[c * S + s] * log(1 - zeta_c_cur) - log(temp1_cur);
                loglik_pro = loglik_pro + Z[c * S + s] * log(1 - zeta_c_pro) - log(temp1_pro);
            }
            
            loglik_cur = loglik_cur + (double) S * log(zeta_c_cur);
            loglik_pro = loglik_pro + (double) S * log(zeta_c_pro);
            
            
            logpost_cur = loglik_cur + (a_zeta - 1) * log(zeta_c_cur) + (b_zeta - 1) * log(1 - zeta_c_cur);
            logpost_pro = loglik_pro + (a_zeta - 1) * log(zeta_c_pro) + (b_zeta - 1) * log(1 - zeta_c_pro);
            
            // logpost_cur = logpost_cur / Tmp;
            // logpost_pro = logpost_pro / Tmp;
            
            GetRNGstate();
            u = runif(0.0, 1.0);   
            PutRNGstate();
  
            if(log(u) < logpost_pro - logpost_cur){
                zeta[c] = zeta_c_pro;
            }
        }
    
    } //end for(c)
    
    return;
}







void update_W(int *L, int *Z, double *W, double *Lambda, double *phi, double *psi,
    double *gamma_D, double *nu_D, double *gamma_R, double *nu_R,
    double *l_D0, double *z_D0, double *l_R0, double *z_R0,
    double *N, double *n, double *M, double *m, 
    int S, int T, int C, int G,
    double a_w, double b_w, double d, double d0,
    int *g_fun, double Tmp){
  
    // Update W
    // W: T * (C + 1) matrix. 
    // W_tC means w_t0 in the paper, proportion of the background subclone in sample t.
    // W_t0 means w_t1 in the paper, proportion of normal subclone
    
    
    // only calculate proposed w for the t-th row w_t = (w_t1, ..., w_tC, w_t0)
    double *w_t_pro;
    w_t_pro = (double *)calloc((C + 1), sizeof(double));
    
    double w_tc_cur, w_tc_pro;
    
    int s, t, c, g, c_iter;
    
    double l_D0_num = *l_D0;
    double z_D0_num = *z_D0;
    double l_R0_num = *l_R0;
    double z_R0_num = *z_R0;
    
    double loglik_cur, loglik_pro, logpost_cur, logpost_pro;
    double A_st_cur, A_st_pro, A_tilde_st_cur, A_tilde_st_pro;
    double B_st_cur, B_st_pro, B_tilde_st_cur, B_tilde_st_pro;
    double p_st_cur, p_st_pro, q_st_cur, q_st_pro;
    double lambda_gt0_cur, lambda_gt0_pro;
    double norm_prop;
    double sd_prop;
    double u;
    
    double w_t_sum;
    
    for(t = 0; t < T; t++){
        for(c = 0; c <= C; c++){
            // update w_tc
            
            w_tc_cur = W[c * T + t];
            
            if(c < C){
                sd_prop = 0.08;
            }
            else{
                sd_prop = 0.01;
            }
            
            
            GetRNGstate();
            norm_prop = rnorm(0.0, sd_prop);
            PutRNGstate();
            
            w_tc_pro = w_tc_cur + norm_prop;
            
            w_t_sum = 0;
            // need to change w_{t, -c} as well
            for(c_iter = 0; c_iter <= C; c_iter++){
                if(c_iter != c){
                    w_t_pro[c_iter] = W[c_iter * T + t] * (1 - w_tc_pro) / (1 - w_tc_cur);
                }
                else{
                    w_t_pro[c_iter] = w_tc_pro;
                }
                w_t_sum = w_t_sum + w_t_pro[c_iter];
            }
            
            for(c_iter = 0; c_iter <= C; c_iter++){
                w_t_pro[c_iter] = w_t_pro[c_iter] / w_t_sum;
            }
            
            w_tc_pro = w_t_pro[c];
            
            // otherwise reject
            if((w_tc_pro < 1) && (w_tc_pro > 0) && (w_t_pro[C] < 0.03)){
                loglik_cur = 0.0;
                loglik_pro = 0.0;
            
                for(s = 0; s < S; s++){
                    g = g_fun[s];
                
                    A_st_cur = 0;
                    A_st_pro = 0;
                    A_tilde_st_cur = 0;
                    A_tilde_st_pro = 0;
                
                    B_st_cur = 0;
                    B_st_pro = 0;
                    B_tilde_st_cur = 0;
                    B_tilde_st_pro = 0;
                
                    lambda_gt0_cur = 0;
                    lambda_gt0_pro = 0;
                
                    for(c_iter = 0; c_iter < C; c_iter++){
                        A_st_cur = A_st_cur + W[c_iter * T + t] * L[c_iter * S + s];
                        A_st_pro = A_st_pro + w_t_pro[c_iter] * L[c_iter * S + s];
                        A_tilde_st_cur = A_tilde_st_cur + W[c_iter * T + t] * Z[c_iter * S + s];
                        A_tilde_st_pro = A_tilde_st_pro + w_t_pro[c_iter] * Z[c_iter * S + s];
                    
                        B_st_cur = B_st_cur + W[c_iter * T + t] * L[c_iter * S + s] * Lambda[c_iter * G + g];
                        B_st_pro = B_st_pro + w_t_pro[c_iter] * L[c_iter * S + s] * Lambda[c_iter * G + g];
                        B_tilde_st_cur = B_tilde_st_cur + W[c_iter * T + t] * Z[c_iter * S + s]* Lambda[c_iter * G + g];
                        B_tilde_st_pro = B_tilde_st_pro + w_t_pro[c_iter] * Z[c_iter * S + s]* Lambda[c_iter * G + g];
                    
                        lambda_gt0_cur = lambda_gt0_cur + W[c_iter * T + t] * Lambda[c_iter * G + g];
                        lambda_gt0_pro = lambda_gt0_pro + w_t_pro[c_iter] * Lambda[c_iter * G + g];
                    }
                
                    A_st_cur = A_st_cur + W[C * T + t] * l_D0_num;
                    A_st_pro = A_st_pro + w_t_pro[C] * l_D0_num;
                    A_tilde_st_cur = A_tilde_st_cur + W[C * T + t] * z_D0_num;
                    A_tilde_st_pro = A_tilde_st_pro + w_t_pro[C] * z_D0_num;
                
                    B_st_cur = B_st_cur + W[C * T + t] * l_R0_num * lambda_gt0_cur;
                    B_st_pro = B_st_pro + w_t_pro[C] * l_R0_num * lambda_gt0_pro;
                    B_tilde_st_cur = B_tilde_st_cur + W[C * T + t] * z_R0_num * lambda_gt0_cur;
                    B_tilde_st_pro = B_tilde_st_pro + w_t_pro[C] * z_R0_num * lambda_gt0_pro;
                
                    p_st_cur = A_tilde_st_cur / A_st_cur;
                    p_st_pro = A_tilde_st_pro / A_st_pro;
                    q_st_cur = B_tilde_st_cur / B_st_cur;
                    q_st_pro = B_tilde_st_pro / B_st_pro;
                    
                    // if N_st = 0 and M_st = 0, done
                    loglik_cur = loglik_cur - (1 / gamma_D[t] + N[t * S + s]) * log(1 + gamma_D[t] * phi[t] * A_st_cur / 2) - (1 / gamma_R[t] + M[t * S + s]) * log(1 + gamma_R[t] * psi[t] * B_st_cur / 2);
                    
                    loglik_pro = loglik_pro - (1 / gamma_D[t] + N[t * S + s]) * log(1 + gamma_D[t] * phi[t] * A_st_pro / 2) - (1 / gamma_R[t] + M[t * S + s]) * log(1 + gamma_R[t] * psi[t] * B_st_pro / 2);
                    
                    // if N_st > 0
                    if(N[t * S + s] > MIN_COUNT){
                        loglik_cur = loglik_cur + N[t * S + s] * log(gamma_D[t] * phi[t] * A_st_cur / 2);
                        
                        loglik_pro = loglik_pro + N[t * S + s] * log(gamma_D[t] * phi[t] * A_st_pro / 2);
                        // if n_st = 0, done
                            
                        if((N[t * S + s] - n[t * S + s]) > MIN_COUNT){
                            loglik_cur = loglik_cur + lgamma((1 - p_st_cur) / nu_D[t] + N[t * S + s] - n[t * S + s]) - lgamma((1 - p_st_cur) / nu_D[t]);
                            
                            loglik_pro = loglik_pro + lgamma((1 - p_st_pro) / nu_D[t] + N[t * S + s] - n[t * S + s]) - lgamma((1 - p_st_pro) / nu_D[t]);
                        }
                            
                        if(n[t * S + s] > MIN_COUNT){
                            loglik_cur = loglik_cur + lgamma(p_st_cur / nu_D[t] + n[t * S + s]) - lgamma(p_st_cur / nu_D[t]);
                            
                            loglik_pro = loglik_pro + lgamma(p_st_pro / nu_D[t] + n[t * S + s]) - lgamma(p_st_pro / nu_D[t]);
                        }
                    }
                        
                    // if M_st > 0
                    if(M[t * S + s] > MIN_COUNT){
                        loglik_cur = loglik_cur + M[t * S + s] * log(gamma_R[t] * psi[t] * B_st_cur / 2);
                        
                        loglik_pro = loglik_pro + M[t * S + s] * log(gamma_R[t] * psi[t] * B_st_pro / 2);
                        // if m_st = 0, done
                        
                        if((M[t * S + s] - m[t * S + s]) > MIN_COUNT){
                            loglik_cur = loglik_cur + lgamma((1 - q_st_cur) / nu_R[t] + M[t * S + s] - m[t * S + s]) - lgamma((1 - q_st_cur) / nu_R[t]);
                            
                            loglik_pro = loglik_pro + lgamma((1 - q_st_pro) / nu_R[t] + M[t * S + s] - m[t * S + s]) - lgamma((1 - q_st_pro) / nu_R[t]);
                        }
                            
                        if(m[t * S + s] > MIN_COUNT){
                            loglik_cur = loglik_cur + lgamma(q_st_cur / nu_R[t] + m[t * S + s]) - lgamma(q_st_cur / nu_R[t]);
                            
                            loglik_pro = loglik_pro + lgamma(q_st_pro / nu_R[t] + m[t * S + s]) - lgamma(q_st_pro / nu_R[t]);
                        }
                    }
                    // end calculating loglik
                    
                }// end for(s)
            
                logpost_cur = loglik_cur / Tmp;
                logpost_pro = loglik_pro / Tmp;
            
                logpost_cur = logpost_cur + (a_w - 1) * log(W[0 * T + t]) + (b_w - 1) * log(1 - W[0 * T + t]);
                logpost_pro = logpost_pro + (a_w - 1) * log(w_t_pro[0]) + (b_w - 1) * log(1 - w_t_pro[0]);
            
                for(c_iter = 1; c_iter < C; c_iter++){
                    logpost_cur = logpost_cur + (d - 1) * log(W[c_iter * T + t] / (1 - W[0 * T + t]));
                    logpost_pro = logpost_pro + (d - 1) * log(w_t_pro[c_iter] / (1 - w_t_pro[0]));
                }
            
                logpost_cur = logpost_cur + (d0 - 1) * log(W[C * T + t] / (1 - W[0 * T + t]));
                logpost_pro = logpost_pro + (d0 - 1) * log(w_t_pro[C] / (1 - w_t_pro[0]));
            
                // logpost_cur = logpost_cur / Tmp;
                // logpost_pro = logpost_pro / Tmp;
                
                GetRNGstate();
                u = runif(0.0, 1.0);
                PutRNGstate();
  
                if(log(u) < logpost_pro - logpost_cur){
                    for(c_iter = 0; c_iter <= C; c_iter++){
                        W[c_iter * T + t] = w_t_pro[c_iter];
                    }
                }
            }// end if(w_tc > 0)
        }// end for(c)
    }// end for(t)
  
  free(w_t_pro);
  return;
}





void update_Lambda(int *L, int *Z, double *W, double *Lambda, double *psi,
    double *gamma_R, double *nu_R,
    double *l_R0, double *z_R0,
    double *M, double *m, 
    int S, int T, int C, int G,
    double a_lambda, double b_lambda,
    int *g_fun, double Tmp){
    
    //update Lambda under iid Gamma prior
    
    int g, t, c, s, c_iter, g_iter;
    
    double l_R0_num = *l_R0;
    double z_R0_num = *z_R0;
    
    double lambda_gc_cur, lambda_gc_pro;
    double loglik_cur, loglik_pro, logpost_cur, logpost_pro, logJacobian;
    double B_st_cur, B_tilde_st_cur, lambda_gt0_cur, B_st_pro, B_tilde_st_pro, lambda_gt0_pro;
    double q_st_cur, q_st_pro;
    
    double randomUnif, norm_prop;
    
    
    // sd of M-H normal proposal
    double sd_prop_lambda = 0.1;
    
    
    for(g = 0; g < G; g++){
        for(c = 0; c < C; c++){
            
            //////////////////////////////////////////////////////////////////
            // Update lambda_gc
            //////////////////////////////////////////////////////////////////
            
            lambda_gc_cur = Lambda[c * G + g];
            GetRNGstate();
            norm_prop = rnorm(0.0, sd_prop_lambda);
            PutRNGstate();
            lambda_gc_pro = lambda_gc_cur * exp(norm_prop); 
            
            loglik_cur = 0.0;
            loglik_pro = 0.0;
            
            for(t = 0; t < T; t++){
                for(s = 0; s < S; s++){
                    g_iter = g_fun[s];
                        
                    if(g_iter == g){
                        
                        B_st_cur = 0.0;
                        B_tilde_st_cur = 0.0;
                        lambda_gt0_cur = 0.0;
                        B_st_pro = 0.0;
                        B_tilde_st_pro = 0.0;
                        lambda_gt0_pro = 0.0;
                        
                        for(c_iter = 0; c_iter < C; c_iter++){
                        
                            if(c_iter != c){
                                B_st_cur = B_st_cur + W[c_iter * T + t] * L[c_iter * S + s] * Lambda[c_iter * G + g];
                                B_tilde_st_cur = B_tilde_st_cur + W[c_iter * T + t] * Z[c_iter * S + s] * Lambda[c_iter * G + g];
                                lambda_gt0_cur = lambda_gt0_cur + W[c_iter * T + t] * Lambda[c_iter * G + g];
                            }
                        }// end for(c_iter)
                        
                        B_st_pro = B_st_cur;
                        B_tilde_st_pro = B_tilde_st_cur;
                        lambda_gt0_pro = lambda_gt0_cur;
                            
                        lambda_gt0_cur = lambda_gt0_cur + W[c * T + t] * lambda_gc_cur;
                        lambda_gt0_pro = lambda_gt0_pro + W[c * T + t] * lambda_gc_pro;
                            
                        B_st_cur = B_st_cur + W[c * T + t] * L[c * S + s] * lambda_gc_cur;
                        B_st_pro = B_st_pro + W[c * T + t] * L[c * S + s] * lambda_gc_pro;
                            
                        B_st_cur = B_st_cur + W[C * T + t] * l_R0_num * lambda_gt0_cur;
                        B_st_pro = B_st_pro + W[C * T + t] * l_R0_num * lambda_gt0_pro;
                            
                        B_tilde_st_cur = B_tilde_st_cur + W[c * T + t] * Z[c * S + s] * lambda_gc_cur;
                        B_tilde_st_pro = B_tilde_st_pro + W[c * T + t] * Z[c * S + s] * lambda_gc_pro;
                            
                        B_tilde_st_cur = B_tilde_st_cur + W[C * T + t] * z_R0_num * lambda_gt0_cur;
                        B_tilde_st_pro = B_tilde_st_pro + W[C * T + t] * z_R0_num * lambda_gt0_pro;
                        
                        q_st_cur = B_tilde_st_cur / B_st_cur;
                        q_st_pro = B_tilde_st_pro / B_st_pro;
                        
                        // calculating loglik
                        loglik_cur = loglik_cur - (1 / gamma_R[t] + M[t * S + s]) * log(1 + gamma_R[t] * psi[t] * B_st_cur / 2);
                        
                        loglik_pro = loglik_pro - (1 / gamma_R[t] + M[t * S + s]) * log(1 + gamma_R[t] * psi[t] * B_st_pro / 2);
                        
                        // if M_st > 0
                        if(M[t * S + s] > MIN_COUNT){
                            loglik_cur = loglik_cur + M[t * S + s] * log(gamma_R[t] * psi[t] * B_st_cur / 2);
                        
                            loglik_pro = loglik_pro + M[t * S + s] * log(gamma_R[t] * psi[t] * B_st_pro / 2);
                            // if m_st = 0, done
                        
                            if((M[t * S + s] - m[t * S + s]) > MIN_COUNT){
                                loglik_cur = loglik_cur + lgamma((1 - q_st_cur) / nu_R[t] + M[t * S + s] - m[t * S + s]) - lgamma((1 - q_st_cur) / nu_R[t]);
                            
                                loglik_pro = loglik_pro + lgamma((1 - q_st_pro) / nu_R[t] + M[t * S + s] - m[t * S + s]) - lgamma((1 - q_st_pro) / nu_R[t]);
                            }
                            
                            if(m[t * S + s] > MIN_COUNT){
                                loglik_cur = loglik_cur + lgamma(q_st_cur / nu_R[t] + m[t * S + s]) - lgamma(q_st_cur / nu_R[t]);
                            
                                loglik_pro = loglik_pro + lgamma(q_st_pro / nu_R[t] + m[t * S + s]) - lgamma(q_st_pro / nu_R[t]);
                            }
                        }
                        // end calculating loglik
                        
                    }// end if(g_iter == g)
                }// end for(s)
            }// end for(t)
                
            logpost_cur = loglik_cur / Tmp + (a_lambda - 1) * log(lambda_gc_cur) - b_lambda * lambda_gc_cur;
            logpost_pro = loglik_pro / Tmp + (a_lambda - 1) * log(lambda_gc_pro) - b_lambda * lambda_gc_pro;
            
            // logpost_cur = logpost_cur / Tmp;
            // logpost_pro = logpost_pro / Tmp;
            
            logJacobian = log(lambda_gc_pro) - log(lambda_gc_cur);
            
            GetRNGstate();
            randomUnif = runif(0.0, 1.0);
            PutRNGstate();
            
            if(log(randomUnif) < (logpost_pro - logpost_cur + logJacobian)){
                Lambda[c * G + g] = lambda_gc_pro;
            }
                
                
                
        }// end for(c)
    }// end for(g)
    
    return;
}





// calculate A, B, A_tilde, B_tilde, p, q, lambda0. All (S * T)
void calc_A_B_p_q(double *A, double *B, double *A_tilde, double *B_tilde, double *p, double *q, double *lambda0,
    int *L, int *Z, double *W, double *Lambda, double *l_D0, double *z_D0, double *l_R0, double *z_R0,
    int S, int T, int C, int G, int *g_fun){
    
    int s, t, g, c_iter;
    
    double l_D0_num = *l_D0;
    double z_D0_num = *z_D0;
    double l_R0_num = *l_R0;
    double z_R0_num = *z_R0;
    
    for(s = 0; s < S; s++){
        for(t = 0; t < T; t++){
            g = g_fun[s];
            
            A[t * S + s] = 0.0;
            B[t * S + s] = 0.0;
            
            A_tilde[t * S + s] = 0.0;
            B_tilde[t * S + s] = 0.0;
            
            lambda0[t * S + s] = 0.0;
            
            for(c_iter = 0; c_iter < C; c_iter++){
                A[t * S + s] = A[t * S + s] + W[c_iter * T + t] * L[c_iter * S + s];
                B[t * S + s] = B[t * S + s] + W[c_iter * T + t] * L[c_iter * S + s] * Lambda[c_iter * G + g];
                
                A_tilde[t * S + s] = A_tilde[t * S + s] + W[c_iter * T + t] * Z[c_iter * S + s];
                B_tilde[t * S + s] = B_tilde[t * S + s] + W[c_iter * T + t] * Z[c_iter * S + s] * Lambda[c_iter * G + g];
                
                lambda0[t * S + s] = lambda0[t * S + s] + W[c_iter * T + t] * Lambda[c_iter * G + g];
            }
            
            A[t * S + s] = A[t * S + s] + W[C * T + t] * l_D0_num;
            B[t * S + s] = B[t * S + s] + W[C * T + t] * l_R0_num * lambda0[t * S + s];
            
            A_tilde[t * S + s] = A_tilde[t * S + s] + W[C * T + t] * z_D0_num;
            B_tilde[t * S + s] = B_tilde[t * S + s] + W[C * T + t] * z_R0_num * lambda0[t * S + s];
            
            p[t * S + s] = A_tilde[t * S + s] / A[t * S + s];
            q[t * S + s] = B_tilde[t * S + s] / B[t * S + s];
        }
    }// end for(s)
    // end calculating A, B, p, q
    
    return;
}




void update_phi_psi_gamma_nu_l0_z0(int *L, int *Z, double *W, double *Lambda, double *phi, double *psi,
    double *gamma_D, double *nu_D, double *gamma_R, double *nu_R,
    double *l_D0, double *z_D0, double *l_R0, double *z_R0,
    double *N, double *n, double *M, double *m, 
    int S, int T, int C, int K_min, int K_max, int G,
    double *a_phi, double *b_phi, double *a_psi, double *b_psi,
    double a_gamma_D, double b_gamma_D, double a_gamma_R, double b_gamma_R,
    double a_nu_D, double b_nu_D, double a_nu_R, double b_nu_R,
    int *g_fun, double Tmp){
    
    // update phi_t, psi_t, gamma_D_t, nu_D_t, gamma_R_t, nu_R_t in this one function
    // can avoid repeated calculation of A_st, B_st, p_st, q_st
    
    int s, t;
    
    double phi_t_cur, phi_t_pro;
    double gamma_D_t_cur, gamma_D_t_pro;
    double psi_t_cur, psi_t_pro;
    double gamma_R_t_cur, gamma_R_t_pro;
    double nu_D_t_cur, nu_D_t_pro;
    double nu_R_t_cur, nu_R_t_pro;
    
    double l_D0_cur, l_D0_pro;
    double z_D0_cur, z_D0_pro;
    double l_R0_cur, l_R0_pro;
    double z_R0_cur, z_R0_pro;
    
    double loglik_cur, loglik_pro, logpost_cur, logpost_pro, logJacobian; 
    double norm_prop, randomUnif;
    
    // standard deviation for the normal proporsals
    double sd_prop_phi = 0.3;
    double sd_prop_gamma_D = 0.1;
    double sd_prop_psi = 0.3;
    double sd_prop_gamma_R = 0.1;
    double sd_prop_nu_D = 0.1;
    double sd_prop_nu_R = 0.1;
    double sd_prop_l_D0 = 0.1;
    double sd_prop_z_D0 = 0.1;
    double sd_prop_l_R0 = 0.1;
    double sd_prop_z_R0 = 0.1;
    
    double A_st_pro, A_tilde_st_pro, p_st_pro, B_st_pro, B_tilde_st_pro, q_st_pro;
    
    double *A;
    A = (double *)calloc(S*T, sizeof(double));
    double *B;
    B = (double *)calloc(S*T, sizeof(double));
    double *A_tilde;
    A_tilde = (double *)calloc(S*T, sizeof(double));
    double *B_tilde;
    B_tilde = (double *)calloc(S*T, sizeof(double));
    double *p;
    p = (double *)calloc(S*T, sizeof(double));
    double *q;
    q = (double *)calloc(S*T, sizeof(double));
    
    // should be G*T, but ok for S*T with some same rows
    double *lambda0;
    lambda0 = (double *)calloc(S*T, sizeof(double));
    
    // calculate A, B, p, q. 
    // These will not change for updating phi, psi, gamma and nu
    calc_A_B_p_q(A, B, A_tilde, B_tilde, p, q, lambda0, L, Z, W, Lambda, l_D0, z_D0, l_R0, z_R0, S, T, C, G, g_fun);
    
    
    for(t = 0; t < T; t++){
        
        /////////////////////////////////////////////////////////
        // update phi_t
        /////////////////////////////////////////////////////////
        
        phi_t_cur = phi[t];
        GetRNGstate();
        norm_prop = rnorm(0.0, sd_prop_phi);
        PutRNGstate();
        phi_t_pro = phi_t_cur * exp(norm_prop);
        
        loglik_cur = 0.0;
        loglik_pro = 0.0;
            
        for(s = 0; s < S; s++){    
            loglik_cur = loglik_cur - (1 / gamma_D[t] + N[t * S + s]) * log(1 + gamma_D[t] * phi_t_cur * A[t * S + s] / 2);
            loglik_pro = loglik_pro - (1 / gamma_D[t] + N[t * S + s]) * log(1 + gamma_D[t] * phi_t_pro * A[t * S + s] / 2);
            
            if(N[t * S + s] > MIN_COUNT){
                loglik_cur = loglik_cur + lgamma(N[t * S + s] + 1 / gamma_D[t]) - lgamma(1 / gamma_D[t]) + N[t * S + s] * log(gamma_D[t] * phi_t_cur * A[t * S + s] / 2); 
                
                loglik_pro = loglik_pro + lgamma(N[t * S + s] + 1 / gamma_D[t]) - lgamma(1 / gamma_D[t]) + N[t * S + s] * log(gamma_D[t] * phi_t_pro * A[t * S + s] / 2);
            }
            
        }// end for(s)
        
        logpost_cur = loglik_cur / Tmp + (a_phi[t] - 1) * log(phi_t_cur) - b_phi[t] * phi_t_cur;
        logpost_pro = loglik_pro / Tmp + (a_phi[t] - 1) * log(phi_t_pro) - b_phi[t] * phi_t_pro;
        
        logJacobian = log(phi_t_pro) - log(phi_t_cur);
        
        GetRNGstate();
        randomUnif = runif(0.0, 1.0);
        PutRNGstate();
        
        if(log(randomUnif) < (logpost_pro - logpost_cur + logJacobian)){
            phi[t] = phi_t_pro;
            loglik_cur = loglik_pro;
        }
        
        /////////////////////////////////////////////////////////
        // update gamma_D_t
        /////////////////////////////////////////////////////////
        gamma_D_t_cur = gamma_D[t];
        GetRNGstate();
        norm_prop = rnorm(0.0, sd_prop_gamma_D);
        PutRNGstate();
        gamma_D_t_pro = gamma_D_t_cur * exp(norm_prop);
        
        loglik_pro = 0.0;
        
        for(s = 0; s < S; s++){    
            loglik_pro = loglik_pro - (1 / gamma_D_t_pro + N[t * S + s]) * log(1 + gamma_D_t_pro * phi[t] * A[t * S + s] / 2);
            
            if(N[t * S + s] > MIN_COUNT){
                loglik_pro = loglik_pro + lgamma(N[t * S + s] + 1 / gamma_D_t_pro) - lgamma(1 / gamma_D_t_pro) + N[t * S + s] * log(gamma_D_t_pro * phi[t] * A[t * S + s] / 2);
            }
            
        }// end for(s)
        
        logpost_cur = loglik_cur / Tmp + (a_gamma_D - 1) * log(gamma_D_t_cur) - b_gamma_D * gamma_D_t_cur;
        logpost_pro = loglik_pro / Tmp + (a_gamma_D - 1) * log(gamma_D_t_pro) - b_gamma_D * gamma_D_t_pro;
        
        logJacobian = log(gamma_D_t_pro) - log(gamma_D_t_cur);
        
        GetRNGstate();
        randomUnif = runif(0.0, 1.0);
        PutRNGstate();
        
        if(log(randomUnif) < (logpost_pro - logpost_cur + logJacobian)){
            gamma_D[t] = gamma_D_t_pro;
        }
        
        /////////////////////////////////////////////////////////
        // update psi_t
        /////////////////////////////////////////////////////////
        psi_t_cur = psi[t];
        GetRNGstate();
        norm_prop = rnorm(0.0, sd_prop_psi);
        PutRNGstate();
        psi_t_pro = psi_t_cur * exp(norm_prop);
        
        loglik_cur = 0.0;
        loglik_pro = 0.0;
        
        for(s = 0; s < S; s++){
            loglik_cur = loglik_cur - (1 / gamma_R[t] + M[t * S + s]) * log(1 + gamma_R[t] * psi_t_cur * B[t * S + s] / 2);
            loglik_pro = loglik_pro - (1 / gamma_R[t] + M[t * S + s]) * log(1 + gamma_R[t] * psi_t_pro * B[t * S + s] / 2);
            
            if(M[t * S + s] > MIN_COUNT){
                loglik_cur = loglik_cur + lgamma(M[t * S + s] + 1 / gamma_R[t]) - lgamma(1 / gamma_R[t]) + M[t * S + s] * log(gamma_R[t] * psi_t_cur * B[t * S + s] / 2); 
                
                loglik_pro = loglik_pro + lgamma(M[t * S + s] + 1 / gamma_R[t]) - lgamma(1 / gamma_R[t]) + M[t * S + s] * log(gamma_R[t] * psi_t_pro * B[t * S + s] / 2);
            }
            
        }
        
        logpost_cur = loglik_cur / Tmp + (a_psi[t] - 1) * log(psi_t_cur) - b_psi[t] * psi_t_cur;
        logpost_pro = loglik_pro / Tmp + (a_psi[t] - 1) * log(psi_t_pro) - b_psi[t] * psi_t_pro;
        
        logJacobian = log(psi_t_pro) - log(psi_t_cur);
        
        GetRNGstate();
        randomUnif = runif(0.0, 1.0);
        PutRNGstate();
        
        if(log(randomUnif) < (logpost_pro - logpost_cur + logJacobian)){
            psi[t] = psi_t_pro;
            loglik_cur = loglik_pro;
        }
        
        /////////////////////////////////////////////////////////
        // update gamma_R_t
        /////////////////////////////////////////////////////////
        gamma_R_t_cur = gamma_R[t];
        GetRNGstate();
        norm_prop = rnorm(0.0, sd_prop_gamma_R);
        PutRNGstate();
        gamma_R_t_pro = gamma_R_t_cur * exp(norm_prop);
        
        loglik_pro = 0.0;
        
        for(s = 0; s < S; s++){    
            loglik_pro = loglik_pro - (1 / gamma_R_t_pro + M[t * S + s]) * log(1 + gamma_R_t_pro * psi[t] * B[t * S + s] / 2);
            
            if(M[t * S + s] > MIN_COUNT){
                loglik_pro = loglik_pro + lgamma(M[t * S + s] + 1 / gamma_R_t_pro) - lgamma(1 / gamma_R_t_pro) + M[t * S + s] * log(gamma_R_t_pro * psi[t] * B[t * S + s] / 2);
            }
            
        }// end for(s)
        
        logpost_cur = loglik_cur / Tmp + (a_gamma_R - 1) * log(gamma_R_t_cur) - b_gamma_R * gamma_R_t_cur;
        logpost_pro = loglik_pro / Tmp + (a_gamma_R - 1) * log(gamma_R_t_pro) - b_gamma_R * gamma_R_t_pro;
        
        logJacobian = log(gamma_R_t_pro) - log(gamma_R_t_cur);
        
        // check if this produces the same random uniform every time
        GetRNGstate();
        randomUnif = runif(0.0, 1.0);
        PutRNGstate();
        
        if(log(randomUnif) < (logpost_pro - logpost_cur + logJacobian)){
            gamma_R[t] = gamma_R_t_pro;
        }
        
        
        /////////////////////////////////////////////////////////
        // update nu_D_t
        /////////////////////////////////////////////////////////
        nu_D_t_cur = nu_D[t];
        GetRNGstate();
        norm_prop = rnorm(0.0, sd_prop_nu_D);
        PutRNGstate();
        nu_D_t_pro = nu_D_t_cur * exp(norm_prop);
        
        loglik_cur = 0.0;
        loglik_pro = 0.0;
        
        for(s = 0; s < S; s++){
        
            // if N_st > 0
            if(N[t * S + s] > MIN_COUNT){
                loglik_cur = loglik_cur + lgamma(1 / nu_D_t_cur) - lgamma(1 / nu_D_t_cur + N[t * S + s]);
                        
                loglik_pro = loglik_pro + lgamma(1 / nu_D_t_pro) - lgamma(1 / nu_D_t_pro + N[t * S + s]);
                // if n_st = 0, done
                            
                if((N[t * S + s] - n[t * S + s]) > MIN_COUNT){
                    loglik_cur = loglik_cur + lgamma((1 - p[t * S + s]) / nu_D_t_cur + N[t * S + s] - n[t * S + s]) - lgamma((1 - p[t * S + s]) / nu_D_t_cur);
                            
                    loglik_pro = loglik_pro + lgamma((1 - p[t * S + s]) / nu_D_t_pro + N[t * S + s] - n[t * S + s]) - lgamma((1 - p[t * S + s]) / nu_D_t_pro);
                }
                            
                if(n[t * S + s] > MIN_COUNT){
                    loglik_cur = loglik_cur + lgamma(p[t * S + s] / nu_D_t_cur + n[t * S + s]) - lgamma(p[t * S + s] / nu_D_t_cur);
                            
                    loglik_pro = loglik_pro + lgamma(p[t * S + s] / nu_D_t_pro + n[t * S + s]) - lgamma(p[t * S + s] / nu_D_t_pro);
                }
            }
            
        }// end for(s)
        
        logpost_cur = loglik_cur / Tmp + (a_nu_D - 1) * log(nu_D_t_cur) - b_nu_D * nu_D_t_cur;
        logpost_pro = loglik_pro / Tmp + (a_nu_D - 1) * log(nu_D_t_pro) - b_nu_D * nu_D_t_pro;
        
        logJacobian = log(nu_D_t_pro) - log(nu_D_t_cur);
        
        // check if this produces the same random uniform every time
        GetRNGstate();
        randomUnif = runif(0.0, 1.0);
        PutRNGstate();
        
        if(log(randomUnif) < (logpost_pro - logpost_cur + logJacobian)){
            nu_D[t] = nu_D_t_pro;
        }
        
        
        /////////////////////////////////////////////////////////
        // update nu_R_t
        /////////////////////////////////////////////////////////
        nu_R_t_cur = nu_R[t];
        GetRNGstate();
        norm_prop = rnorm(0.0, sd_prop_nu_R);
        PutRNGstate();
        nu_R_t_pro = nu_R_t_cur * exp(norm_prop);
        
        loglik_cur = 0.0;
        loglik_pro = 0.0;
        
        for(s = 0; s < S; s++){
            
            // if M_st > 0
            if(M[t * S + s] > MIN_COUNT){
                loglik_cur = loglik_cur + lgamma(1 / nu_R_t_cur) - lgamma(1 / nu_R_t_cur + M[t * S + s]);
                        
                loglik_pro = loglik_pro + lgamma(1 / nu_R_t_pro) - lgamma(1 / nu_R_t_pro + M[t * S + s]);
                // if m_st = 0, done
                            
                if((M[t * S + s] - m[t * S + s]) > MIN_COUNT){
                    loglik_cur = loglik_cur + lgamma((1 - q[t * S + s]) / nu_R_t_cur + M[t * S + s] - m[t * S + s]) - lgamma((1 - q[t * S + s]) / nu_R_t_cur);
                            
                    loglik_pro = loglik_pro + lgamma((1 - q[t * S + s]) / nu_R_t_pro + M[t * S + s] - m[t * S + s]) - lgamma((1 - q[t * S + s]) / nu_R_t_pro);
                }
                            
                if(m[t * S + s] > MIN_COUNT){
                    loglik_cur = loglik_cur + lgamma(q[t * S + s] / nu_R_t_cur + m[t * S + s]) - lgamma(q[t * S + s] / nu_R_t_cur);
                            
                    loglik_pro = loglik_pro + lgamma(q[t * S + s] / nu_R_t_pro + m[t * S + s]) - lgamma(q[t * S + s] / nu_R_t_pro);
                }
            }
            
        }
        
        logpost_cur = loglik_cur / Tmp + (a_nu_R - 1) * log(nu_R_t_cur) - b_nu_R * nu_R_t_cur;
        logpost_pro = loglik_pro / Tmp + (a_nu_R - 1) * log(nu_R_t_pro) - b_nu_R * nu_R_t_pro;
        
        logJacobian = log(nu_R_t_pro) - log(nu_R_t_cur);
        
        // check if this produces the same random uniform every time
        GetRNGstate();
        randomUnif = runif(0.0, 1.0);
        PutRNGstate();
        
        if(log(randomUnif) < (logpost_pro - logpost_cur + logJacobian)){
            nu_R[t] = nu_R_t_pro;
        }
        
    }// end for(t)
    
    /////////////////////////////////////////////////////////
    // update l_D0, z_D0, l_R0, z_R0
    /////////////////////////////////////////////////////////
    // update all together so that can save some computational time
    // even if acceptance probability is low, not important parameters anyway..
    
    l_D0_cur = *l_D0;
    GetRNGstate();
    norm_prop = rnorm(0.0, sd_prop_l_D0);
    PutRNGstate();
    l_D0_pro = l_D0_cur * exp(norm_prop);
    
    z_D0_cur = *z_D0;
    GetRNGstate();
    norm_prop = rnorm(0.0, sd_prop_z_D0);
    PutRNGstate();
    z_D0_pro = z_D0_cur * exp(norm_prop);
    
    l_R0_cur = *l_R0;
    GetRNGstate();
    norm_prop = rnorm(0.0, sd_prop_l_R0);
    PutRNGstate();
    l_R0_pro = l_R0_cur * exp(norm_prop);
    
    z_R0_cur = *z_R0;
    GetRNGstate();
    norm_prop = rnorm(0.0, sd_prop_z_R0);
    PutRNGstate();
    z_R0_pro = z_R0_cur * exp(norm_prop);
    
    // otherwise reject directly
    if((l_D0_pro > (double) K_min) && (l_D0_pro < (double) K_max) && (z_D0_pro < l_D0_pro) && (l_R0_pro > (double) K_min) && (l_R0_pro < (double) K_max) && (z_R0_pro < l_R0_pro)){
        loglik_cur = 0.0;
        loglik_pro = 0.0;
        
        for(s = 0; s < S; s++){
            for(t = 0; t < T; t++){
            
                A_st_pro = A[t * S + s] - W[C * T + t] * l_D0_cur + W[C * T + t] * l_D0_pro;
                A_tilde_st_pro = A_tilde[t * S + s] - W[C * T + t] * z_D0_cur + W[C * T + t] * z_D0_pro;
                p_st_pro = A_tilde_st_pro / A_st_pro;
            
                B_st_pro = B[t * S + s] - W[C * T + t] * l_R0_cur * lambda0[t * S + s] + W[C * T + t] * l_R0_pro * lambda0[t * S + s];
                B_tilde_st_pro = B_tilde[t * S + s] - W[C * T + t] * z_R0_cur * lambda0[t * S + s] + W[C * T + t] * z_R0_pro * lambda0[t * S + s]; 
                q_st_pro = B_tilde_st_pro / B_st_pro;
                
                loglik_cur = loglik_cur - (1 / gamma_D[t] + N[t * S + s]) * log(1 + gamma_D[t] * phi[t] * A[t * S + s] / 2) - (1 / gamma_R[t] + M[t * S + s]) * log(1 + gamma_R[t] * psi[t] * B[t * S + s] / 2);
                
                loglik_pro = loglik_pro - (1 / gamma_D[t] + N[t * S + s]) * log(1 + gamma_D[t] * phi[t] * A_st_pro / 2) - (1 / gamma_R[t] + M[t * S + s]) * log(1 + gamma_R[t] * psi[t] * B_st_pro / 2);
                
                
                // if N_st > 0
                if(N[t * S + s] > MIN_COUNT){
                    loglik_cur = loglik_cur + N[t * S + s] * log(gamma_D[t] * phi[t] * A[t * S + s] / 2);
                    
                    loglik_pro = loglik_pro + N[t * S + s] * log(gamma_D[t] * phi[t] * A_st_pro / 2);
                    // if n_st = 0, done
                            
                    if((N[t * S + s] - n[t * S + s]) > MIN_COUNT){
                        loglik_cur = loglik_cur + lgamma((1 - p[t * S + s]) / nu_D[t] + N[t * S + s] - n[t * S + s]) - lgamma((1 - p[t * S + s]) / nu_D[t]);
                        
                        loglik_pro = loglik_pro + lgamma((1 - p_st_pro) / nu_D[t] + N[t * S + s] - n[t * S + s]) - lgamma((1 - p_st_pro) / nu_D[t]);
                    }
                            
                    if(n[t * S + s] > MIN_COUNT){
                        loglik_cur = loglik_cur + lgamma(p[t * S + s] / nu_D[t] + n[t * S + s]) - lgamma(p[t * S + s] / nu_D[t]);
                        
                        loglik_pro = loglik_pro + lgamma(p_st_pro / nu_D[t] + n[t * S + s]) - lgamma(p_st_pro / nu_D[t]);
                    }
                }
                        
                // if M_st > 0
                if(M[t * S + s] > MIN_COUNT){
                    loglik_cur = loglik_cur + M[t * S + s] * log(gamma_R[t] * psi[t] * B[t * S + s] / 2);
                    
                    loglik_pro = loglik_pro + M[t * S + s] * log(gamma_R[t] * psi[t] * B_st_pro / 2);
                    // if m_st = 0, done
                    
                    if((M[t * S + s] - m[t * S + s]) > MIN_COUNT){
                        loglik_cur = loglik_cur + lgamma((1 - q[t * S + s]) / nu_R[t] + M[t * S + s] - m[t * S + s]) - lgamma((1 - q[t * S + s]) / nu_R[t]);
                        
                        loglik_pro = loglik_pro + lgamma((1 - q_st_pro) / nu_R[t] + M[t * S + s] - m[t * S + s]) - lgamma((1 - q_st_pro) / nu_R[t]);
                    }
                        
                    if(m[t * S + s] > MIN_COUNT){
                        loglik_cur = loglik_cur + lgamma(q[t * S + s] / nu_R[t] + m[t * S + s])  - lgamma(q[t * S + s] / nu_R[t]);
                        
                        loglik_pro = loglik_pro + lgamma(q_st_pro / nu_R[t] + m[t * S + s]) - lgamma(q_st_pro / nu_R[t]);
                    }
                }
                // end calculating loglik
            }
        }
        
        logpost_cur = loglik_cur / Tmp - log(l_D0_cur) - log(l_R0_cur);
        logpost_pro = loglik_pro / Tmp - log(l_D0_pro) - log(l_R0_pro);
        
        logJacobian = log(l_D0_pro) + log(z_D0_pro) + log(l_R0_pro) + log(z_R0_pro) - log(l_D0_cur) - log(z_D0_cur) - log(l_R0_cur) - log(z_R0_cur);
        
        GetRNGstate();
        randomUnif = runif(0.0, 1.0);
        PutRNGstate();
        
        if(log(randomUnif) < (logpost_pro - logpost_cur + logJacobian)){
            *l_D0 = l_D0_pro;
            *z_D0 = z_D0_pro;
            *l_R0 = l_R0_pro;
            *z_R0 = z_R0_pro;
        }
        
    } // end if(l_D0 ...)
    // finish updating l_D0, z_D0, l_R0, z_R0

    free(A);
    free(B);
    free(A_tilde);
    free(B_tilde);
    free(p);
    free(q);
    free(lambda0);
    return;
}



//////////////////////////////////////////////////
// Calculate log-likelihood 
//////////////////////////////////////////////////
void calc_loglik(double *loglik, int *L, int *Z, double *W, double *Lambda, 
    double *phi, double *psi,
    double *gamma_D, double *nu_D, double *gamma_R, double *nu_R,
    double *l_D0, double *z_D0, double *l_R0, double *z_R0,
    double *N, double *n, double *M, double *m, 
    int *S, int *T, int *C, int *G,
    int *g_fun){
    
    *loglik = 0.0;
    
    int s, t, c, g;
    
    double A_st, B_st, A_tilde_st, B_tilde_st, p_st, q_st, lambda_gt0;
    
    for(t = 0; t < (*T); t++){
        for(s = 0; s < (*S); s++){
            g = g_fun[s];
            A_st = 0;
            B_st = 0;
            A_tilde_st = 0;
            B_tilde_st = 0;
            lambda_gt0 = 0;
                        
            for(c = 0; c < (*C); c++){
                A_st = A_st + W[c * (*T) + t] * L[c * (*S) + s];
                B_st = B_st + W[c * (*T) + t] * L[c * (*S) + s] * Lambda[c * (*G) + g];
                A_tilde_st = A_tilde_st + W[c * (*T) + t] * Z[c * (*S) + s];
                B_tilde_st = B_tilde_st + W[c * (*T) + t] * Z[c * (*S) + s]* Lambda[c * (*G) + g];
                lambda_gt0 = lambda_gt0 + W[c * (*T) + t] * Lambda[c * (*G) + g];
            }
            
            A_st = A_st + W[(*C) * (*T) + t] * (*l_D0);
            A_tilde_st = A_tilde_st + W[(*C) * (*T) + t] * (*z_D0);
            p_st = A_tilde_st / A_st;
            B_st = B_st + W[(*C) * (*T) + t] * (*l_R0) * lambda_gt0;
            B_tilde_st = B_tilde_st + W[(*C) * (*T) + t] * (*z_R0) * lambda_gt0;
            q_st = B_tilde_st / B_st;
            
            // calculating loglik[s, t]
            (*loglik) = (*loglik) - (1 / gamma_D[t] + N[t * (*S) + s]) * log(1 + gamma_D[t] * phi[t] * A_st / 2) - (1 / gamma_R[t] + M[t * (*S) + s]) * log(1 + gamma_R[t] * psi[t] * B_st / 2);
            
            // if N_st > 0
            if(N[t * (*S) + s] > MIN_COUNT){
                (*loglik) = (*loglik) + lgamma(N[t * (*S) + s] + 1 / gamma_D[t]) - lgamma(1 / gamma_D[t]) + N[t * (*S) + s] * log(gamma_D[t] * phi[t] * A_st / 2) + lgamma(1 / nu_D[t]) - lgamma(1 / nu_D[t] + N[t * (*S) + s]);
                
                // if n_st = 0, done
                        
                if((N[t * (*S) + s] - n[t * (*S) + s]) > MIN_COUNT){
                    (*loglik) = (*loglik) + lgamma((1 - p_st) / nu_D[t] + N[t * (*S) + s] - n[t * (*S) + s]) - lgamma((1 - p_st) / nu_D[t]);
                }
                
                if(n[t * (*S) + s] > MIN_COUNT){
                    (*loglik) = (*loglik) + lgamma(p_st / nu_D[t] + n[t * (*S) + s]) - lgamma(p_st / nu_D[t]);
                }
            }
                        
            // if M_st > 0
            if(M[t * (*S) + s] > MIN_COUNT){
                (*loglik) = (*loglik) + lgamma(M[t * (*S) + s] + 1 / gamma_R[t]) - lgamma(1 / gamma_R[t]) + M[t * (*S) + s] * log(gamma_R[t] * psi[t] * B_st / 2) + lgamma(1 / nu_R[t]) - lgamma(1 / nu_R[t] + M[t * (*S) + s]);
                
                // if m_st = 0, done
                
                if((M[t * (*S) + s] - m[t * (*S) + s]) > MIN_COUNT){
                    (*loglik) = (*loglik) + lgamma((1 - q_st) / nu_R[t] + M[t * (*S) + s] - m[t * (*S) + s]) - lgamma((1 - q_st) / nu_R[t]);
                }
                
                if(m[t * (*S) + s] > MIN_COUNT){
                    (*loglik) = (*loglik) + lgamma(q_st / nu_R[t] + m[t * (*S) + s]) - lgamma(q_st / nu_R[t]);
                }
            }
            // end calculating loglik[s, t]
            
        }
    }
    
    if(isnan((*loglik))){
        (*loglik) = -INFINITY;
    }
    
    return;
}




void calc_logpost(double *loglik, double *logpost, 
    int *L, int *Z, double *W, double *Lambda, 
    double *pai, double *zeta, double *phi, double *psi,
    double *gamma_D, double *nu_D, double *gamma_R, double *nu_R,
    double *l_D0, double *z_D0, double *l_R0, double *z_R0,
    double *N, double *n, double *M, double *m, 
    int *S, int *T, int *C, int *G, int *K_min, int *K_max,
    double *a_w, double *b_w, double *d, double *d0,
    double *a_lambda, double *b_lambda,
    double *a_pai, double *b_pai, double *a_zeta, double *b_zeta,
    double *a_phi, double *b_phi, double *a_psi, double *b_psi,
    double *a_gamma_D, double *b_gamma_D, double *a_gamma_R, double *b_gamma_R,
    double *a_nu_D, double *b_nu_D, double *a_nu_R, double *b_nu_R,
    int *g_fun){
    
    
    int s, t, c, g, k1;
    
    double logprior = 0.0;
    
    
    double temp1, temp2;
    
    // likelihood
    calc_loglik(loglik, L, Z, W, Lambda, phi, psi, gamma_D, nu_D, gamma_R, nu_R, l_D0, z_D0, l_R0, z_R0, N, n, M, m, S, T, C, G, g_fun);
    
    
    // prior for L, Z, pai and zeta
    for(c = 1; c < (*C); c++){
        for(s = 0; s < (*S); s++){
            logprior = logprior + abs(L[c * (*S) + s] - 2) * log(1 - pai[c]);
            logprior = logprior + Z[c * (*S) + s] * log(1 - zeta[c]);
            temp2 = 0;
            for(k1 = 0; k1 <= L[c * (*S) + s]; k1++){
                temp2 = temp2 + zeta[c] * pow(1 - zeta[c], k1);
            }
            logprior = logprior - log(temp2);
        }
        
        temp1 = 0;
        for(k1 = (*K_min); k1 <= (*K_max); k1++){
            temp1 = temp1 + pai[c] * pow(1 - pai[c], abs(k1 - 2));
        }
        logprior = logprior + (*S) * log(pai[c]) - (*S) * log(temp1);
        logprior = logprior + (*S) * log(zeta[c]);
        
        logprior = logprior + ((*a_pai) - 1) * log(pai[c]) + ((*b_pai) - 1) * log(1 - pai[c]);
        logprior = logprior + ((*a_zeta) - 1) * log(zeta[c]) + ((*b_zeta) - 1) * log(1 - zeta[c]);
    }
    
    
    // prior for W, phi, psi, gamma_D, nu_D, gamma_R, nu_R
    for(t = 0; t < (*T); t++){
        logprior = logprior + ((*a_w) - 1) * log(W[0 * (*T) + t]) + ((*b_w) - 1 - ((*C) - 1) * ((*d) - 1) - ((*d0) - 1)) * log(1 - W[0 * (*T) + t]);
        for(c = 1; c < (*C); c++){
            logprior = logprior + ((*d) - 1) * log(W[c * (*T) + t]);
        }
        logprior = logprior + ((*d0) - 1) * log(W[(*C) * (*T) + t]);
        
        logprior = logprior + (a_phi[t] - 1) * log(phi[t]) - b_phi[t] * phi[t];
        logprior = logprior + (a_psi[t] - 1) * log(psi[t]) - b_psi[t] * psi[t];
        logprior = logprior + ((*a_gamma_D) - 1) * log(gamma_D[t]) - (*b_gamma_D) * gamma_D[t];
        logprior = logprior + ((*a_nu_D) - 1) * log(nu_D[t]) - (*b_nu_D) * nu_D[t];
        logprior = logprior + ((*a_gamma_R) - 1) * log(gamma_R[t]) - (*b_gamma_R) * gamma_R[t];
        logprior = logprior + ((*a_nu_R) - 1) * log(nu_R[t]) - (*b_nu_R) * nu_R[t];
        
    }
    
    // prior for l_D0, z_D0, l_R0, z_R0 
    logprior = logprior - log(*l_D0) - log(*l_R0);
    
    // prior for Lambda
    for(g = 0; g < (*G); g++){
        for(c = 0; c < (*C); c++){
            logprior = logprior + ((*a_lambda) - 1) * log(Lambda[c * (*G) + g]) - (*b_lambda) * Lambda[c * (*G) + g];
        }
    }
    
    (*logpost) = (*loglik) + logprior;
    
    if(isnan((*logpost))){
        (*logpost) = -INFINITY;
    }
    
    return;
}





// swap the state of two Markov chains
void PT_swap(int **L1, int **Z1, double **W1, double **Lambda1, double **pai1, double **zeta1,
    double **phi1, double **psi1, double **gamma_D1, double **nu_D1, double **gamma_R1, double **nu_R1,
    double *l_D01, double *z_D01, double *l_R01, double *z_R01, double *loglik1, double *lpost1, 
    int **L2, int **Z2, double **W2, double **Lambda2, double **pai2, double **zeta2,
    double **phi2, double **psi2, double **gamma_D2, double **nu_D2, double **gamma_R2, double **nu_R2,
    double *l_D02, double *z_D02, double *l_R02, double *z_R02, double *loglik2, double *lpost2){
    
    
    int *temp_int;
    double *temp_double;
    double temp_double2;
    
    temp_int = *L1;
    *L1 = *L2;
    *L2 = temp_int;
    
    temp_int = *Z1;
    *Z1 = *Z2;
    *Z2 = temp_int;
    
    temp_double = *W1;
    *W1 = *W2;
    *W2 = temp_double;
    
    temp_double = *Lambda1;
    *Lambda1 = *Lambda2;
    *Lambda2 = temp_double;
    
    temp_double = *pai1;
    *pai1 = *pai2;
    *pai2 = temp_double;
    
    temp_double = *zeta1;
    *zeta1 = *zeta2;
    *zeta2 = temp_double;
    
    temp_double = *phi1;
    *phi1 = *phi2;
    *phi2 = temp_double;
    
    temp_double = *psi1;
    *psi1 = *psi2;
    *psi2 = temp_double;
    
    temp_double = *gamma_D1;
    *gamma_D1 = *gamma_D2;
    *gamma_D2 = temp_double;
    
    temp_double = *nu_D1;
    *nu_D1 = *nu_D2;
    *nu_D2 = temp_double;
    
    temp_double = *gamma_R1;
    *gamma_R1 = *gamma_R2;
    *gamma_R2 = temp_double;
    
    temp_double = *nu_R1;
    *nu_R1 = *nu_R2;
    *nu_R2 = temp_double;
    
    temp_double2 = *l_D01;
    *l_D01 = *l_D02;
    *l_D02 = temp_double2;
    
    temp_double2 = *z_D01;
    *z_D01 = *z_D02;
    *z_D02 = temp_double2;
    
    temp_double2 = *l_R01;
    *l_R01 = *l_R02;
    *l_R02 = temp_double2;
    
    temp_double2 = *z_R01;
    *z_R01 = *z_R02;
    *z_R02 = temp_double2;
    
    temp_double2 = *loglik1;
    *loglik1 = *loglik2;
    *loglik2 = temp_double2;
    
    temp_double2 = *lpost1;
    *lpost1 = *lpost2;
    *lpost2 = temp_double2;
    
    return;
}




////////////////////////////////////////////////////////////////////
// Main MCMC function for R
////////////////////////////////////////////////////////////////////
    
void RNDClone_MCMC(int *L_PTR, int *Z_PTR, double *W_PTR, double *Lambda_PTR, 
    double *pai_PTR, double *zeta_PTR, double *phi_PTR, double *psi_PTR,
    double *gamma_D_PTR, double *nu_D_PTR, double *gamma_R_PTR, double *nu_R_PTR,
    double *l_D0_PTR, double *z_D0_PTR, double *l_R0_PTR, double *z_R0_PTR,
    double *loglik_PTR, double *logpost_PTR,
    double *N, double *n, double *M, double *m,
    int *S, int *T, int *C, int *G, int *K_min, int *K_max,
    double *a_w, double *b_w, double *d, double *d0,
    double *a_lambda, double *b_lambda,
    double *a_pai, double *b_pai, double *a_zeta, double *b_zeta,
    double *a_phi, double *b_phi, double *a_psi, double *b_psi,
    double *a_gamma_D, double *b_gamma_D, double *a_gamma_R, double *b_gamma_R,
    double *a_nu_D, double *b_nu_D, double *a_nu_R, double *b_nu_R,
    int *niter, int *g_fun, double *Delta, int *nThread){
    
    // main MCMC function
    // L_PTR, Z_PTR, W_PTR, Lambda_PTR, ...., z_R0_PTR: parameters
    // N, n, M, m: data
    // S, T, C, G, K: constants, # of SNVs, # of samples, # of subclones (random), # of genes, maximum copy #
    // a_w, ..., b_nu_R: hyperparameters (fixed)
    // niter: number of MCMC iterations 
    // g_fun: a vector storing the corresponding gene g of SNV s, g_fun(SNV s) = gene g
    // Delta: nThread-dimensional vector of tempratures.
    // nThread: number of temperatures to run parallel tempering (for the PT version, it's equal to number of threads)
    
    
    int s, t, c, g;
    // i: iteration, from 1 to niter
    int i, rank, rank_partner;
	double u, lalpha;
    
    int S_num = *S;
    int T_num = *T;
    int C_num = *C;
    int G_num = *G;
    int K_min_num = *K_min;
    int K_max_num = *K_max;
    
    double a_w_num = *a_w;
    double b_w_num = *b_w;
    double d_num = *d;
    double d0_num = *d0;
    double a_lambda_num = *a_lambda;
    double b_lambda_num = *b_lambda;
    double a_pai_num = *a_pai;
    double b_pai_num = *b_pai;
    double a_zeta_num = *a_zeta;
    double b_zeta_num = *b_zeta;
    double a_gamma_D_num = *a_gamma_D;
    double b_gamma_D_num = *b_gamma_D;
    double a_gamma_R_num = *a_gamma_R;
    double b_gamma_R_num = *b_gamma_R;
    double a_nu_D_num = *a_nu_D;
    double b_nu_D_num = *b_nu_D;
    double a_nu_R_num = *a_nu_R;
    double b_nu_R_num = *b_nu_R;
    
    int niter_num = *niter;
    int nThread_num = *nThread;
    ////////////////////////////////////////////////////////////////////
    // Create PT arrays for saving parameter values for each temprature
    // Memory allocation
    ////////////////////////////////////////////////////////////////////
    
    int *L_PT[nThread_num]; // S * C
    int *Z_PT[nThread_num]; // S * C
    double *W_PT[nThread_num]; // T * (C+1)
    double *Lambda_PT[nThread_num]; // G * C
    
    double *pai_PT[nThread_num]; // C
    double *zeta_PT[nThread_num]; // C
    
    double *phi_PT[nThread_num]; // T
    double *psi_PT[nThread_num]; // T
    double *gamma_D_PT[nThread_num]; // T
    double *nu_D_PT[nThread_num]; // T
    double *gamma_R_PT[nThread_num]; // T
    double *nu_R_PT[nThread_num]; // T
    
    double l_D0_PT[nThread_num]; // 1
    double z_D0_PT[nThread_num]; // 1
    double l_R0_PT[nThread_num]; // 1
    double z_R0_PT[nThread_num]; // 1
    
    double loglik_PT[nThread_num];
    double logpost_PT[nThread_num];	// 1
    
	
	for(rank = 0; rank < nThread_num; rank++){
        
        L_PT[rank] = (int*)(calloc(S_num*C_num, sizeof(int)));
        Z_PT[rank] = (int*)(calloc(S_num*C_num, sizeof(int)));
		W_PT[rank] = (double*)(calloc(T_num*(C_num+1), sizeof(double)));
        Lambda_PT[rank] = (double*)(calloc(G_num*C_num, sizeof(double)));
        pai_PT[rank] = (double*)(calloc(C_num, sizeof(double)));
        zeta_PT[rank] = (double*)(calloc(C_num, sizeof(double)));
        
        phi_PT[rank] = (double*)(calloc(T_num, sizeof(double)));
        psi_PT[rank] = (double*)(calloc(T_num, sizeof(double)));
        gamma_D_PT[rank] = (double*)(calloc(T_num, sizeof(double)));
        nu_D_PT[rank] = (double*)(calloc(T_num, sizeof(double)));
        gamma_R_PT[rank] = (double*)(calloc(T_num, sizeof(double)));
        nu_R_PT[rank] = (double*)(calloc(T_num, sizeof(double)));
        
	}  
	
    ////////////////////////////////////////////////////////////////////
    // Setting initial values of the PT arrays = input values
    ////////////////////////////////////////////////////////////////////
    
    for(rank = 0; rank < nThread_num; rank++){
		
        // Setting L and Z
        for(s = 0; s < S_num; s++){
            for(c = 0; c < C_num; c++){
                L_PT[rank][c * S_num + s] = L_PTR[rank * S_num * C_num + c * S_num + s];
                Z_PT[rank][c * S_num + s] = Z_PTR[rank * S_num * C_num + c * S_num + s];
            }
        }
        
        // Setting W
        for(t = 0; t < T_num; t++){
            for(c = 0; c <= C_num; c++){
                W_PT[rank][c * T_num + t] = W_PTR[rank * T_num * (C_num + 1) + c * T_num + t];
            }
        }
        
        // Setting Lambda
        for(g = 0; g < G_num; g++){
            for(c = 0; c < C_num; c++){
                Lambda_PT[rank][c * G_num + g] = Lambda_PTR[rank * G_num * C_num + c * G_num + g];
            }
        }
        
        // Setting pai and zeta
        for(c = 0; c < C_num; c++){
            pai_PT[rank][c] = pai_PTR[rank * C_num + c];
            zeta_PT[rank][c] = zeta_PTR[rank * C_num + c];
        }
        
        // Setting phi, psi, gamma_D, nu_D, gamma_R, nu_R
        for(t = 0; t < T_num; t++){
            phi_PT[rank][t] = phi_PTR[rank * T_num + t];
            psi_PT[rank][t] = psi_PTR[rank * T_num + t];
            gamma_D_PT[rank][t] = gamma_D_PTR[rank * T_num + t];
            nu_D_PT[rank][t] = nu_D_PTR[rank * T_num + t];
            gamma_R_PT[rank][t] = gamma_R_PTR[rank * T_num + t];
            nu_R_PT[rank][t] = nu_R_PTR[rank * T_num + t];
        }
        
        // Setting l_D0, z_D0, l_R0, z_R0
        l_D0_PT[rank] = l_D0_PTR[rank];
        z_D0_PT[rank] = z_D0_PTR[rank];
        l_R0_PT[rank] = l_R0_PTR[rank];
        z_R0_PT[rank] = z_R0_PTR[rank];
        
        loglik_PT[rank] = loglik_PTR[rank];
	    logpost_PT[rank] = logpost_PTR[rank];
	
	} // end for(rank)
    
    ///////////////////////////////////////////////////////////////////////////////
	// MCMC update.. Using for loop for parallel tempering
	///////////////////////////////////////////////////////////////////////////////
	
	for(i = 0; i < niter_num; i++){
	    for(rank = 0; rank < nThread_num; rank++){
	        
	        update_L(L_PT[rank], Z_PT[rank], W_PT[rank], Lambda_PT[rank], pai_PT[rank], zeta_PT[rank], phi_PT[rank], psi_PT[rank], gamma_D_PT[rank], nu_D_PT[rank], gamma_R_PT[rank], nu_R_PT[rank], &l_D0_PT[rank], &z_D0_PT[rank], &l_R0_PT[rank], &z_R0_PT[rank], N, n, M, m, S_num, T_num, C_num, K_min_num, K_max_num, G_num, g_fun, Delta[rank]);
    
            update_pai(L_PT[rank], pai_PT[rank], S_num, C_num, K_min_num, K_max_num, a_pai_num, b_pai_num);
    
            update_Z(L_PT[rank], Z_PT[rank], W_PT[rank], Lambda_PT[rank], zeta_PT[rank], nu_D_PT[rank], nu_R_PT[rank], &l_D0_PT[rank], &z_D0_PT[rank], &l_R0_PT[rank], &z_R0_PT[rank], N, n, M, m, S_num, T_num, C_num, K_min_num, K_max_num, G_num, g_fun, Delta[rank]);
            
            update_zeta(L_PT[rank], Z_PT[rank], zeta_PT[rank], S_num, C_num, a_zeta_num, b_zeta_num);
    
            update_W(L_PT[rank], Z_PT[rank], W_PT[rank], Lambda_PT[rank], phi_PT[rank], psi_PT[rank], gamma_D_PT[rank], nu_D_PT[rank], gamma_R_PT[rank], nu_R_PT[rank], &l_D0_PT[rank], &z_D0_PT[rank], &l_R0_PT[rank], &z_R0_PT[rank], N, n, M, m, S_num, T_num, C_num, G_num, a_w_num, b_w_num, d_num, d0_num, g_fun, Delta[rank]);
            
            update_Lambda(L_PT[rank], Z_PT[rank], W_PT[rank], Lambda_PT[rank], psi_PT[rank], gamma_R_PT[rank], nu_R_PT[rank], &l_R0_PT[rank], &z_R0_PT[rank], M, m, S_num, T_num, C_num, G_num, a_lambda_num, b_lambda_num, g_fun, Delta[rank]);
            
            update_phi_psi_gamma_nu_l0_z0(L_PT[rank], Z_PT[rank], W_PT[rank], Lambda_PT[rank], phi_PT[rank], psi_PT[rank], gamma_D_PT[rank], nu_D_PT[rank], gamma_R_PT[rank], nu_R_PT[rank], &l_D0_PT[rank], &z_D0_PT[rank], &l_R0_PT[rank], &z_R0_PT[rank], N, n, M, m, S_num, T_num, C_num, K_min_num, K_max_num, G_num, a_phi, b_phi, a_psi, b_psi, a_gamma_D_num, b_gamma_D_num, a_gamma_R_num, b_gamma_R_num, a_nu_D_num, b_nu_D_num, a_nu_R_num, b_nu_R_num, g_fun, Delta[rank]);
            
            calc_logpost(&loglik_PT[rank], &logpost_PT[rank], L_PT[rank], Z_PT[rank], W_PT[rank], Lambda_PT[rank], pai_PT[rank], zeta_PT[rank], phi_PT[rank], psi_PT[rank], gamma_D_PT[rank], nu_D_PT[rank], gamma_R_PT[rank], nu_R_PT[rank], &l_D0_PT[rank], &z_D0_PT[rank], &l_R0_PT[rank], &z_R0_PT[rank], N, n, M, m, S, T, C, G, K_min, K_max, a_w, b_w, d, d0, a_lambda, b_lambda, a_pai, b_pai, a_zeta, b_zeta, a_phi, b_phi, a_psi, b_psi, a_gamma_D, b_gamma_D, a_gamma_R, b_gamma_R, a_nu_D, b_nu_D, a_nu_R, b_nu_R, g_fun);
            
            
	    }
	    
	    for(rank = 0; rank < (nThread_num - 1); rank++){
	    
	        GetRNGstate();
            u = runif(0.0, 1.0);
            PutRNGstate();
	        
	        rank_partner = rank + 1;
	        lalpha = (1/Delta[rank] - 1/Delta[rank_partner]) * (loglik_PT[rank_partner] - loglik_PT[rank]);
	        
	        if(log(u) < lalpha){
                PT_swap(&L_PT[rank], &Z_PT[rank], &W_PT[rank], &Lambda_PT[rank], &pai_PT[rank], &zeta_PT[rank], &phi_PT[rank], &psi_PT[rank], &gamma_D_PT[rank], &nu_D_PT[rank], &gamma_R_PT[rank], &nu_R_PT[rank], &l_D0_PT[rank], &z_D0_PT[rank], &l_R0_PT[rank], &z_R0_PT[rank], &loglik_PT[rank], &logpost_PT[rank], &L_PT[rank_partner], &Z_PT[rank_partner], &W_PT[rank_partner], &Lambda_PT[rank_partner], &pai_PT[rank_partner], &zeta_PT[rank_partner], &phi_PT[rank_partner], &psi_PT[rank_partner], &gamma_D_PT[rank_partner], &nu_D_PT[rank_partner], &gamma_R_PT[rank_partner], &nu_R_PT[rank_partner], &l_D0_PT[rank_partner], &z_D0_PT[rank_partner], &l_R0_PT[rank_partner], &z_R0_PT[rank_partner], &loglik_PT[rank_partner], &logpost_PT[rank_partner]);
            }
	    }
	
	} // end for(i)
	
    ///////////////////////////////////////////////////////////////////////////////
	// End MCMC update
	///////////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////
    // Setting updated parameter values = current PT array values
    ////////////////////////////////////////////////////////////////////
    
    
    for(rank = 0; rank < nThread_num; rank++){
		
        // Setting L and Z
        for(s = 0; s < S_num; s++){
            for(c = 0; c < C_num; c++){
                L_PTR[rank * S_num * C_num + c * S_num + s] = L_PT[rank][c * S_num + s];
                Z_PTR[rank * S_num * C_num + c * S_num + s] = Z_PT[rank][c * S_num + s];
            }
        }
        
        // Setting W
        for(t = 0; t < T_num; t++){
            for(c = 0; c <= C_num; c++){
                W_PTR[rank * T_num * (C_num + 1) + c * T_num + t] = W_PT[rank][c * T_num + t];
            }
        }
        
        // Setting Lambda
        for(g = 0; g < G_num; g++){
            for(c = 0; c < C_num; c++){
                Lambda_PTR[rank * G_num * C_num + c * G_num + g] = Lambda_PT[rank][c * G_num + g];
            }
        }
        
        // Setting pai and zeta
        for(c = 0; c < C_num; c++){
            pai_PTR[rank * C_num + c] = pai_PT[rank][c];
            zeta_PTR[rank * C_num + c] = zeta_PT[rank][c];
        }
        
        // Setting phi, psi, gamma_D, nu_D, gamma_R, nu_R
        for(t = 0; t < T_num; t++){
            phi_PTR[rank * T_num + t] = phi_PT[rank][t];
            psi_PTR[rank * T_num + t] = psi_PT[rank][t];
            gamma_D_PTR[rank * T_num + t] = gamma_D_PT[rank][t];
            nu_D_PTR[rank * T_num + t] = nu_D_PT[rank][t];
            gamma_R_PTR[rank * T_num + t] = gamma_R_PT[rank][t];
            nu_R_PTR[rank * T_num + t] = nu_R_PT[rank][t];
        }
        
        // Setting l_D0, z_D0, l_R0, z_R0
        l_D0_PTR[rank] = l_D0_PT[rank];
        z_D0_PTR[rank] = z_D0_PT[rank];
        l_R0_PTR[rank] = l_R0_PT[rank];
        z_R0_PTR[rank] = z_R0_PT[rank];
        
        loglik_PTR[rank] = loglik_PT[rank];
	    logpost_PTR[rank] = logpost_PT[rank];
	
	} // end for(rank)
	
	
	////////////////////////////////////////////////////////////////////
    // Memory free
    ////////////////////////////////////////////////////////////////////
	
    
	for(rank = 0; rank < nThread_num; rank++){
	    free(L_PT[rank]);
		free(Z_PT[rank]);
		free(W_PT[rank]);
		free(Lambda_PT[rank]);
		free(pai_PT[rank]);
		free(zeta_PT[rank]);
		
		free(phi_PT[rank]);
		free(psi_PT[rank]);
		free(gamma_D_PT[rank]);
		free(nu_D_PT[rank]);
		free(gamma_R_PT[rank]);
		free(nu_R_PT[rank]);
	}
    
    return;
}






/////////////////////////////////////////////////////////////////////////////////////////////////////////
// DNA-based subclone reconstruction. 
// Only update L, Z, W
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void update_L_DNA(int *L, int *Z, double *W, double *pai, double *zeta, double *phi,
    double *gamma_D, double *nu_D,
    double *l_D0, double *z_D0,
    double *N, double *n,
    int S, int T, int C, int K_min, int K_max,
    double Tmp){
    
    // element-wise update L, sequentially update l_sc
    
    int s, t, c, k, c_iter, z_sc, k1, max_z_sc_K_min;
    
    double l_D0_num = *l_D0;
    double z_D0_num = *z_D0;
    
    double temp1;
    double A_st, A_tilde_st, p_st;
    double loglik, logpost;
    double max_log_prob_l_sc = 0.0, sum_prob_l_sc;
    double randomUnif;
    
    
    double *prob_l_sc; 
    prob_l_sc = (double *)calloc((K_max + 1), sizeof(double));
    
    for(s = 0; s < S; s++){
        for(c = 1; c < C; c++){
            // update l_sc (c >= 1, c = 0 corresponds to the normal subclone)
            z_sc = Z[c * S + s];
            max_z_sc_K_min = max(z_sc, K_min);
            
            if(max_z_sc_K_min == K_max){
                L[c * S + s] = K_max;
            }
            else{
                for(k = max_z_sc_K_min; k <= K_max; k++){
                    // calculate p(l_sc = k | ...)
                    loglik = 0;
                    for(t = 0; t < T; t++){
                        A_st = 0;
                        A_tilde_st = 0;
                        
                        for(c_iter = 0; c_iter < C; c_iter++){
                            if(c_iter != c){
                                A_st = A_st + W[c_iter * T + t] * L[c_iter * S + s];
                            }
                            else{
                                A_st = A_st + W[c_iter * T + t] * k;
                            }
                            
                            A_tilde_st = A_tilde_st + W[c_iter * T + t] * Z[c_iter * S + s];
                        }
                    
                        A_st = A_st + W[C * T + t] * l_D0_num;
                        A_tilde_st = A_tilde_st + W[C * T + t] * z_D0_num;
                        p_st = A_tilde_st / A_st;
                        
                        // if N_st = 0 done
                        loglik = loglik - (1 / gamma_D[t] + N[t * S + s]) * log(1 + gamma_D[t] * phi[t] * A_st / 2);
                        
                        // if N_st > 0
                        if(N[t * S + s] > MIN_COUNT){
                            loglik = loglik + N[t * S + s] * log(gamma_D[t] * phi[t] * A_st / 2);
                            // if n_st = 0, done
                            
                            if((N[t * S + s] - n[t * S + s]) > MIN_COUNT){
                                loglik = loglik + lgamma((1 - p_st) / nu_D[t] + N[t * S + s] - n[t * S + s]) - lgamma((1 - p_st) / nu_D[t]);
                            }
                            
                            if(n[t * S + s] > MIN_COUNT){
                                loglik = loglik + lgamma(p_st / nu_D[t] + n[t * S + s]) - lgamma(p_st / nu_D[t]);
                            }
                        }
                        
                    }
                    
                    temp1 = 0;
                    for(k1 = 0; k1 <= k; k1++){
                      temp1 = temp1 + pow(1 - zeta[c], k1);
                    }
                    
                    // setting pai[c, 0, .. K_min - 1] = 0
                    logpost = loglik / Tmp + (double) abs(k - 2) * log(1 - pai[c]) - log(temp1);
                    
                    prob_l_sc[k] = logpost;
                    if(k == max_z_sc_K_min){
                        max_log_prob_l_sc = logpost;
                    }
                    else{
                        if(max_log_prob_l_sc < logpost) max_log_prob_l_sc = logpost;
                    }
                } // end for(k). Now we have p(l_sc = k | ...) for all k, stored in prob_l_sc
                
                sum_prob_l_sc = 0;
                
                for(k = max_z_sc_K_min; k <= K_max; k++){
                    prob_l_sc[k] = prob_l_sc[k] - max_log_prob_l_sc;
                    prob_l_sc[k] = exp(prob_l_sc[k]);
                    sum_prob_l_sc = sum_prob_l_sc + prob_l_sc[k];
                }
                
                for(k = 0; k < max_z_sc_K_min; k++){
                    prob_l_sc[k] = 0;
                }
                
                GetRNGstate();
                randomUnif = runif(0.0, 1.0);
                PutRNGstate();
  
                // l_sc ~ Discrete(0:K, prob_l_sc). Sample l_sc.
                k = max_z_sc_K_min;
                prob_l_sc[k] = prob_l_sc[k] / sum_prob_l_sc;
                while(randomUnif > prob_l_sc[k]){
  	                k++;
  	                prob_l_sc[k] = prob_l_sc[k] / sum_prob_l_sc + prob_l_sc[k - 1];
                }
  
                L[c * S + s] = k;
            
            }// end if(max_z_sc_K_min == K_max)
        }// end for(c)
    }// end for(s)
    
    free(prob_l_sc);
    return;
}





void update_Z_DNA(int *L, int *Z, double *W, double *zeta,
    double *nu_D,
    double *l_D0, double *z_D0,
    double *N, double *n,
    int S, int T, int C, int K_min, int K_max,
    double Tmp){
    
    int s, t, c, k, c_iter, l_sc;
    
    double l_D0_num = *l_D0;
    double z_D0_num = *z_D0;
    
    double A_st, A_tilde_st, p_st;
    double loglik, logpost;
    double max_log_prob_z_sc, sum_prob_z_sc;
    double randomUnif;
    
    double *prob_z_sc; 
    prob_z_sc = (double *)calloc((K_max + 1), sizeof(double));
    
    for(s = 0; s < S; s++){
        for(c = 1; c < C; c++){
            // update z_sc (c >= 1, c = 0 corresponds to the normal subclone)
            
            l_sc = L[c * S + s];
            // min(l_sc, K_max) is always l_sc so no need
            
            if(l_sc == 0){
                Z[c * S + s] = 0;
            }
            else{
                for(k = 0; k <= l_sc; k++){
                    // calculate p(z_sc = k | ...)
                    loglik = 0;
                    for(t = 0; t < T; t++){
                        A_st = 0;
                        A_tilde_st = 0;
                        
                        for(c_iter = 0; c_iter < C; c_iter++){
                            if(c_iter != c){
                                A_tilde_st = A_tilde_st + W[c_iter * T + t] * Z[c_iter * S + s];
                            }
                            else{
                                A_tilde_st = A_tilde_st + W[c_iter * T + t] * k;
                            }
                            A_st = A_st + W[c_iter * T + t] * L[c_iter * S + s];
                        }
                    
                        A_st = A_st + W[C * T + t] * l_D0_num;
                        A_tilde_st = A_tilde_st + W[C * T + t] * z_D0_num;
                        p_st = A_tilde_st / A_st;
                        
                        // if N_st > 0
                        if(N[t * S + s] > MIN_COUNT){
                            if((N[t * S + s] - n[t * S + s]) > MIN_COUNT){
                                loglik = loglik + lgamma((1 - p_st) / nu_D[t] + N[t * S + s] - n[t * S + s]) - lgamma((1 - p_st) / nu_D[t]);
                            }
                            
                            if(n[t * S + s] > MIN_COUNT){
                                loglik = loglik + lgamma(p_st / nu_D[t] + n[t * S + s]) - lgamma(p_st / nu_D[t]);
                            }
                        }
                        
                    }
                
                    logpost = loglik / Tmp + (double) k * log(1 - zeta[c]);
                
                    prob_z_sc[k] = logpost;
                    if(k == 0){
                        max_log_prob_z_sc = logpost;
                    }
                    else{
                        if(max_log_prob_z_sc < logpost) max_log_prob_z_sc = logpost;
                    }
                } // end for(k). Now we have p(z_sc = k | ...) for all k, stored in prob_z_sc
                
                sum_prob_z_sc = 0;
                
                for(k = 0; k <= l_sc; k++){
                    prob_z_sc[k] = prob_z_sc[k] - max_log_prob_z_sc;
                    prob_z_sc[k] = exp(prob_z_sc[k]);
                    sum_prob_z_sc = sum_prob_z_sc + prob_z_sc[k];
                }
                
                for(k = l_sc + 1; k <= K_max; k++){
                    prob_z_sc[k] = 0;
                }
                
                GetRNGstate();
                randomUnif = runif(0.0, 1.0);
                PutRNGstate();
  
                // z_sc ~ Discrete(0:K, prob_z_sc). Sample z_sc.
                k = 0;
                prob_z_sc[k] = prob_z_sc[k] / sum_prob_z_sc;
                while(randomUnif > prob_z_sc[k]){
  	                k++;
  	                prob_z_sc[k] = prob_z_sc[k] / sum_prob_z_sc + prob_z_sc[k - 1];
                }
                
                Z[c * S + s] = k;
            }// end if(l_sc == 0)
        }// end for(c)
    }// end for(s)
    
    free(prob_z_sc);
    return;
}



void update_W_DNA(int *L, int *Z, double *W, double *phi,
    double *gamma_D, double *nu_D,
    double *l_D0, double *z_D0,
    double *N, double *n,
    int S, int T, int C,
    double a_w, double b_w, double d, double d0,
    double Tmp){
  
    // Update W
    // W: T * (C + 1) matrix. 
    // W_tC means w_t0 in the paper, proportion of the background subclone in sample t.
    // W_t0 means w_t1 in the paper, proportion of normal subclone
    
    
    // only calculate proposed w for the t-th row w_t = (w_t1, ..., w_tC, w_t0)
    double *w_t_pro;
    w_t_pro = (double *)calloc((C + 1), sizeof(double));
    
    double w_tc_cur, w_tc_pro;
    
    int s, t, c, c_iter;
    
    double l_D0_num = *l_D0;
    double z_D0_num = *z_D0;
    
    double loglik_cur, loglik_pro, logpost_cur, logpost_pro;
    double A_st_cur, A_st_pro, A_tilde_st_cur, A_tilde_st_pro;
    double p_st_cur, p_st_pro;
    double norm_prop;
    double sd_prop;
    double u;
    
    double w_t_sum;
    
    for(t = 0; t < T; t++){
        for(c = 0; c <= C; c++){
            // update w_tc
            
            w_tc_cur = W[c * T + t];
            
            if(c < C){
                sd_prop = 0.08;
            }
            else{
                sd_prop = 0.01;
            }
            
            GetRNGstate();
            norm_prop = rnorm(0.0, sd_prop);
            PutRNGstate();
            
            w_tc_pro = w_tc_cur + norm_prop;
            
            w_t_sum = 0;
            // need to change w_{t, -c} as well
            for(c_iter = 0; c_iter <= C; c_iter ++){
                if(c_iter != c){
                    w_t_pro[c_iter] = W[c_iter * T + t] * (1 - w_tc_pro) / (1 - w_tc_cur);
                }
                else{
                    w_t_pro[c_iter] = w_tc_pro;
                }
                w_t_sum = w_t_sum + w_t_pro[c_iter];
            }
            
            for(c_iter = 0; c_iter <= C; c_iter++){
                w_t_pro[c_iter] = w_t_pro[c_iter] / w_t_sum;
            }
            
            w_tc_pro = w_t_pro[c];
            
            // otherwise reject
            if((w_tc_pro < 1) && (w_tc_pro > 0) && (w_t_pro[C] < 0.03)){
                loglik_cur = 0.0;
                loglik_pro = 0.0;
            
                for(s = 0; s < S; s++){
                
                    A_st_cur = 0;
                    A_st_pro = 0;
                    A_tilde_st_cur = 0;
                    A_tilde_st_pro = 0;
                
                    for(c_iter = 0; c_iter < C; c_iter++){
                        A_st_cur = A_st_cur + W[c_iter * T + t] * L[c_iter * S + s];
                        A_st_pro = A_st_pro + w_t_pro[c_iter] * L[c_iter * S + s];
                        A_tilde_st_cur = A_tilde_st_cur + W[c_iter * T + t] * Z[c_iter * S + s];
                        A_tilde_st_pro = A_tilde_st_pro + w_t_pro[c_iter] * Z[c_iter * S + s];
                    }
                
                    A_st_cur = A_st_cur + W[C * T + t] * l_D0_num;
                    A_st_pro = A_st_pro + w_t_pro[C] * l_D0_num;
                    A_tilde_st_cur = A_tilde_st_cur + W[C * T + t] * z_D0_num;
                    A_tilde_st_pro = A_tilde_st_pro + w_t_pro[C] * z_D0_num;
                
                    p_st_cur = A_tilde_st_cur / A_st_cur;
                    p_st_pro = A_tilde_st_pro / A_st_pro;
                    
                    // if N_st = 0 and M_st = 0, done
                    loglik_cur = loglik_cur - (1 / gamma_D[t] + N[t * S + s]) * log(1 + gamma_D[t] * phi[t] * A_st_cur / 2);
                    
                    loglik_pro = loglik_pro - (1 / gamma_D[t] + N[t * S + s]) * log(1 + gamma_D[t] * phi[t] * A_st_pro / 2);
                    
                    // if N_st > 0
                    if(N[t * S + s] > MIN_COUNT){
                        loglik_cur = loglik_cur + N[t * S + s] * log(gamma_D[t] * phi[t] * A_st_cur / 2);
                        
                        loglik_pro = loglik_pro + N[t * S + s] * log(gamma_D[t] * phi[t] * A_st_pro / 2);
                        // if n_st = 0, done
                            
                        if((N[t * S + s] - n[t * S + s]) > MIN_COUNT){
                            loglik_cur = loglik_cur + lgamma((1 - p_st_cur) / nu_D[t] + N[t * S + s] - n[t * S + s]) - lgamma((1 - p_st_cur) / nu_D[t]);
                            
                            loglik_pro = loglik_pro + lgamma((1 - p_st_pro) / nu_D[t] + N[t * S + s] - n[t * S + s]) - lgamma((1 - p_st_pro) / nu_D[t]);
                        }
                            
                        if(n[t * S + s] > MIN_COUNT){
                            loglik_cur = loglik_cur + lgamma(p_st_cur / nu_D[t] + n[t * S + s]) - lgamma(p_st_cur / nu_D[t]);
                            
                            loglik_pro = loglik_pro + lgamma(p_st_pro / nu_D[t] + n[t * S + s]) - lgamma(p_st_pro / nu_D[t]);
                        }
                    }
                    // end calculating loglik
                
                }// end for(s)
            
                logpost_cur = loglik_cur / Tmp;
                logpost_pro = loglik_pro / Tmp;
            
                logpost_cur = logpost_cur + (a_w - 1) * log(W[0 * T + t]) + (b_w - 1) * log(1 - W[0 * T + t]);
                logpost_pro = logpost_pro + (a_w - 1) * log(w_t_pro[0]) + (b_w - 1) * log(1 - w_t_pro[0]);
            
                for(c_iter = 1; c_iter < C; c_iter++){
                    logpost_cur = logpost_cur + (d - 1) * log(W[c_iter * T + t] / (1 - W[0 * T + t]));
                    logpost_pro = logpost_pro + (d - 1) * log(w_t_pro[c_iter] / (1 - w_t_pro[0]));
                }
            
                logpost_cur = logpost_cur + (d0 - 1) * log(W[C * T + t] / (1 - W[0 * T + t]));
                logpost_pro = logpost_pro + (d0 - 1) * log(w_t_pro[C] / (1 - w_t_pro[0]));
                
                
                GetRNGstate();
                u = runif(0.0, 1.0);
                PutRNGstate();
  
                if(log(u) < logpost_pro - logpost_cur){
                    for(c_iter = 0; c_iter <= C; c_iter++){
                        W[c_iter * T + t] = w_t_pro[c_iter];
                    }
                }
            }// end if(w_tc > 0)
        }// end for(c)
    }// end for(t)
  
  free(w_t_pro);
  return;
}




// For DNA related quantities only. calculate A, A_tilde, p. All (S * T)
void calc_A_p(double *A, double *A_tilde, double *p,
    int *L, int *Z, double *W, double *l_D0, double *z_D0,
    int S, int T, int C){
    
    int s, t, c_iter;
    
    double l_D0_num = *l_D0;
    double z_D0_num = *z_D0;
    
    for(s = 0; s < S; s++){
        for(t = 0; t < T; t++){
            
            A[t * S + s] = 0.0;
            A_tilde[t * S + s] = 0.0;
            
            for(c_iter = 0; c_iter < C; c_iter++){
                A[t * S + s] = A[t * S + s] + W[c_iter * T + t] * L[c_iter * S + s];
                A_tilde[t * S + s] = A_tilde[t * S + s] + W[c_iter * T + t] * Z[c_iter * S + s];
            }
            
            A[t * S + s] = A[t * S + s] + W[C * T + t] * l_D0_num;
            A_tilde[t * S + s] = A_tilde[t * S + s] + W[C * T + t] * z_D0_num;
            
            p[t * S + s] = A_tilde[t * S + s] / A[t * S + s];
        }
    }// end for(s)
    // end calculating A, B, p, q
    
    return;
}




void update_phi_gamma_nu_l0_z0_DNA(int *L, int *Z, double *W, double *phi,
    double *gamma_D, double *nu_D,
    double *l_D0, double *z_D0,
    double *N, double *n,
    int S, int T, int C, int K_min, int K_max,
    double *a_phi, double *b_phi,
    double a_gamma_D, double b_gamma_D,
    double a_nu_D, double b_nu_D,
    double Tmp){
    
    // update phi_t, gamma_D_t, nu_D_t, l_D0, z_D0 in this one function
    // can avoid repeated calculation of A_st, B_st, p_st, q_st
    
    int s, t;
    
    double phi_t_cur, phi_t_pro;
    double gamma_D_t_cur, gamma_D_t_pro;
    double nu_D_t_cur, nu_D_t_pro;
    
    double l_D0_cur, l_D0_pro;
    double z_D0_cur, z_D0_pro;
    
    double loglik_cur, loglik_pro, logpost_cur, logpost_pro, logJacobian; 
    double norm_prop, randomUnif;
    
    // standard deviation for the normal proporsals
    double sd_prop_phi = 0.3;
    double sd_prop_gamma_D = 0.1;
    double sd_prop_nu_D = 0.1;
    double sd_prop_l_D0 = 0.1;
    double sd_prop_z_D0 = 0.1;
    
    double A_st_pro, A_tilde_st_pro, p_st_pro;
    
    double *A;
    A = (double *)calloc(S*T, sizeof(double));
    double *A_tilde;
    A_tilde = (double *)calloc(S*T, sizeof(double));
    double *p;
    p = (double *)calloc(S*T, sizeof(double));
    
    
    // calculate A, p. 
    // These will not change for updating phi, gamma_D and nu_D
    calc_A_p(A, A_tilde, p, L, Z, W, l_D0, z_D0, S, T, C);
    
    for(t = 0; t < T; t++){
        
        /////////////////////////////////////////////////////////
        // update phi_t
        /////////////////////////////////////////////////////////
        
        phi_t_cur = phi[t];
        GetRNGstate();
        norm_prop = rnorm(0.0, sd_prop_phi);
        PutRNGstate();
        phi_t_pro = phi_t_cur * exp(norm_prop);
        
        loglik_cur = 0.0;
        loglik_pro = 0.0;
        
        for(s = 0; s < S; s++){    
            loglik_cur = loglik_cur - (1 / gamma_D[t] + N[t * S + s]) * log(1 + gamma_D[t] * phi_t_cur * A[t * S + s] / 2);
            loglik_pro = loglik_pro - (1 / gamma_D[t] + N[t * S + s]) * log(1 + gamma_D[t] * phi_t_pro * A[t * S + s] / 2);
            
            if(N[t * S + s] > MIN_COUNT){
                loglik_cur = loglik_cur + lgamma(N[t * S + s] + 1 / gamma_D[t]) - lgamma(1 / gamma_D[t]) + N[t * S + s] * log(gamma_D[t] * phi_t_cur * A[t * S + s] / 2); 
                
                loglik_pro = loglik_pro + lgamma(N[t * S + s] + 1 / gamma_D[t]) - lgamma(1 / gamma_D[t]) + N[t * S + s] * log(gamma_D[t] * phi_t_pro * A[t * S + s] / 2);
            }
            
        }// end for(s)
        
        logpost_cur = loglik_cur / Tmp + (a_phi[t] - 1) * log(phi_t_cur) - b_phi[t] * phi_t_cur;
        logpost_pro = loglik_pro / Tmp + (a_phi[t] - 1) * log(phi_t_pro) - b_phi[t] * phi_t_pro;
        
        logJacobian = log(phi_t_pro) - log(phi_t_cur);
        
        GetRNGstate();
        randomUnif = runif(0.0, 1.0);
        PutRNGstate();
        
        if(log(randomUnif) < (logpost_pro - logpost_cur + logJacobian)){
            phi[t] = phi_t_pro;
            loglik_cur = loglik_pro;
        }
        
        /////////////////////////////////////////////////////////
        // update gamma_D_t
        /////////////////////////////////////////////////////////
        gamma_D_t_cur = gamma_D[t];
        
        GetRNGstate();
        norm_prop = rnorm(0.0, sd_prop_gamma_D);
        PutRNGstate();
        
        gamma_D_t_pro = gamma_D_t_cur * exp(norm_prop);
        
        loglik_pro = 0.0;
        
        for(s = 0; s < S; s++){    
            loglik_pro = loglik_pro - (1 / gamma_D_t_pro + N[t * S + s]) * log(1 + gamma_D_t_pro * phi[t] * A[t * S + s] / 2);
            
            if(N[t * S + s] > MIN_COUNT){
                loglik_pro = loglik_pro + lgamma(N[t * S + s] + 1 / gamma_D_t_pro) - lgamma(1 / gamma_D_t_pro) + N[t * S + s] * log(gamma_D_t_pro * phi[t] * A[t * S + s] / 2);
            }
            
        }// end for(s)
        
        logpost_cur = loglik_cur / Tmp + (a_gamma_D - 1) * log(gamma_D_t_cur) - b_gamma_D * gamma_D_t_cur;
        logpost_pro = loglik_pro / Tmp + (a_gamma_D - 1) * log(gamma_D_t_pro) - b_gamma_D * gamma_D_t_pro;
        
        logJacobian = log(gamma_D_t_pro) - log(gamma_D_t_cur);
        
        GetRNGstate();
        randomUnif = runif(0.0, 1.0);
        PutRNGstate();
        
        if(log(randomUnif) < (logpost_pro - logpost_cur + logJacobian)){
            gamma_D[t] = gamma_D_t_pro;
        }
        
        /////////////////////////////////////////////////////////
        // update nu_D_t
        /////////////////////////////////////////////////////////
        nu_D_t_cur = nu_D[t];
        GetRNGstate();
        norm_prop = rnorm(0.0, sd_prop_nu_D);
        PutRNGstate();
        nu_D_t_pro = nu_D_t_cur * exp(norm_prop);
        
        loglik_cur = 0.0;
        loglik_pro = 0.0;
        
        for(s = 0; s < S; s++){
        
            // if N_st > 0
            if(N[t * S + s] > MIN_COUNT){
                loglik_cur = loglik_cur + lgamma(1 / nu_D_t_cur) - lgamma(1 / nu_D_t_cur + N[t * S + s]);
                        
                loglik_pro = loglik_pro + lgamma(1 / nu_D_t_pro) - lgamma(1 / nu_D_t_pro + N[t * S + s]);
                // if n_st = 0, done
                            
                if((N[t * S + s] - n[t * S + s]) > MIN_COUNT){
                    loglik_cur = loglik_cur + lgamma((1 - p[t * S + s]) / nu_D_t_cur + N[t * S + s] - n[t * S + s]) - lgamma((1 - p[t * S + s]) / nu_D_t_cur);
                            
                    loglik_pro = loglik_pro + lgamma((1 - p[t * S + s]) / nu_D_t_pro + N[t * S + s] - n[t * S + s]) - lgamma((1 - p[t * S + s]) / nu_D_t_pro);
                }
                            
                if(n[t * S + s] > MIN_COUNT){
                    loglik_cur = loglik_cur + lgamma(p[t * S + s] / nu_D_t_cur + n[t * S + s]) - lgamma(p[t * S + s] / nu_D_t_cur);
                            
                    loglik_pro = loglik_pro + lgamma(p[t * S + s] / nu_D_t_pro + n[t * S + s]) - lgamma(p[t * S + s] / nu_D_t_pro);
                }
            }
            
        }// end for(s)
        
        logpost_cur = loglik_cur / Tmp + (a_nu_D - 1) * log(nu_D_t_cur) - b_nu_D * nu_D_t_cur;
        logpost_pro = loglik_pro / Tmp + (a_nu_D - 1) * log(nu_D_t_pro) - b_nu_D * nu_D_t_pro;
        
        logJacobian = log(nu_D_t_pro) - log(nu_D_t_cur);
        
        // check if this produces the same random uniform every time
        GetRNGstate();
        randomUnif = runif(0.0, 1.0);
        PutRNGstate();
        
        if(log(randomUnif) < (logpost_pro - logpost_cur + logJacobian)){
            nu_D[t] = nu_D_t_pro;
        }
        
    }// end for(t)
    
    /////////////////////////////////////////////////////////
    // update l_D0, z_D0
    /////////////////////////////////////////////////////////
    // update all together so that can save some computational time
    // even if acceptance probability is low, not important parameters anyway..
    
    l_D0_cur = *l_D0;
    GetRNGstate();
    norm_prop = rnorm(0.0, sd_prop_l_D0);
    PutRNGstate();
    l_D0_pro = l_D0_cur * exp(norm_prop);
    
    z_D0_cur = *z_D0;
    GetRNGstate();
    norm_prop = rnorm(0.0, sd_prop_z_D0);
    PutRNGstate();
    z_D0_pro = z_D0_cur * exp(norm_prop);
    
    // otherwise reject directly
    if((l_D0_pro > (double) K_min) && (l_D0_pro < (double) K_max) && (z_D0_pro < l_D0_pro)){
        loglik_cur = 0.0;
        loglik_pro = 0.0;
        
        for(s = 0; s < S; s++){
            for(t = 0; t < T; t++){
            
                A_st_pro = A[t * S + s] - W[C * T + t] * l_D0_cur + W[C * T + t] * l_D0_pro;
                A_tilde_st_pro = A_tilde[t * S + s] - W[C * T + t] * z_D0_cur + W[C * T + t] * z_D0_pro;
                p_st_pro = A_tilde_st_pro / A_st_pro;
                
                
                loglik_cur = loglik_cur - (1 / gamma_D[t] + N[t * S + s]) * log(1 + gamma_D[t] * phi[t] * A[t * S + s] / 2);
                
                loglik_pro = loglik_pro - (1 / gamma_D[t] + N[t * S + s]) * log(1 + gamma_D[t] * phi[t] * A_st_pro / 2);
                
                // if N_st > 0
                if(N[t * S + s] > MIN_COUNT){
                    loglik_cur = loglik_cur + N[t * S + s] * log(gamma_D[t] * phi[t] * A[t * S + s] / 2);
                    
                    loglik_pro = loglik_pro + N[t * S + s] * log(gamma_D[t] * phi[t] * A_st_pro / 2);
                    // if n_st = 0, done
                            
                    if((N[t * S + s] - n[t * S + s]) > MIN_COUNT){
                        loglik_cur = loglik_cur + lgamma((1 - p[t * S + s]) / nu_D[t] + N[t * S + s] - n[t * S + s]) - lgamma((1 - p[t * S + s]) / nu_D[t]);
                        
                        loglik_pro = loglik_pro + lgamma((1 - p_st_pro) / nu_D[t] + N[t * S + s] - n[t * S + s]) - lgamma((1 - p_st_pro) / nu_D[t]);
                    }
                            
                    if(n[t * S + s] > MIN_COUNT){
                        loglik_cur = loglik_cur + lgamma(p[t * S + s] / nu_D[t] + n[t * S + s]) - lgamma(p[t * S + s] / nu_D[t]);
                        
                        loglik_pro = loglik_pro + lgamma(p_st_pro / nu_D[t] + n[t * S + s]) - lgamma(p_st_pro / nu_D[t]);
                    }
                }
                // end calculating loglik
                
            }
        }// end for(s)
        
        logpost_cur = loglik_cur / Tmp - log(l_D0_cur);
        logpost_pro = loglik_pro / Tmp - log(l_D0_pro);
        
        
        logJacobian = log(l_D0_pro) + log(z_D0_pro) - log(l_D0_cur) - log(z_D0_cur);
        
        GetRNGstate();
        randomUnif = runif(0.0, 1.0);
        PutRNGstate();
        
        if(log(randomUnif) < (logpost_pro - logpost_cur + logJacobian)){
            *l_D0 = l_D0_pro;
            *z_D0 = z_D0_pro;
        }
        
    } // end if(l_D0 ...)
    // finish updating l_D0, z_D0, l_R0, z_R0

    free(A);
    free(A_tilde);
    free(p);
    return;
}






void calc_loglik_DNA(double *loglik, int *L, int *Z, double *W,
    double *phi,
    double *gamma_D, double *nu_D,
    double *l_D0, double *z_D0,
    double *N, double *n,
    int *S, int *T, int *C){
    
    *loglik = 0.0;
    
    int s, t, c;
    
    double A_st, A_tilde_st, p_st;
    
    for(t = 0; t < (*T); t++){
        for(s = 0; s < (*S); s++){
            A_st = 0;
            A_tilde_st = 0;
                        
            for(c = 0; c < (*C); c++){
                A_st = A_st + W[c * (*T) + t] * L[c * (*S) + s];
                A_tilde_st = A_tilde_st + W[c * (*T) + t] * Z[c * (*S) + s];
            }
            
            A_st = A_st + W[(*C) * (*T) + t] * (*l_D0);
            A_tilde_st = A_tilde_st + W[(*C) * (*T) + t] * (*z_D0);
            p_st = A_tilde_st / A_st;
            
            // calculating loglik[s, t]
            (*loglik) = (*loglik) - (1 / gamma_D[t] + N[t * (*S) + s]) * log(1 + gamma_D[t] * phi[t] * A_st / 2);
            
            // if N_st > 0
            if(N[t * (*S) + s] > MIN_COUNT){
                (*loglik) = (*loglik) + lgamma(N[t * (*S) + s] + 1 / gamma_D[t]) - lgamma(1 / gamma_D[t]) + N[t * (*S) + s] * log(gamma_D[t] * phi[t] * A_st / 2) + lgamma(1 / nu_D[t]) - lgamma(1 / nu_D[t] + N[t * (*S) + s]);
                
                // if n_st = 0, done
                        
                if((N[t * (*S) + s] - n[t * (*S) + s]) > MIN_COUNT){
                    (*loglik) = (*loglik) + lgamma((1 - p_st) / nu_D[t] + N[t * (*S) + s] - n[t * (*S) + s]) - lgamma((1 - p_st) / nu_D[t]);
                }
                
                if(n[t * (*S) + s] > MIN_COUNT){
                    (*loglik) = (*loglik) + lgamma(p_st / nu_D[t] + n[t * (*S) + s]) - lgamma(p_st / nu_D[t]);
                }
            }
            // end calculating loglik[s, t]
            
        }
    }
    
    if(isnan((*loglik))){
        (*loglik) = -INFINITY;
    }
    
    return;
}




void calc_logpost_DNA(double *loglik, double *logpost, 
    int *L, int *Z, double *W,
    double *pai, double *zeta, double *phi,
    double *gamma_D, double *nu_D,
    double *l_D0, double *z_D0,
    double *N, double *n,
    int *S, int *T, int *C, int *K_min, int *K_max,
    double *a_w, double *b_w, double *d, double *d0,
    double *a_pai, double *b_pai, double *a_zeta, double *b_zeta,
    double *a_phi, double *b_phi,
    double *a_gamma_D, double *b_gamma_D,
    double *a_nu_D, double *b_nu_D){
    
    
    int s, t, c, k1;
    
    double logprior = 0.0;
    
    double temp1, temp2;
    
    // likelihood
    calc_loglik_DNA(loglik, L, Z, W, phi, gamma_D, nu_D, l_D0, z_D0, N, n, S, T, C);
    
    // prior for L, Z, pai and zeta
    for(c = 1; c < (*C); c++){
        for(s = 0; s < (*S); s++){
            logprior = logprior + abs(L[c * (*S) + s] - 2) * log(1 - pai[c]);
            logprior = logprior + Z[c * (*S) + s] * log(1 - zeta[c]);
            temp2 = 0;
            for(k1 = 0; k1 <= L[c * (*S) + s]; k1++){
                temp2 = temp2 + zeta[c] * pow(1 - zeta[c], k1);
            }
            logprior = logprior - log(temp2);
        }
        
        temp1 = 0;
        for(k1 = (*K_min); k1 <= (*K_max); k1++){
            temp1 = temp1 + pai[c] * pow(1 - pai[c], abs(k1 - 2));
        }
        logprior = logprior + (*S) * log(pai[c]) - (*S) * log(temp1);
        logprior = logprior + (*S) * log(zeta[c]);
        
        logprior = logprior + ((*a_pai) - 1) * log(pai[c]) + ((*b_pai) - 1) * log(1 - pai[c]);
        logprior = logprior + ((*a_zeta) - 1) * log(zeta[c]) + ((*b_zeta) - 1) * log(1 - zeta[c]);
    }
    
    
    // prior for W, phi, gamma_D, nu_D
    for(t = 0; t < (*T); t++){
        logprior = logprior + ((*a_w) - 1) * log(W[0 * (*T) + t]) + ((*b_w) - 1 - ((*C) - 1) * ((*d) - 1) - ((*d0) - 1)) * log(1 - W[0 * (*T) + t]);
        for(c = 1; c < (*C); c++){
            logprior = logprior + ((*d) - 1) * log(W[c * (*T) + t]);
        }
        logprior = logprior + ((*d0) - 1) * log(W[(*C) * (*T) + t]);
        
        logprior = logprior + (a_phi[t] - 1) * log(phi[t]) - b_phi[t] * phi[t];
        logprior = logprior + ((*a_gamma_D) - 1) * log(gamma_D[t]) - (*b_gamma_D) * gamma_D[t];
        logprior = logprior + ((*a_nu_D) - 1) * log(nu_D[t]) - (*b_nu_D) * nu_D[t];
        
    }
    
    // prior for l_D0, z_D0
    logprior = logprior - log(*l_D0);
    
    
    (*logpost) = (*loglik) + logprior;
    
    if(isnan((*logpost))){
        (*logpost) = -INFINITY;
    }
    
    return;
}







void PT_swap_DNA(int **L1, int **Z1, double **W1, double **pai1, double **zeta1,
    double **phi1, double **gamma_D1, double **nu_D1,
    double *l_D01, double *z_D01, double *loglik1, double *lpost1, 
    int **L2, int **Z2, double **W2, double **pai2, double **zeta2,
    double **phi2, double **gamma_D2, double **nu_D2,
    double *l_D02, double *z_D02, double *loglik2, double *lpost2){
    
    
    int *temp_int;
    double *temp_double;
    double temp_double2;
    
    temp_int = *L1;
    *L1 = *L2;
    *L2 = temp_int;
    
    temp_int = *Z1;
    *Z1 = *Z2;
    *Z2 = temp_int;
    
    temp_double = *W1;
    *W1 = *W2;
    *W2 = temp_double;
    
    temp_double = *pai1;
    *pai1 = *pai2;
    *pai2 = temp_double;
    
    temp_double = *zeta1;
    *zeta1 = *zeta2;
    *zeta2 = temp_double;
    
    temp_double = *phi1;
    *phi1 = *phi2;
    *phi2 = temp_double;
    
    temp_double = *gamma_D1;
    *gamma_D1 = *gamma_D2;
    *gamma_D2 = temp_double;
    
    temp_double = *nu_D1;
    *nu_D1 = *nu_D2;
    *nu_D2 = temp_double;
    
    temp_double2 = *l_D01;
    *l_D01 = *l_D02;
    *l_D02 = temp_double2;
    
    temp_double2 = *z_D01;
    *z_D01 = *z_D02;
    *z_D02 = temp_double2;
    
    temp_double2 = *loglik1;
    *loglik1 = *loglik2;
    *loglik2 = temp_double2;
    
    temp_double2 = *lpost1;
    *lpost1 = *lpost2;
    *lpost2 = temp_double2;
    
    return;
}



void DClone_MCMC(int *L_PTR, int *Z_PTR, double *W_PTR,
    double *pai_PTR, double *zeta_PTR, double *phi_PTR,
    double *gamma_D_PTR, double *nu_D_PTR,
    double *l_D0_PTR, double *z_D0_PTR,
    double *loglik_PTR, double *logpost_PTR,
    double *N, double *n,
    int *S, int *T, int *C, int *K_min, int *K_max,
    double *a_w, double *b_w, double *d, double *d0,
    double *a_pai, double *b_pai, double *a_zeta, double *b_zeta,
    double *a_phi, double *b_phi,
    double *a_gamma_D, double *b_gamma_D,
    double *a_nu_D, double *b_nu_D,
    int *niter, double *Delta, int *nThread){
    
    // inferring L, Z, W using DNA data only
    // DClone!!
    
    
    int s, t, c;
    // i: iteration, from 1 to niter
    int i, rank, rank_partner;
	double u, lalpha;
    
    int S_num = *S;
    int T_num = *T;
    int C_num = *C;
    int K_min_num = *K_min;
    int K_max_num = *K_max;
    
    double a_w_num = *a_w;
    double b_w_num = *b_w;
    double d_num = *d;
    double d0_num = *d0;
    double a_pai_num = *a_pai;
    double b_pai_num = *b_pai;
    double a_zeta_num = *a_zeta;
    double b_zeta_num = *b_zeta;
    double a_gamma_D_num = *a_gamma_D;
    double b_gamma_D_num = *b_gamma_D;
    double a_nu_D_num = *a_nu_D;
    double b_nu_D_num = *b_nu_D;
    
    int niter_num = *niter;
    int nThread_num = *nThread;
    
    ////////////////////////////////////////////////////////////////////
    // Create PT arrays for saving parameter values for each temprature
    // Memory allocation
    ////////////////////////////////////////////////////////////////////
    
    int *L_PT[nThread_num]; // S * C
    int *Z_PT[nThread_num]; // S * C
    double *W_PT[nThread_num]; // T * (C+1)
    
    double *pai_PT[nThread_num]; // C
    double *zeta_PT[nThread_num]; // C
    
    double *phi_PT[nThread_num]; // T
    double *gamma_D_PT[nThread_num]; // T
    double *nu_D_PT[nThread_num]; // T
    
    double l_D0_PT[nThread_num]; // 1
    double z_D0_PT[nThread_num]; // 1
    
    double loglik_PT[nThread_num];
    double logpost_PT[nThread_num];	// 1
    
	
	for(rank = 0; rank < nThread_num; rank++){
        
        L_PT[rank] = (int*)(calloc(S_num*C_num, sizeof(int)));
        Z_PT[rank] = (int*)(calloc(S_num*C_num, sizeof(int)));
		W_PT[rank] = (double*)(calloc(T_num*(C_num+1), sizeof(double)));
		
        pai_PT[rank] = (double*)(calloc(C_num, sizeof(double)));
        zeta_PT[rank] = (double*)(calloc(C_num, sizeof(double)));
        
        phi_PT[rank] = (double*)(calloc(T_num, sizeof(double)));
        gamma_D_PT[rank] = (double*)(calloc(T_num, sizeof(double)));
        nu_D_PT[rank] = (double*)(calloc(T_num, sizeof(double)));
        
	}  
	
    ////////////////////////////////////////////////////////////////////
    // Setting initial values of the PT arrays = input values
    ////////////////////////////////////////////////////////////////////
    
    for(rank = 0; rank < nThread_num; rank++){
		
        // Setting L and Z
        for(s = 0; s < S_num; s++){
            for(c = 0; c < C_num; c++){
                L_PT[rank][c * S_num + s] = L_PTR[rank * S_num * C_num + c * S_num + s];
                Z_PT[rank][c * S_num + s] = Z_PTR[rank * S_num * C_num + c * S_num + s];
            }
        }
        
        // Setting W
        for(t = 0; t < T_num; t++){
            for(c = 0; c <= C_num; c++){
                W_PT[rank][c * T_num + t] = W_PTR[rank * T_num * (C_num + 1) + c * T_num + t];
            }
        }
        
        
        // Setting pai and zeta
        for(c = 0; c < C_num; c++){
            pai_PT[rank][c] = pai_PTR[rank * C_num + c];
            zeta_PT[rank][c] = zeta_PTR[rank * C_num + c];
        }
        
        // Setting phi, psi, gamma_D, nu_D, gamma_R, nu_R
        for(t = 0; t < T_num; t++){
            phi_PT[rank][t] = phi_PTR[rank * T_num + t];
            gamma_D_PT[rank][t] = gamma_D_PTR[rank * T_num + t];
            nu_D_PT[rank][t] = nu_D_PTR[rank * T_num + t];
        }
        
        // Setting l_D0, z_D0, l_R0, z_R0
        l_D0_PT[rank] = l_D0_PTR[rank];
        z_D0_PT[rank] = z_D0_PTR[rank];
        
        loglik_PT[rank] = loglik_PTR[rank];
	    logpost_PT[rank] = logpost_PTR[rank];
	
	} // end for(rank)
    
    ///////////////////////////////////////////////////////////////////////////////
	// MCMC update.. Using for loop for parallel tempering
	///////////////////////////////////////////////////////////////////////////////
	
	for(i = 0; i < niter_num; i++){
	    for(rank = 0; rank < nThread_num; rank++){
	        
	        update_L_DNA(L_PT[rank], Z_PT[rank], W_PT[rank], pai_PT[rank], zeta_PT[rank], phi_PT[rank], gamma_D_PT[rank], nu_D_PT[rank], &l_D0_PT[rank], &z_D0_PT[rank], N, n, S_num, T_num, C_num, K_min_num, K_max_num, Delta[rank]);
    
            update_pai(L_PT[rank], pai_PT[rank], S_num, C_num, K_min_num, K_max_num, a_pai_num, b_pai_num);
    
            update_Z_DNA(L_PT[rank], Z_PT[rank], W_PT[rank], zeta_PT[rank], nu_D_PT[rank], &l_D0_PT[rank], &z_D0_PT[rank], N, n, S_num, T_num, C_num, K_min_num, K_max_num, Delta[rank]);
            
            update_zeta(L_PT[rank], Z_PT[rank], zeta_PT[rank], S_num, C_num, a_zeta_num, b_zeta_num);
    
            update_W_DNA(L_PT[rank], Z_PT[rank], W_PT[rank], phi_PT[rank], gamma_D_PT[rank], nu_D_PT[rank], &l_D0_PT[rank], &z_D0_PT[rank], N, n, S_num, T_num, C_num, a_w_num, b_w_num, d_num, d0_num, Delta[rank]);
            
            update_phi_gamma_nu_l0_z0_DNA(L_PT[rank], Z_PT[rank], W_PT[rank], phi_PT[rank], gamma_D_PT[rank], nu_D_PT[rank], &l_D0_PT[rank], &z_D0_PT[rank], N, n, S_num, T_num, C_num, K_min_num, K_max_num, a_phi, b_phi, a_gamma_D_num, b_gamma_D_num, a_nu_D_num, b_nu_D_num, Delta[rank]);
            
            calc_logpost_DNA(&loglik_PT[rank], &logpost_PT[rank], L_PT[rank], Z_PT[rank], W_PT[rank], pai_PT[rank], zeta_PT[rank], phi_PT[rank], gamma_D_PT[rank], nu_D_PT[rank], &l_D0_PT[rank], &z_D0_PT[rank], N, n, S, T, C, K_min, K_max, a_w, b_w, d, d0, a_pai, b_pai, a_zeta, b_zeta, a_phi, b_phi, a_gamma_D, b_gamma_D, a_nu_D, b_nu_D);
            
            
	    }
	    
	    for(rank = 0; rank < (nThread_num - 1); rank++){
	    
	        GetRNGstate();
            u = runif(0.0, 1.0);
            PutRNGstate();
	        
	        rank_partner = rank + 1;
	        lalpha = (1/Delta[rank] - 1/Delta[rank_partner]) * (loglik_PT[rank_partner] - loglik_PT[rank]);
	        
	        if(log(u) < lalpha){
                PT_swap_DNA(&L_PT[rank], &Z_PT[rank], &W_PT[rank], &pai_PT[rank], &zeta_PT[rank], &phi_PT[rank], &gamma_D_PT[rank], &nu_D_PT[rank], &l_D0_PT[rank], &z_D0_PT[rank], &loglik_PT[rank], &logpost_PT[rank], &L_PT[rank_partner], &Z_PT[rank_partner], &W_PT[rank_partner], &pai_PT[rank_partner], &zeta_PT[rank_partner], &phi_PT[rank_partner], &gamma_D_PT[rank_partner], &nu_D_PT[rank_partner], &l_D0_PT[rank_partner], &z_D0_PT[rank_partner], &loglik_PT[rank_partner], &logpost_PT[rank_partner]);
            }
	    }
	
	} // end for(i)
	
    ///////////////////////////////////////////////////////////////////////////////
	// End MCMC update
	///////////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////
    // Setting updated parameter values = current PT array values
    ////////////////////////////////////////////////////////////////////
    
    
    for(rank = 0; rank < nThread_num; rank++){
		
        // Setting L and Z
        for(s = 0; s < S_num; s++){
            for(c = 0; c < C_num; c++){
                L_PTR[rank * S_num * C_num + c * S_num + s] = L_PT[rank][c * S_num + s];
                Z_PTR[rank * S_num * C_num + c * S_num + s] = Z_PT[rank][c * S_num + s];
            }
        }
        
        // Setting W
        for(t = 0; t < T_num; t++){
            for(c = 0; c <= C_num; c++){
                W_PTR[rank * T_num * (C_num + 1) + c * T_num + t] = W_PT[rank][c * T_num + t];
            }
        }
        
        
        // Setting pai and zeta
        for(c = 0; c < C_num; c++){
            pai_PTR[rank * C_num + c] = pai_PT[rank][c];
            zeta_PTR[rank * C_num + c] = zeta_PT[rank][c];
        }
        
        // Setting phi, psi, gamma_D, nu_D, gamma_R, nu_R
        for(t = 0; t < T_num; t++){
            phi_PTR[rank * T_num + t] = phi_PT[rank][t];
            gamma_D_PTR[rank * T_num + t] = gamma_D_PT[rank][t];
            nu_D_PTR[rank * T_num + t] = nu_D_PT[rank][t];
        }
        
        // Setting l_D0, z_D0, l_R0, z_R0
        l_D0_PTR[rank] = l_D0_PT[rank];
        z_D0_PTR[rank] = z_D0_PT[rank];
        
        loglik_PTR[rank] = loglik_PT[rank];
	    logpost_PTR[rank] = logpost_PT[rank];
	
	} // end for(rank)
	
	
	////////////////////////////////////////////////////////////////////
    // Memory free
    ////////////////////////////////////////////////////////////////////
	
    
	for(rank = 0; rank < nThread_num; rank++){
	    free(L_PT[rank]);
		free(Z_PT[rank]);
		free(W_PT[rank]);
		
		free(pai_PT[rank]);
		free(zeta_PT[rank]);
		
		free(phi_PT[rank]);
		free(gamma_D_PT[rank]);
		free(nu_D_PT[rank]);
	}
    
    return;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// End DNA only..
/////////////////////////////////////////////////////////////////////////////////////////////////////////





/////////////////////////////////////////////////////////////////////////////////////////////////////////
// RNA-based subclone reconstruction. 
// Only update W and Lambda
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// Update W using RNA data only
void update_W_RNA(int *L, int *Z, double *W, double *Lambda, double *psi,
    double *gamma_R, double *nu_R,
    double *l_R0, double *z_R0,
    double *M, double *m, 
    int S, int T, int C, int G,
    double a_w, double b_w, double d, double d0,
    int *g_fun, double Tmp){
  
    // Update W
    // W: T * (C + 1) matrix. 
    // W_tC means w_t0 in the paper, proportion of the background subclone in sample t.
    // W_t0 means w_t1 in the paper, proportion of normal subclone
    
    
    // only calculate proposed w for the t-th row w_t = (w_t1, ..., w_tC, w_t0)
    double *w_t_pro;
    w_t_pro = (double *)calloc((C + 1), sizeof(double));
    
    double w_tc_cur, w_tc_pro;
    
    int s, t, c, g, c_iter;
    
    double l_R0_num = *l_R0;
    double z_R0_num = *z_R0;
    
    double loglik_cur, loglik_pro, logpost_cur, logpost_pro;
    double B_st_cur, B_st_pro, B_tilde_st_cur, B_tilde_st_pro;
    double q_st_cur, q_st_pro;
    double lambda_gt0_cur, lambda_gt0_pro;
    double norm_prop;
    double sd_prop;
    double u;
    
    double w_t_sum;
    
    for(t = 0; t < T; t++){
        for(c = 0; c <= C; c++){
            // update w_tc
            
            w_tc_cur = W[c * T + t];
            
            if(c < C){
                sd_prop = 0.08;
            }
            else{
                sd_prop = 0.01;
            }
            
            
            GetRNGstate();
            norm_prop = rnorm(0.0, sd_prop);
            PutRNGstate();
            
            w_tc_pro = w_tc_cur + norm_prop;
            
            w_t_sum = 0;
            // need to change w_{t, -c} as well
            for(c_iter = 0; c_iter <= C; c_iter++){
                if(c_iter != c){
                    w_t_pro[c_iter] = W[c_iter * T + t] * (1 - w_tc_pro) / (1 - w_tc_cur);
                }
                else{
                    w_t_pro[c_iter] = w_tc_pro;
                }
                w_t_sum = w_t_sum + w_t_pro[c_iter];
            }
            
            for(c_iter = 0; c_iter <= C; c_iter++){
                w_t_pro[c_iter] = w_t_pro[c_iter] / w_t_sum;
            }
            
            w_tc_pro = w_t_pro[c];
            
            // otherwise reject
            if((w_tc_pro < 1) && (w_tc_pro > 0) && (w_t_pro[C] < 0.03)){
                loglik_cur = 0.0;
                loglik_pro = 0.0;
            
                for(s = 0; s < S; s++){
                    g = g_fun[s];
                
                    B_st_cur = 0;
                    B_st_pro = 0;
                    B_tilde_st_cur = 0;
                    B_tilde_st_pro = 0;
                
                    lambda_gt0_cur = 0;
                    lambda_gt0_pro = 0;
                
                    for(c_iter = 0; c_iter < C; c_iter++){

                        B_st_cur = B_st_cur + W[c_iter * T + t] * L[c_iter * S + s] * Lambda[c_iter * G + g];
                        B_st_pro = B_st_pro + w_t_pro[c_iter] * L[c_iter * S + s] * Lambda[c_iter * G + g];
                        B_tilde_st_cur = B_tilde_st_cur + W[c_iter * T + t] * Z[c_iter * S + s]* Lambda[c_iter * G + g];
                        B_tilde_st_pro = B_tilde_st_pro + w_t_pro[c_iter] * Z[c_iter * S + s]* Lambda[c_iter * G + g];
                    
                        lambda_gt0_cur = lambda_gt0_cur + W[c_iter * T + t] * Lambda[c_iter * G + g];
                        lambda_gt0_pro = lambda_gt0_pro + w_t_pro[c_iter] * Lambda[c_iter * G + g];
                    }
                
                    B_st_cur = B_st_cur + W[C * T + t] * l_R0_num * lambda_gt0_cur;
                    B_st_pro = B_st_pro + w_t_pro[C] * l_R0_num * lambda_gt0_pro;
                    B_tilde_st_cur = B_tilde_st_cur + W[C * T + t] * z_R0_num * lambda_gt0_cur;
                    B_tilde_st_pro = B_tilde_st_pro + w_t_pro[C] * z_R0_num * lambda_gt0_pro;
                    
                    q_st_cur = B_tilde_st_cur / B_st_cur;
                    q_st_pro = B_tilde_st_pro / B_st_pro;
                    
                    // if N_st = 0 and M_st = 0, done
                    loglik_cur = loglik_cur - (1 / gamma_R[t] + M[t * S + s]) * log(1 + gamma_R[t] * psi[t] * B_st_cur / 2);
                    
                    loglik_pro = loglik_pro - (1 / gamma_R[t] + M[t * S + s]) * log(1 + gamma_R[t] * psi[t] * B_st_pro / 2);
                        
                    // if M_st > 0
                    if(M[t * S + s] > MIN_COUNT){
                        loglik_cur = loglik_cur + M[t * S + s] * log(gamma_R[t] * psi[t] * B_st_cur / 2);
                        
                        loglik_pro = loglik_pro + M[t * S + s] * log(gamma_R[t] * psi[t] * B_st_pro / 2);
                        // if m_st = 0, done
                        
                        if((M[t * S + s] - m[t * S + s]) > MIN_COUNT){
                            loglik_cur = loglik_cur + lgamma((1 - q_st_cur) / nu_R[t] + M[t * S + s] - m[t * S + s]) - lgamma((1 - q_st_cur) / nu_R[t]);
                            
                            loglik_pro = loglik_pro + lgamma((1 - q_st_pro) / nu_R[t] + M[t * S + s] - m[t * S + s]) - lgamma((1 - q_st_pro) / nu_R[t]);
                        }
                            
                        if(m[t * S + s] > MIN_COUNT){
                            loglik_cur = loglik_cur + lgamma(q_st_cur / nu_R[t] + m[t * S + s]) - lgamma(q_st_cur / nu_R[t]);
                            
                            loglik_pro = loglik_pro + lgamma(q_st_pro / nu_R[t] + m[t * S + s]) - lgamma(q_st_pro / nu_R[t]);
                        }
                    }
                    // end calculating loglik
                    
                }// end for(s)
            
                logpost_cur = loglik_cur / Tmp;
                logpost_pro = loglik_pro / Tmp;
            
                logpost_cur = logpost_cur + (a_w - 1) * log(W[0 * T + t]) + (b_w - 1) * log(1 - W[0 * T + t]);
                logpost_pro = logpost_pro + (a_w - 1) * log(w_t_pro[0]) + (b_w - 1) * log(1 - w_t_pro[0]);
            
                for(c_iter = 1; c_iter < C; c_iter++){
                    logpost_cur = logpost_cur + (d - 1) * log(W[c_iter * T + t] / (1 - W[0 * T + t]));
                    logpost_pro = logpost_pro + (d - 1) * log(w_t_pro[c_iter] / (1 - w_t_pro[0]));
                }
            
                logpost_cur = logpost_cur + (d0 - 1) * log(W[C * T + t] / (1 - W[0 * T + t]));
                logpost_pro = logpost_pro + (d0 - 1) * log(w_t_pro[C] / (1 - w_t_pro[0]));
            
                // logpost_cur = logpost_cur / Tmp;
                // logpost_pro = logpost_pro / Tmp;
                
                GetRNGstate();
                u = runif(0.0, 1.0);
                PutRNGstate();
  
                if(log(u) < logpost_pro - logpost_cur){
                    for(c_iter = 0; c_iter <= C; c_iter++){
                        W[c_iter * T + t] = w_t_pro[c_iter];
                    }
                }
            }// end if(w_tc > 0)
        }// end for(c)
    }// end for(t)
  
  free(w_t_pro);
  return;
}


// For RNA related quantities only. calculate B, B_tilde, q and lambda0. All (S * T)
void calc_B_q(double *B, double *B_tilde, double *q, double *lambda0,
    int *L, int *Z, double *W, double *Lambda, double *l_R0, double *z_R0,
    int S, int T, int C, int G, int *g_fun){
    
    int s, t, g, c_iter;
    
    double l_R0_num = *l_R0;
    double z_R0_num = *z_R0;
    
    for(s = 0; s < S; s++){
        for(t = 0; t < T; t++){
            g = g_fun[s];
            
            B[t * S + s] = 0.0;
            B_tilde[t * S + s] = 0.0;
            
            lambda0[t * S + s] = 0.0;
            
            for(c_iter = 0; c_iter < C; c_iter++){
                B[t * S + s] = B[t * S + s] + W[c_iter * T + t] * L[c_iter * S + s] * Lambda[c_iter * G + g];
                B_tilde[t * S + s] = B_tilde[t * S + s] + W[c_iter * T + t] * Z[c_iter * S + s] * Lambda[c_iter * G + g];
                
                lambda0[t * S + s] = lambda0[t * S + s] + W[c_iter * T + t] * Lambda[c_iter * G + g];
            }
            
            B[t * S + s] = B[t * S + s] + W[C * T + t] * l_R0_num * lambda0[t * S + s];
            B_tilde[t * S + s] = B_tilde[t * S + s] + W[C * T + t] * z_R0_num * lambda0[t * S + s];
            q[t * S + s] = B_tilde[t * S + s] / B[t * S + s];
        }
    }// end for(s)
    // end calculating A, B, p, q
    
    return;
}




void update_psi_gamma_nu_l0_z0_RNA(int *L, int *Z, double *W, double *Lambda, double *psi,
    double *gamma_R, double *nu_R,
    double *l_R0, double *z_R0,
    double *M, double *m,  
    int S, int T, int C, int K_min, int K_max, int G,
    double *a_psi, double *b_psi,
    double a_gamma_R, double b_gamma_R,
    double a_nu_R, double b_nu_R,
    int *g_fun, double Tmp){
    
    // update phi_t, psi_t, gamma_D_t, nu_D_t, gamma_R_t, nu_R_t in this one function
    // can avoid repeated calculation of A_st, B_st, p_st, q_st
    
    int s, t;
    
    double psi_t_cur, psi_t_pro;
    double gamma_R_t_cur, gamma_R_t_pro;
    double nu_R_t_cur, nu_R_t_pro;
    
    double l_R0_cur, l_R0_pro;
    double z_R0_cur, z_R0_pro;
    
    double loglik_cur, loglik_pro, logpost_cur, logpost_pro, logJacobian; 
    double norm_prop, randomUnif;
    
    // standard deviation for the normal proporsals
    double sd_prop_psi = 0.3;
    double sd_prop_gamma_R = 0.1;
    double sd_prop_nu_R = 0.1;
    double sd_prop_l_R0 = 0.1;
    double sd_prop_z_R0 = 0.1;
    
    double B_st_pro, B_tilde_st_pro, q_st_pro;
    
    double *B;
    B = (double *)calloc(S*T, sizeof(double));
    double *B_tilde;
    B_tilde = (double *)calloc(S*T, sizeof(double));
    double *q;
    q = (double *)calloc(S*T, sizeof(double));
    
    // should be G*T, but ok for S*T with some same rows
    double *lambda0;
    lambda0 = (double *)calloc(S*T, sizeof(double));
    
    // calculate B, q. 
    // These will not change for updating psi, gamma_R and nu_R
    calc_B_q(B, B_tilde, q, lambda0, L, Z, W, Lambda, l_R0, z_R0, S, T, C, G, g_fun);
    
    for(t = 0; t < T; t++){
    
        /////////////////////////////////////////////////////////
        // update psi_t
        /////////////////////////////////////////////////////////
        psi_t_cur = psi[t];
        GetRNGstate();
        norm_prop = rnorm(0.0, sd_prop_psi);
        PutRNGstate();
        psi_t_pro = psi_t_cur * exp(norm_prop);
        
        loglik_cur = 0.0;
        loglik_pro = 0.0;
        
        for(s = 0; s < S; s++){
            loglik_cur = loglik_cur - (1 / gamma_R[t] + M[t * S + s]) * log(1 + gamma_R[t] * psi_t_cur * B[t * S + s] / 2);
            loglik_pro = loglik_pro - (1 / gamma_R[t] + M[t * S + s]) * log(1 + gamma_R[t] * psi_t_pro * B[t * S + s] / 2);
            
            if(M[t * S + s] > MIN_COUNT){
                loglik_cur = loglik_cur + lgamma(M[t * S + s] + 1 / gamma_R[t]) - lgamma(1 / gamma_R[t]) + M[t * S + s] * log(gamma_R[t] * psi_t_cur * B[t * S + s] / 2); 
                
                loglik_pro = loglik_pro + lgamma(M[t * S + s] + 1 / gamma_R[t]) - lgamma(1 / gamma_R[t]) + M[t * S + s] * log(gamma_R[t] * psi_t_pro * B[t * S + s] / 2);
            }
            
        }
        
        logpost_cur = loglik_cur / Tmp + (a_psi[t] - 1) * log(psi_t_cur) - b_psi[t] * psi_t_cur;
        logpost_pro = loglik_pro / Tmp + (a_psi[t] - 1) * log(psi_t_pro) - b_psi[t] * psi_t_pro;
        
        logJacobian = log(psi_t_pro) - log(psi_t_cur);
        
        GetRNGstate();
        randomUnif = runif(0.0, 1.0);
        PutRNGstate();
        
        if(log(randomUnif) < (logpost_pro - logpost_cur + logJacobian)){
            psi[t] = psi_t_pro;
            loglik_cur = loglik_pro;
        }
        
        /////////////////////////////////////////////////////////
        // update gamma_R_t
        /////////////////////////////////////////////////////////
        gamma_R_t_cur = gamma_R[t];
        GetRNGstate();
        norm_prop = rnorm(0.0, sd_prop_gamma_R);
        PutRNGstate();
        gamma_R_t_pro = gamma_R_t_cur * exp(norm_prop);
        
        loglik_pro = 0.0;
        
        for(s = 0; s < S; s++){    
            loglik_pro = loglik_pro - (1 / gamma_R_t_pro + M[t * S + s]) * log(1 + gamma_R_t_pro * psi[t] * B[t * S + s] / 2);
            
            if(M[t * S + s] > MIN_COUNT){
                loglik_pro = loglik_pro + lgamma(M[t * S + s] + 1 / gamma_R_t_pro) - lgamma(1 / gamma_R_t_pro) + M[t * S + s] * log(gamma_R_t_pro * psi[t] * B[t * S + s] / 2);
            }
            
        }// end for(s)
        
        logpost_cur = loglik_cur / Tmp + (a_gamma_R - 1) * log(gamma_R_t_cur) - b_gamma_R * gamma_R_t_cur;
        logpost_pro = loglik_pro / Tmp + (a_gamma_R - 1) * log(gamma_R_t_pro) - b_gamma_R * gamma_R_t_pro;
        
        
        logJacobian = log(gamma_R_t_pro) - log(gamma_R_t_cur);
        
        // check if this produces the same random uniform every time
        GetRNGstate();
        randomUnif = runif(0.0, 1.0);
        PutRNGstate();
        
        if(log(randomUnif) < (logpost_pro - logpost_cur + logJacobian)){
            gamma_R[t] = gamma_R_t_pro;
        }
        
        
        /////////////////////////////////////////////////////////
        // update nu_R_t
        /////////////////////////////////////////////////////////
        nu_R_t_cur = nu_R[t];
        GetRNGstate();
        norm_prop = rnorm(0.0, sd_prop_nu_R);
        PutRNGstate();
        nu_R_t_pro = nu_R_t_cur * exp(norm_prop);
        
        loglik_cur = 0.0;
        loglik_pro = 0.0;
        
        for(s = 0; s < S; s++){
            
            // if M_st > 0
            if(M[t * S + s] > MIN_COUNT){
                loglik_cur = loglik_cur + lgamma(1 / nu_R_t_cur) - lgamma(1 / nu_R_t_cur + M[t * S + s]);
                        
                loglik_pro = loglik_pro + lgamma(1 / nu_R_t_pro) - lgamma(1 / nu_R_t_pro + M[t * S + s]);
                // if m_st = 0, done
                            
                if((M[t * S + s] - m[t * S + s]) > MIN_COUNT){
                    loglik_cur = loglik_cur + lgamma((1 - q[t * S + s]) / nu_R_t_cur + M[t * S + s] - m[t * S + s]) - lgamma((1 - q[t * S + s]) / nu_R_t_cur);
                            
                    loglik_pro = loglik_pro + lgamma((1 - q[t * S + s]) / nu_R_t_pro + M[t * S + s] - m[t * S + s]) - lgamma((1 - q[t * S + s]) / nu_R_t_pro);
                }
                            
                if(m[t * S + s] > MIN_COUNT){
                    loglik_cur = loglik_cur + lgamma(q[t * S + s] / nu_R_t_cur + m[t * S + s]) - lgamma(q[t * S + s] / nu_R_t_cur);
                            
                    loglik_pro = loglik_pro + lgamma(q[t * S + s] / nu_R_t_pro + m[t * S + s]) - lgamma(q[t * S + s] / nu_R_t_pro);
                }
            }
            
        }
        
        logpost_cur = loglik_cur / Tmp + (a_nu_R - 1) * log(nu_R_t_cur) - b_nu_R * nu_R_t_cur;
        logpost_pro = loglik_pro / Tmp + (a_nu_R - 1) * log(nu_R_t_pro) - b_nu_R * nu_R_t_pro;
        
        logJacobian = log(nu_R_t_pro) - log(nu_R_t_cur);
        
        // check if this produces the same random uniform every time
        GetRNGstate();
        randomUnif = runif(0.0, 1.0);
        PutRNGstate();
        
        if(log(randomUnif) < (logpost_pro - logpost_cur + logJacobian)){
            nu_R[t] = nu_R_t_pro;
        }
        
    }// end for(t)
    
    /////////////////////////////////////////////////////////
    // update l_R0, z_R0
    /////////////////////////////////////////////////////////
    
    l_R0_cur = *l_R0;
    GetRNGstate();
    norm_prop = rnorm(0.0, sd_prop_l_R0);
    PutRNGstate();
    l_R0_pro = l_R0_cur * exp(norm_prop);
    
    z_R0_cur = *z_R0;
    GetRNGstate();
    norm_prop = rnorm(0.0, sd_prop_z_R0);
    PutRNGstate();
    z_R0_pro = z_R0_cur * exp(norm_prop);
    
    // otherwise reject directly
    if((l_R0_pro > (double) K_min) && (l_R0_pro < (double) K_max) && (z_R0_pro < l_R0_pro)){
        loglik_cur = 0.0;
        loglik_pro = 0.0;
        
        for(s = 0; s < S; s++){
            for(t = 0; t < T; t++){
            
                B_st_pro = B[t * S + s] - W[C * T + t] * l_R0_cur * lambda0[t * S + s] + W[C * T + t] * l_R0_pro * lambda0[t * S + s];
                B_tilde_st_pro = B_tilde[t * S + s] - W[C * T + t] * z_R0_cur * lambda0[t * S + s] + W[C * T + t] * z_R0_pro * lambda0[t * S + s]; 
                q_st_pro = B_tilde_st_pro / B_st_pro;
                
                loglik_cur = loglik_cur - (1 / gamma_R[t] + M[t * S + s]) * log(1 + gamma_R[t] * psi[t] * B[t * S + s] / 2);
                
                loglik_pro = loglik_pro - (1 / gamma_R[t] + M[t * S + s]) * log(1 + gamma_R[t] * psi[t] * B_st_pro / 2);
                
                // if M_st > 0
                if(M[t * S + s] > MIN_COUNT){
                    loglik_cur = loglik_cur + M[t * S + s] * log(gamma_R[t] * psi[t] * B[t * S + s] / 2);
                    
                    loglik_pro = loglik_pro + M[t * S + s] * log(gamma_R[t] * psi[t] * B_st_pro / 2);
                    // if m_st = 0, done
                    
                    if((M[t * S + s] - m[t * S + s]) > MIN_COUNT){
                        loglik_cur = loglik_cur + lgamma((1 - q[t * S + s]) / nu_R[t] + M[t * S + s] - m[t * S + s]) - lgamma((1 - q[t * S + s]) / nu_R[t]);
                        
                        loglik_pro = loglik_pro + lgamma((1 - q_st_pro) / nu_R[t] + M[t * S + s] - m[t * S + s]) - lgamma((1 - q_st_pro) / nu_R[t]);
                    }
                        
                    if(m[t * S + s] > MIN_COUNT){
                        loglik_cur = loglik_cur + lgamma(q[t * S + s] / nu_R[t] + m[t * S + s])  - lgamma(q[t * S + s] / nu_R[t]);
                        
                        loglik_pro = loglik_pro + lgamma(q_st_pro / nu_R[t] + m[t * S + s]) - lgamma(q_st_pro / nu_R[t]);
                    }
                }
                // end calculating loglik
            }
        }
        
        logpost_cur = loglik_cur / Tmp - log(l_R0_cur);
        logpost_pro = loglik_pro / Tmp - log(l_R0_pro);
        
        logJacobian = log(l_R0_pro) + log(z_R0_pro) - log(l_R0_cur) - log(z_R0_cur);
        
        GetRNGstate();
        randomUnif = runif(0.0, 1.0);
        PutRNGstate();
        
        if(log(randomUnif) < (logpost_pro - logpost_cur + logJacobian)){
            *l_R0 = l_R0_pro;
            *z_R0 = z_R0_pro;
        }
        
    } // end if(l_D0 ...)
    // finish updating l_D0, z_D0, l_R0, z_R0

    free(B);
    free(B_tilde);
    free(q);
    free(lambda0);
    return;
}






void calc_loglik_RNA(double *loglik, 
    int *L, int *Z, double *W, double *Lambda, 
    double *psi,
    double *gamma_R, double *nu_R,
    double *l_R0, double *z_R0,
    double *M, double *m, 
    int *S, int *T, int *C, int *G,
    int *g_fun){
    
    *loglik = 0.0;
    
    int s, t, c, g;
    
    double B_st, B_tilde_st, q_st, lambda_gt0;
    
    for(t = 0; t < (*T); t++){
        for(s = 0; s < (*S); s++){
            g = g_fun[s];
            
            B_st = 0;
            B_tilde_st = 0;
            lambda_gt0 = 0;
                        
            for(c = 0; c < (*C); c++){
                B_st = B_st + W[c * (*T) + t] * L[c * (*S) + s] * Lambda[c * (*G) + g];
                B_tilde_st = B_tilde_st + W[c * (*T) + t] * Z[c * (*S) + s]* Lambda[c * (*G) + g];
                lambda_gt0 = lambda_gt0 + W[c * (*T) + t] * Lambda[c * (*G) + g];
            }
            
            B_st = B_st + W[(*C) * (*T) + t] * (*l_R0) * lambda_gt0;
            B_tilde_st = B_tilde_st + W[(*C) * (*T) + t] * (*z_R0) * lambda_gt0;
            q_st = B_tilde_st / B_st;
            
            // calculating loglik[s, t]
            (*loglik) = (*loglik) - (1 / gamma_R[t] + M[t * (*S) + s]) * log(1 + gamma_R[t] * psi[t] * B_st / 2);
                        
            // if M_st > 0
            if(M[t * (*S) + s] > MIN_COUNT){
                (*loglik) = (*loglik) + lgamma(M[t * (*S) + s] + 1 / gamma_R[t]) - lgamma(1 / gamma_R[t]) + M[t * (*S) + s] * log(gamma_R[t] * psi[t] * B_st / 2) + lgamma(1 / nu_R[t]) - lgamma(1 / nu_R[t] + M[t * (*S) + s]);
                
                // if m_st = 0, done
                
                if((M[t * (*S) + s] - m[t * (*S) + s]) > MIN_COUNT){
                    (*loglik) = (*loglik) + lgamma((1 - q_st) / nu_R[t] + M[t * (*S) + s] - m[t * (*S) + s]) - lgamma((1 - q_st) / nu_R[t]);
                }
                
                if(m[t * (*S) + s] > MIN_COUNT){
                    (*loglik) = (*loglik) + lgamma(q_st / nu_R[t] + m[t * (*S) + s]) - lgamma(q_st / nu_R[t]);
                }
            }
            // end calculating loglik[s, t]
            
        }
    }
    
    if(isnan((*loglik))){
        (*loglik) = -INFINITY;
    }
    
    return;
}




void calc_logpost_RNA(double *loglik, double *logpost, 
    int *L, int *Z, double *W, double *Lambda, 
    double *psi,
    double *gamma_R, double *nu_R,
    double *l_R0, double *z_R0,
    double *M, double *m, 
    int *S, int *T, int *C, int *G,
    double *a_w, double *b_w, double *d, double *d0,
    double *a_lambda, double *b_lambda,
    double *a_psi, double *b_psi,
    double *a_gamma_R, double *b_gamma_R,
    double *a_nu_R, double *b_nu_R,
    int *g_fun){
    
    
    int t, c, g;
    
    double logprior = 0.0;
    
    
    // likelihood
    calc_loglik_RNA(loglik, L, Z, W, Lambda, psi, gamma_R, nu_R, l_R0, z_R0, M, m, S, T, C, G, g_fun);
    
    // prior for W, psi, gamma_R, nu_R
    for(t = 0; t < (*T); t++){
        
        logprior = logprior + ((*a_w) - 1) * log(W[0 * (*T) + t]) + ((*b_w) - 1 - ((*C) - 1) * ((*d) - 1) - ((*d0) - 1)) * log(1 - W[0 * (*T) + t]);
        for(c = 1; c < (*C); c++){
            logprior = logprior + ((*d) - 1) * log(W[c * (*T) + t]);
        }
        logprior = logprior + ((*d0) - 1) * log(W[(*C) * (*T) + t]);

        logprior = logprior + (a_psi[t] - 1) * log(psi[t]) - b_psi[t] * psi[t];
        logprior = logprior + ((*a_gamma_R) - 1) * log(gamma_R[t]) - (*b_gamma_R) * gamma_R[t];
        logprior = logprior + ((*a_nu_R) - 1) * log(nu_R[t]) - (*b_nu_R) * nu_R[t];
        
    }
    
    // prior for l_R0, z_R0 
    logprior = logprior - log(*l_R0);
    
    // prior for Lambda
    for(g = 0; g < (*G); g++){
        for(c = 0; c < (*C); c++){
            logprior = logprior + ((*a_lambda) - 1) * log(Lambda[c * (*G) + g]) - (*b_lambda) * Lambda[c * (*G) + g];
        }
    }
    
    (*logpost) = (*loglik) + logprior;
    
    if(isnan((*logpost))){
        (*logpost) = -INFINITY;
    }
    
    return;
}





// swap the state of two Markov chains
void PT_swap_RNA(double **W1, double **Lambda1,
    double **psi1, double **gamma_R1, double **nu_R1,
    double *l_R01, double *z_R01, double *loglik1, double *lpost1, 
    double **W2, double **Lambda2,
    double **psi2, double **gamma_R2, double **nu_R2,
    double *l_R02, double *z_R02, double *loglik2, double *lpost2){
    
    double *temp_double;
    double temp_double2;
    
    temp_double = *W1;
    *W1 = *W2;
    *W2 = temp_double;

    temp_double = *Lambda1;
    *Lambda1 = *Lambda2;
    *Lambda2 = temp_double;
    
    temp_double = *psi1;
    *psi1 = *psi2;
    *psi2 = temp_double;
    
    temp_double = *gamma_R1;
    *gamma_R1 = *gamma_R2;
    *gamma_R2 = temp_double;
    
    temp_double = *nu_R1;
    *nu_R1 = *nu_R2;
    *nu_R2 = temp_double;
    
    temp_double2 = *l_R01;
    *l_R01 = *l_R02;
    *l_R02 = temp_double2;
    
    temp_double2 = *z_R01;
    *z_R01 = *z_R02;
    *z_R02 = temp_double2;
    
    temp_double2 = *loglik1;
    *loglik1 = *loglik2;
    *loglik2 = temp_double2;
    
    temp_double2 = *lpost1;
    *lpost1 = *lpost2;
    *lpost2 = temp_double2;
    
    return;
}



void RClone_MCMC(int *L, int *Z, 
    double *W_PTR, double *Lambda_PTR, 
    double *psi_PTR,
    double *gamma_R_PTR, double *nu_R_PTR,
    double *l_R0_PTR, double *z_R0_PTR,
    double *loglik_PTR, double *logpost_PTR,
    double *M, double *m,
    int *S, int *T, int *C, int *G, int *K_min, int *K_max,
    double *a_w, double *b_w, double *d, double *d0,
    double *a_lambda, double *b_lambda,
    double *a_psi, double *b_psi,
    double *a_gamma_R, double *b_gamma_R,
    double *a_nu_R, double *b_nu_R,
    int *niter, int *g_fun, double *Delta, int *nThread){
    
    // main MCMC function
    // L, Z, W_PTR, Lambda_PTR, ...., z_R0_PTR: parameters
    // N, n, M, m: data
    // S, T, C, G, K: constants, # of SNVs, # of samples, # of subclones (random), # of genes, maximum copy #
    // a_w, ..., b_nu_R: hyperparameters (fixed)
    // niter: number of MCMC iterations 
    // g_fun: a vector storing the corresponding gene g of SNV s, g_fun(SNV s) = gene g
    // Delta: nThread-dimensional vector of tempratures.
    // nThread: number of temperatures to run parallel tempering (for the PT version, it's equal to number of threads)
    
    
    int t, c, g;
    // i: iteration, from 1 to niter
    int i, rank, rank_partner;
	double u, lalpha;
    
    int S_num = *S;
    int T_num = *T;
    int C_num = *C;
    int G_num = *G;
    int K_min_num = *K_min;
    int K_max_num = *K_max;
    
    double a_w_num = *a_w;
    double b_w_num = *b_w;
    double d_num = *d;
    double d0_num = *d0;
    double a_lambda_num = *a_lambda;
    double b_lambda_num = *b_lambda;
    double a_gamma_R_num = *a_gamma_R;
    double b_gamma_R_num = *b_gamma_R;
    double a_nu_R_num = *a_nu_R;
    double b_nu_R_num = *b_nu_R;
    
    int niter_num = *niter;
    int nThread_num = *nThread;
    ////////////////////////////////////////////////////////////////////
    // Create PT arrays for saving parameter values for each temprature
    // Memory allocation
    ////////////////////////////////////////////////////////////////////
    
    double *W_PT[nThread_num]; // T * (C+1)
    double *Lambda_PT[nThread_num]; // G * C
    
    double *psi_PT[nThread_num]; // T
    double *gamma_R_PT[nThread_num]; // T
    double *nu_R_PT[nThread_num]; // T
    
    double l_R0_PT[nThread_num]; // 1
    double z_R0_PT[nThread_num]; // 1
    
    double loglik_PT[nThread_num];
    double logpost_PT[nThread_num];	// 1
    
	
	for(rank = 0; rank < nThread_num; rank++){

        W_PT[rank] = (double*)(calloc(T_num*(C_num+1), sizeof(double)));
        Lambda_PT[rank] = (double*)(calloc(G_num*C_num, sizeof(double)));
        
        psi_PT[rank] = (double*)(calloc(T_num, sizeof(double)));
        gamma_R_PT[rank] = (double*)(calloc(T_num, sizeof(double)));
        nu_R_PT[rank] = (double*)(calloc(T_num, sizeof(double)));
        
	}  
	
    ////////////////////////////////////////////////////////////////////
    // Setting initial values of the PT arrays = input values
    ////////////////////////////////////////////////////////////////////
    
    for(rank = 0; rank < nThread_num; rank++){
        
        // Setting W
        for(t = 0; t < T_num; t++){
            for(c = 0; c <= C_num; c++){
                W_PT[rank][c * T_num + t] = W_PTR[rank * T_num * (C_num + 1) + c * T_num + t];
            }
        }

        // Setting Lambda
        for(g = 0; g < G_num; g++){
            for(c = 0; c < C_num; c++){
                Lambda_PT[rank][c * G_num + g] = Lambda_PTR[rank * G_num * C_num + c * G_num + g];
            }
        }
        
        
        // Setting phi, psi, gamma_D, nu_D, gamma_R, nu_R
        for(t = 0; t < T_num; t++){
            psi_PT[rank][t] = psi_PTR[rank * T_num + t];
            gamma_R_PT[rank][t] = gamma_R_PTR[rank * T_num + t];
            nu_R_PT[rank][t] = nu_R_PTR[rank * T_num + t];
        }
        
        // Setting l_D0, z_D0, l_R0, z_R0
        l_R0_PT[rank] = l_R0_PTR[rank];
        z_R0_PT[rank] = z_R0_PTR[rank];
        
        loglik_PT[rank] = loglik_PTR[rank];
	    logpost_PT[rank] = logpost_PTR[rank];
	
	} // end for(rank)
    
    ///////////////////////////////////////////////////////////////////////////////
	// MCMC update.. Using for loop for parallel tempering
	///////////////////////////////////////////////////////////////////////////////
	
	for(i = 0; i < niter_num; i++){
	    for(rank = 0; rank < nThread_num; rank++){
            
            update_W_RNA(L, Z, W_PT[rank], Lambda_PT[rank], psi_PT[rank], gamma_R_PT[rank], nu_R_PT[rank], &l_R0_PT[rank], &z_R0_PT[rank], M, m, S_num, T_num, C_num, G_num, a_w_num, b_w_num, d_num, d0_num, g_fun, Delta[rank]);

            update_Lambda(L, Z, W_PT[rank], Lambda_PT[rank], psi_PT[rank], gamma_R_PT[rank], nu_R_PT[rank], &l_R0_PT[rank], &z_R0_PT[rank], M, m, S_num, T_num, C_num, G_num, a_lambda_num, b_lambda_num, g_fun, Delta[rank]);
            
            update_psi_gamma_nu_l0_z0_RNA(L, Z, W_PT[rank], Lambda_PT[rank], psi_PT[rank], gamma_R_PT[rank], nu_R_PT[rank], &l_R0_PT[rank], &z_R0_PT[rank], M, m, S_num, T_num, C_num, K_min_num, K_max_num, G_num, a_psi, b_psi, a_gamma_R_num, b_gamma_R_num, a_nu_R_num, b_nu_R_num, g_fun, Delta[rank]);
            
            calc_logpost_RNA(&loglik_PT[rank], &logpost_PT[rank], 
                L, Z, W_PT[rank], Lambda_PT[rank], psi_PT[rank], 
                gamma_R_PT[rank], nu_R_PT[rank], &l_R0_PT[rank], &z_R0_PT[rank], 
                M, m, S, T, C, G, a_w, b_w, d, d0,
                a_lambda, b_lambda, a_psi, b_psi, 
                a_gamma_R, b_gamma_R, a_nu_R, b_nu_R, g_fun);
            
	    }
	    
	    for(rank = 0; rank < (nThread_num - 1); rank++){
	    
	        GetRNGstate();
            u = runif(0.0, 1.0);
            PutRNGstate();
	        
	        rank_partner = rank + 1;
	        lalpha = (1/Delta[rank] - 1/Delta[rank_partner]) * (loglik_PT[rank_partner] - loglik_PT[rank]);
	        
	        if(log(u) < lalpha){
                PT_swap_RNA(&W_PT[rank], &Lambda_PT[rank], 
                    &psi_PT[rank], &gamma_R_PT[rank], 
                    &nu_R_PT[rank], &l_R0_PT[rank], &z_R0_PT[rank], 
                    &loglik_PT[rank], &logpost_PT[rank], 
                    &W_PT[rank_partner], &Lambda_PT[rank_partner], 
                    &psi_PT[rank_partner], &gamma_R_PT[rank_partner], 
                    &nu_R_PT[rank_partner], &l_R0_PT[rank_partner], &z_R0_PT[rank_partner], 
                    &loglik_PT[rank_partner], &logpost_PT[rank_partner]);
            }
	    }
	
	} // end for(i)
	
    ///////////////////////////////////////////////////////////////////////////////
	// End MCMC update
	///////////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////
    // Setting updated parameter values = current PT array values
    ////////////////////////////////////////////////////////////////////
    
    
    for(rank = 0; rank < nThread_num; rank++){
        
        // Setting W
        for(t = 0; t < T_num; t++){
            for(c = 0; c <= C_num; c++){
                W_PTR[rank * T_num * (C_num + 1) + c * T_num + t] = W_PT[rank][c * T_num + t];
            }
        }

        // Setting Lambda
        for(g = 0; g < G_num; g++){
            for(c = 0; c < C_num; c++){
                Lambda_PTR[rank * G_num * C_num + c * G_num + g] = Lambda_PT[rank][c * G_num + g];
            }
        }
        
        
        // Setting psi, gamma_R, nu_R
        for(t = 0; t < T_num; t++){
            psi_PTR[rank * T_num + t] = psi_PT[rank][t];
            gamma_R_PTR[rank * T_num + t] = gamma_R_PT[rank][t];
            nu_R_PTR[rank * T_num + t] = nu_R_PT[rank][t];
        }
        
        // Setting l_R0, z_R0
        l_R0_PTR[rank] = l_R0_PT[rank];
        z_R0_PTR[rank] = z_R0_PT[rank];
        
        loglik_PTR[rank] = loglik_PT[rank];
	    logpost_PTR[rank] = logpost_PT[rank];
	
	} // end for(rank)
	
	
	////////////////////////////////////////////////////////////////////
    // Memory free
    ////////////////////////////////////////////////////////////////////
	
    
	for(rank = 0; rank < nThread_num; rank++){
        free(W_PT[rank]);
		free(Lambda_PT[rank]);
		
		free(psi_PT[rank]);
		free(gamma_R_PT[rank]);
		free(nu_R_PT[rank]);
	}
    
    return;
}



















/////////////////////////////////////////////////////////////////////////////////////////////////////////
// RNA-based subclone reconstruction with fixed W. 
// Only update Lambda
/////////////////////////////////////////////////////////////////////////////////////////////////////////

void calc_logpost_RNA_fixed_W(double *loglik, double *logpost, 
    int *L, int *Z, double *W, double *Lambda, 
    double *psi,
    double *gamma_R, double *nu_R,
    double *l_R0, double *z_R0,
    double *M, double *m, 
    int *S, int *T, int *C, int *G,
    double *a_lambda, double *b_lambda,
    double *a_psi, double *b_psi,
    double *a_gamma_R, double *b_gamma_R,
    double *a_nu_R, double *b_nu_R,
    int *g_fun){
    
    
    int t, c, g;
    
    double logprior = 0.0;
    
    
    // likelihood
    calc_loglik_RNA(loglik, L, Z, W, Lambda, psi, gamma_R, nu_R, l_R0, z_R0, M, m, S, T, C, G, g_fun);
    
    // prior for psi, gamma_R, nu_R
    for(t = 0; t < (*T); t++){
    
        logprior = logprior + (a_psi[t] - 1) * log(psi[t]) - b_psi[t] * psi[t];
        logprior = logprior + ((*a_gamma_R) - 1) * log(gamma_R[t]) - (*b_gamma_R) * gamma_R[t];
        logprior = logprior + ((*a_nu_R) - 1) * log(nu_R[t]) - (*b_nu_R) * nu_R[t];
        
    }
    
    // prior for l_R0, z_R0 
    logprior = logprior - log(*l_R0);
    
    // prior for Lambda
    for(g = 0; g < (*G); g++){
        for(c = 0; c < (*C); c++){
            logprior = logprior + ((*a_lambda) - 1) * log(Lambda[c * (*G) + g]) - (*b_lambda) * Lambda[c * (*G) + g];
        }
    }
    
    (*logpost) = (*loglik) + logprior;
    
    if(isnan((*logpost))){
        (*logpost) = -INFINITY;
    }
    
    return;
}





// swap the state of two Markov chains
void PT_swap_RNA_fixed_W(double **Lambda1,
    double **psi1, double **gamma_R1, double **nu_R1,
    double *l_R01, double *z_R01, double *loglik1, double *lpost1, 
    double **Lambda2,
    double **psi2, double **gamma_R2, double **nu_R2,
    double *l_R02, double *z_R02, double *loglik2, double *lpost2){
    
    double *temp_double;
    double temp_double2;
    
    temp_double = *Lambda1;
    *Lambda1 = *Lambda2;
    *Lambda2 = temp_double;
    
    temp_double = *psi1;
    *psi1 = *psi2;
    *psi2 = temp_double;
    
    temp_double = *gamma_R1;
    *gamma_R1 = *gamma_R2;
    *gamma_R2 = temp_double;
    
    temp_double = *nu_R1;
    *nu_R1 = *nu_R2;
    *nu_R2 = temp_double;
    
    temp_double2 = *l_R01;
    *l_R01 = *l_R02;
    *l_R02 = temp_double2;
    
    temp_double2 = *z_R01;
    *z_R01 = *z_R02;
    *z_R02 = temp_double2;
    
    temp_double2 = *loglik1;
    *loglik1 = *loglik2;
    *loglik2 = temp_double2;
    
    temp_double2 = *lpost1;
    *lpost1 = *lpost2;
    *lpost2 = temp_double2;
    
    return;
}




void RClone_MCMC_fixed_W(int *L, int *Z, double *W, double *Lambda_PTR, 
    double *psi_PTR,
    double *gamma_R_PTR, double *nu_R_PTR,
    double *l_R0_PTR, double *z_R0_PTR,
    double *loglik_PTR, double *logpost_PTR,
    double *M, double *m,
    int *S, int *T, int *C, int *G, int *K_min, int *K_max,
    double *a_lambda, double *b_lambda,
    double *a_psi, double *b_psi,
    double *a_gamma_R, double *b_gamma_R,
    double *a_nu_R, double *b_nu_R,
    int *niter, int *g_fun, double *Delta, int *nThread){
    
    // main MCMC function
    // L_PTR, Z_PTR, W_PTR, Lambda_PTR, ...., z_R0_PTR: parameters
    // N, n, M, m: data
    // S, T, C, G, K: constants, # of SNVs, # of samples, # of subclones (random), # of genes, maximum copy #
    // a_w, ..., b_nu_R: hyperparameters (fixed)
    // niter: number of MCMC iterations 
    // g_fun: a vector storing the corresponding gene g of SNV s, g_fun(SNV s) = gene g
    // Delta: nThread-dimensional vector of tempratures.
    // nThread: number of temperatures to run parallel tempering (for the PT version, it's equal to number of threads)
    
    
    int t, c, g;
    // i: iteration, from 1 to niter
    int i, rank, rank_partner;
	double u, lalpha;
    
    int S_num = *S;
    int T_num = *T;
    int C_num = *C;
    int G_num = *G;
    int K_min_num = *K_min;
    int K_max_num = *K_max;
    
    double a_lambda_num = *a_lambda;
    double b_lambda_num = *b_lambda;
    double a_gamma_R_num = *a_gamma_R;
    double b_gamma_R_num = *b_gamma_R;
    double a_nu_R_num = *a_nu_R;
    double b_nu_R_num = *b_nu_R;
    
    int niter_num = *niter;
    int nThread_num = *nThread;
    ////////////////////////////////////////////////////////////////////
    // Create PT arrays for saving parameter values for each temprature
    // Memory allocation
    ////////////////////////////////////////////////////////////////////
    
    double *Lambda_PT[nThread_num]; // G * C
    
    double *psi_PT[nThread_num]; // T
    double *gamma_R_PT[nThread_num]; // T
    double *nu_R_PT[nThread_num]; // T
    
    double l_R0_PT[nThread_num]; // 1
    double z_R0_PT[nThread_num]; // 1
    
    double loglik_PT[nThread_num];
    double logpost_PT[nThread_num];	// 1
    
	
	for(rank = 0; rank < nThread_num; rank++){
        Lambda_PT[rank] = (double*)(calloc(G_num*C_num, sizeof(double)));
        
        psi_PT[rank] = (double*)(calloc(T_num, sizeof(double)));
        gamma_R_PT[rank] = (double*)(calloc(T_num, sizeof(double)));
        nu_R_PT[rank] = (double*)(calloc(T_num, sizeof(double)));
        
	}  
	
    ////////////////////////////////////////////////////////////////////
    // Setting initial values of the PT arrays = input values
    ////////////////////////////////////////////////////////////////////
    
    for(rank = 0; rank < nThread_num; rank++){
    
        // Setting Lambda
        for(g = 0; g < G_num; g++){
            for(c = 0; c < C_num; c++){
                Lambda_PT[rank][c * G_num + g] = Lambda_PTR[rank * G_num * C_num + c * G_num + g];
            }
        }
        
        
        // Setting phi, psi, gamma_D, nu_D, gamma_R, nu_R
        for(t = 0; t < T_num; t++){
            psi_PT[rank][t] = psi_PTR[rank * T_num + t];
            gamma_R_PT[rank][t] = gamma_R_PTR[rank * T_num + t];
            nu_R_PT[rank][t] = nu_R_PTR[rank * T_num + t];
        }
        
        // Setting l_D0, z_D0, l_R0, z_R0
        l_R0_PT[rank] = l_R0_PTR[rank];
        z_R0_PT[rank] = z_R0_PTR[rank];
        
        loglik_PT[rank] = loglik_PTR[rank];
	    logpost_PT[rank] = logpost_PTR[rank];
	
	} // end for(rank)
    
    ///////////////////////////////////////////////////////////////////////////////
	// MCMC update.. Using for loop for parallel tempering
	///////////////////////////////////////////////////////////////////////////////
	
	for(i = 0; i < niter_num; i++){
	    for(rank = 0; rank < nThread_num; rank++){
            
            update_Lambda(L, Z, W, Lambda_PT[rank], psi_PT[rank], gamma_R_PT[rank], nu_R_PT[rank], &l_R0_PT[rank], &z_R0_PT[rank], M, m, S_num, T_num, C_num, G_num, a_lambda_num, b_lambda_num, g_fun, Delta[rank]);
            
            update_psi_gamma_nu_l0_z0_RNA(L, Z, W, Lambda_PT[rank], psi_PT[rank], gamma_R_PT[rank], nu_R_PT[rank], &l_R0_PT[rank], &z_R0_PT[rank], M, m, S_num, T_num, C_num, K_min_num, K_max_num, G_num, a_psi, b_psi, a_gamma_R_num, b_gamma_R_num, a_nu_R_num, b_nu_R_num, g_fun, Delta[rank]);
            
            calc_logpost_RNA_fixed_W(&loglik_PT[rank], &logpost_PT[rank], L, Z, W, Lambda_PT[rank], psi_PT[rank], gamma_R_PT[rank], nu_R_PT[rank], &l_R0_PT[rank], &z_R0_PT[rank], M, m, S, T, C, G, a_lambda, b_lambda, a_psi, b_psi, a_gamma_R, b_gamma_R, a_nu_R, b_nu_R, g_fun);
            
            
	    }
	    
	    for(rank = 0; rank < (nThread_num - 1); rank++){
	    
	        GetRNGstate();
            u = runif(0.0, 1.0);
            PutRNGstate();
	        
	        rank_partner = rank + 1;
	        lalpha = (1/Delta[rank] - 1/Delta[rank_partner]) * (loglik_PT[rank_partner] - loglik_PT[rank]);
	        
	        if(log(u) < lalpha){
                PT_swap_RNA_fixed_W(&Lambda_PT[rank], &psi_PT[rank], &gamma_R_PT[rank], &nu_R_PT[rank], &l_R0_PT[rank], &z_R0_PT[rank], &loglik_PT[rank], &logpost_PT[rank], &Lambda_PT[rank_partner], &psi_PT[rank_partner], &gamma_R_PT[rank_partner], &nu_R_PT[rank_partner], &l_R0_PT[rank_partner], &z_R0_PT[rank_partner], &loglik_PT[rank_partner], &logpost_PT[rank_partner]);
            }
	    }
	
	} // end for(i)
	
    ///////////////////////////////////////////////////////////////////////////////
	// End MCMC update
	///////////////////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////
    // Setting updated parameter values = current PT array values
    ////////////////////////////////////////////////////////////////////
    
    
    for(rank = 0; rank < nThread_num; rank++){
        
        // Setting Lambda
        for(g = 0; g < G_num; g++){
            for(c = 0; c < C_num; c++){
                Lambda_PTR[rank * G_num * C_num + c * G_num + g] = Lambda_PT[rank][c * G_num + g];
            }
        }
        
        
        // Setting psi, gamma_R, nu_R
        for(t = 0; t < T_num; t++){
            psi_PTR[rank * T_num + t] = psi_PT[rank][t];
            gamma_R_PTR[rank * T_num + t] = gamma_R_PT[rank][t];
            nu_R_PTR[rank * T_num + t] = nu_R_PT[rank][t];
        }
        
        // Setting l_R0, z_R0
        l_R0_PTR[rank] = l_R0_PT[rank];
        z_R0_PTR[rank] = z_R0_PT[rank];
        
        loglik_PTR[rank] = loglik_PT[rank];
	    logpost_PTR[rank] = logpost_PT[rank];
	
	} // end for(rank)
	
	
	////////////////////////////////////////////////////////////////////
    // Memory free
    ////////////////////////////////////////////////////////////////////
	
    
	for(rank = 0; rank < nThread_num; rank++){
		free(Lambda_PT[rank]);
		
		free(psi_PT[rank]);
		free(gamma_R_PT[rank]);
		free(nu_R_PT[rank]);
	}
    
    return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// End RNA only..
/////////////////////////////////////////////////////////////////////////////////////////////////////////

} // closing extern "C"
