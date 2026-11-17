functions{
    row_vector lpmf_actual_cs_for_given_index(
        int N_cs_actual,
        real index_cases_local,
        real r0,
        real vmr_minus_one
    ) {
        row_vector[N_cs_actual] cs_actual_lpmf;
        real kappa = r0 / vmr_minus_one;
        real log_r0_div_kappa = log(vmr_minus_one);
        real log_vmr = log(1 + vmr_minus_one);
        vector[6] ones = rep_vector(1., 6);
        matrix[N_cs_actual, 6] tmp;

        for (i in 1:N_cs_actual) {
            tmp[i, 1] = (1.*index_cases_local) / (1.*(index_cases_local + (i - 1)));
            tmp[i, 2] = kappa * index_cases_local + kappa * (i - 1) + (i - 1);
            tmp[i, 3] = kappa * index_cases_local + kappa * (i - 1);
            tmp[i, 4] = (i - 1) + 1;
            tmp[i, 5] = (i - 1);
            tmp[i, 6] = kappa * index_cases_local + kappa * (i - 1) + (i - 1);
        }
        tmp[, 1] = log(tmp[, 1]);
        tmp[, 2] = lgamma(tmp[, 2]);
        tmp[, 3] = -lgamma(tmp[, 3]);
        tmp[, 4] = -lgamma(tmp[, 4]);
        tmp[, 5] *= log_r0_div_kappa;
        tmp[, 6] *= -log_vmr;

        cs_actual_lpmf = (tmp * ones)';
        return cs_actual_lpmf;
    }
    
   row_vector lpmf_obs_cs(
        int N_cs_actual,
        int N_jcases_sbt,
        row_vector cs_actual_lpmf,
        matrix bin_lpmf,
        int c_local
    ) {
        row_vector[N_jcases_sbt] cs_obs_lpmf;
        
        cs_obs_lpmf[1] = log_sum_exp(cs_actual_lpmf[1:N_cs_actual] + bin_lpmf[1, 1:N_cs_actual]);
        for (i in 1:(N_jcases_sbt - 1)) {
            cs_obs_lpmf[i + 1] = log_sum_exp(cs_actual_lpmf[(i + 1):N_cs_actual] + bin_lpmf[(i + 1), (i + 1):N_cs_actual]);
        }
        return cs_obs_lpmf;
    }
}


data {
    int<lower=1> N_sbts;                    
    int<lower=1> N_cs_obs;				    
    array[N_sbts] int<lower=1> N_cs_actual;	
    int<lower=1> N_icases;					
    array[N_sbts] int<lower=1> N_jcases;
    int<lower=1> sampling_n;                
    int<lower=1, upper=sampling_n> sampling_k;
    array[N_icases, max(N_jcases), N_sbts] int<lower=0> cs_obs; 
    int<lower=0> N_origins;					// number of geographic origins
    int<lower=0> N_sgs_m;					// number of obs subgraphs for index case m: replaced by N_sgs_m
    matrix<lower=0>[N_sbts,N_origins] pr_origins;		  // proportion of subgraphs from each origin
    int<lower=1, upper=N_icases> INDEX_CASE_IDX;
    array[max(N_cs_actual)] int<lower=1> z;
}

transformed data {
    int<lower=1> max_N_cs_actual = max(N_cs_actual);      
    vector[max_N_cs_actual + 1] bin_successes;   
    for (i in 0:max_N_cs_actual) {
        bin_successes[i + 1] = i;
    }
}

parameters {
    real<lower=-23, upper=0> log_r0_overall;
    real<lower=0> log_r0_sd;
    vector<lower=-23, upper=0>[N_sbts] log_r0_sbts;
    real<lower=0> inv_vmr_minus_one_overall;
    vector<lower=0>[N_sbts] vmr_minus_one_sbts;
    real<lower=0, upper=1> rho;
    matrix[max_N_cs_actual+1, max_N_cs_actual+1] bin_lpmf;
}

generated quantities {
    array[N_sbts] real<lower=0> N_sgs_pre_st;    // the total number of subgraphs by subtype
    array[N_sbts] int<lower=0> N_sgs_unobs;      // the estimated number of unobserved subgraphs by subtype
    array[N_sbts] real<lower=0> N_sgs_e;         // the estimated number of subgraphs
    array[N_sbts] int<lower=0> N_inf_x;
    array[N_sbts] int<lower=0> N_inf_e;
    array[N_sbts,N_sgs_m] int actual_cs;
    array[N_sbts,5000] int actual_cs_unobs;
    matrix[N_sbts,N_sgs_m] obs_cs_pre;           // may not be necessary
    matrix[N_sbts,N_sgs_m] origins_subgraphs_x;
    matrix[N_sbts,5000] origins_subgraphs_e;
    vector[max_N_cs_actual] unsampled_ch;
    array[N_sbts] real<lower=0> pr_unsampled;

    {
        int c = INDEX_CASE_IDX;
        matrix[N_sbts,max_N_cs_actual] cs_actual_lpmf;
        matrix[N_sbts,max_N_cs_actual] cs_actual_pmf;
        vector[max_N_cs_actual+1] ones = rep_vector(1., max_N_cs_actual+1);
        vector[N_sbts] r0 = exp( log_r0_sbts );
        matrix[N_sbts,max_N_cs_actual] cs_obs_pre_lpmf;
        matrix[N_sbts,max_N_cs_actual] cs_obs_pre_pmf;

        for(s in 1:N_sbts){
            N_sgs_pre_st[s] = sum(cs_obs[c,,s]);

            cs_actual_lpmf[s,] = lpmf_actual_cs_for_given_index(
                    max_N_cs_actual, c, r0[s], vmr_minus_one_sbts[s]);
          
            cs_actual_lpmf[s,] -=  log( exp(cs_actual_lpmf[s,]) * ones[1:max_N_cs_actual]);

            cs_actual_pmf[s,] = exp(cs_actual_lpmf[s,]);


          	// prob a chain is unsampled
           	for(i in 1:max_N_cs_actual){
           	  	unsampled_ch[i] = ((1-rho)^i / i*1.0) / c*1.0;
           	}
            pr_unsampled[s] = unsampled_ch' * cs_actual_pmf[s,]'; 
            
            // number of unobserved clusters    
            if(N_sgs_pre_st[s] > 0){
              	N_sgs_unobs[s] = neg_binomial_rng(N_sgs_pre_st[s], (1-pr_unsampled[s]) / pr_unsampled[s]);
            }else{
               	N_sgs_unobs[s] = 0;
            }
            
            // may not be necessary                
            N_inf_x[s] = 0;
            for(i in 1:N_sgs_m)
            {
                if(i<=N_sgs_pre_st[s]){
                    actual_cs[s,i] =  categorical_rng(to_vector(cs_actual_pmf[s,]));
                    N_inf_x[s] = N_inf_x[s] + (actual_cs[s,i] - 1);
                }else{
                    actual_cs[s,i] =  0;
                }
            }
            
            // sum of observed and unobserved clusters (actual number of clusters by subtype: l_s)
            N_sgs_e[s] = N_sgs_pre_st[s] + N_sgs_unobs[s];
            
            N_inf_e[s] = 0;
            for(i in 1:5000)
            {
                if(i<=N_sgs_e[s]){
                    actual_cs_unobs[s,i] =  categorical_rng(to_vector(cs_actual_pmf[s,]));
                    N_inf_e[s] = N_inf_e[s] + (actual_cs_unobs[s,i]);
                }else{
                    actual_cs_unobs[s,i] =  0;
                }
            }
                                   
            cs_obs_pre_lpmf[s,1:(max_N_cs_actual)] =
                lpmf_obs_cs(
                    max_N_cs_actual,
                    max_N_cs_actual,
                    cs_actual_lpmf[s,],
                    bin_lpmf,
                    c);

            cs_obs_pre_lpmf[s,1:(max_N_cs_actual)] -=  log( exp(cs_obs_pre_lpmf[s,1:(max_N_cs_actual)]) * ones[1:(max_N_cs_actual)]);
			cs_obs_pre_pmf[s,] = exp(cs_obs_pre_lpmf[s,]);

			// Not necessary
            for(i in 1:N_sgs_m)
            {
              if(i<=N_sgs_pre_st[s]){
                obs_cs_pre[s,i] =  categorical_rng(to_vector(cs_obs_pre_pmf[s,]));
              }else{
                obs_cs_pre[s,i] =  positive_infinity();
              }
            }
            
            for(i in 1:N_sgs_m)
            {
            	if(i<=N_sgs_pre_st[s]){
                origins_subgraphs_x[s,i] = categorical_rng(to_vector(pr_origins[s,]));
            	}else{
            		origins_subgraphs_x[s,i] =  positive_infinity();
            	}
            }

            for(i in 1:5000)
            {
                if(i<=N_sgs_e[s]){
                    origins_subgraphs_e[s,i] = categorical_rng(to_vector(pr_origins[s,]));
                }else{
                    origins_subgraphs_e[s,i] = positive_infinity();
                }
            }               
        }
    }
}
