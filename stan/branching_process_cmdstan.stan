functions {
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

    matrix calculate_bin_lpmf_matrix(
        int max_N_cs_actual,
        int N_cs_obs,
        vector bin_successes,
        real rho
    ) {
        matrix[N_cs_obs + 1, max_N_cs_actual + 1] bin_lpmf;
        real log_rho = log(rho);
        real log_one_minus_rho = log(1 - rho);
        int tmp_int;

        for (j in 1:(max_N_cs_actual + 1)) {
            tmp_int = j;
            if(j > N_cs_obs + 1){
            	tmp_int = N_cs_obs + 1;
            }
            
            bin_lpmf[1:tmp_int, j] = bin_successes[1:tmp_int] * log_rho;
            bin_lpmf[1:tmp_int, j] += ((rep_vector(j - 1, tmp_int) - bin_successes[1:tmp_int]) * log_one_minus_rho);
            for (i in 1:tmp_int) {
                bin_lpmf[i, j] += lchoose(1.*(j - 1), bin_successes[i]);
            }
            if (tmp_int < (N_cs_obs + 1)) {
                bin_lpmf[(tmp_int + 1):(N_cs_obs + 1), j] = rep_vector(negative_infinity(), N_cs_obs + 1 - tmp_int);
            }
        }
        return bin_lpmf;
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

    real log_dens_for_given_index(
        array[,] int index_cases,
        int start,
        int end,
        int N_sbts,
        int max_N_cs_actual,
        array[] int N_jcases,
        array[,,] int cs_obs_pre,
        vector log_r0,
        vector vmr_minus_one,
        matrix bin_lpmf
    ) {
        real lpmf = 0.0;
        int C_slice = end - start + 1;
        vector[N_sbts] r0 = exp(log_r0);

        for (c_slice in 1:C_slice) {
            int c = c_slice + start - 1;
            for (s in 1:N_sbts) {
                vector[max_N_cs_actual] ones = rep_vector(1., max_N_cs_actual);
                
                if (sum(cs_obs_pre[c, , s]) == 0) continue;
                row_vector[max_N_cs_actual] cs_actual_lpmf = lpmf_actual_cs_for_given_index(
                    max_N_cs_actual, c, r0[s], vmr_minus_one[s]
                );
                cs_actual_lpmf[1:max_N_cs_actual] -=  log( exp(cs_actual_lpmf[1:max_N_cs_actual]) * ones[1:max_N_cs_actual]);
                
                row_vector[max_N_cs_actual] cs_obs_pre_lpmf = lpmf_obs_cs(
                    max_N_cs_actual, max_N_cs_actual, cs_actual_lpmf, bin_lpmf, c
                );
                
                lpmf += to_row_vector(cs_obs_pre[c,1:N_jcases[s],s]) * cs_obs_pre_lpmf[1:N_jcases[s] ]';
            }
        }
        return lpmf;
    }
}

data {
    int<lower=1> N_sbts;
    int<lower=1> N_cs_obs;
    array[N_sbts] int<lower=1> N_cs_actual;
    int<lower=1> N_icases;
    array[N_sbts] int<lower=1> N_jcases;
    array[N_icases, N_sbts] int<lower=-1> index_cases;
    array[N_icases, max(N_jcases), N_sbts] int<lower=0> cs_obs;
    int<lower=1> sampling_n;
    int<lower=1, upper=sampling_n> sampling_k;
}

transformed data {
    int<lower=1> max_N_cs_actual = max(N_cs_actual);
    array[N_icases, max(N_jcases), N_sbts] int<lower=0> cs_obs_pre = cs_obs;
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
}

transformed parameters {
    matrix[max_N_cs_actual + 1, max_N_cs_actual + 1] bin_lpmf = calculate_bin_lpmf_matrix(
        max_N_cs_actual, max_N_cs_actual, bin_successes, rho
    );
}

model {
    target += normal_lpdf(log_r0_overall | -0.55751, 0.5); // maybe need to adjust these priors
    target += exponential_lpdf(log_r0_sd | 0.1);
    target += normal_lpdf(log_r0_sbts | log_r0_overall, log_r0_sd);
    target += exponential_lpdf(inv_vmr_minus_one_overall | 1);
    target += exponential_lpdf(vmr_minus_one_sbts | inv_vmr_minus_one_overall);
    target += beta_lpdf(rho | sampling_k + 0.5, sampling_n - sampling_k + 0.5);
    target += reduce_sum(
    	log_dens_for_given_index, 
    	index_cases, 
    	1,
        N_sbts,
        max_N_cs_actual,
        N_jcases,
        cs_obs_pre,
        log_r0_sbts,
        vmr_minus_one_sbts,
        bin_lpmf
    );
}
