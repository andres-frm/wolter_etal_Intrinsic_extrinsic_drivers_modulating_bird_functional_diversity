sapply(c('dplyr', 'ggplot2', 'patchwork', 
         'rethinking', 'cmdstanr', 'magrittr', 
         'forcats', 'cowplot', 'readxl', 'sf'), library, character.only = T)

extrafont::loadfonts('win')

source('functions_mod_diagnostics.r')

bayesian_R2 <- 
  function(fit, y_obs, mu = 'mu', boot = 2e3) {
    y_pred <- fit$draws(mu, format = 'matrix')
    
    r2 <- 
      apply(y_pred, 1, FUN = 
              function(x) {
                y_p <- x
                
                var_residual <- var(y_obs - y_p)
                var_total <- var(y_obs)
                
                1 - (var_residual / var_total)
                
              })
    
    sapply(1:boot, FUN = 
             function(x) {
               mean(sample(r2, 1e3, replace = F))
             })
  }

d <- read_xlsx('full_data_for_model_03_apr_25.xlsx', sheet = 1, 
               col_names = T, na = 'NA')
str(d)

cords <- d[, c("Individual", "Lon", "Lat")]

# latitude/longitude coordinates to UTM
rs_sf <- st_as_sf(cords, coords = c("Lon", "Lat"), crs = 4326)

rs_utm <- st_transform(rs_sf, crs = 31982)

utm_coords <- st_coordinates(rs_utm) |> 
  as_tibble() |> 
  rename(x_utm = X, y_utm = Y)

# distance matrix
dist_inds <- 
  sapply(seq_along(utm_coords$x_utm), FUN = 
           function(i) {
             x <- utm_coords$x_utm[i]
             y <- utm_coords$y_utm[i]
             
             dist <- sqrt((x - utm_coords$x_utm)^2 + (y - utm_coords$y_utm)^2)
             
           })

colnames(dist_inds) <- d$Individual
rownames(dist_inds) <- d$Individual

responses <- d[, grep('obs_22', colnames(d))]
predictors <- d[, c('S', "Height", "Crop", "Duration", 
                    "Fruit", "Tree", "Infrastructure")]

apply(predictors, 2, function(x) sum(is.na(x)))

predictors <- apply(predictors, 2, function(x) as.vector(scale(x)))

d$ind_id <- as.numeric(gsub('^(I)([a-z]*)([0-9]*)$', '\\3', d$Individual))

dat_1 <- cbind(d[, "ind_id"], responses, predictors) |> as_tibble()

dat_1 <- 
  dat_1 |> 
  rename(S_std = S)

dat_1$S <- sapply(d$S, function(x) if (x == 0) 0.1 else x)

colnames(dat_1) <- tolower(colnames(dat_1))

summary(dat_1)

# ===== model 1: effects on species richness
# 

dat <- lapply(dat_1, function(x) x)

dat$N <- length(dat$ind_id)
dat$N_ind <- max(dat$ind_id)
dat$dist_mat <- dist_inds
names(dat)[2:4] <- c('fric', 'fdiv', 'freg')
names(dat)

dat$fric <- as.vector(scale(dat$fric))

# ========== Model FRic ======

cat(file = 'model_FRic.stan', 
    '
    functions {
      matrix cov_GPL2(matrix x, 
                      real eta,
                      real rho,
                      real delta) {
                      
                      int N = dims(x)[1];
                      matrix[N, N] K;
                      matrix[N, N] L_K;
    
                      for (i in 1:(N-1)) {
                        K[i, i] = eta + delta;
                        for (j in (i+1):N) {
                          K[i, j] = eta * exp(-rho * square(x[i, j]));
                          K[j, i] = K[i, j];
                        }
                      }
                      
                      K[N, N] = eta + delta;
                      L_K = cholesky_decompose(K);
                      return L_K;
                      }
    }
    
    data {
      int N;
      int N_ind;
      matrix[N_ind, N_ind] dist_mat;
      array[N] int ind_id;
      vector[N] fric;
      vector[N] s_std;
      vector[N] height;
      vector[N] crop;
      vector[N] duration;
      vector[N] fruit;
      vector[N] tree;
      vector[N] infrastructure;
    }
    
    parameters {
    
      
      // model FRic
      real alpha_FRic;
      real<lower = 0> sigma_FRic;
      real beta_S;
      real beta_forest_FRic;
      real beta_height_FRic;
      real beta_crop_FRic;
      real beta_duration_FRic;
      real beta_fruit_FRic;
      real beta_urban_FRic;
      vector[N_ind] z_ind_FRic;
      real<lower = 0> eta_FRic;
      real<lower = 0> rho_FRic;
    }
    
    transformed parameters {
      
    
      // pars FRic
      vector[N_ind] ind_FRic;
      matrix[N_ind, N_ind] L_K_FRic;
      L_K_FRic = cov_GPL2(dist_mat, eta_FRic, rho_FRic, 0.001);
      ind_FRic = L_K_FRic * z_ind_FRic;
    }
    
    model {
    
      // prior FRic
      alpha_FRic ~ normal(0, 1);
      sigma_FRic ~ exponential(1);
      z_ind_FRic ~ normal(0, 0.25);
      eta_FRic ~ exponential(4);
      rho_FRic ~ exponential(1);
      beta_forest_FRic ~ normal(0, 0.25);
      beta_height_FRic ~ normal(0, 0.25);
      beta_crop_FRic ~ normal(0, 0.25);
      beta_duration_FRic ~ normal(0, 0.25);
      beta_fruit_FRic ~ normal(0, 0.25);
      beta_urban_FRic ~ normal(0, 0.25);
      beta_S ~ normal(0, 1);
    
      // model FRic
      fric ~ student_t(7, 
                        alpha_FRic + beta_forest_FRic * tree +
                        beta_height_FRic * height +
                        beta_crop_FRic * crop +
                        beta_duration_FRic * duration +
                        beta_fruit_FRic * fruit +
                        beta_urban_FRic * infrastructure +
                        ind_FRic[ind_id], sigma_FRic);
    }
    
    generated quantities {
      
      array[N] real ppcheck_FRic;
      vector[N] mu_FRic;
        
      mu_FRic = alpha_FRic + beta_forest_FRic * tree +
                beta_height_FRic * height +
                beta_crop_FRic * crop +
                beta_duration_FRic * duration +
                beta_fruit_FRic * fruit +
                beta_urban_FRic * infrastructure +
                ind_FRic[ind_id];
    
      ppcheck_FRic = student_t_rng(7, mu_FRic, sigma_FRic);
    
    }
    ')

cat(file = 'model_FRic.stan', 
    '
    functions {
      matrix cov_GPL2(matrix x, 
                      real eta,
                      real rho,
                      real delta) {
                      
                      int N = dims(x)[1];
                      matrix[N, N] K;
                      matrix[N, N] L_K;
    
                      for (i in 1:(N-1)) {
                        K[i, i] = eta + delta;
                        for (j in (i+1):N) {
                          K[i, j] = eta * exp(-rho * square(x[i, j]));
                          K[j, i] = K[i, j];
                        }
                      }
                      
                      K[N, N] = eta + delta;
                      L_K = cholesky_decompose(K);
                      return L_K;
                      }
    }
    
    data {
      int N;
      int N_ind;
      matrix[N_ind, N_ind] dist_mat;
      array[N] int ind_id;
      vector[N] fric;
      vector[N] s_std;
      vector[N] height;
      vector[N] crop;
      vector[N] duration;
      vector[N] fruit;
      vector[N] tree;
      vector[N] infrastructure;
    }
    
    parameters {
    
      // model S
      real alpha_S;
      real<lower = 0> sigma_S;
      real beta_forest_S;
      real beta_height_S;
      real beta_crop_S;
      real beta_duration_S;
      real beta_fruit_S;
      real beta_urban_S;
      vector[N_ind] z_ind_S;
      real<lower = 0> eta_S;
      real<lower = 0> rho_S;
    
      // model FRic
      real alpha_FRic;
      real<lower = 0> sigma_FRic;
      real beta_S;
      real beta_forest_FRic;
      real beta_height_FRic;
      real beta_crop_FRic;
      real beta_duration_FRic;
      real beta_fruit_FRic;
      real beta_urban_FRic;
      vector[N_ind] z_ind_FRic;
      real<lower = 0> eta_FRic;
      real<lower = 0> rho_FRic;
    }
    
    transformed parameters {
      // pars S
      vector[N_ind] ind_S;
      matrix[N_ind, N_ind] L_K_S;
      L_K_S = cov_GPL2(dist_mat, eta_S, rho_S, 0.001);
      ind_S = L_K_S * z_ind_S;
    
      // pars FRic
      vector[N_ind] ind_FRic;
      matrix[N_ind, N_ind] L_K_FRic;
      L_K_FRic = cov_GPL2(dist_mat, eta_FRic, rho_FRic, 0.001);
      ind_FRic = L_K_FRic * z_ind_FRic;
    }
    
    model {
    
      // prior S
      alpha_S ~ normal(0, 1);
      sigma_S ~ exponential(1);
      z_ind_S ~ normal(0, 0.25);
      eta_S ~ exponential(4);
      rho_S ~ exponential(1);
      beta_forest_S ~ normal(0, 0.25);
      beta_height_S ~ normal(0, 0.25);
      beta_crop_S ~ normal(0, 0.25);
      beta_duration_S ~ normal(0, 0.25);
      beta_fruit_S ~ normal(0, 0.25);
      beta_urban_S ~ normal(0, 0.25);
    
      // prior FRic
      alpha_FRic ~ normal(0, 1);
      sigma_FRic ~ exponential(1);
      z_ind_FRic ~ normal(0, 0.25);
      eta_FRic ~ exponential(4);
      rho_FRic ~ exponential(1);
      beta_forest_FRic ~ normal(0, 0.25);
      beta_height_FRic ~ normal(0, 0.25);
      beta_crop_FRic ~ normal(0, 0.25);
      beta_duration_FRic ~ normal(0, 0.25);
      beta_fruit_FRic ~ normal(0, 0.25);
      beta_urban_FRic ~ normal(0, 0.25);
      beta_S ~ normal(0, 1);
    
      // model S
      s_std ~ student_t(7, 
                        alpha_S + 
                        beta_forest_S * tree +
                        beta_height_S * height +
                        beta_crop_S * crop +
                        beta_duration_S * duration +
                        beta_fruit_S * fruit +
                        beta_urban_S * infrastructure +
                        ind_S[ind_id], sigma_S);
    
      // model FRic
      fric ~ student_t(7, 
                        alpha_FRic + 
                        beta_S * s_std +
                        beta_forest_FRic * tree +
                        beta_height_FRic * height +
                        beta_crop_FRic * crop +
                        beta_duration_FRic * duration +
                        beta_fruit_FRic * fruit +
                        beta_urban_FRic * infrastructure +
                        ind_FRic[ind_id], sigma_FRic);
    }
    
    generated quantities {
      array[N] real ppcheck_S;
      array[N] real ppcheck_FRic;
      vector[N] mu_S;
      vector[N] mu_FRic;
    
      mu_S = alpha_S + 
             beta_forest_S * tree +
             beta_height_S * height +
             beta_crop_S * crop +
             beta_duration_S * duration +
             beta_fruit_S * fruit +
             beta_urban_S * infrastructure +
             ind_S[ind_id];
    
      ppcheck_S = student_t_rng(7, mu_S, sigma_S);
    
        
      mu_FRic = alpha_FRic + 
                beta_S * s_std +
                beta_forest_FRic * tree +
                beta_height_FRic * height +
                beta_crop_FRic * crop +
                beta_duration_FRic * duration +
                beta_fruit_FRic * fruit +
                beta_urban_FRic * infrastructure +
                ind_FRic[ind_id];
    
      ppcheck_FRic = student_t_rng(7, mu_FRic, sigma_FRic);
    
    }
    ')

file <- paste0(getwd(), '/model_FRic.stan')
fit_FRic <- cmdstan_model(file, compile = T)

indx_na_vars <- 
  unlist(lapply(dat, function(x) sum(is.na(x)) == 0), 
         use.names = F)

mod_FRic <- 
  fit_FRic$sample(
    data = dat[indx_na_vars], 
    iter_warmup = 500, 
    iter_sampling = 4e3,
    chains = 3, 
    parallel_chains = 3,
    thin = 3,
    seed = 5
  )

sum_FRic <- mod_FRic$summary()

mod_diagnostics(mod_FRic, sum_FRic)

sum_FRic |> print(n = 20)

indx_pars <- 
  c(grep('^beta', sum_FRic$variable), 
    grep('^sigma', sum_FRic$variable), 
    grep('^alpha', sum_FRic$variable), 
    grep('^ind', sum_FRic$variable))

length(indx_pars)

par(mfrow = c(4, 5), mar = c(4, 4, 1, 1))
for (i in 1:17) {
  trace_plot(mod_FRic, sum_FRic$variable[indx_pars[i]], 3)
}
par(mfrow = c(1, 1))

ppcheck_S <- mod_FRic$draws('ppcheck_S', format = 'matrix')
ppcheck_FRic <- mod_FRic$draws('ppcheck_FRic', format = 'matrix')

par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
plot(density(dat$s_std), xlim = c(-4, 5), 
     ylim = c(0, 0.7), main = '', xlab = 'N species (z-scores)')
for (i in 1:500) lines(density(ppcheck_S[i, ]), lwd = 0.1)
lines(density(dat$s_std), col = 'red', lwd = 1.5)

plot(density(dat$fric), xlim = c(-4, 4), 
     ylim = c(0, 0.7), main = '', xlab = 'FRic (z-scores)')
for (i in 1:500) lines(density(ppcheck_FRic[i, ]), lwd = 0.1)
lines(density(dat$fric), col = 'red', lwd = 1.5)
par(mfrow = c(1, 1))

par(mfrow = c(2, 2))
plot(density(apply(ppcheck_S, 1, mean)), xlab = 'AVG. S (estimated)', 
     main = '')
abline(v = mean(dat$s_std), col = 'red')
plot(density(apply(ppcheck_S, 1, sd)), xlab = 'SD. S (estimated)', 
     main = '')
abline(v = sd(dat$s_std), col = 'red')

plot(density(apply(ppcheck_FRic, 1, mean)), xlab = 'AVG. FRic (estimated)', 
     main = '')
abline(v = mean(dat$fric), col = 'red')
plot(density(apply(ppcheck_FRic, 1, sd)), xlab = 'SD. FRic (estimated)', 
     main = '')
abline(v = sd(dat$fric), col = 'red')
par(mfrow = c(1, 1))

R2_S <- bayesian_R2(mod_FRic, dat$s_std, 'mu_S')
quantile(R2_S)

R2_FRic <- bayesian_R2(mod_FRic, dat$s_std, 'mu_FRic')
quantile(R2_FRic)

plot(density(R2_S), col = 'lightblue', main = '', lwd = 3)

betas_S <- mod_FRic$draws(sum_FRic$variable[grep('^bet(.*)(_S)$', sum_FRic$variable)], 
                       format = 'df')
betas_S <- betas_S[, grep('^beta', colnames(betas_S))]

apply(betas_S, 2, function(x) mean(x > 0))
apply(betas_S, 2, function(x) mean(x))
apply(betas_S, 2, function(x) sd(x))
apply(betas_S, 2, 
      function(x) {
        tibble(li = quantile(x, 0.025), 
               ls = quantile(x, 0.975))
      })

betas_FRic <- mod_FRic$draws(sum_FRic$variable[grep('^bet(.*)(_FRic)$', 
                                                    sum_FRic$variable)], 
                          format = 'df')
betas_FRic <- betas_FRic[, grep('^beta', colnames(betas_FRic))]

betas_FRic$beta_S <- betas_S$beta_S
apply(betas_FRic, 2, function(x) mean(x > 0))
apply(betas_FRic, 2, function(x) mean(x))
apply(betas_FRic, 2, function(x) sd(x))
apply(betas_FRic, 2, 
      function(x) {
        tibble(li = quantile(x, 0.025), 
               ls = quantile(x, 0.975))
      })

# ========= conditional effects =========

sum_FRic[indx_pars, ]
grep('^bet(.*)(_FRic)$', 
     sum_FRic$variable[indx_pars])
post_FRic <- mod_FRic$draws(sum_FRic[indx_pars, ]$variable, format = 'df')

post_FRic <- post_FRic[, sum_FRic[indx_pars, ]$variable]

apply(post_FRic, 2, function(x) mean(x > 0))
apply(post_FRic[, grep('^beta', colnames(post_FRic))], 2, mean)
apply(post_FRic[, grep('^beta', colnames(post_FRic))], 2, 
      function(x) mean(x > 0))

# conditional effects 

post_FRic 

conditional_effect_FRic <- 
  function(par, x_var, N = 500, FRic = T) {
    
    x <- seq(min(x_var), max(x_var), length.out = N)
    
    ind_FRic <- grep('^ind_FRic', colnames(post_FRic))
    ind_S <- grep('^ind_S', colnames(post_FRic))
    
    est <- 
      lapply(x, FUN = 
               function(i) {
                 if (FRic) {
                   mu_ <- 
                     post_FRic$alpha_FRic +
                     post_FRic[[par]] * i +
                     apply(post_FRic[, ind_FRic], 1, mean)
                   
                   mu_pred <- 
                     rstudent(N, 7, mu_, post_FRic$sigma_FRic)
                 } else {
                   mu_ <- 
                     post_FRic$alpha_S +
                     post_FRic[[par]] * i +
                     apply(post_FRic[, ind_S], 1, mean)
                   
                   mu_pred <- 
                     rstudent(N, 7, mu_, post_FRic$sigma_S)
                 }
                 
                 tibble(mu = mean(mu_), 
                        li = quantile(mu_, 0.025), 
                        ls = quantile(mu_, 0.975), 
                        li_pred = quantile(mu_pred, 0.025), 
                        ls_pred = quantile(mu_pred, 0.975))
               })
    est <- do.call('rbind', est)
    est$x <- x
    est
    
  }

# ==== effect height =====

mean(post_FRic$beta_height_S > post_FRic$beta_height_FRic)

rel_imp_height_FRic <- 
  post_FRic$beta_height_FRic / post_FRic$beta_height_S

median(rel_imp_height_FRic)

est_height_FRic <- conditional_effect_FRic('beta_height_FRic', dat$height)
est_height_S <- conditional_effect_FRic('beta_height_S', dat$height)

plot_cond_heigth_FRic <- 
  ggplot() +
  geom_ribbon(data = est_height_FRic, 
              aes(x, mu, ymin = li_pred, ymax = ls_pred), 
              alpha = 0.5, fill = 'lightblue2') +
  geom_ribbon(data = est_height_FRic, 
              aes(x, mu, ymin = li, ymax = ls), 
              alpha = 0.7, fill = 'lightblue') +
  geom_line(data = est_height_FRic, 
            aes(x, mu), color = 'tan1', linewidth = 1, alpha = 0.7) +
  geom_point(data = tibble(x = dat$height, 
                           y = dat$fric), 
             aes(x, y), color = 'tan1', alpha = 0.5) +
  theme_classic() +
  labs(x = 'Tree height (z-scores)', 
       y = 'FRic (z-scores)') +
  theme(legend.position = c(0.2, 0.9),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = 'Times New Roman', size = 10), 
        axis.line = element_line(linewidth = 0.25))

plot_cond_heigth_S <- 
  ggplot() +
  geom_ribbon(data = est_height_S, 
              aes(x, mu, ymin = li_pred, ymax = ls_pred), 
              alpha = 0.5, fill = 'lightblue2') +
  geom_ribbon(data = est_height_S, 
              aes(x, mu, ymin = li, ymax = ls), 
              alpha = 0.7, fill = 'lightblue') +
  geom_line(data = est_height_S, 
            aes(x, mu), color = 'tan1', linewidth = 1, alpha = 0.7) +
  geom_point(data = tibble(x = dat$height, 
                           y = dat$s_std), 
             aes(x, y), color = 'tan1', alpha = 0.5) +
  theme_classic() +
  labs(x = 'Tree height (z-scores)', 
       y = 'S (z-scores)') +
  theme(legend.position = c(0.2, 0.9),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = 'Times New Roman', size = 10), 
        axis.line = element_line(linewidth = 0.25))


causal_effect_FRic <- 
  function(par_EQ1, par_EQ2, x_var, x_label, total_effect = T, negative = F) {
    ind_s <- grep('^ind_S', colnames(post_FRic))
    ind_FRic <- grep('^ind_FRic', colnames(post_FRic))
    mu_S <- mean(d$S)
    SD_S <- sd(d$S)
    
    mu_FRic <- mean(d$fric_obs_22)
    SD_FRic <- sd(d$fric_obs_22)
    
    estimated <- 
      lapply(c(min, max), FUN = 
               function(.fun) {
                 
                 est_S <- 
                   post_FRic$alpha_S +
                   post_FRic[[par_EQ1]] * .fun(x_var) +
                   apply(post_FRic[, ind_s], 1, mean)
                 
                 if (total_effect) {
                   message('Calculating total effect')
                 } else {
                   message('Calculating conditional effect')
                   est_S <- min(dat$s_std)
                 }
                 
                 est_FRic <- 
                   post_FRic$alpha_FRic +
                   post_FRic$beta_S * est_S +
                   post_FRic[[par_EQ2]] * .fun(x_var) +
                   apply(post_FRic[, ind_FRic], 1, mean)
                 
                 est_FRic <- mu_FRic + est_FRic * SD_FRic
                 est_S <- mu_S + est_S * SD_S
                 
                 list(S = tibble(est = est_S, 
                                 Y = 'S', 
                                 X = x_label), 
                      FRic = tibble(est = est_FRic, 
                                    Y = 'FRic', 
                                    X = x_label))
                 
               })
    
    names(estimated) <- c('Min', 'Max')
    
    estimated <- 
      lapply(names(estimated), FUN = 
               function(x) {
                 d <- estimated[[x]]
                 d$S$intervention <- x
                 d$FRic$intervention <- x
                 d
               })
    names(estimated) <- c('Min', 'Max')
    
    estimated <- 
      lapply(estimated, FUN = 
               function(x) {
                 do.call('rbind', x)
               })
    estimated <- do.call('rbind', estimated)
    estimated <- split(estimated, estimated$Y)
    lapply(estimated, FUN =
             function(x) {
               
               if (negative) {
                 x2 <- x[x$intervention == 'Max', ]
                 x2$est <-
                   x[x$intervention == 'Min', ]$est -
                   x[x$intervention == 'Max', ]$est
                 
                 x2$intervention <- 'Contrast'
                 z <- (x2$est * 100)/x[x$intervention == 'Min', ]$est
                 x$porcentaje <- NA
                 x2$porcentaje <- z
                 x2 <- rbind(x, x2)
                 x2$intervention <- as.factor(x2$intervention)
                 x2$intervention <- factor(x2$intervention, 
                                           levels = c('Min', 'Max', 
                                                      'Contrast'))
                 x2
               } else {
                 x2 <- x[x$intervention == 'Max', ]
                 x2$est <-
                   x[x$intervention == 'Max', ]$est -
                   x[x$intervention == 'Min', ]$est
                 
                 x2$intervention <- 'Contrast'
                 z <- (x2$est * 100)/x[x$intervention == 'Max', ]$est
                 x$porcentaje <- NA
                 x2$porcentaje <- z
                 x2 <- rbind(x, x2)
                 x2$intervention <- as.factor(x2$intervention)
                 x2$intervention <- factor(x2$intervention, 
                                           levels = c('Min', 'Max', 
                                                      'Contrast'))
                 x2
               }
               
             })
  }

effect_height <- causal_effect_FRic('beta_height_S', 
                                    'beta_height_FRic', 
                                    dat$height, 'Height', T, negative = F)

effect_height$FRic$total_effect <- 'Total effect'

conditional_height_FRic <- causal_effect_FRic('beta_height_S', 
                                              'beta_height_FRic', 
                                              dat$height, 'Height', F, negative = F)[[1]]

conditional_height_FRic$total_effect <- 'Conditional effect'

plot_CE_height_FRic <- 
  rbind(effect_height$FRic, 
      conditional_height_FRic) |> 
  group_by(intervention, total_effect) |> 
  transmute(mu = mean(est), 
            li = quantile(est, 0.025), 
            ls = quantile(est, 0.975)) |> 
  unique() |> 
  ggplot() +
  geom_errorbar(aes(intervention, ymin = li, ymax = ls, 
                    color = total_effect), width = 0, 
                position = position_dodge(width = 0.3)) +
  geom_point(aes(intervention, mu, color = total_effect), 
             position = position_dodge(width = 0.3)) +
  scale_color_manual(values = c('tan1', 'lightblue')) +
  scale_fill_manual(values = c('tan1', 'lightblue')) +
  lims(y = c(-25, 50)) +
  theme_classic() +
  labs(x = 'Intervention of tree height', 
       y = 'FRic') +
  geom_hline(yintercept = 0, linetype = 3) +
  theme(legend.position = c(0.37, 0.15),
        legend.background = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = 'Times New Roman', size = 10), 
        axis.line = element_line(linewidth = 0.25)) 

# figure height FRic

plot_heigh_FRic <- 
  plot_grid(plot_cond_heigth_S, 
          plot_cond_heigth_FRic, 
          plot_CE_height_FRic, ncol = 3, 
          labels = c('(i)', '(ii)', '(iii)'), 
          label_fontfamily = 'Times New Roman', 
          label_fontface = 'plain', 
          label_size = 10, 
          label_x = 0.175)


cont_height_FRic <- na.omit(effect_height$FRic)

conditional_height_FRic <- na.omit(conditional_height_FRic)

diff_height_FRic <- cont_height_FRic$est - conditional_height_FRic$est

mean(effect_height$FRic[effect_height$FRic$intervention != 'Min', ]$est >
       effect_height$FRic[effect_height$FRic$intervention == 'Min', ]$est) 

mean(effect_height$S[effect_height$S$intervention != 'Min', ]$est >
       effect_height$S[effect_height$S$intervention == 'Min', ]$est) 

tribble(~ value,               ~estimation, 
        'P(total effect > 0)', mean(cont_height_FRic$est > 0), 
        'mu total effect', mean(cont_height_FRic$est), 
        'sd total effect', sd(cont_height_FRic$est), 
        'AVG change relative to maximum (total) (%)', mean(cont_height_FRic$porcentaje), 
        'SD change relative to maximum (total) (%)', sd(cont_height_FRic$porcentaje), 
        'P(conditional effect > 0)', mean(conditional_height_FRic$est > 0), 
        'mu conditional effect', mean(conditional_height_FRic$est), 
        'sd conditional effect', sd(conditional_height_FRic$est), 
        'AVG change relative to maximum (conditional) (%)', mean(conditional_height_FRic$porcentaje), 
        'SD change relative to maximum (conditional) (%)', sd(conditional_height_FRic$porcentaje),
        'P(total > conditional)', mean(cont_height_FRic$est > conditional_height_FRic$est), 
        'AVG diff total-conditional (%)', (100-(mean(diff_height_FRic)*100) / mean(cont_height_FRic$est)))

plot(density(diff_height_FRic))
# ==== effect urbanization =====

apply(post_FRic[, grep('^beta', colnames(post_FRic))], 2, 
      function(x) mean(x))
apply(post_FRic[, grep('^beta', colnames(post_FRic))], 2, 
      function(x) sd(x))
apply(post_FRic[, grep('^beta', colnames(post_FRic))], 2, 
      function(x) mean(x < 0))

((mean(post_FRic$beta_urban_S) - mean(post_FRic$beta_urban_FRic)) /
  mean(post_FRic$beta_urban_FRic)) * 100

est_urban_FRic <- conditional_effect_FRic('beta_urban_FRic', dat$infrastructure)
est_urban_S <- conditional_effect_FRic('beta_urban_S', dat$infrastructure)

plot_cond_urban_FRic <- 
  ggplot() +
  geom_ribbon(data = est_urban_FRic, 
              aes(x, mu, ymin = li_pred, ymax = ls_pred), 
              alpha = 0.5, fill = 'lightblue2') +
  geom_ribbon(data = est_urban_FRic, 
              aes(x, mu, ymin = li, ymax = ls), 
              alpha = 0.7, fill = 'lightblue') +
  geom_line(data = est_urban_FRic, 
            aes(x, mu), color = 'tan1', linewidth = 1, alpha = 0.7) +
  geom_point(data = tibble(x = dat$infrastructure, 
                           y = dat$fric), 
             aes(x, y), color = 'tan1', alpha = 0.5) +
  theme_classic() +
  labs(x = 'Infrastructure cover (z-scores)', 
       y = 'FRic (z-scores)') +
  theme(legend.position = c(0.2, 0.9),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = 'Times New Roman', size = 10), 
        axis.line = element_line(linewidth = 0.25))

plot_cond_urban_S <- 
  ggplot() +
    geom_ribbon(data = est_urban_S, 
                aes(x, mu, ymin = li_pred, ymax = ls_pred), 
                alpha = 0.5, fill = 'lightblue2') +
    geom_ribbon(data = est_urban_S, 
                aes(x, mu, ymin = li, ymax = ls), 
                alpha = 0.7, fill = 'lightblue') +
    geom_line(data = est_urban_S, 
              aes(x, mu), color = 'tan1', linewidth = 1, alpha = 0.7) +
    geom_point(data = tibble(x = dat$infrastructure, 
                             y = dat$s_std), 
               aes(x, y), color = 'tan1', alpha = 0.5) +
    theme_classic() +
    labs(x = 'Infrastructure cover (z-scores)', 
         y = 'S (z-scores)') +
    theme(legend.position = c(0.2, 0.9),
          legend.title = element_blank(),
          panel.grid = element_blank(),
          text = element_text(family = 'Times New Roman', size = 10), 
          axis.line = element_line(linewidth = 0.25))

effect_urban <- causal_effect_FRic('beta_urban_S', 
                                    'beta_urban_FRic', 
                                    dat$infrastructure, 'Urban', T)

effect_urban$FRic$total_effect <- 'Total effect'

conditional_urban_FRic <- causal_effect_FRic('beta_urban_S', 
                                              'beta_urban_FRic', 
                                              dat$infrastructure, 'Urban', F)[[1]]

conditional_urban_FRic$total_effect <- 'Conditional effect'

plot_CE_urban_FRic <- 
  rbind(effect_urban$FRic, 
        conditional_urban_FRic) |> 
  group_by(intervention, total_effect) |> 
  transmute(mu = mean(est), 
            li = quantile(est, 0.025), 
            ls = quantile(est, 0.975)) |> 
  unique() |> 
  ggplot() +
  geom_errorbar(aes(intervention, ymin = li, ymax = ls, 
                    color = total_effect), width = 0, 
                position = position_dodge(width = 0.3)) +
  geom_point(aes(intervention, mu, color = total_effect), 
             position = position_dodge(width = 0.3)) +
  scale_color_manual(values = c('tan1', 'lightblue')) +
  scale_fill_manual(values = c('tan1', 'lightblue')) +
  theme_classic() +
  labs(x = 'Intervention of\n infrastructure cover', 
       y = 'FRic') +
  lims(y = c(-25, 50)) +
  geom_hline(yintercept = 0, linetype = 3) +
  theme(legend.position = 'none',
        legend.title = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = 'Times New Roman', size = 10), 
        axis.line = element_line(linewidth = 0.25)) 

# figure urban FRic

plot_urban_FRic <- 
  plot_grid(plot_cond_urban_S, 
            plot_cond_urban_FRic, 
            plot_CE_urban_FRic, ncol = 3, 
            labels = c('(i)', '(ii)', '(iii)'), 
            label_fontfamily = 'Times New Roman', 
            label_fontface = 'plain', 
            label_size = 10, 
            label_x = 0.175)


cont_urban_FRic <- na.omit(effect_urban$FRic)

conditional_urban_FRic <- na.omit(conditional_urban_FRic)

diff_urban_FRic <- cont_urban_FRic$est - conditional_urban_FRic$est

tribble(~ value,               ~estimation, 
        'P(total effect < 0)', mean(cont_urban_FRic$est > 0), 
        'mu total effect', mean(cont_urban_FRic$est), 
        'sd total effect', sd(cont_urban_FRic$est), 
        'AVG change relative to maximum (total) (%)', mean(cont_urban_FRic$porcentaje), 
        'SD change relative to maximum (total) (%)', sd(cont_urban_FRic$porcentaje), 
        'P(conditional effect > 0)', mean(conditional_urban_FRic$est > 0), 
        'mu conditional effect', mean(conditional_urban_FRic$est), 
        'sd conditional effect', sd(conditional_urban_FRic$est), 
        'AVG change relative to maximum (conditional) (%)', mean(conditional_urban_FRic$porcentaje), 
        'SD change relative to maximum (conditional) (%)', sd(conditional_urban_FRic$porcentaje),
        'P(total < conditional)', mean(cont_urban_FRic$est > conditional_urban_FRic$est), 
        'AVG diff total-conditional (%)', (mean(diff_urban_FRic)*100 / mean(cont_urban_FRic$est)-100))

# ==== effect forest ======

apply(post_FRic[, grep('^beta', colnames(post_FRic))], 2, 
      function(x) mean(x))
apply(post_FRic[, grep('^beta', colnames(post_FRic))], 2, 
      function(x) sd(x))
apply(post_FRic[, grep('^beta', colnames(post_FRic))], 2, 
      function(x) mean(x < 0))

((mean(post_FRic$beta_forest_FRic) - mean(post_FRic$beta_forest_S)) /
    mean(post_FRic$beta_forest_S)) * 100

est_forest_FRic <- conditional_effect_FRic('beta_forest_FRic', dat$tree)
est_forest_S <- conditional_effect_FRic('beta_forest_S', dat$tree)

plot_cond_forest_FRic <- 
  ggplot() +
  geom_ribbon(data = est_forest_FRic, 
              aes(x, mu, ymin = li_pred, ymax = ls_pred), 
              alpha = 0.5, fill = 'lightblue2') +
  geom_ribbon(data = est_forest_FRic, 
              aes(x, mu, ymin = li, ymax = ls), 
              alpha = 0.7, fill = 'lightblue') +
  geom_line(data = est_forest_FRic, 
            aes(x, mu), color = 'tan1', linewidth = 1, alpha = 0.7) +
  geom_point(data = tibble(x = dat$tree, 
                           y = dat$fric), 
             aes(x, y), color = 'tan1', alpha = 0.5) +
  theme_classic() +
  labs(x = 'Tree cover (z-scores)', 
       y = 'FRic (z-scores)') +
  theme(legend.position = c(0.2, 0.9),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = 'Times New Roman', size = 10), 
        axis.line = element_line(linewidth = 0.25))

plot_cond_forest_S <- 
  ggplot() +
  geom_ribbon(data = est_forest_S, 
              aes(x, mu, ymin = li_pred, ymax = ls_pred), 
              alpha = 0.5, fill = 'lightblue2') +
  geom_ribbon(data = est_forest_S, 
              aes(x, mu, ymin = li, ymax = ls), 
              alpha = 0.7, fill = 'lightblue') +
  geom_line(data = est_forest_S, 
            aes(x, mu), color = 'tan1', linewidth = 1, alpha = 0.7) +
  geom_point(data = tibble(x = dat$infrastructure, 
                           y = dat$s_std), 
             aes(x, y), color = 'tan1', alpha = 0.5) +
  theme_classic() +
  labs(x = 'Tree cover (z-scores)', 
       y = 'S (z-scores)') +
  theme(legend.position = c(0.2, 0.9),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = 'Times New Roman', size = 10), 
        axis.line = element_line(linewidth = 0.25))

effect_forest <- causal_effect_FRic('beta_forest_S', 
                                   'beta_forest_FRic', 
                                   dat$tree, 'forest', T, negative = T)

effect_forest$FRic$total_effect <- 'Total effect'

conditional_forest_FRic <- causal_effect_FRic('beta_forest_S', 
                                             'beta_forest_FRic', 
                                             dat$tree, 'forest', F, negative = T)[[1]]

conditional_forest_FRic$total_effect <- 'Conditional effect'

plot_CE_forest_FRic <- 
  rbind(effect_forest$FRic, 
        conditional_forest_FRic) |> 
  group_by(intervention, total_effect) |> 
  transmute(mu = mean(est), 
            li = quantile(est, 0.025), 
            ls = quantile(est, 0.975)) |> 
  unique() |> 
  ggplot() +
  geom_errorbar(aes(intervention, ymin = li, ymax = ls, 
                    color = total_effect), width = 0, 
                position = position_dodge(width = 0.3)) +
  geom_point(aes(intervention, mu, color = total_effect), 
             position = position_dodge(width = 0.3)) +
  scale_color_manual(values = c('tan1', 'lightblue')) +
  scale_fill_manual(values = c('tan1', 'lightblue')) +
  lims(y = c(-25, 50)) +
  theme_classic() +
  labs(x = 'Intervention of\n tree cover', 
       y = 'FRic') +
  geom_hline(yintercept = 0, linetype = 3) +
  theme(legend.position = 'none',
        legend.title = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = 'Times New Roman', size = 10), 
        axis.line = element_line(linewidth = 0.25)) 

# figure forest FRic

plot_forest_FRic <- 
  plot_grid(plot_cond_forest_S, 
            plot_cond_forest_FRic, 
            plot_CE_forest_FRic, ncol = 3, 
            labels = c('(i)', '(ii)', '(iii)'), 
            label_fontfamily = 'Times New Roman', 
            label_fontface = 'plain', 
            label_size = 10, 
            label_x = 0.175)

cont_forest_FRic <- na.omit(effect_forest$FRic)

conditional_forest_FRic <- na.omit(conditional_forest_FRic)

diff_forest_FRic <- cont_forest_FRic$est - conditional_forest_FRic$est

tribble(~value,               ~estimation, 
        'P(total effect < 0)', mean(cont_forest_FRic$est > 0), 
        'mu total effect', mean(cont_forest_FRic$est), 
        'sd total effect', sd(cont_forest_FRic$est), 
        'AVG change relative to maximum (total) (%)', mean(cont_forest_FRic$porcentaje), 
        'SD change relative to maximum (total) (%)', sd(cont_forest_FRic$porcentaje), 
        'P(conditional effect > 0)', mean(conditional_forest_FRic$est > 0), 
        'mu conditional effect', mean(conditional_forest_FRic$est), 
        'sd conditional effect', sd(conditional_forest_FRic$est), 
        'AVG change relative to maximum (conditional) (%)', median(conditional_forest_FRic$porcentaje), 
        'SD change relative to maximum (conditional) (%)', sd(conditional_forest_FRic$porcentaje),
        'P(total < conditional)', mean(cont_forest_FRic$est > conditional_forest_FRic$est), 
        'AVG diff total-conditional (%)', (mean(diff_forest_FRic)*100 / mean(cont_forest_FRic$est)))


# ===== plots FRic mod ===

plot_grid(plot_heigh_FRic, 
          plot_urban_FRic, 
          plot_forest_FRic,
          ncol = 1, 
          labels = c('(a)', '(b)', '(c)'), 
          label_fontfamily = 'Times New Roman', 
          label_size = 12, 
          label_fontface = 'plain', 
          label_x = -0.007)

ggsave('FRic_mod.jpg', width = 20, height = 18, units = 'cm', dpi = 1000)

apply(post_FRic[, grep('^beta', colnames(post_FRic))], 2, 
      function(x) mean(x))

apply(post_FRic[, grep('^beta', colnames(post_FRic))], 2, 
           function(x) mean(x))

sapply(c('_S', '_FRic'), FUN = 
         function(x) {
           
           v <- 
           apply(post_FRic[, grep('^beta', colnames(post_FRic))], 2, 
                 function(x) mean(x))
           
           v <- v[grep(x, names(v))]
           
           t <- 
             c(grep('urban', names(v)), 
               grep('forest', names(v)))
           
           z <- 
             sapply(v[t], FUN = 
                    function(i) {
                      i <- abs(i)
                      ((v[grep('height', names(v))] - i) / i) * 100
                    })
           
           tibble(mu = mean(z), 
                  sd = sd(z))
           
         })


# ======== model FDiv =====

dat2 <- dat_1

indx_na_FDiv <- dat2$ind_id[which(is.na(dat2$fdiv_obs_22))]

dist_inds_FDiv <- dist_inds[-indx_na_FDiv, -indx_na_FDiv]

dat2 <- dat2[-indx_na_FDiv ,]

dat2$ind_id <- 1:nrow(dat2)

colnames(dist_inds_FDiv) <- paste0('Ind', dat2$ind_id)
rownames(dist_inds_FDiv) <- paste0('Ind', dat2$ind_id)

dat2 <- lapply(dat2, function(x) x)

dat2$N <- length(dat2$ind_id)
dat2$N_ind <- max(dat2$ind_id)
dat2$dist_mat <- dist_inds_FDiv
names(dat2)[2:4] <- c('fric', 'fdiv', 'freg')
names(dat2)

summary(dat2$fdiv)
plot(density(dat2$fdiv))
dat2$fdiv <- as.vector(scale(log(dat2$fdiv)))

cat(file = 'model_FDiv.stan', 
    '
    functions {
      matrix cov_GPL2(matrix x, 
                      real eta,
                      real rho,
                      real delta) {
                      
                      int N = dims(x)[1];
                      matrix[N, N] K;
                      matrix[N, N] L_K;
    
                      for (i in 1:(N-1)) {
                        K[i, i] = eta + delta;
                        for (j in (i+1):N) {
                          K[i, j] = eta * exp(-rho * square(x[i, j]));
                          K[j, i] = K[i, j];
                        }
                      }
                      
                      K[N, N] = eta + delta;
                      L_K = cholesky_decompose(K);
                      return L_K;
                      }
    }
    
    data {
      int N;
      int N_ind;
      matrix[N_ind, N_ind] dist_mat;
      array[N] int ind_id;
      vector[N] fdiv;
      vector[N] s_std;
      vector[N] height;
      vector[N] crop;
      vector[N] duration;
      vector[N] fruit;
      vector[N] tree;
      vector[N] infrastructure;
    }
    
    parameters {
    
      // model S
      real alpha_S;
      real<lower = 0> sigma_S;
      real beta_forest_S;
      real beta_height_S;
      real beta_crop_S;
      real beta_duration_S;
      real beta_fruit_S;
      real beta_urban_S;
      vector[N_ind] z_ind_S;
      real<lower = 0> eta_S;
      real<lower = 0> rho_S;
    
      // model FDiv
      real alpha_FDiv;
      real<lower = 0> sigma_FDiv;
      real beta_S;
      real beta_forest_FDiv;
      real beta_height_FDiv;
      real beta_crop_FDiv;
      real beta_duration_FDiv;
      real beta_fruit_FDiv;
      real beta_urban_FDiv;
      vector[N_ind] z_ind_FDiv;
      real<lower = 0> eta_FDiv;
      real<lower = 0> rho_FDiv;
    }
    
    transformed parameters {
      // pars S
      vector[N_ind] ind_S;
      matrix[N_ind, N_ind] L_K_S;
      L_K_S = cov_GPL2(dist_mat, eta_S, rho_S, 0.001);
      ind_S = L_K_S * z_ind_S;
    
      // pars FDiv
      vector[N_ind] ind_FDiv;
      matrix[N_ind, N_ind] L_K_FDiv;
      L_K_FDiv = cov_GPL2(dist_mat, eta_FDiv, rho_FDiv, 0.001);
      ind_FDiv = L_K_FDiv * z_ind_FDiv;
    }
    
    model {
    
      // prior S
      alpha_S ~ normal(0, 1);
      sigma_S ~ exponential(1);
      z_ind_S ~ normal(0, 0.25);
      eta_S ~ exponential(4);
      rho_S ~ exponential(1);
      beta_forest_S ~ normal(0, 0.25);
      beta_height_S ~ normal(0, 0.25);
      beta_crop_S ~ normal(0, 0.25);
      beta_duration_S ~ normal(0, 0.25);
      beta_fruit_S ~ normal(0, 0.25);
      beta_urban_S ~ normal(0, 0.25);
    
      // prior FDiv
      alpha_FDiv ~ normal(0, 1);
      sigma_FDiv ~ exponential(1);
      z_ind_FDiv ~ normal(0, 0.25);
      eta_FDiv ~ exponential(4);
      rho_FDiv ~ exponential(1);
      beta_forest_FDiv ~ normal(0, 0.25);
      beta_height_FDiv ~ normal(0, 0.25);
      beta_crop_FDiv ~ normal(0, 0.25);
      beta_duration_FDiv ~ normal(0, 0.25);
      beta_fruit_FDiv ~ normal(0, 0.25);
      beta_urban_FDiv ~ normal(0, 0.25);
      beta_S ~ normal(0, 1);
    
      // model S
      s_std ~ student_t(7, 
                        alpha_S + 
                        beta_forest_S * tree +
                        beta_height_S * height +
                        beta_crop_S * crop +
                        beta_duration_S * duration +
                        beta_fruit_S * fruit +
                        beta_urban_S * infrastructure +
                        ind_S[ind_id], sigma_S);
    
      // model FDiv
      fdiv ~ student_t(7, 
                        alpha_FDiv + 
                        beta_S * s_std +
                        beta_forest_FDiv * tree +
                        beta_height_FDiv * height +
                        beta_crop_FDiv * crop +
                        beta_duration_FDiv * duration +
                        beta_fruit_FDiv * fruit +
                        beta_urban_FDiv * infrastructure +
                        ind_FDiv[ind_id], sigma_FDiv);
    }
    
    generated quantities {
      array[N] real ppcheck_S;
      array[N] real ppcheck_FDiv;
      vector[N] mu_S;
      vector[N] mu_FDiv;
    
      mu_S = alpha_S + 
             beta_forest_S * tree +
             beta_height_S * height +
             beta_crop_S * crop +
             beta_duration_S * duration +
             beta_fruit_S * fruit +
             beta_urban_S * infrastructure +
             ind_S[ind_id];
    
      ppcheck_S = student_t_rng(7, mu_S, sigma_S);
    
        
      mu_FDiv = alpha_FDiv + 
                beta_S * s_std +
                beta_forest_FDiv * tree +
                beta_height_FDiv * height +
                beta_crop_FDiv * crop +
                beta_duration_FDiv * duration +
                beta_fruit_FDiv * fruit +
                beta_urban_FDiv * infrastructure +
                ind_FDiv[ind_id];
    
      ppcheck_FDiv = student_t_rng(7, mu_FDiv, sigma_FDiv);
    
    }
    ')

file <- paste0(getwd(), '/model_FDiv.stan')
fit_FDiv <- cmdstan_model(file, compile = T)

indx_na_vars <- 
  unlist(lapply(dat2, function(x) sum(is.na(x)) == 0), 
         use.names = F)

mod_FDiv <- 
  fit_FDiv$sample(
    data = dat2[indx_na_vars], 
    iter_warmup = 500, 
    iter_sampling = 4e3,
    chains = 3, 
    parallel_chains = 3,
    thin = 3,
    seed = 5
  )

sum_FDiv <- mod_FDiv$summary()

mod_diagnostics(mod_FDiv, sum_FDiv)

sum_FDiv |> print(n = 20)

indx_pars <- 
  c(grep('^beta', sum_FDiv$variable), 
    grep('^sigma', sum_FDiv$variable), 
    grep('^alpha', sum_FDiv$variable), 
    grep('^ind', sum_FDiv$variable))

length(indx_pars)

par(mfrow = c(4, 5), mar = c(4, 4, 1, 1))
for (i in 1:17) {
  trace_plot(mod_FDiv, sum_FDiv$variable[indx_pars[i]], 3)
}
par(mfrow = c(1, 1))

ppcheck_S <- mod_FDiv$draws('ppcheck_S', format = 'matrix')
ppcheck_FDiv <- mod_FDiv$draws('ppcheck_FDiv', format = 'matrix')

par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
plot(density(dat2$s_std), xlim = c(-4, 5), 
     ylim = c(0, 0.7), main = '', xlab = 'N species (z-scores)')
for (i in 1:500) lines(density(ppcheck_S[i, ]), lwd = 0.1)
lines(density(dat2$s_std), col = 'red', lwd = 1.5)

plot(density(dat2$fdiv), xlim = c(-4, 5), 
     ylim = c(0, 0.7), main = '', xlab = 'FDiv (z-scores)')
for (i in 1:500) lines(density(ppcheck_FDiv[i, ]), lwd = 0.1)
lines(density(dat2$fdiv), col = 'red', lwd = 1.5)
par(mfrow = c(1, 1))

par(mfrow = c(2, 2))
plot(density(apply(ppcheck_S, 1, mean)), xlab = 'AVG. S (estimated)', 
     main = '')
abline(v = mean(dat2$s_std), col = 'red')
plot(density(apply(ppcheck_S, 1, sd)), xlab = 'SD. S (estimated)', 
     main = '')
abline(v = sd(dat2$s_std), col = 'red')

plot(density(apply(ppcheck_FDiv, 1, mean)), xlab = 'AVG. FDiv (estimated)', 
     main = '')
abline(v = mean(dat2$fdiv), col = 'red')
plot(density(apply(ppcheck_FDiv, 1, sd)), xlab = 'SD. FDiv (estimated)', 
     main = '')
abline(v = sd(dat2$fdiv), col = 'red')
par(mfrow = c(1, 1))

R2_S <- bayesian_R2(mod_FDiv, dat2$s_std, 'mu_S')
quantile(R2_S)

R2_FDiv <- bayesian_R2(mod_FDiv, dat2$fdiv, 'mu_FDiv')
quantile(R2_FDiv)

plot(density(R2_S), col = 'lightblue', main = '', lwd = 3)

betas_S <- mod_FDiv$draws(sum_FDiv$variable[grep('^bet(.*)(_S)$', sum_FDiv$variable)], 
                          format = 'df')
betas_S <- betas_S[, grep('^beta', colnames(betas_S))]

apply(betas_S, 2, function(x) mean(x > 0))
apply(betas_S, 2, function(x) mean(x))
apply(betas_S, 2, function(x) sd(x)/sqrt(length(dat2$ind_id)))
apply(betas_S, 2, 
      function(x) {
        tibble(li = quantile(x, 0.025), 
               ls = quantile(x, 0.975))
      })

betas_FDiv <- mod_FDiv$draws(sum_FDiv$variable[grep('^bet(.*)(_FDiv)$', 
                                                    sum_FDiv$variable)], 
                             format = 'df')
betas_FDiv <- betas_FDiv[, grep('^beta', colnames(betas_FDiv))]

betas_FDiv$beta_S <- betas_S$beta_S

apply(betas_FDiv, 2, function(x) mean(x))
apply(betas_FDiv, 2, function(x) sd(x))
apply(betas_FDiv, 2, function(x) mean(x > 0))
apply(betas_FDiv, 2, 
      function(x) {
        tibble(li = quantile(x, 0.025), 
               ls = quantile(x, 0.975))
      })

# ========= conditional effects =========

sum_FDiv[indx_pars, ]

post_FDiv <- mod_FDiv$draws(sum_FDiv[indx_pars, ]$variable, format = 'df')

post_FDiv <- post_FDiv[, sum_FDiv[indx_pars, ]$variable]

apply(post_FDiv, 2, function(x) mean(x > 0))

apply(post_FDiv[, grep('^beta', colnames(post_FDiv))], 2, mean)
apply(post_FDiv[, grep('^beta', colnames(post_FDiv))], 2, sd)
apply(post_FDiv[, grep('^beta', colnames(post_FDiv))], 2, 
      function(x) mean(x > 0))

# conditional effects 

post_FDiv 

conditional_effect_FDiv <- 
  function(par, x_var, N = 500, FDiv = T) {
    
    x <- seq(min(x_var), max(x_var), length.out = N)
    
    ind_FDiv <- grep('^ind_FDiv', colnames(post_FDiv))
    ind_S <- grep('^ind_S', colnames(post_FDiv))
    
    est <- 
      lapply(x, FUN = 
               function(i) {
                 if (FDiv) {
                   mu_ <- 
                     post_FDiv$alpha_FDiv +
                     post_FDiv[[par]] * i +
                     apply(post_FDiv[, ind_FDiv], 1, mean)
                   
                   mu_pred <- 
                     rstudent(N, 7, mu_, post_FDiv$sigma_FDiv)
                 } else {
                   mu_ <- 
                     post_FDiv$alpha_S +
                     post_FDiv[[par]] * i +
                     apply(post_FDiv[, ind_S], 1, mean)
                   
                   mu_pred <- 
                     rstudent(N, 7, mu_, post_FDiv$sigma_S)
                 }
                 
                 tibble(mu = mean(mu_), 
                        li = quantile(mu_, 0.025), 
                        ls = quantile(mu_, 0.975), 
                        li_pred = quantile(mu_pred, 0.025), 
                        ls_pred = quantile(mu_pred, 0.975))
               })
    est <- do.call('rbind', est)
    est$x <- x
    est
    
  }

# ==== effect height ======

apply(betas_S, 2, function(x) mean(x > 0))
apply(betas_FDiv, 2, function(x) mean(x > 0))
apply(betas_FDiv, 2, function(x) mean(x))
apply(betas_S, 2, function(x) mean(x))

est_height_FDiv <- conditional_effect_FDiv('beta_height_FDiv', 
                                           dat2$height, FDiv = T)
est_height_S_FDiv <- conditional_effect_FDiv('beta_height_S', 
                                             dat2$height, FDiv = F)

plot_cond_heigth_FDiv <- 
  ggplot() +
  geom_ribbon(data = est_height_FDiv, 
              aes(x, mu, ymin = li_pred, ymax = ls_pred), 
              alpha = 0.5, fill = 'lightblue2') +
  geom_ribbon(data = est_height_FDiv, 
              aes(x, mu, ymin = li, ymax = ls), 
              alpha = 0.7, fill = 'lightblue') +
  geom_line(data = est_height_FDiv, 
            aes(x, mu), color = 'tan1', linewidth = 1, alpha = 0.7) +
  geom_point(data = tibble(x = dat2$height, 
                           y = dat2$fdiv), 
             aes(x, y), color = 'tan1', alpha = 0.5) +
  theme_classic() +
  labs(x = 'Tree height (z-scores)', 
       y = 'FDiv (z-scores)') +
  theme(legend.position = c(0.2, 0.9),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = 'Times New Roman', size = 10), 
        axis.line = element_line(linewidth = 0.25))

plot_cond_heigth_S_FDiv <- 
  ggplot() +
  geom_ribbon(data = est_height_S_FDiv, 
              aes(x, mu, ymin = li_pred, ymax = ls_pred), 
              alpha = 0.5, fill = 'lightblue2') +
  geom_ribbon(data = est_height_S_FDiv, 
              aes(x, mu, ymin = li, ymax = ls), 
              alpha = 0.7, fill = 'lightblue') +
  geom_line(data = est_height_S_FDiv, 
            aes(x, mu), color = 'tan1', linewidth = 1, alpha = 0.7) +
  geom_point(data = tibble(x = dat2$height, 
                           y = dat2$s_std), 
             aes(x, y), color = 'tan1', alpha = 0.5) +
  theme_classic() +
  labs(x = 'Tree height (z-scores)', 
       y = 'S (z-scores)') +
  theme(legend.position = c(0.2, 0.9),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = 'Times New Roman', size = 10), 
        axis.line = element_line(linewidth = 0.25))



causal_effect_FDiv <- 
  function(par_EQ1, par_EQ2, x_var, x_label, total_effect = T, negative = F) {
    ind_s <- grep('^ind_S', colnames(post_FDiv))
    ind_FDiv <- grep('^ind_FDiv', colnames(post_FDiv))
    mu_S <- mean(d$S)
    SD_S <- sd(d$S)
    
    FDiv <- dat_1$fdiv_obs_22[-which(is.na(dat_1$fdiv_obs_22))]
    
    mu_FDiv <- mean(FDiv)
    SD_FDiv <- sd(FDiv)
    
    estimated <- 
      lapply(c(min, max), FUN = 
               function(.fun) {
                 
                 est_S <- 
                   post_FDiv$alpha_S +
                   post_FDiv[[par_EQ1]] * .fun(x_var) +
                   apply(post_FDiv[, ind_s], 1, mean)
                 
                 if (total_effect) {
                   message('Calculating total effect')
                 } else {
                   message('Calculating conditional effect')
                   est_S <- min(dat2$s_std)
                 }
                 
                 est_FDiv <- 
                   post_FDiv$alpha_FDiv +
                   post_FDiv$beta_S * est_S +
                   post_FDiv[[par_EQ2]] * .fun(x_var) +
                   apply(post_FDiv[, ind_FDiv], 1, mean)
                 
                 est_FDiv <- mu_FDiv + est_FDiv * SD_FDiv
                 est_S <- mu_S + est_S * SD_S
                 
                 list(S = tibble(est = est_S, 
                                 Y = 'S', 
                                 X = x_label), 
                      FDiv = tibble(est = est_FDiv, 
                                    Y = 'FDiv', 
                                    X = x_label))
                 
               })
    
    names(estimated) <- c('Min', 'Max')
    
    estimated <- 
      lapply(names(estimated), FUN = 
               function(x) {
                 d <- estimated[[x]]
                 d$S$intervention <- x
                 d$FDiv$intervention <- x
                 d
               })
    names(estimated) <- c('Min', 'Max')
    
    estimated <- 
      lapply(estimated, FUN = 
               function(x) {
                 do.call('rbind', x)
               })
    estimated <- do.call('rbind', estimated)
    estimated <- split(estimated, estimated$Y)
    lapply(estimated, FUN =
             function(x) {
               if (negative) {
                 
                 x2 <- x[x$intervention == 'Max', ]
                 x2$est <-
                   x[x$intervention == 'Min', ]$est -
                   x[x$intervention == 'Max', ]$est
                 
                 x2$intervention <- 'Contrast'
                 z <- (x2$est * 100)/x[x$intervention == 'Min', ]$est
                 x$porcentaje <- NA
                 x2$porcentaje <- z
                 x2 <- rbind(x, x2)
                 x2$intervention <- as.factor(x2$intervention)
                 x2$intervention <- factor(x2$intervention, 
                                           levels = c('Min', 'Max', 
                                                      'Contrast'))
                 x2
                 
               } else {
                 x2 <- x[x$intervention == 'Max', ]
                 x2$est <-
                   x[x$intervention == 'Max', ]$est -
                   x[x$intervention == 'Min', ]$est
                 
                 x2$intervention <- 'Contrast'
                 z <- (x2$est * 100)/x[x$intervention == 'Max', ]$est
                 x$porcentaje <- NA
                 x2$porcentaje <- z
                 x2 <- rbind(x, x2)
                 x2$intervention <- as.factor(x2$intervention)
                 x2$intervention <- factor(x2$intervention, 
                                           levels = c('Min', 'Max', 
                                                      'Contrast'))
                 x2
               }
             })
  }

effect_height_FDiv <- causal_effect_FDiv('beta_height_S', 
                                         'beta_height_FDiv', 
                                         dat2$height, 'Height', T, negative = F)

effect_height_FDiv$FDiv$total_effect <- 'Total effect'

conditional_height_FDiv <- causal_effect_FDiv('beta_height_S', 
                                              'beta_height_FDiv', 
                                              dat$height, 'Height', F, 
                                              negative = F)[[1]]

conditional_height_FDiv$total_effect <- 'Conditional effect'

plot_CE_height_FDiv <- 
  rbind(effect_height_FDiv$FDiv, 
        conditional_height_FDiv) |> 
  group_by(intervention, total_effect) |> 
  transmute(mu = mean(est), 
            li = quantile(est, 0.025), 
            ls = quantile(est, 0.975)) |> 
  unique() |> 
  ggplot() +
  geom_errorbar(aes(intervention, ymin = li, ymax = ls, 
                    color = total_effect), width = 0, 
                position = position_dodge(width = 0.3)) +
  geom_point(aes(intervention, mu, color = total_effect), 
             position = position_dodge(width = 0.3)) +
  scale_color_manual(values = c('tan1', 'lightblue')) +
  scale_fill_manual(values = c('tan1', 'lightblue')) +
  lims(y = c(-1, 2)) +
  theme_classic() +
  labs(x = 'Intervention of tree height', 
       y = 'FDiv') +
  geom_hline(yintercept = 0, linetype = 3) +
  theme(legend.position = c(0.37, 0.17),
        legend.background = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = 'Times New Roman', size = 10), 
        axis.line = element_line(linewidth = 0.25)) 

# figure height FDiv

plot_heigh_FDiv <- 
  plot_grid(plot_cond_heigth_S_FDiv, 
            plot_cond_heigth_FDiv, 
            plot_CE_height_FDiv, ncol = 3, 
            labels = c('(i)', '(ii)', '(iii)'), 
            label_fontfamily = 'Times New Roman', 
            label_fontface = 'plain', 
            label_size = 10, 
            label_x = 0.175)

cont_height_FDiv <- na.omit(effect_height_FDiv$FDiv)

conditional_height_FDiv <- na.omit(conditional_height_FDiv)

diff_height_FDiv <- cont_height_FDiv$est - conditional_height_FDiv$est


tribble(~ value,               ~estimation, 
        'P(total effect > 0)', mean(cont_height_FDiv$est > 0), 
        'mu total effect', mean(cont_height_FDiv$est), 
        'sd total effect', sd(cont_height_FDiv$est), 
        'AVG change relative to maximum (total) (%)', mean(cont_height_FDiv$porcentaje), 
        'SD change relative to maximum (total) (%)', sd(cont_height_FDiv$porcentaje), 
        'P(conditional effect > 0)', mean(conditional_height_FDiv$est > 0), 
        'mu conditional effect', mean(conditional_height_FDiv$est), 
        'sd conditional effect', sd(conditional_height_FDiv$est), 
        'AVG change relative to maximum (conditional) (%)', mean(conditional_height_FDiv$porcentaje), 
        'SD change relative to maximum (conditional) (%)', sd(conditional_height_FDiv$porcentaje),
        'P(total > conditional)', mean(cont_height_FDiv$est > conditional_height_FDiv$est), 
        'AVG diff total-conditional (%)', mean(diff_height_FDiv)*100 / mean(cont_height_FDiv$est))

# ==== effect urbanization ======

apply(post_FDiv[, grep('^beta', colnames(post_FDiv))], 2, 
      function(x) mean(x < 0))
apply(post_FDiv[, grep('^beta', colnames(post_FDiv))], 2, 
      function(x) mean(x))

est_urban_FDiv <- conditional_effect_FDiv('beta_urban_FDiv', 
                                          dat2$infrastructure, 
                                          FDiv = T)
est_urban_S_FDiv <- conditional_effect_FDiv('beta_urban_S', 
                                            dat2$infrastructure, 
                                            FDiv = F)

plot_cond_urban_FDiv <- 
  ggplot() +
  geom_ribbon(data = est_urban_FDiv, 
              aes(x, mu, ymin = li_pred, ymax = ls_pred), 
              alpha = 0.5, fill = 'lightblue2') +
  geom_ribbon(data = est_urban_FDiv, 
              aes(x, mu, ymin = li, ymax = ls), 
              alpha = 0.7, fill = 'lightblue') +
  geom_line(data = est_urban_FDiv, 
            aes(x, mu), color = 'tan1', linewidth = 1, alpha = 0.7) +
  geom_point(data = tibble(x = dat2$infrastructure, 
                           y = dat2$fdiv), 
             aes(x, y), color = 'tan1', alpha = 0.5) +
  theme_classic() +
  labs(x = 'Infrastructure cover (z-scores)', 
       y = 'FDiv (z-scores)') +
  theme(legend.position = c(0.2, 0.9),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = 'Times New Roman', size = 10), 
        axis.line = element_line(linewidth = 0.25))

plot_cond_urban_S_FDiv <- 
  ggplot() +
  geom_ribbon(data = est_urban_S_FDiv, 
              aes(x, mu, ymin = li_pred, ymax = ls_pred), 
              alpha = 0.5, fill = 'lightblue2') +
  geom_ribbon(data = est_urban_S_FDiv, 
              aes(x, mu, ymin = li, ymax = ls), 
              alpha = 0.7, fill = 'lightblue') +
  geom_line(data = est_urban_S_FDiv, 
            aes(x, mu), color = 'tan1', linewidth = 1, alpha = 0.7) +
  geom_point(data = tibble(x = dat2$infrastructure, 
                           y = dat2$s_std), 
             aes(x, y), color = 'tan1', alpha = 0.5) +
  theme_classic() +
  labs(x = 'Infrastructure cover (z-scores)', 
       y = 'S (z-scores)') +
  theme(legend.position = c(0.2, 0.9),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = 'Times New Roman', size = 10), 
        axis.line = element_line(linewidth = 0.25))

effect_urban_FDiv <- causal_effect_FDiv('beta_urban_S', 
                                        'beta_urban_FDiv', 
                                         dat2$infrastructure, 
                                         'Urban', T, negative = F)

effect_urban_FDiv$FDiv$total_effect <- 'Total effect'

conditional_urban_FDiv <- causal_effect_FDiv('beta_urban_S', 
                                             'beta_urban_FDiv', 
                                             dat2$infrastructure, 
                                             'Urban', F, 
                                             negative = F)[[1]]

conditional_urban_FDiv$total_effect <- 'Conditional effect'

plot_CE_urban_FDiv <- 
  rbind(effect_urban_FDiv$FDiv, 
        conditional_urban_FDiv) |> 
  group_by(intervention, total_effect) |> 
  transmute(mu = mean(est), 
            li = quantile(est, 0.025), 
            ls = quantile(est, 0.975)) |> 
  unique() |> 
  ggplot() +
  geom_errorbar(aes(intervention, ymin = li, ymax = ls, 
                    color = total_effect), width = 0, 
                position = position_dodge(width = 0.3)) +
  geom_point(aes(intervention, mu, color = total_effect), 
             position = position_dodge(width = 0.3)) +
  scale_color_manual(values = c('tan1', 'lightblue')) +
  scale_fill_manual(values = c('tan1', 'lightblue')) +
  theme_classic() +
  labs(x = 'Intervention of\n infrastructure cover', 
       y = 'FDiv') +
  lims(y = c(-0.7, 2)) +
  geom_hline(yintercept = 0, linetype = 3) +
  theme(legend.position = 'none',
        legend.title = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = 'Times New Roman', size = 10), 
        axis.line = element_line(linewidth = 0.25)) 

# figure urban FDiv

plot_urban_FDiv <- 
  plot_grid(plot_cond_urban_S_FDiv, 
            plot_cond_urban_FDiv, 
            plot_CE_urban_FDiv, ncol = 3, 
            labels = c('(i)', '(ii)', '(iii)'), 
            label_fontfamily = 'Times New Roman', 
            label_fontface = 'plain', 
            label_size = 10, 
            label_x = 0.175)

cont_urban_FDiv <- na.omit(effect_urban_FDiv$FDiv)

conditional_urban_FDiv <- na.omit(conditional_urban_FDiv)

diff_urban_FDiv <- cont_urban_FDiv$est - conditional_urban_FDiv$est

tribble(~ value,               ~estimation, 
        'P(total effect < 0)', mean(cont_urban_FDiv$est < 0), 
        'mu total effect', mean(cont_urban_FDiv$est), 
        'sd total effect', sd(cont_urban_FDiv$est), 
        'AVG change relative to maximum (total) (%)', mean(cont_urban_FDiv$porcentaje), 
        'SD change relative to maximum (total) (%)', sd(cont_urban_FDiv$porcentaje), 
        'P(conditional effect < 0)', mean(conditional_urban_FDiv$est < 0), 
        'mu conditional effect', mean(conditional_urban_FDiv$est), 
        'sd conditional effect', sd(conditional_urban_FDiv$est), 
        'AVG change relative to maximum (conditional) (%)', mean(conditional_urban_FDiv$porcentaje), 
        'SD change relative to maximum (conditional) (%)', sd(conditional_urban_FDiv$porcentaje),
        'P(total < conditional)', mean(cont_urban_FDiv$est < conditional_urban_FDiv$est), 
        'AVG diff total-conditional (%)', mean(diff_urban_FDiv)*100 / mean(cont_urban_FDiv$est))

# ==== effect crop size =====

apply(post_FDiv[, grep('^beta', colnames(post_FDiv))], 2, 
      function(x) mean(x > 0))

apply(post_FDiv[, grep('^beta', colnames(post_FDiv))], 2, 
      function(x) mean(x))

((0.190306201 - 0.143893985)/0.143893985)*100

est_crop_FDiv <- conditional_effect_FDiv('beta_crop_FDiv', 
                                           dat2$crop, 
                                           FDiv = T)
est_crop_S_FDiv <- conditional_effect_FDiv('beta_crop_S', 
                                             dat2$crop, 
                                             FDiv = F)

plot_cond_crop_FDiv <- 
  ggplot() +
  geom_ribbon(data = est_crop_FDiv, 
              aes(x, mu, ymin = li_pred, ymax = ls_pred), 
              alpha = 0.5, fill = 'lightblue2') +
  geom_ribbon(data = est_crop_FDiv, 
              aes(x, mu, ymin = li, ymax = ls), 
              alpha = 0.7, fill = 'lightblue') +
  geom_line(data = est_crop_FDiv, 
            aes(x, mu), color = 'tan1', linewidth = 1, alpha = 0.7) +
  geom_point(data = tibble(x = dat2$crop, 
                           y = dat2$fdiv), 
             aes(x, y), color = 'tan1', alpha = 0.5) +
  theme_classic() +
  labs(x = 'Crop size (z-scores)', 
       y = 'FDiv (z-scores)') +
  theme(legend.position = c(0.2, 0.9),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = 'Times New Roman', size = 10), 
        axis.line = element_line(linewidth = 0.25))

plot_cond_crop_S_FDiv <- 
  ggplot() +
  geom_ribbon(data = est_crop_S_FDiv, 
              aes(x, mu, ymin = li_pred, ymax = ls_pred), 
              alpha = 0.5, fill = 'lightblue2') +
  geom_ribbon(data = est_crop_S_FDiv, 
              aes(x, mu, ymin = li, ymax = ls), 
              alpha = 0.7, fill = 'lightblue') +
  geom_line(data = est_crop_S_FDiv, 
            aes(x, mu), color = 'tan1', linewidth = 1, alpha = 0.7) +
  geom_point(data = tibble(x = dat2$crop, 
                           y = dat2$s_std), 
             aes(x, y), color = 'tan1', alpha = 0.5) +
  theme_classic() +
  labs(x = 'Crop size (z-scores)', 
       y = 'S (z-scores)') +
  theme(legend.position = c(0.2, 0.9),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = 'Times New Roman', size = 10), 
        axis.line = element_line(linewidth = 0.25))

effect_crop_FDiv <- causal_effect_FDiv('beta_crop_S', 
                                    'beta_crop_FDiv', 
                                    dat2$crop, 'crop', T, 
                                    negative = F)

effect_crop_FDiv$FDiv$total_effect <- 'Total effect'

conditional_crop_FDiv <- causal_effect_FDiv('beta_crop_S', 
                                              'beta_crop_FDiv', 
                                              dat2$crop, 'crop', F, 
                                            negative = F)[[1]]

conditional_crop_FDiv$total_effect <- 'Conditional effect'

plot_CE_crop_FDiv <- 
  rbind(effect_crop_FDiv$FDiv, 
        conditional_crop_FDiv) |> 
  group_by(intervention, total_effect) |> 
  transmute(mu = mean(est), 
            li = quantile(est, 0.025), 
            ls = quantile(est, 0.975)) |> 
  unique() |> 
  ggplot() +
  geom_errorbar(aes(intervention, ymin = li, ymax = ls, 
                    color = total_effect), width = 0, 
                position = position_dodge(width = 0.3)) +
  geom_point(aes(intervention, mu, color = total_effect), 
             position = position_dodge(width = 0.3)) +
  scale_color_manual(values = c('tan1', 'lightblue')) +
  scale_fill_manual(values = c('tan1', 'lightblue')) +
  lims(y = c(-0.7, 2.5)) +
  theme_classic() +
  labs(x = 'Intervention of\n crop size', 
       y = 'FDiv') +
  geom_hline(yintercept = 0, linetype = 3) +
  theme(legend.position = 'none',
        legend.title = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = 'Times New Roman', size = 10), 
        axis.line = element_line(linewidth = 0.25)) 

# figure crop FDiv

plot_crop_FDiv <- 
  plot_grid(plot_cond_crop_S_FDiv, 
            plot_cond_crop_FDiv, 
            plot_CE_crop_FDiv, ncol = 3, 
            labels = c('(i)', '(ii)', '(iii)'), 
            label_fontfamily = 'Times New Roman', 
            label_fontface = 'plain', 
            label_size = 10, 
            label_x = 0.175)

cont_crop_FDiv <- na.omit(effect_crop_FDiv$FDiv)

conditional_crop_FDiv <- na.omit(conditional_crop_FDiv)

diff_crop_FDiv <- cont_crop_FDiv$est - conditional_crop_FDiv$est

tribble(~ value,               ~estimation, 
        'P(total effect > 0)', mean(cont_crop_FDiv$est > 0), 
        'mu total effect', mean(cont_crop_FDiv$est), 
        'sd total effect', sd(cont_crop_FDiv$est), 
        'AVG change relative to maximum (total) (%)', mean(cont_crop_FDiv$porcentaje), 
        'SD change relative to maximum (total) (%)', sd(cont_crop_FDiv$porcentaje), 
        'P(conditional effect > 0)', mean(conditional_crop_FDiv$est > 0), 
        'mu conditional effect', mean(conditional_crop_FDiv$est), 
        'sd conditional effect', sd(conditional_crop_FDiv$est), 
        'AVG change relative to maximum (conditional) (%)', mean(conditional_crop_FDiv$porcentaje), 
        'SD change relative to maximum (conditional) (%)', sd(conditional_crop_FDiv$porcentaje),
        'P(total > conditional)', mean(cont_crop_FDiv$est > conditional_crop_FDiv$est), 
        'AVG diff total-conditional (%)', mean(diff_crop_FDiv)*100 / mean(cont_crop_FDiv$est))


# ===== plots FDiv mod ===

plot_grid(plot_heigh_FDiv, 
          plot_urban_FDiv, 
          plot_crop_FDiv,
          ncol = 1, 
          labels = c('(a)', '(b)', '(c)'), 
          label_fontfamily = 'Times New Roman', 
          label_size = 12, 
          label_fontface = 'plain', 
          label_x = -0.007)

ggsave('FDiv_mod.jpg', width = 20, height = 18, units = 'cm', dpi = 1000)



sapply(c('_S', '_FDiv'), FUN = 
        function(x) {
          
          v <- 
            apply(post_FDiv[, grep('^beta', colnames(post_FDiv))], 2, 
                  function(x) mean(x))
          
          v <- v[grep(x, names(v))]
          
          t <- 
            c(grep('height', names(v)), 
              grep('crop', names(v)))
          
          z <- 
            sapply(v[t], FUN = 
                     function(i) {
                       i <- abs(i)
                       ((abs(v[grep('urban', names(v))]) - i) / i) * 100
                     })
          
          tibble(mu = mean(z), 
                 sd = sd(z))
          
        })


# ======== model FReg =====

dat3 <- dat_1

indx_na_FReg <- dat3$ind_id[which(is.na(dat3$freg_obs_22))]

dist_inds_FReg <- dist_inds[-indx_na_FReg, -indx_na_FReg]

dat3 <- dat3[-indx_na_FReg ,]

dat3$ind_id <- 1:nrow(dat3)

colnames(dist_inds_FReg) <- paste0('Ind', dat3$ind_id)
rownames(dist_inds_FReg) <- paste0('Ind', dat3$ind_id)

dat3 <- lapply(dat3, function(x) x)

dat3$N <- length(dat3$ind_id)
dat3$N_ind <- max(dat3$ind_id)
dat3$dist_mat <- dist_inds_FReg
names(dat3)[2:4] <- c('fric', 'fdiv', 'freg')
names(dat3)

summary(dat3$freg)
plot(density(dat3$freg))
dat3$freg <- as.vector(scale(dat3$freg))

cat(file = 'model_FReg.stan', 
    '
    functions {
      matrix cov_GPL2(matrix x, 
                      real eta,
                      real rho,
                      real delta) {
                      
                      int N = dims(x)[1];
                      matrix[N, N] K;
                      matrix[N, N] L_K;
    
                      for (i in 1:(N-1)) {
                        K[i, i] = eta + delta;
                        for (j in (i+1):N) {
                          K[i, j] = eta * exp(-rho * square(x[i, j]));
                          K[j, i] = K[i, j];
                        }
                      }
                      
                      K[N, N] = eta + delta;
                      L_K = cholesky_decompose(K);
                      return L_K;
                      }
    }
    
    data {
      int N;
      int N_ind;
      matrix[N_ind, N_ind] dist_mat;
      array[N] int ind_id;
      vector[N] freg;
      vector[N] s_std;
      vector[N] height;
      vector[N] crop;
      vector[N] duration;
      vector[N] fruit;
      vector[N] tree;
      vector[N] infrastructure;
    }
    
    parameters {
    
      // model S
      real alpha_S;
      real<lower = 0> sigma_S;
      real beta_forest_S;
      real beta_height_S;
      real beta_crop_S;
      real beta_duration_S;
      real beta_fruit_S;
      real beta_urban_S;
      vector[N_ind] z_ind_S;
      real<lower = 0> eta_S;
      real<lower = 0> rho_S;
    
      // model FReg
      real alpha_FReg;
      real<lower = 0> sigma_FReg;
      real beta_S;
      real beta_forest_FReg;
      real beta_height_FReg;
      real beta_crop_FReg;
      real beta_duration_FReg;
      real beta_fruit_FReg;
      real beta_urban_FReg;
      vector[N_ind] z_ind_FReg;
      real<lower = 0> eta_FReg;
      real<lower = 0> rho_FReg;
    }
    
    transformed parameters {
      // pars S
      vector[N_ind] ind_S;
      matrix[N_ind, N_ind] L_K_S;
      L_K_S = cov_GPL2(dist_mat, eta_S, rho_S, 0.001);
      ind_S = L_K_S * z_ind_S;
    
      // pars FReg
      vector[N_ind] ind_FReg;
      matrix[N_ind, N_ind] L_K_FReg;
      L_K_FReg = cov_GPL2(dist_mat, eta_FReg, rho_FReg, 0.001);
      ind_FReg = L_K_FReg * z_ind_FReg;
    }
    
    model {
    
      // prior S
      alpha_S ~ normal(0, 1);
      sigma_S ~ exponential(1);
      z_ind_S ~ normal(0, 0.25);
      eta_S ~ exponential(4);
      rho_S ~ exponential(1);
      beta_forest_S ~ normal(0, 0.25);
      beta_height_S ~ normal(0, 0.25);
      beta_crop_S ~ normal(0, 0.25);
      beta_duration_S ~ normal(0, 0.25);
      beta_fruit_S ~ normal(0, 0.25);
      beta_urban_S ~ normal(0, 0.25);
    
      // prior FReg
      alpha_FReg ~ normal(0, 1);
      sigma_FReg ~ exponential(1);
      z_ind_FReg ~ normal(0, 0.25);
      eta_FReg ~ exponential(4);
      rho_FReg ~ exponential(1);
      beta_forest_FReg ~ normal(0, 0.25);
      beta_height_FReg ~ normal(0, 0.25);
      beta_crop_FReg ~ normal(0, 0.25);
      beta_duration_FReg ~ normal(0, 0.25);
      beta_fruit_FReg ~ normal(0, 0.25);
      beta_urban_FReg ~ normal(0, 0.25);
      beta_S ~ normal(0, 1);
    
      // model S
      s_std ~ student_t(7, 
                        alpha_S + 
                        beta_forest_S * tree +
                        beta_height_S * height +
                        beta_crop_S * crop +
                        beta_duration_S * duration +
                        beta_fruit_S * fruit +
                        beta_urban_S * infrastructure +
                        ind_S[ind_id], sigma_S);
    
      // model FReg
      freg ~ student_t(7, 
                        alpha_FReg + 
                        beta_S * s_std +
                        beta_forest_FReg * tree +
                        beta_height_FReg * height +
                        beta_crop_FReg * crop +
                        beta_duration_FReg * duration +
                        beta_fruit_FReg * fruit +
                        beta_urban_FReg * infrastructure +
                        ind_FReg[ind_id], sigma_FReg);
    }
    
    generated quantities {
      array[N] real ppcheck_S;
      array[N] real ppcheck_FReg;
      vector[N] mu_S;
      vector[N] mu_FReg;
    
      mu_S = alpha_S + 
             beta_forest_S * tree +
             beta_height_S * height +
             beta_crop_S * crop +
             beta_duration_S * duration +
             beta_fruit_S * fruit +
             beta_urban_S * infrastructure +
             ind_S[ind_id];
    
      ppcheck_S = student_t_rng(7, mu_S, sigma_S);
    
        
      mu_FReg = alpha_FReg + 
                beta_S * s_std +
                beta_forest_FReg * tree +
                beta_height_FReg * height +
                beta_crop_FReg * crop +
                beta_duration_FReg * duration +
                beta_fruit_FReg * fruit +
                beta_urban_FReg * infrastructure +
                ind_FReg[ind_id];
    
      ppcheck_FReg = student_t_rng(7, mu_FReg, sigma_FReg);
    
    }
    ')

file <- paste0(getwd(), '/model_FReg.stan')
fit_FReg <- cmdstan_model(file, compile = T)

indx_na_vars <- 
  unlist(lapply(dat3, function(x) sum(is.na(x)) == 0), 
         use.names = F)

mod_FReg <- 
  fit_FReg$sample(
    data = dat3[indx_na_vars], 
    iter_warmup = 500, 
    iter_sampling = 4e3,
    chains = 3, 
    parallel_chains = 3,
    thin = 3,
    seed = 5
  )

sum_FReg <- mod_FReg$summary()

mod_diagnostics(mod_FReg, sum_FReg)

sum_FReg |> print(n = 20)

indx_pars <- 
  c(grep('^beta', sum_FReg$variable), 
    grep('^sigma', sum_FReg$variable), 
    grep('^alpha', sum_FReg$variable), 
    grep('^ind', sum_FReg$variable))

length(indx_pars)

par(mfrow = c(4, 5), mar = c(4, 4, 1, 1))
for (i in 1:17) {
  trace_plot(mod_FReg, sum_FReg$variable[indx_pars[i]], 3)
}
par(mfrow = c(1, 1))

ppcheck_S <- mod_FReg$draws('ppcheck_S', format = 'matrix')
ppcheck_FReg <- mod_FReg$draws('ppcheck_FReg', format = 'matrix')

par(mfrow = c(1, 2), mar = c(4, 4, 1, 1))
plot(density(dat3$s_std), xlim = c(-4, 5), 
     ylim = c(0, 0.7), main = '', xlab = 'N species (z-scores)')
for (i in 1:500) lines(density(ppcheck_S[i, ]), lwd = 0.1)
lines(density(dat3$s_std), col = 'red', lwd = 1.5)

plot(density(dat3$freg), xlim = c(-5, 5), 
     ylim = c(0, 0.7), main = '', xlab = 'FReg (z-scores)')
for (i in 1:500) lines(density(ppcheck_FReg[i, ]), lwd = 0.1)
lines(density(dat3$freg), col = 'red', lwd = 1.5)
par(mfrow = c(1, 1))

par(mfrow = c(2, 2))
plot(density(apply(ppcheck_S, 1, mean)), xlab = 'AVG. S (estimated)', 
     main = '')
abline(v = mean(dat3$s_std), col = 'red')
plot(density(apply(ppcheck_S, 1, sd)), xlab = 'SD. S (estimated)', 
     main = '')
abline(v = sd(dat3$s_std), col = 'red')

plot(density(apply(ppcheck_FReg, 1, mean)), xlab = 'AVG. FReg (estimated)', 
     main = '')
abline(v = mean(dat3$freg), col = 'red')
plot(density(apply(ppcheck_FReg, 1, sd)), xlab = 'SD. FReg (estimated)', 
     main = '')
abline(v = sd(dat3$freg), col = 'red')
par(mfrow = c(1, 1))

R2_S <- bayesian_R2(mod_FReg, dat3$s_std, 'mu_S')
quantile(R2_S)

R2_FReg <- bayesian_R2(mod_FReg, dat3$freg, 'mu_FReg')
quantile(R2_FReg)

plot(density(R2_S), col = 'lightblue', main = '', lwd = 3)

betas_S <- mod_FReg$draws(sum_FReg$variable[grep('^bet(.*)(_S)$', sum_FReg$variable)], 
                          format = 'df')
betas_S <- betas_S[, grep('^beta', colnames(betas_S))]

apply(betas_S, 2, function(x) mean(x > 0))
apply(betas_S, 2, function(x) mean(x))
apply(betas_S, 2, function(x) sd(x))
apply(betas_S, 2, 
      function(x) {
        tibble(li = quantile(x, 0.025), 
               ls = quantile(x, 0.975))
      })

betas_FReg <- mod_FReg$draws(sum_FReg$variable[grep('^bet(.*)(_FReg)$', 
                                                    sum_FReg$variable)], 
                             format = 'df')
betas_FReg <- betas_FReg[, grep('^beta', colnames(betas_FReg))]

betas_FReg$beta_S <- betas_S$beta_S

apply(betas_FReg, 2, function(x) mean(x > 0))
apply(betas_FReg, 2, function(x) mean(x))
apply(betas_FReg, 2, function(x) sd(x))
apply(betas_FReg, 2, 
      function(x) {
        tibble(li = quantile(x, 0.025), 
               ls = quantile(x, 0.975))
      })

# ========= conditional effects =========

sum_FReg[indx_pars, ]

post_FReg <- mod_FReg$draws(sum_FReg[indx_pars, ]$variable, format = 'df')

post_FReg <- post_FReg[, sum_FReg[indx_pars, ]$variable]

apply(post_FReg, 2, function(x) mean(x > 0))

apply(post_FReg[, grep('^beta', colnames(post_FReg))], 2, mean)
apply(post_FReg[, grep('^beta', colnames(post_FReg))], 2, sd)
apply(post_FReg[, grep('^beta', colnames(post_FReg))], 2, 
      function(x) mean(x > 0))

# conditional effects 

post_FReg 

conditional_effect_FReg <- 
  function(par, x_var, N = 500, FReg = T) {
    
    x <- seq(min(x_var), max(x_var), length.out = N)
    
    ind_FReg <- grep('^ind_FReg', colnames(post_FReg))
    ind_S <- grep('^ind_S', colnames(post_FReg))
    
    est <- 
      lapply(x, FUN = 
               function(i) {
                 if (FReg) {
                   mu_ <- 
                     post_FReg$alpha_FReg +
                     post_FReg[[par]] * i +
                     apply(post_FReg[, ind_FReg], 1, mean)
                   
                   mu_pred <- 
                     rstudent(N, 7, mu_, post_FReg$sigma_FReg)
                 } else {
                   mu_ <- 
                     post_FReg$alpha_S +
                     post_FReg[[par]] * i +
                     apply(post_FReg[, ind_S], 1, mean)
                   
                   mu_pred <- 
                     rstudent(N, 7, mu_, post_FReg$sigma_S)
                 }
                 
                 tibble(mu = mean(mu_), 
                        li = quantile(mu_, 0.025), 
                        ls = quantile(mu_, 0.975), 
                        li_pred = quantile(mu_pred, 0.025), 
                        ls_pred = quantile(mu_pred, 0.975))
               })
    est <- do.call('rbind', est)
    est$x <- x
    est
    
  }

# ==== effect crop ======

apply(betas_S, 2, function(x) mean(x > 0))
apply(betas_FReg, 2, function(x) mean(x > 0))
apply(betas_FReg, 2, function(x) mean(x))
apply(betas_S, 2, function(x) mean(x))

est_crop_FReg <- conditional_effect_FReg('beta_crop_FReg', 
                                          dat3$crop, FReg = T)
est_crop_S_FReg <- conditional_effect_FReg('beta_crop_S', 
                                            dat3$crop, FReg = F)

plot_cond_crop_FReg <- 
  ggplot() +
  geom_ribbon(data = est_crop_FReg, 
              aes(x, mu, ymin = li_pred, ymax = ls_pred), 
              alpha = 0.5, fill = 'lightblue2') +
  geom_ribbon(data = est_crop_FReg, 
              aes(x, mu, ymin = li, ymax = ls), 
              alpha = 0.7, fill = 'lightblue') +
  geom_line(data = est_crop_FReg, 
            aes(x, mu), color = 'tan1', linewidth = 1, alpha = 0.7) +
  geom_point(data = tibble(x = dat3$crop, 
                           y = dat3$freg), 
             aes(x, y), color = 'tan1', alpha = 0.5) +
  theme_classic() +
  labs(x = 'Crop size (z-scores)', 
       y = 'FReg (z-scores)') +
  theme(legend.position = c(0.2, 0.9),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = 'Times New Roman', size = 10), 
        axis.line = element_line(linewidth = 0.25))

plot_cond_crop_S_FReg <- 
  ggplot() +
  geom_ribbon(data = est_crop_S_FReg, 
              aes(x, mu, ymin = li_pred, ymax = ls_pred), 
              alpha = 0.5, fill = 'lightblue2') +
  geom_ribbon(data = est_crop_S_FReg, 
              aes(x, mu, ymin = li, ymax = ls), 
              alpha = 0.7, fill = 'lightblue') +
  geom_line(data = est_crop_S_FReg, 
            aes(x, mu), color = 'tan1', linewidth = 1, alpha = 0.7) +
  geom_point(data = tibble(x = dat3$crop, 
                           y = dat3$s_std), 
             aes(x, y), color = 'tan1', alpha = 0.5) +
  theme_classic() +
  labs(x = 'Crop size (z-scores)', 
       y = 'S (z-scores)') +
  theme(legend.position = c(0.2, 0.9),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = 'Times New Roman', size = 10), 
        axis.line = element_line(linewidth = 0.25))



causal_effect_FReg <- 
  function(par_EQ1, par_EQ2, x_var, x_label, total_effect = T) {
    ind_s <- grep('^ind_S', colnames(post_FReg))
    ind_FReg <- grep('^ind_FReg', colnames(post_FReg))
    mu_S <- mean(d$S)
    SD_S <- sd(d$S)
    
    FReg <- dat_1$freg_obs_22[-which(is.na(dat_1$freg_obs_22))]
    
    mu_FReg <- mean(FReg)
    SD_FReg <- sd(FReg)
    
    estimated <- 
      lapply(c(min, max), FUN = 
               function(.fun) {
                 
                 est_S <- 
                   post_FReg$alpha_S +
                   post_FReg[[par_EQ1]] * .fun(x_var) +
                   apply(post_FReg[, ind_s], 1, mean)
                 
                 if (total_effect) {
                   message('Calculating total effect')
                 } else {
                   message('Calculating conditional effect')
                   est_S <- min(dat3$s_std)
                 }
                 
                 est_FReg <- 
                   post_FReg$alpha_FReg +
                   post_FReg$beta_S * est_S +
                   post_FReg[[par_EQ2]] * .fun(x_var) +
                   apply(post_FReg[, ind_FReg], 1, mean)
                 
                 est_FReg <- mu_FReg + est_FReg * SD_FReg
                 est_S <- mu_S + est_S * SD_S
                 
                 list(S = tibble(est = est_S, 
                                 Y = 'S', 
                                 X = x_label), 
                      FReg = tibble(est = est_FReg, 
                                    Y = 'FReg', 
                                    X = x_label))
                 
               })
    
    names(estimated) <- c('Min', 'Max')
    
    estimated <- 
      lapply(names(estimated), FUN = 
               function(x) {
                 d <- estimated[[x]]
                 d$S$intervention <- x
                 d$FReg$intervention <- x
                 d
               })
    names(estimated) <- c('Min', 'Max')
    
    estimated <- 
      lapply(estimated, FUN = 
               function(x) {
                 do.call('rbind', x)
               })
    estimated <- do.call('rbind', estimated)
    estimated <- split(estimated, estimated$Y)
    lapply(estimated, FUN =
             function(x) {
               x2 <- x[x$intervention == 'Max', ]
               x2$est <-
                 x[x$intervention == 'Max', ]$est -
                 x[x$intervention == 'Min', ]$est
               
               x2$intervention <- 'Contrast'
               z <- (x2$est * 100)/x[x$intervention == 'Max', ]$est
               x$porcentaje <- NA
               x2$porcentaje <- z
               x2 <- rbind(x, x2)
               x2$intervention <- as.factor(x2$intervention)
               x2$intervention <- factor(x2$intervention, 
                                         levels = c('Min', 'Max', 
                                                    'Contrast'))
               x2
             })
  }

effect_crop_FReg <- causal_effect_FReg('beta_crop_S', 
                                        'beta_crop_FReg', 
                                        dat3$crop, 'crop', T)

effect_crop_FReg$FReg$total_effect <- 'Total effect'

conditional_crop_FReg <- causal_effect_FReg('beta_crop_S', 
                                             'beta_crop_FReg', 
                                             dat3$crop, 'crop', F)[[1]]

conditional_crop_FReg$total_effect <- 'Conditional effect'

plot_CE_crop_FReg <- 
  rbind(effect_crop_FReg$FReg, 
        conditional_crop_FReg) |> 
  group_by(intervention, total_effect) |> 
  transmute(mu = mean(est), 
            li = quantile(est, 0.025), 
            ls = quantile(est, 0.975)) |> 
  unique() |> 
  ggplot() +
  geom_errorbar(aes(intervention, ymin = li, ymax = ls, 
                    color = total_effect), width = 0, 
                position = position_dodge(width = 0.3)) +
  geom_point(aes(intervention, mu, color = total_effect), 
             position = position_dodge(width = 0.3)) +
  scale_color_manual(values = c('tan1', 'lightblue')) +
  scale_fill_manual(values = c('tan1', 'lightblue')) +
  lims(y = c(-0.5, 1.5)) +
  theme_classic() +
  labs(x = 'Intervention of crop size', 
       y = 'FReg') +
  geom_hline(yintercept = 0, linetype = 3) +
  theme(legend.position = c(0.37, 0.875),
        legend.background = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = 'Times New Roman', size = 10), 
        axis.line = element_line(linewidth = 0.25)) 

# figure crop FReg

plot_crop_FReg <- 
  plot_grid(plot_cond_crop_S_FReg, 
            plot_cond_crop_FReg, 
            plot_CE_crop_FReg, ncol = 3, 
            labels = c('(a)', '(b)', '(c)'), 
            label_fontfamily = 'Times New Roman', 
            label_fontface = 'plain', 
            label_size = 10, 
            label_x = 0.175)

plot_crop_FReg

ggsave('FReg_mod.jpg', width = 20, height = 7, units = 'cm', dpi = 1000)

cont_crop_FReg <- na.omit(effect_crop_FReg$FReg)

conditional_crop_FReg <- na.omit(conditional_crop_FReg)

diff_crop_FReg <- conditional_crop_FReg$est - cont_crop_FReg$est

tribble(~ value,               ~estimation, 
        'P(total effect > 0)', mean(cont_crop_FReg$est > 0), 
        'mu total effect', mean(cont_crop_FReg$est), 
        'sd total effect', sd(cont_crop_FReg$est), 
        'AVG change relative to maximum (total) (%)', mean(cont_crop_FReg$porcentaje), 
        'SD change relative to maximum (total) (%)', sd(cont_crop_FReg$porcentaje), 
        'P(conditional effect > 0)', mean(conditional_crop_FReg$est > 0), 
        'mu conditional effect', mean(conditional_crop_FReg$est), 
        'sd conditional effect', sd(conditional_crop_FReg$est), 
        'AVG change relative to maximum (conditional) (%)', mean(conditional_crop_FReg$porcentaje), 
        'SD change relative to maximum (conditional) (%)', sd(conditional_crop_FReg$porcentaje),
        'P(total > conditional)', mean(cont_crop_FReg$est < conditional_crop_FReg$est), 
        'AVG diff total-conditional (%)', 100-(mean(diff_crop_FReg)*100 / mean(conditional_crop_FReg$est)))



#========== effect fruiting duration ========


est_duration_FReg <- conditional_effect_FReg('beta_duration_FReg', 
                                         dat3$duration, FReg = T)
est_duration_S_FReg <- conditional_effect_FReg('beta_duration_S', 
                                           dat3$duration, FReg = F)

plot_cond_duration_FReg <- 
  ggplot() +
  geom_ribbon(data = est_duration_FReg, 
              aes(x, mu, ymin = li_pred, ymax = ls_pred), 
              alpha = 0.5, fill = 'lightblue2') +
  geom_ribbon(data = est_duration_FReg, 
              aes(x, mu, ymin = li, ymax = ls), 
              alpha = 0.7, fill = 'lightblue') +
  geom_line(data = est_duration_FReg, 
            aes(x, mu), color = 'tan1', linewidth = 1, alpha = 0.7) +
  geom_point(data = tibble(x = dat3$duration, 
                           y = dat3$freg), 
             aes(x, y), color = 'tan1', alpha = 0.5) +
  theme_classic() +
  labs(x = 'Fruit duration (z-scores)', 
       y = 'FReg (z-scores)') +
  theme(legend.position = c(0.2, 0.9),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = 'Times New Roman', size = 10), 
        axis.line = element_line(linewidth = 0.25))

plot_cond_duration_S_FReg <- 
  ggplot() +
  geom_ribbon(data = est_duration_S_FReg, 
              aes(x, mu, ymin = li_pred, ymax = ls_pred), 
              alpha = 0.5, fill = 'lightblue2') +
  geom_ribbon(data = est_duration_S_FReg, 
              aes(x, mu, ymin = li, ymax = ls), 
              alpha = 0.7, fill = 'lightblue') +
  geom_line(data = est_duration_S_FReg, 
            aes(x, mu), color = 'tan1', linewidth = 1, alpha = 0.7) +
  geom_point(data = tibble(x = dat3$duration, 
                           y = dat3$s_std), 
             aes(x, y), color = 'tan1', alpha = 0.5) +
  theme_classic() +
  labs(x = 'Fruit duration (z-scores)', 
       y = 'S (z-scores)') +
  theme(legend.position = c(0.2, 0.9),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = 'Times New Roman', size = 10), 
        axis.line = element_line(linewidth = 0.25))



effect_duration_FReg <- causal_effect_FReg('beta_duration_S', 
                                       'beta_duration_FReg', 
                                       dat3$duration, 'duration', T)

effect_duration_FReg$FReg$total_effect <- 'Total effect'

conditional_duration_FReg <- causal_effect_FReg('beta_duration_S', 
                                            'beta_duration_FReg', 
                                            dat3$duration, 'duration', F)[[1]]

conditional_duration_FReg$total_effect <- 'Conditional effect'

plot_CE_duration_FReg <- 
  rbind(effect_duration_FReg$FReg, 
        conditional_duration_FReg) |> 
  group_by(intervention, total_effect) |> 
  transmute(mu = mean(est), 
            li = quantile(est, 0.025), 
            ls = quantile(est, 0.975)) |> 
  unique() |> 
  ggplot() +
  geom_errorbar(aes(intervention, ymin = li, ymax = ls, 
                    color = total_effect), width = 0, 
                position = position_dodge(width = 0.3)) +
  geom_point(aes(intervention, mu, color = total_effect), 
             position = position_dodge(width = 0.3)) +
  scale_color_manual(values = c('tan1', 'lightblue')) +
  scale_fill_manual(values = c('tan1', 'lightblue')) +
  lims(y = c(-0.5, 1.5)) +
  theme_classic() +
  labs(x = 'Intervention of duration size', 
       y = 'FReg') +
  geom_hline(yintercept = 0, linetype = 3) +
  theme(legend.position = c(0.37, 0.875),
        legend.background = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        text = element_text(family = 'Times New Roman', size = 10), 
        axis.line = element_line(linewidth = 0.25)) 

# figure duration FReg

plot_duration_FReg <- 
  plot_grid(plot_cond_duration_S_FReg, 
            plot_cond_duration_FReg, 
            plot_CE_duration_FReg, ncol = 3, 
            labels = c('(a)', '(b)', '(c)'), 
            label_fontfamily = 'Times New Roman', 
            label_fontface = 'plain', 
            label_size = 10, 
            label_x = 0.175)

plot_duration_FReg

ggsave('FReg_mod.jpg', width = 20, height = 7, units = 'cm', dpi = 1000)

cont_duration_FReg <- na.omit(effect_duration_FReg$FReg)

conditional_duration_FReg <- na.omit(conditional_duration_FReg)

diff_duration_FReg <- mean(conditional_duration_FReg$est) - mean(cont_duration_FReg$est)

tribble(~ value,               ~estimation, 
        'P(total effect > 0)', mean(cont_duration_FReg$est > 0), 
        'mu total effect', mean(cont_duration_FReg$est), 
        'sd total effect', sd(cont_duration_FReg$est)/sqrt(length(cont_duration_FReg$est)), 
        'AVG change relative to maximum (total) (%)', mean(cont_duration_FReg$porcentaje), 
        'SD change relative to maximum (total) (%)', sd(cont_duration_FReg$porcentaje)/sqrt(length(cont_duration_FReg$porcentaje)), 
        'P(conditional effect > 0)', mean(conditional_duration_FReg$est > 0), 
        'mu conditional effect', mean(conditional_duration_FReg$est), 
        'sd conditional effect', sd(conditional_duration_FReg$est)/sqrt(length(conditional_duration_FReg$est)), 
        'AVG change relative to maximum (conditional) (%)', mean(conditional_duration_FReg$porcentaje), 
        'SD change relative to maximum (conditional) (%)', sd(conditional_duration_FReg$porcentaje)/sqrt(length(conditional_duration_FReg$porcentaje)),
        'P(total > conditional)', mean(cont_duration_FReg$est > conditional_duration_FReg$est), 
        'AVG diff total-conditional (%)', 100-(mean(diff_duration_FReg)*100 / mean(conditional_duration_FReg$est)))



