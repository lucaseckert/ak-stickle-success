#### ASSORTMENT TESTS ####

## pacakges
library(tidyverse)
library(purrr)
library(ggh4x)

## constants
source(here::here("02_tidy-data/constants.R"))

## freq data
freq_data<-read.csv("02_tidy-data/assortment/freq-data.csv") %>%
  select(lake, genotype, count_exp=count_exp_raw, count_obs)

#### HELPER FUNCTIONS ####
split_geno <- function(g) {
  c(substr(g, 1, 2), substr(g, 3, 4))
}

make_geno <- function(a1, a2) {
  paste(sort(c(a1, a2)), collapse = "")
}

offspring_prob <- function(g_off, g_m, g_f) {
  off <- split_geno(g_off)
  m <- split_geno(g_m)
  f <- split_geno(g_f)
  
  probs <- 0
  for (gm in m) {
    for (gf in f) {
      g <- make_geno(gm, gf)
      if (g == g_off) probs <- probs + 0.25
    }
  }
  probs
}

#### ASSORTMENT MODEL ####
fit_assortment_counts <- function(genotype, count_obs) {
  
  ## 1. Setup
  genos <- sort(unique(genotype))
  obs <- setNames(count_obs[match(genos, genotype)], genos)
  n_total <- sum(obs)
  geno_split <- lapply(genos, split_geno)
  
  ## 2. Frequencies
  all_alleles <- unlist(geno_split)
  weights <- rep(obs, each = 2)
  p_counts <- tapply(weights, all_alleles, sum)
  p <- p_counts / sum(p_counts)
  
  ## 3. Probability & Likelihood
  calc_expected <- function(a) {
    probs <- sapply(geno_split, function(g) {
      a1 <- g[1]; a2 <- g[2]
      if (a1 == a2) { p[a1]^2 + a * p[a1] * (1 - p[a1]) } 
      else { 2 * p[a1] * p[a2] * (1 - a) }
    })
    return(probs)
  }
  
  loglik <- function(a) {
    exp_p <- calc_expected(a)
    if (any(exp_p <= 1e-10) || any(is.na(exp_p))) return(-1e12)
    sum(obs * log(exp_p))
  }
  
  ## 4. Fit the MLE
  fit <- optimize(loglik, interval = c(-0.999, 0.999), maximum = TRUE)
  a_hat <- fit$maximum
  max_ll <- fit$objective
  
  ## 5. NEW: FIND LR CONFIDENCE INTERVALS
  # We look for 'a' where LogLik(a) = MaxLL - 1.92
  target_val <- max_ll - 1.92
  
  # Helper to find where (LogLik - Target) = 0
  root_fn <- function(a) loglik(a) - target_val
  
  # Find Lower Bound (search from -0.999 to a_hat)
  lower_ci <- tryCatch(
    uniroot(root_fn, interval = c(-0.999, a_hat))$root,
    error = function(e) -0.999 # Hits boundary
  )
  
  # Find Upper Bound (search from a_hat to 0.999)
  upper_ci <- tryCatch(
    uniroot(root_fn, interval = c(a_hat, 0.999))$root,
    error = function(e) 0.999 # Hits boundary
  )
  
  ## 6. Diagnostics & Output
  final_probs <- calc_expected(a_hat)
  diag_df <- data.frame(genotype=genos, obs=obs, exp=round(final_probs*n_total, 2))
  
  list(
    a_hat = a_hat,
    lower_ci = lower_ci,
    upper_ci = upper_ci,
    logLik = max_ll,
    LR = 2 * (max_ll - loglik(0)),
    p_value = pchisq(2 * (max_ll - loglik(0)), df = 1, lower.tail = FALSE),
    diagnostics = diag_df,
    n = n_total
  )
}

## assortment results
results <- freq_data %>%
  group_by(lake) %>%
  group_split() %>%
  setNames(unique(freq_data$lake)) %>%
  lapply(function(df) {
    fit_assortment_counts(
      genotype = df$genotype,
      count_obs = df$count_obs
    )
  })

## assortment summary
summary_df <- lapply(names(results), function(lk) {
  res <- results[[lk]]
  data.frame(
    lake     = lk,
    a_hat    = res$a_hat,
    low_ci   = res$lower_ci,
    high_ci  = res$upper_ci,
    p_value  = res$p_value,
    n        = res$n
  )
}) %>% 
  bind_rows() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

## export
summary_df %>% 
  mutate(across(c(a_hat, low_ci, high_ci), round, 3),
         p_value = if_else(p_value < 0.001, "<0.001", as.character(round(p_value, 3))),
         p_adj = if_else(p_adj < 0.001, "<0.001", as.character(round(p_adj, 3)))) %>% 
  mutate(lake = factor(lake, levels = c(limnetic_lakes, benthic_lakes, both_lakes))) %>%
  mutate(ci=paste0("[", low_ci, ", ", high_ci, "]")) %>%
  arrange(lake) %>%
  flextable::flextable() %>%
  print(preview = "docx")

## residuals
all_diagnostics <- results %>%
  map_dfr(~ .x$diagnostics, .id = "lake")
rownames(all_diagnostics) <- NULL

write.csv(all_diagnostics, "02_tidy-data/assortment/assortment-residuals.csv", row.names = F, quote = F)

## export
all_diagnostics %>%
  filter(lake %in% c("RL","GL","LP")) %>% 
  mutate(lake = factor(lake, levels = c("RL","LP","GL")),
         class = if_else(substr(genotype,1,2) == substr(genotype,3,4), "pure", "admixed"),
         residual=obs-exp) %>%
  select(lake, class, genotype, obs, exp, residual) %>%
  arrange(lake, class, genotype) %>%
  flextable::flextable() %>%
  print(preview = "docx")
  
all_diagnostics %>%
  filter(lake %in% c("RL","GL","LP")) %>% 
  mutate(lake = factor(lake, levels = c("RL","LP","GL")),
         class = if_else(substr(genotype,1,2) == substr(genotype,3,4), "pure", "admixed"),
         residual=obs-exp) %>%
  group_by(lake, class) %>%
  summarize(mean_residual=mean(residual), .groups = "drop")

##### LEAVE ONE OUT FUNCTION ####
leave_one_allele_out_simple <- function(df) {
  
  # --- Fit full model ---
  full_fit <- fit_assortment_counts(
    genotype = df$genotype,
    count_obs = df$count_obs
  )
  
  if(!is.list(full_fit)) return(NULL)
  
  a_full <- full_fit$a_hat
  ll_full <- full_fit$logLik
  n_full <- full_fit$n
  
  # --- Extract alleles ---
  split_geno <- strsplit(as.character(df$genotype), "(?<=..)", perl=TRUE)
  alleles <- unique(unlist(split_geno))
  
  # --- Loop over alleles ---
  results <- lapply(alleles, function(al) {
    
    # Remove genotypes containing allele 'al'
    keep <- !sapply(split_geno, function(g) al %in% g)
    df_minus <- df[keep, ]
    
    if(nrow(df_minus) < 2) return(NULL)
    
    fit_minus <- try(
      fit_assortment_counts(
        genotype = df_minus$genotype,
        count_obs = df_minus$count_obs
      ),
      silent = TRUE
    )
    
    if(!is.list(fit_minus)) return(NULL)
    
    n_minus <- fit_minus$n
    
    data.frame(
      allele = al,
      a_full = a_full,
      a_minus = fit_minus$a_hat,
      delta_a = a_full - fit_minus$a_hat,
      logLik_full = ll_full,
      logLik_minus = fit_minus$logLik,
      delta_logLik = ll_full - fit_minus$logLik,
      n_full = n_full,
      n_minus = n_minus,
      delta_n = n_full - n_minus,
      prop_removed = (n_full - n_minus) / n_full
    )
  })
  
  bind_rows(results)
}


## leave one out results
allele_influence_a <- imap_dfr(
  split(freq_data, freq_data$lake),
  function(df, lake_name) {
    
    res <- leave_one_allele_out_simple(df)
    
    if(!is.null(res)) {
      res$lake <- lake_name
      return(res)
    }
  }
)

#### ALLELIC FIS ####

## scaled allelic f
allelic_f_unbiased_scaled <- function(df) {
  # 1. Basic counts
  n_ind <- sum(df$count_obs) 
  n_alleles <- 2 * n_ind      
  
  # Split genotypes
  g_split <- strsplit(as.character(df$genotype), "(?<=..)", perl=TRUE)
  
  # 2. Allele Frequencies
  weights <- rep(df$count_obs, each = 2)
  all_alleles <- unlist(g_split)
  p_counts <- tapply(weights, all_alleles, sum)
  p <- p_counts / n_alleles
  
  alleles <- names(p)
  
  # 3. Calculate Unbiased Expected Heterozygosity for the whole locus (Hs)
  # This is the "Total Diversity" we use for scaling
  Hs_unbiased_locus <- (n_alleles / (n_alleles - 1)) * (1 - sum(p^2))
  
  # 4. Per-Allele Stats
  results <- lapply(alleles, function(al) {
    p_i <- p[al]
    
    # Unbiased Expected Heterozygosity for this specific allele
    exp_het_i <- (n_alleles / (n_alleles - 1)) * (2 * p_i * (1 - p_i)) * n_ind
    
    # Observed Heterozygotes for this allele
    obs_het_i <- sum(df$count_obs[sapply(g_split, function(x) sum(x == al) == 1)])
    
    F_i <- 1 - (obs_het_i / exp_het_i)
    
    data.frame(
      allele = al,
      p = p_i,
      obs_het = obs_het_i,
      exp_het = exp_het_i,
      F_i = F_i
    )
  })
  
  res <- bind_rows(results)
  
  # 5. THE SCALING STEP
  # We weight each allele by its portion of the total locus diversity
  # sum(res$weight) will equal 1
  res$weight <- res$exp_het / (Hs_unbiased_locus * n_ind * 2) 
  
  # Scaled contribution
  res$F_i_scaled <- res$F_i * res$weight
  
  # 6. Global FIS (Sum of scaled contributions)
  # This will now match hierfstat's logic
  global_FIS <- sum(res$F_i_scaled)
  res$global_FIS <- global_FIS
  
  return(res)
}

## scaled allelic F results
allelic_f_summary <- freq_data %>%
  group_by(lake) %>%
  do(allelic_f_unbiased_scaled(.))

allelic_f_means<-allelic_f_summary %>% 
  group_by(allele) %>% 
  summarize(mean_F = mean(F_i), .groups = "drop") %>% 
  mutate(allele = factor(allele, levels = sources_morph))

allele_annos<- allelic_f_summary %>% 
  mutate(lake=factor(lake, levels = c(limnetic_lakes, benthic_lakes, both_lakes))) %>%
  select(lake, global_FIS) %>% 
  distinct() %>% 
  filter(global_FIS< -0.099) %>% 
  mutate(x=1, 
         y=0.11,
         label = paste0("F[IS] == ", round(global_FIS, 2)))

allelic_f_sources<-allelic_f_summary %>%
  mutate(allele = factor(allele, levels = sources_morph),
         fac = "Source Population") %>%
  ggplot() +
  geom_hline(yintercept = 0, color = "grey50", linetype = "dashed") +
  geom_point(aes(x = allele, y = F_i, fill = allele), shape = 21, size = 4, alpha = 1) +
  geom_point(data=allelic_f_means, aes(x=allele, y=mean_F), color="black", size=4, shape=3, stroke=1) +
  scale_fill_manual(values=source_cols_alt) +
  labs(y = expression(F[IS]),
       x="Source Population")+
  theme(legend.position = "none",
        axis.text.x = element_text(size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(size=12, face="bold"))+
  facet_wrap(~fac)

allelic_f_lakes<-allelic_f_summary %>%
  mutate(allele = factor(allele, levels = sources_morph),
         lake=factor(lake, levels = c(limnetic_lakes, benthic_lakes, both_lakes))) %>%
  ggplot(aes(x=1)) +
  geom_hline(yintercept = 0, color = "grey50", linetype = "dashed") +
  geom_point(aes(y = F_i, fill = allele), shape = 21, size = 4, alpha = 1) +
  geom_point(data = mutate(allelic_f_summary,
                           lake=factor(lake, levels = c(limnetic_lakes, benthic_lakes, both_lakes))),
             aes(x=1, y=global_FIS), color="black", size=4, shape=3, stroke=1) +
  geom_text(
    data = allele_annos,
    aes(x, y, label = label),
    parse = T,
    inherit.aes = FALSE,
    size=3.5)+
  scale_fill_manual(values=source_cols_alt) +
  labs(y = expression(F[italic(i)]))+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size=12, face="bold"),
        axis.title.y = element_text(size=14))+
  guides(fill = guide_legend(nrow = 1))+
  facet_wrap2(~lake, ncol = 9,
              strip = strip_themed(
                background_x = elem_list_rect(fill = all_lake_cols[-1])))

allelic_f_plot<-plot_grid(allelic_f_lakes, allelic_f_sources,
                          rel_widths = c(0.7, 0.3),
                          align="h", axis="b")

save_fig(
  plot = allelic_f_plot,
  filename = "03_figures/assortment/allele-f.tiff",
  width = 1, height=0.3)


#### HEIRCHAL ASSORTMENT MODEL ####

fit_hierarchical_assortment <- function(genotype, count_obs, eco_map) {
  # --- 1. Setup & Frequencies ---
  force(genotype); force(count_obs); force(eco_map)
  genos <- sort(unique(genotype))
  obs <- setNames(count_obs[match(genos, genotype)], genos)
  n_total <- sum(obs)
  geno_split <- lapply(genos, split_geno)
  
  all_alleles <- unlist(geno_split)
  weights <- rep(obs, each = 2)
  p_counts <- tapply(weights, all_alleles, sum)
  p_env <- p_counts / sum(p_counts)
  force(p_env)
  
  # --- 2. Probability & Likelihood ---
  calc_expected_dual <- function(a, e) {
    p <- p_env
    probs <- sapply(geno_split, function(g) {
      a1 <- g[1]; a2 <- g[2]; p1 <- p[a1]; p2 <- p[a2]
      prob <- if(a1 == a2) p1^2 else 2 * p1 * p2
      if(a1 == a2) { prob <- prob + a * p1 * (1 - p1) } 
      else { prob <- prob * (1 - a) }
      if(eco_map[a1] == eco_map[a2]) { prob <- prob * (1 + e) }
      return(max(1e-12, prob))
    })
    return(probs / sum(probs))
  }
  
  loglik <- function(par) {
    if(par[1] <= -0.999 || par[1] >= 0.999) return(-1e12)
    if(par[2] <= -0.999 || par[2] >= 10) return(-1e12)
    sum(obs * log(calc_expected_dual(par[1], par[2])))
  }
  
  # --- 3. Full Model Fit ---
  fit <- optim(par = c(0, 0), fn = loglik, control = list(fnscale = -1))
  a_hat <- fit$par[1]; e_hat <- fit$par[2]; ll_full <- fit$value
  
  # --- 4. Profile Likelihood CI Logic ---
  # Threshold: MaxLL - 1.92
  target <- ll_full - 1.92
  
  # Function to find the Profile Likelihood for one parameter
  # It fixes one param and optimizes the other
  get_profile_ci <- function(param_idx) {
    other_idx <- if(param_idx == 1) 2 else 1
    
    # Root function: find where Profile LL - target = 0
    profile_root <- function(val) {
      # Optimize the 'other' parameter while current one is fixed at 'val'
      fixed_fit <- optimize(function(other_val) {
        p_vec <- numeric(2)
        p_vec[param_idx] <- val
        p_vec[other_idx] <- other_val
        loglik(p_vec)
      }, interval = if(other_idx == 1) c(-0.99, 0.99) else c(-0.99, 5), maximum = TRUE)
      
      return(fixed_fit$objective - target)
    }
    
    # Find Lower/Upper bounds
    lower_limit <- if(param_idx == 1) -0.99 else -0.99
    upper_limit <- if(param_idx == 1) 0.99 else 5
    
    low_ci <- tryCatch(uniroot(profile_root, c(lower_limit, fit$par[param_idx]))$root, error = function(e) NA)
    upp_ci <- tryCatch(uniroot(profile_root, c(fit$par[param_idx], upper_limit))$root, error = function(e) NA)
    return(c(low_ci, upp_ci))
  }
  
  ci_a <- get_profile_ci(1)
  ci_e <- get_profile_ci(2)
  
  # --- 5. Nested Null Model (for LR_eco_effect) ---
  fit_null_e <- optimize(function(a) loglik(c(a, 0)), interval=c(-0.99, 0.99), maximum=TRUE)
  
  list(
    a_hat = a_hat, a_low = ci_a[1], a_high = ci_a[2],
    e_hat = e_hat, e_low = ci_e[1], e_high = ci_e[2],
    ll_full = ll_full,
    LR_eco_effect = 2 * (ll_full - fit_null_e$objective),
    p_val_eco = pchisq(max(0, 2 * (ll_full - fit_null_e$objective)), 1, lower.tail = FALSE),
    n = n_total
  )
}

#### ASSORTMENT BY ECOTYPE ####

## ecotype map
eco_map<- c(
  "FG" = "BB",
  "TL" = "BB",
  "WB" = "BB",
  "WT" = "BB",
  "LG" = "LL",
  "SL" = "LL",
  "SR" = "LL",
  "WK" = "LL"
)

## ecotype data
freq_data_both<-filter(freq_data, lake %in% both_lakes)

## ecotype results
results_list_eco <- lapply(split(freq_data_both, freq_data_both$lake), function(df) {
  try(fit_hierarchical_assortment(df$genotype, df$count_obs, eco_map))
})

## ecotype summary
summary_hierarchical_eco <- imap_dfr(results_list_eco, function(res, name) {
  if(is.list(res)) {
    data.frame(
      lake = name,
      a_hat = round(res$a_hat, 3),
      a_low = round(res$a_low, 3),
      a_high = round(res$a_high, 3),
      e_hat = round(res$e_hat, 3),
      e_low = round(res$e_low, 3),
      e_high = round(res$e_high, 3),
      logLik = round(res$ll_full, 2),
      LR = round(res$LR_eco_effect, 2),
      p_value = round(res$p_val_eco, 4)
    )
  }
}) %>% mutate(p_adj = p.adjust(p_value, method = "BH"))

summary_hierarchical_eco %>% 
  mutate(across(c(a_hat, a_low, a_high, e_hat, e_low, e_high), round, 3),
         p_val_eco = if_else(p_value < 0.001, "<0.001", as.character(round(p_value, 3))),
         p_adj_eco = if_else(p_adj < 0.001, "<0.001", as.character(round(p_adj, 3)))) %>% 
  mutate(lake = factor(lake, levels = both_lakes)) %>%
  mutate(e_ci = paste0("[", e_low, ", ", e_high, "]"),
         a_ci = paste0("[", a_low, ", ", a_high, "]")) %>%
  select(lake, e_hat, e_ci, a_hat, a_ci, logLik, LR, p_value, p_adj) %>%
  arrange(lake) %>%
  flextable::flextable() %>%
  print(preview = "docx")

#### ASSORTMENT BY REGION ####

## region map
region_map<-c(
  "FG" = "MS",
  "TL" = "KN",
  "WB" = "MS",
  "WT" = "KN",
  "LG" = "MS",
  "SL" = "KN",
  "SR" = "MS",
  "WK" = "KN"
)

## region results
results_list_region <- lapply(split(freq_data, freq_data$lake), function(df) {
  try(fit_hierarchical_assortment(df$genotype, df$count_obs, region_map))
})

## region summary
summary_hierarchical_region <- imap_dfr(results_list_region, function(res, name) {
  if(is.list(res)) {
    data.frame(
      lake = name,
      a_hat = round(res$a_hat, 3),
      a_low = round(res$a_low, 3),
      a_high = round(res$a_high, 3),
      r_hat = round(res$e_hat, 3),
      r_low = round(res$e_low, 3),
      r_high = round(res$e_high, 3),
      logLik = round(res$ll_full, 2),
      LR = round(res$LR_eco_effect, 2),
      p_value = round(res$p_val_eco, 4)
    )
  }
}) %>% mutate(p_adj = p.adjust(p_value, method = "BH"))

summary_hierarchical_region %>% 
  mutate(across(c(a_hat, a_low, a_high, r_hat, r_low, r_high), round, 3),
         p_val_eco = if_else(p_value < 0.001, "<0.001", as.character(round(p_value, 3))),
         p_adj_eco = if_else(p_adj < 0.001, "<0.001", as.character(round(p_adj, 3)))) %>% 
  mutate(lake = factor(lake, levels = c(limnetic_lakes, benthic_lakes, both_lakes))) %>%
  mutate(r_ci = paste0("[", r_low, ", ", r_high, "]"),
         a_ci = paste0("[", a_low, ", ", a_high, "]")) %>%
  select(lake, r_hat, r_ci, a_hat, a_ci, logLik, LR, p_value, p_adj) %>%
  arrange(lake) %>%
  flextable::flextable() %>%
  print(preview = "docx")

## all results
results_summary<-bind_rows(
  summary_df %>% mutate(group = "source") %>% select(lake, group, estimate = a_hat, p_value, p_adj),
  summary_hierarchical_region %>% mutate(group = "region") %>% select(lake, group, estimate = r_hat, p_value, p_adj),
  summary_hierarchical_eco %>% mutate(group = "ecotype") %>% select(lake, group, estimate = e_hat, p_value, p_adj)
)

## saving
write.csv(results_summary, file = "02_tidy-data/assortment/assortment-test-results.csv",
          row.names = F, quote = F)

## CIs
param_labs<-c(paste0("Source: hat(a)[s]"), paste0("Region: hat(a)[r]"), paste0("Ecotype: hat(a)[e]"))
all_results<-bind_rows(summary_df %>% mutate(group = "source") %>% rename(a_low=low_ci, a_high=high_ci),
                       summary_hierarchical_region %>% mutate(group = "region"),
                       summary_hierarchical_eco %>% mutate(group = "ecotype")) %>% 
  select(lake, group, ends_with("hat"), ends_with("low"), ends_with("high")) %>% 
  pivot_longer(cols = -c(lake, group), names_to = c("param", "stat"), names_sep = "_") %>% 
  pivot_wider(names_from = stat, values_from = value) %>% 
  filter(!is.na(hat)) %>% 
  mutate(lake = factor(lake, levels = c(limnetic_lakes, benthic_lakes, both_lakes)),
         group = factor(group, levels = c("source", "region", "ecotype"), labels = c("Source", "Region", "Ecotype")),
         param = factor(param, levels = c("a", "r", "e"), labels = param_labs))

ci_plot<-ggplot(all_results, aes(x=group, y=hat, color=param))+
  geom_hline(yintercept = 0, color = "grey50", linetype = "dashed") +
  geom_point(position = position_dodge(width = 0.5), size=3)+
  geom_errorbar(aes(ymin=low, ymax=high), width=0.5, position = position_dodge(width = 0.5), linewidth = 0.5)+
  labs(y = "Estimate", x="Model", color="Parameter:")+
  scale_color_manual(labels = function(x) parse(text = x), values=c("#619CFF","#F8766D","#00BA38"))+
  theme(legend.position = "top",
        legend.title = element_text(size=12),
        axis.text.x = element_text(angle = 45, vjust = 0.75, hjust = 0.6))+
  facet_grid2(
    cols = vars(lake),
    scales = "free_x",
    strip = strip_themed(
      background_x = elem_list_rect(fill = all_lake_cols[-1])))
  
save_fig(
  plot = ci_plot,
  filename = "03_figures/assortment/ci-fig.tiff",
  width = 1, height=0.33)

#### BAD TESTS ####

reg_data<-freq_data %>%
  mutate(a1 = region_map[substr(genotype, 1, 2)],
         a2 = region_map[substr(genotype, 3, 4)]) %>%
  mutate(genotype = map2_chr(a1, a2, ~paste(sort(c(.x, .y)), collapse = ""))) %>%
  group_by(lake, genotype) %>%
  summarize(count_obs = sum(count_obs))

## results
reg_results <- reg_data %>%
  group_by(lake) %>%
  group_split() %>%
  setNames(unique(reg_data$lake)) %>%
  lapply(function(df) {
    fit_assortment_counts(
      genotype = df$genotype,
      count_obs = df$count_obs
    )
  })

## assortment summary
reg_summary <- lapply(names(reg_results), function(lk) {
  res <- reg_results[[lk]]
  data.frame(
    lake     = lk,
    a_hat    = res$a_hat,
    low_ci   = res$lower_ci,
    high_ci  = res$upper_ci,
    logLik   = res$logLik,
    LR       = res$LR,
    p_value  = res$p_value,
    n        = res$n
  )
}) %>% 
  bind_rows() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

eco_data<-freq_data %>%
  filter(lake %in% both_lakes) %>% 
  mutate(a1 = eco_map[substr(genotype, 1, 2)],
         a2 = eco_map[substr(genotype, 3, 4)]) %>%
  mutate(genotype = map2_chr(a1, a2, ~paste(sort(c(.x, .y)), collapse = ""))) %>%
  group_by(lake, genotype) %>%
  summarize(count_obs = sum(count_obs))

## results
eco_results <- eco_data %>%
  group_by(lake) %>%
  group_split() %>%
  setNames(unique(eco_data$lake)) %>%
  lapply(function(df) {
    fit_assortment_counts(
      genotype = df$genotype,
      count_obs = df$count_obs
    )
  })

## assortment summary
eco_summary <- lapply(names(eco_results), function(lk) {
  res <- eco_results[[lk]]
  data.frame(
    lake     = lk,
    a_hat    = res$a_hat,
    low_ci   = res$lower_ci,
    high_ci  = res$upper_ci,
    logLik   = res$logLik,
    LR       = res$LR,
    p_value  = res$p_value,
    n        = res$n
  )
}) %>% 
  bind_rows() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

both_summary<-bind_rows(
  reg_summary %>% mutate(group = "region"),
  eco_summary %>% mutate(group = "ecotype")) %>% 
  mutate(across(c(a_hat, low_ci, high_ci, logLik, LR), round, 3),
       p_value = if_else(p_value < 0.001, "<0.001", as.character(round(p_value, 3))),
       p_adj = if_else(p_adj < 0.001, "<0.001", as.character(round(p_adj, 3)))) %>%
  mutate(ci=paste0("[", low_ci, ", ", high_ci, "]")) %>%
  mutate(lake = factor(lake, levels = c(limnetic_lakes, benthic_lakes, both_lakes))) %>%
  arrange(group,lake) %>%
  select(group,lake,a_hat,ci,logLik,LR,p_value,p_adj,n) %>%
  flextable::flextable() %>%
  print(preview = "docx")

#### ALL TOGETHER ####

fit_triple_assortment <- function(genotype, count_obs, eco_map, region_map) {
  
  ## 1. Setup
  force(genotype); force(count_obs); force(eco_map); force(region_map)
  genos <- sort(unique(genotype))
  obs <- setNames(count_obs[match(genos, genotype)], genos)
  n_total <- sum(obs)
  geno_split <- lapply(genos, split_geno)
  
  ## 2. Frequencies
  all_alleles <- unlist(geno_split)
  weights <- rep(obs, each = 2)
  p_counts <- tapply(weights, all_alleles, sum)
  p_env <- p_counts / sum(p_counts)
  
  ## 3. The Triple Probability Function
  calc_expected_triple <- function(a, e, r) {
    p <- p_env
    probs <- sapply(geno_split, function(g) {
      a1 <- g[1]; a2 <- g[2]; p1 <- p[a1]; p2 <- p[a2]
      
      # Level 1: Allele Assortment (Baseline)
      prob <- if(a1 == a2) p1^2 + a*p1*(1-p1) else 2*p1*p2*(1-a)
      
      # Level 2: Region Bonus (Apply if in same region)
      if(region_map[a1] == region_map[a2]) {
        prob <- prob * (1 + r)
      }
      
      # Level 3: Ecotype Bonus (Apply if same ecotype)
      if(eco_map[a1] == eco_map[a2]) {
        prob <- prob * (1 + e)
      }
      
      return(max(1e-12, prob))
    })
    return(probs / sum(probs))
  }
  
  ## 4. Likelihood
  loglik <- function(par) {
    # par = c(a, e, r)
    if(par[1] <= -0.99 || par[1] >= 0.99) return(-1e12)
    if(par[2] <= -0.99 || par[2] >= 10) return(-1e12)
    if(par[3] <= -0.99 || par[3] >= 10) return(-1e12)
    
    sum(obs * log(calc_expected_triple(par[1], par[2], par[3])))
  }
  
  ## 5. Fit
  fit <- optim(par = c(0, 0, 0), fn = loglik, control = list(fnscale = -1))
  max_ll <- fit$value
  
  ## 6. Profile CI Logic (Simplified for 3D)
  target <- max_ll - 1.92
  
  get_3d_profile_ci <- function(param_idx) {
    other_indices <- setdiff(1:3, param_idx)
    
    profile_root <- function(val) {
      # Optimize the TWO other parameters while param_idx is fixed
      fixed_fit <- optim(par = c(0, 0), fn = function(others) {
        p_vec <- numeric(3)
        p_vec[param_idx] <- val
        p_vec[other_indices] <- others
        loglik(p_vec)
      }, control = list(fnscale = -1))
      return(fixed_fit$value - target)
    }
    
    low <- tryCatch(uniroot(profile_root, c(-0.9, fit$par[param_idx]))$root, error = function(e) NA)
    upp <- tryCatch(uniroot(profile_root, c(fit$par[param_idx], 5))$root, error = function(e) NA)
    return(c(low, upp))
  }
  
  # Note: Profiling in 3D is computationally heavier
  ci_a <- get_3d_profile_ci(1)
  ci_e <- get_3d_profile_ci(2)
  ci_r <- get_3d_profile_ci(3)
  
  list(
    a_hat = fit$par[1], a_ci = ci_a,
    e_hat = fit$par[2], e_ci = ci_e,
    r_hat = fit$par[3], r_ci = ci_r,
    ll_full = max_ll,
    n = n_total
  )
}

# 1. Run the analysis
results_list <- freq_data_both %>%
  group_split(lake) %>%
  setNames(unique(freq_data_both$lake)) %>%
  map(~ try(fit_triple_assortment(.x$genotype, .x$count_obs, eco_map, region_map)))

# 2. Extract into a summary table
summary_triple <- imap_dfr(results_list, function(res, name) {
  if(is.list(res)) {
    data.frame(
      lake   = name,
      a_hat  = res$a_hat,
      a_low  = res$a_ci[1],
      a_high = res$a_ci[2],
      e_hat  = res$e_hat,
      e_low  = res$e_ci[1],
      e_high = res$e_ci[2],
      r_hat  = res$r_hat,
      r_low  = res$r_ci[1],
      r_high = res$r_ci[2],
      logLik = res$ll_full,
      n      = res$n
    )
  }
})


