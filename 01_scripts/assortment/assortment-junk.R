##### Assortment JUNK #####

eco_data<-freq_data %>%
  filter(lake %in% both_lakes) %>% 
  mutate(a1 = eco_map[substr(genotype, 1, 2)],
         a2 = eco_map[substr(genotype, 3, 4)]) %>%
  mutate(genotype = map2_chr(a1, a2, ~paste(sort(c(.x, .y)), collapse = ""))) %>%
  group_by(lake, genotype) %>%
  summarize(count_obs = sum(count_obs))

eco_results <- eco_data %>%
  group_by(lake) %>%
  group_split() %>%
  setNames(unique(eco_data$lake)) %>%
  lapply(function(df) {
    fit_assortment_counts(
      genotype = df$genotype, # This will now be "BB", "LL", or "BL"
      count_obs = df$count_obs # This is the total count of offspring for each ecotype combination in this lake
    )
  })

## clean up results
eco_summary_df <- bind_rows( 
  lapply(names(eco_results), function(lk) {
    cbind(lake = lk, as.data.frame(eco_results[[lk]])) # This will give us a summary of the ecotype-level assortment test for each lake
  })
) %>% 
  mutate(p_adj = p.adjust(p_value, method = "BH"))

## diags
all_diagnostics <- results %>%
  map_dfr(~ .x$diagnostics, .id = "lake")
rownames(all_diagnostics) <- NULL

all_diagnostics %>% group_by(genotype) %>% 
  summarise(mean_residual = mean(residual), .groups = "drop") %>% view()

allele_df <- all_diagnostics %>%
  mutate(a1 = substr(genotype, 1, 2),
         a2 = substr(genotype, 3, 4)) %>%
  pivot_longer(cols = c(a1, a2), names_to = "pos", values_to = "allele") %>% 
  mutate(type=if_else(substr(genotype, 1, 2) == substr(genotype, 3, 4), "homozygote", "heterozygote"))

allele_trends<-allele_df %>%
  group_by(allele) %>%
  summarize(mean_residual = mean(residual),
            freq_excess = sum(residual > 0),
            total_obs = n())

geno_trends<-allele_df %>%
  group_by(genotype) %>%
  summarize(mean_residual = mean(residual),
            freq_excess = sum(residual > 0),
            total_obs = n())

allele_order<-allele_df %>% 
  group_by(allele) %>% 
  summarize(mean_exp = mean(exp)) %>% 
  arrange(-mean_exp) %>% 
  pull(allele)

allele_means<-allele_df %>% 
  group_by(allele, type) %>% 
  summarize(mean_residual = mean(residual), .groups = "drop")

allele_means_sub<-allele_df %>% 
  filter(lake %in% c("LP","RL","FL","GL")) %>% 
  group_by(allele, type) %>% 
  summarize(mean_residual = mean(residual), .groups = "drop")

allele_plot_df<-bind_rows(
  mutate(filter(allele_df, lake %in% c("LP","RL","FL","GL")), group="Non-Random Lakes"),
  mutate(allele_df, group="All Lakes"))

allele_plot_means<-bind_rows(allele_means_sub %>% mutate(group="Non-Random Lakes"),
                             allele_means %>% mutate(group="All Lakes"))

allele_plot_df %>% mutate(allele=factor(allele, levels=allele_order)) %>%
  ggplot(aes(x = type, y = residual, fill = type))+
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed", linewidth = 0.5)+
  geom_point(shape=21, position = position_dodge(width = 1), size = 2.5, alpha = 0.5)+
  geom_line(
    data = mutate(allele_plot_means, allele=factor(allele, levels=allele_order)),
    aes(x = type, y = mean_residual, group = interaction(allele, group)),
    linewidth = 0.5,
    color = "black",
    linetype = "dashed") +
  geom_point(
    data = mutate(allele_plot_means, allele=factor(allele, levels=allele_order)),
    aes(x = type, y = mean_residual),
    position = position_dodge(width = 1),
    shape = 3, size = 3, stroke = 1,
    show.legend = FALSE) +
  scale_fill_manual(values = geno_cols, labels = geno_labs)+
  labs(y = "Obs. - Exp. Count",
       title = "Highest  ---------------  Mean Source Population Frequency  ---------------  Lowest")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust=0.5, size=12))+
  guides(fill = guide_legend(
    override.aes = list(
      size = 4,
      alpha = 1,
      shape = 21)))+
  facet_grid2(
    cols = vars(allele),
    rows = vars(group),
    strip = strip_themed(
      background_x = elem_list_rect(
        fill=c(limnetic_col,benthic_col,benthic_col,limnetic_col,limnetic_col,limnetic_col,benthic_col,benthic_col))))

## from HWE
raw_allele_df <- freq_data %>%
  mutate(residual = count_obs-count_exp) %>% 
  mutate(a1 = substr(genotype, 1, 2),
         a2 = substr(genotype, 3, 4)) %>%
  pivot_longer(cols = c(a1, a2), names_to = "pos", values_to = "allele") %>% 
  mutate(type=if_else(substr(genotype, 1, 2) == substr(genotype, 3, 4), "homozygote", "heterozygote"))

raw_allele_plot_df<-bind_rows(
  mutate(filter(raw_allele_df, lake %in% c("LP","RL","FL","GL")), group="Non-Random Lakes"),
  mutate(raw_allele_df, group="All Lakes"))

raw_allele_means<-raw_allele_df %>% 
  group_by(allele, type) %>% 
  summarize(mean_residual = mean(residual), .groups = "drop")

raw_allele_means_sub<-raw_allele_df %>% 
  filter(lake %in% c("LP","RL","FL","GL")) %>% 
  group_by(allele, type) %>% 
  summarize(mean_residual = mean(residual), .groups = "drop")

raw_allele_plot_means<-bind_rows(raw_allele_means_sub %>% mutate(group="Non-Random Lakes"),
                                 raw_allele_means %>% mutate(group="All Lakes"))

raw_allele_plot_df %>% mutate(allele=factor(allele, levels=allele_order)) %>%
  ggplot(aes(x = type, y = residual, fill = type))+
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed", linewidth = 0.5)+
  geom_point(shape=21, position = position_dodge(width = 1), size = 2.5, alpha = 0.5)+
  geom_line(
    data = mutate(raw_allele_plot_means, allele=factor(allele, levels=allele_order)),
    aes(x = type, y = mean_residual, group = interaction(allele, group)),
    linewidth = 0.5,
    color = "black",
    linetype = "dashed") +
  geom_point(
    data = mutate(raw_allele_plot_means, allele=factor(allele, levels=allele_order)),
    aes(x = type, y = mean_residual),
    position = position_dodge(width = 1),
    shape = 3, size = 3, stroke = 1,
    show.legend = FALSE) +
  scale_fill_manual(values = geno_cols, labels = geno_labs)+
  labs(y = "Obs. - Exp. Count",
       title = "Highest  ---------------  Mean Source Population Frequency  ---------------  Lowest")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        plot.title = element_text(hjust=0.5, size=12))+
  guides(fill = guide_legend(
    override.aes = list(
      size = 4,
      alpha = 1,
      shape = 21)))+
  facet_grid2(
    cols = vars(allele),
    rows = vars(group),
    strip = strip_themed(
      background_x = elem_list_rect(
        fill=c(limnetic_col,benthic_col,benthic_col,limnetic_col,limnetic_col,limnetic_col,benthic_col,benthic_col))))

allelic_f <- function(df) {
  # 1. Total number of individuals
  n_total <- sum(df$count_obs)
  
  # 2. Calculate Allele Frequencies (p) correctly
  # Split genotypes into individual alleles
  g_split <- strsplit(as.character(df$genotype), "(?<=..)", perl=TRUE)
  all_alleles <- unlist(g_split)
  
  # Repeat alleles based on the observed count of that genotype
  weights <- rep(df$count_obs, each = 2)
  p_counts <- tapply(weights, all_alleles, sum)
  p <- p_counts / (2 * n_total)
  
  unique_alleles <- names(p)
  
  # 3. Calculate F for each allele
  results <- lapply(unique_alleles, function(al) {
    p_i <- p[al]
    
    # EXPECTED: Under HWE, probability of being a het for allele 'al' 
    # is 2 * p_i * (Sum of all other p's). Since (Sum of all others) = (1 - p_i):
    exp_het_hwe <- 2 * p_i * (1 - p_i) * n_total
    
    # OBSERVED: Count fish that have exactly ONE copy of 'al'
    # This avoids counting homozygotes (two copies) or non-carriers (zero copies)
    obs_het <- sum(df$count_obs[
      sapply(g_split, function(x) sum(x == al) == 1)
    ])
    
    data.frame(
      allele = al,
      freq = round(p_i, 3),
      obs_het = obs_het,
      exp_het_hwe = round(exp_het_hwe, 2),
      F_allele = round(1 - (obs_het / exp_het_hwe), 3)
    )
  })
  
  return(bind_rows(results))
}


allelic_f_eco_summary <- freq_data %>%
  filter(lake %in% both_lakes) %>%
  mutate(
    a1 = substr(genotype, 1, 2),
    a2 = substr(genotype, 3, 4),
    eco1 = eco_map[a1],
    eco2 = eco_map[a2]) %>% 
  mutate(genotype=paste0(eco1,eco2)) %>% 
  group_by(lake) %>%
  do(calc_driver_F(.))

calc_ecotype_F <- function(df, eco_map) {
  
  n_total <- sum(df$count_obs)
  
  # Split genotypes
  g_split <- strsplit(as.character(df$genotype), "(?<=..)", perl=TRUE)
  
  # Allele frequencies
  all_alleles <- unlist(g_split)
  weights <- rep(df$count_obs, each = 2)
  p_counts <- tapply(weights, all_alleles, sum)
  p <- p_counts / (2 * n_total)
  
  # Total frequency per ecotype
  eco_of_allele <- eco_map[names(p)]
  P_eco <- tapply(p, eco_of_allele, sum)
  
  # Observed counts per ecotype pairing
  obs_counts <- sapply(g_split, function(g) {
    eco_map[g[1]] == eco_map[g[2]]
  })
  
  obs_within <- sum(df$count_obs[obs_counts])
  
  # Expected under HWE
  exp_within <- sum(P_eco^2) * n_total
  
  F_ecotype <- (obs_within - exp_within) / exp_within
  
  data.frame(
    F_ecotype = F_ecotype,
    obs_within = obs_within,
    exp_within = exp_within
  )
}

eco_f_results<-freq_data %>%
  filter(lake %in% both_lakes) %>%
  group_by(lake) %>%
  do(calc_ecotype_F(., eco_map))


eco_diags<-results_list_eco$GL$diagnostics %>%
  mutate(
    a1 = substr(genotype, 1, 2),
    a2 = substr(genotype, 3, 4),
    eco1 = eco_map[a1],
    eco2 = eco_map[a2],
    pairing_type = case_when(
      eco1 == eco2 & a1 == a2 ~ "Homozygote (Same Allele)",
      eco1 == eco2 & a1 != a2 ~ "Within-Ecotype Heterozygote",
      eco1 != eco2           ~ "Between-Ecotype Heterozygote"
    )
  )

ggplot(eco_diags, aes(x = pairing_type, y = res, fill = pairing_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "GL Diagnostics: Ecotype Assortment Signal",
       y = "Residual (Observed - Expected Count)",
       x = "Mating Pair Type")

## allelic fis by ecotype
calc_allele_ecotype <- function(df, mapping) {
  n_total <- sum(df$count_obs)
  
  # 1. Get Allele Frequencies
  g_split <- strsplit(as.character(df$genotype), "(?<=..)", perl=TRUE)
  all_alleles <- unlist(g_split)
  weights <- rep(df$count_obs, each = 2)
  p_counts <- tapply(weights, all_alleles, sum)
  p <- p_counts / sum(p_counts)
  
  unique_alleles <- names(p)
  
  # 2. Iterate through each allele
  results <- lapply(unique_alleles, function(al) {
    p_i <- p[al]
    eco_i <- mapping[al] # Which ecotype does this allele belong to?
    
    # Identify which other alleles are same vs different ecotype
    others <- names(p)[names(p) != al]
    same_eco_alleles <- others[mapping[others] == eco_i]
    diff_eco_alleles <- others[mapping[others] != eco_i]
    
    # CALCULATE EXPECTED (HWE)
    # Expected same-eco hets: 2 * p_i * sum(p_j in same ecotype)
    exp_same <- 2 * p_i * sum(p[same_eco_alleles]) * n_total
    # Expected diff-eco hets: 2 * p_i * sum(p_k in diff ecotype)
    exp_diff <- 2 * p_i * sum(p[diff_eco_alleles]) * n_total
    
    # CALCULATE OBSERVED
    # Find genotypes that are hets for 'al'
    is_het <- sapply(g_split, function(x) sum(x == al) == 1)
    het_df <- df[is_het, ]
    
    # Of those hets, which are same-eco vs diff-eco?
    obs_same <- sum(het_df$count_obs[sapply(strsplit(het_df$genotype, "(?<=..)", perl=TRUE), 
                                            function(x) mapping[x[x != al]] == eco_i)])
    obs_diff <- sum(het_df$count_obs[sapply(strsplit(het_df$genotype, "(?<=..)", perl=TRUE), 
                                            function(x) mapping[x[x != al]] != eco_i)])
    
    data.frame(
      allele = al,
      ecotype = eco_i,
      type = c("Intra-Ecotype", "Inter-Ecotype"),
      obs = c(obs_same, obs_diff),
      exp = c(exp_same, exp_diff),
      F_val = c(1 - (obs_same/exp_same), 1 - (obs_diff/exp_diff))
    )
  })
  
  bind_rows(results)
}

allele_eco_summary <- freq_data %>%
  filter(lake %in% both_lakes) %>%
  group_by(lake) %>%
  do(calc_allele_ecotype(., eco_map))

## LEAVE ONE OUT
leave_one_allele_out <- function(df, eco_map) {
  
  # --- Full model ---
  full_fit <- fit_hierarchical_assortment(
    df$genotype,
    df$count_obs,
    eco_map
  )
  
  if(!is.list(full_fit)) return(NULL)
  
  e_full <- full_fit$e_hat
  
  # --- Extract alleles ---
  split_geno <- strsplit(as.character(df$genotype), "(?<=..)", perl=TRUE)
  alleles <- unique(unlist(split_geno))
  
  # --- Loop over alleles ---
  results <- lapply(alleles, function(al) {
    
    # Keep genotypes NOT containing this allele
    keep <- !sapply(split_geno, function(g) al %in% g)
    
    df_minus <- df[keep, ]
    
    # Skip if too little data remains
    if(nrow(df_minus) < 2) return(NULL)
    
    fit_minus <- try(
      fit_hierarchical_assortment(
        df_minus$genotype,
        df_minus$count_obs,
        eco_map
      ),
      silent = TRUE
    )
    
    if(!is.list(fit_minus)) return(NULL)
    
    data.frame(
      allele = al,
      e_full = e_full,
      e_minus = fit_minus$e_hat,
      delta_e = e_full - fit_minus$e_hat
    )
  })
  
  bind_rows(results)
}

## results
allele_influence <- imap_dfr(
  split(freq_data_both, freq_data_both$lake),
  function(df, lake_name) {
    
    res <- leave_one_allele_out(df, eco_map)
    
    if(!is.null(res)) {
      res$lake <- lake_name
      return(res)
    }
  }
)


summary_hierarchical_region %>%
  select(-a_hat) %>% 
  mutate(across(is.numeric, round, 3)) %>%
  mutate(lake = factor(lake, levels = c(limnetic_lakes, benthic_lakes, both_lakes))) %>%
  arrange(lake) %>%
  flextable::flextable() %>% 
  print(preview = "docx")
  


