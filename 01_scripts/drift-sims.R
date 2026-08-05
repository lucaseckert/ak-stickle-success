#### DRIFT SIMULATIONS ####

## Description: adapted from Dan Bolnick's code
## Author: Lucas Eckert (lucas.eckert@mail.mcgill.ca)

## packages
library(tidyverse)

## data
stocks<-read.csv("02_tidy-data/other-data/actual-stocks.csv")

# N_replicate_simulations remains the total number of runs you want to do.
N_replicate_simulations <- 1000

# Parameters (These can remain the same or be adjusted)
N_offspring <- 50
Sample2021 <- 100 
N_genotypes <- 8 # Based on the original code, the simulation uses 4 frequencies

### 1. Define Custom Colonists and Founder Numbers ###

# Define your new Colonists matrix
count_matrix<-stocks %>% select(-lake,-total) %>% as.matrix()
My_Colonists<-prop.table(count_matrix, margin = 1)

# Define the corresponding Founder Population Size (N) for EACH ROW/SCENARIO
My_Founder_Ns <- stocks$total
N_scenarios <- nrow(My_Colonists)
lakes<-stocks$lake

# Initialize recordkeeping for ALL scenarios (N_scenarios * N_replicate_simulations rows)
recordkeeping <- data.frame(matrix(nrow = N_replicate_simulations * N_scenarios, ncol = N_genotypes + 1))
colnames(recordkeeping) <- c(paste0("Change_", colnames(count_matrix)), "Lake")

### 2. Modified Simulation Loop ###

# Outer loop iterates over each scenario (row in My_Colonists)
for(s in 1:N_scenarios){
  
  Founder_N <- My_Founder_Ns[s]
  Founder_pool_frequencies <- My_Colonists[s, 1:N_genotypes]
  
  # Inner loop runs the replicate simulations for the current scenario
  for(i in 1:N_replicate_simulations){
    
    # Calculate the starting index for recording results
    record_index <- (s - 1) * N_replicate_simulations + i
    
    # --- Core Simulation Steps (Same as original, but using scenario-specific N) ---
    
    # F0 Generation
    F0  <- rmultinom(n = 1, size = Founder_N, prob = Founder_pool_frequencies)
    F0 <- t(F0)/Founder_N
    
    # F1 Generation
    F1_size <- Founder_N * N_offspring
    F1 <- rmultinom(n = 1, size = F1_size, prob = F0)
    F1 <- t(F1)/F1_size
    
    # F2 Generation
    F2_size <- Founder_N * N_offspring^2
    F2<- rmultinom(n = 1, size = F2_size, prob = F1)
    F2 <- t(F2)/F2_size
    
    # Sample
    Sample <- rmultinom(n = 1, size = Sample2021 , prob = F2)
    Sample <- t(Sample)/( Sample2021 )
    
    # Calculate Change
    Change <- Sample - Founder_pool_frequencies
    
    # --- Record Results for this replicate ---
    
    recordkeeping[record_index, 1:N_genotypes] <- Change
    recordkeeping[record_index, "Lake"] <- lakes[s]
  }
}

recordkeeping_long<-recordkeeping %>% pivot_longer(cols = c(starts_with("Change")), names_to = "source",
                               names_prefix = "Change_")

for(i in 1:nrow(recordkeeping_long)){
  Lake<-recordkeeping_long$Lake[i]
  Source<-recordkeeping_long$source[i]
  if(filter(stocks,lake==Lake)[[Source]]==0){
    recordkeeping_long$value[i]<-NA
  }
}

saveRDS(recordkeeping_long, file = "02_tidy-data/other-data/drift-sims.RDS")
drift_sims<-readRDS("02_tidy-data/other-data/drift-sims.RDS") %>% 
  rename(lake=Lake)

## for factors
sources<-c("LG","SR","WK","SL","WT","WB","FG","TL")
lakes<-c("CC","LL","LP","CL","FL","HL","RL","LO","GL")
source_cols_alt<-c("#08519c","#3182bd","#6baed6","#bdd7e7",
                   "#006d2c","#31a354","#74c476","#bae4b3")

cis<-drift_sims %>% 
  group_by(lake,source) %>% 
  summarize(ci_low=quantile(value, 0.025, na.rm=T),
            ci_high=quantile(value, 0.975, na.rm=T))

success<-read.csv("02_tidy-data/success-data.csv") %>% 
  filter(time=="delta_net") %>% 
  select(lake,source,delta) %>% 
  left_join(cis) %>% 
  mutate(sig=factor(if_else(delta<ci_low | delta>ci_high, 1,0))) %>% 
  mutate(source=factor(source, levels=sources),
         lake=factor(lake, levels=lakes))

drift_sims<-drift_sims %>% mutate(source=factor(source, levels=sources),
                                  lake=factor(lake, levels=lakes))

ggplot(drift_sims, aes(x=value, fill=source))+
  geom_histogram(binwidth = 0.01, alpha=0.75) +
  geom_vline(
    data = success,                                       
    aes(xintercept = delta, color = source, linetype = sig),       
    linewidth=1, alpha=1) +
  scale_linetype_manual(values=c(2,1), guide="none")+
  scale_fill_manual(values=source_cols_alt)+
  scale_color_manual(values=source_cols_alt, guide="none")+
  facet_wrap(~lake, scales="free")+
  labs(x="Change in Proportion", fill="Source\nPopulation")+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        legend.position = "right",
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank(),
        axis.text = element_text(size=12),
        strip.text = element_text(size=12),
        strip.background = element_rect(fill = "grey95"),
        legend.text = element_text(size=12),
        legend.title = element_text(size = 14, hjust=0.5, face="bold"))

ggsave(filename = "03_figures/drift-sims.png",
       device = "png",
       width = 1440, height=1040, units = "px",
       dpi=300, scale=2)

success$P_value <- NA
success$N_replicates <- NA # To store the count of simulations used

## p values
for (j in 1:nrow(success)) {
  
  current_lake <- as.character(success$lake[j])
  current_source <- as.character(success$source[j])
  
  sim_subset <- recordkeeping_long[
    recordkeeping_long$Lake == current_lake & 
      recordkeeping_long$source == current_source, 
  ]
  
  # Skip if no simulation data is found for this combination
  if (nrow(sim_subset) == 0) next
  
  # --- Perform the Two-Tailed Test ---
  observed_delta <- success$delta[j]
  
  # Get the vector of simulated changes for this specific genotype
  simulated_changes <- sim_subset$value
  
  # Count how many simulated changes (absolute value) are >= the observed change (absolute value)
  n_extreme <- sum(abs(simulated_changes) >= abs(observed_delta))
  
  # Calculate the P-value
  P_value <- n_extreme / length(simulated_changes)
  
  # --- Store the Result ---
  success$P_value[j] <- P_value
  success$N_replicates[j] <- length(simulated_changes)
}

## new method
p_val_results <- drift_sims %>%
  group_by(lake, source) %>%
  summarize(
    # Get observed delta for this specific lake/source
    # (We use first() because delta is the same for all sim rows in a group)
    obs_delta = success$delta[success$lake == cur_group()$lake & 
                                success$source == cur_group()$source][1],
    
    # Count extreme values (Two-tailed)
    n_extreme = sum(abs(value) >= abs(obs_delta), na.rm = TRUE),
    n_sim = sum(!is.na(value)),
    
    # Calculate P-value with the +1 correction
    p_rand = (n_extreme + 1) / (n_sim + 1),
    
    # Calculate null distribution metrics for reporting
    null_mean = mean(value, na.rm = TRUE),
    ci_low = quantile(value, 0.025, na.rm = TRUE),
    ci_high = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# 2. Merge back into your success dataframe
success_final <- success %>%
  select(lake, source, delta) %>%
  left_join(p_val_results, by = c("lake", "source")) %>%
  mutate(
    sig = factor(if_else(p_rand <= 0.05, 1, 0)),
    source = factor(source, levels = sources),
    lake = factor(lake, levels = lakes)
  ) %>% 
  mutate(across(where(is.numeric), round, 3)) %>% 
  mutate(ci=paste0(ci_low, " - ", ci_high)) %>% 
  select(lake,source,delta,null_mean,ci,p_rand)

## export
write.csv(success_final, file = "03_figures/tables/drift-sim-results.csv",
          row.names = F, quote = F)

