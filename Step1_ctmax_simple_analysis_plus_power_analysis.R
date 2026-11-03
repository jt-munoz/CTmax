
# Read data
################################################################################

# Clear workspace
rm(list = ls())

# Packages
library(tidyverse)

# Read data
data_orig <- readxl::read_excel("/Users/cfree/Dropbox/Chris/UCSB/workshops/mongolia_2026/ctmax/data/CTMax_Data_Eg_Uur2.xlsx", na="NA")


# Format data
################################################################################

# Format data
colnames(data_orig)
data <- data_orig %>%
  # Simplify
  select(species, "ctmax...2", "site") %>%
  # Rename
  rename("ctmax_c"="ctmax...2") %>%
  # Add lake vs river
  mutate(site_type=ifelse(site !="lake", "river", "lake")) %>%
  # Format CTMAX
  mutate(ctmax_c=as.numeric(ctmax_c)) %>%
  # Remove empty rows
  na.omit()

# Inspect
str(data)

# Sample sizes
count(data, species, site_type)

# Plot quick
g1 <- ggplot(data, aes(x=ctmax_c, y=species, fill=site_type)) +
  geom_boxplot() +
  # Labels
  labs(x=expression("CT"["max"]*" (°C)"), y="", fill="Site type") +
  scale_x_continuous(breaks=seq(24,32,1)) +
  # Theme
  theme_bw() +
  theme(legend.position = "top")
g1

# Environment within species
################################################################################

# Test grayling river vs lake
g_lake <- data %>%
  filter(site_type=="lake" & species=="grayling") %>%
  pull(ctmax_c)
g_river <- data %>%
  filter(site_type=="river" & species=="grayling") %>%
  pull(ctmax_c)
t.test(g_lake, g_river)


# Test lenok river vs lake
l_lake <- data %>%
  filter(site_type=="lake" & species=="lenok") %>%
  pull(ctmax_c)
l_river <- data %>%
  filter(site_type=="river" & species=="lenok") %>%
  pull(ctmax_c)
t.test(l_lake, l_river)




# Lenok power analysis
################################################################################

# Standard deviation from lake
sd_lake <- sd(l_lake)

# Simulation parameters
nsim <- 100
ctmax_vals <- seq(27, 30, 0.5)
nsamples_vals <- seq(10, 40, 1)

# Build simulation dataframe
sim_df <- expand.grid(ctmax_c=ctmax_vals,
                      nsamples=nsamples_vals)

# Loop through simulation dataframe
i <- 1
for(i in 1:nrow(sim_df)){

  # Extract parametrers
  nsample <- sim_df$nsamples[i]
  ctmax_c <- sim_df$ctmax_c[i]

  # Simulate CT max values
  for(j in 1:nsim){

    # Print where we are
    print(paste(i, j, sep="-"))

    # Simulate N CTmax values for river lenok
    ctmax_sim <- rnorm(n=nsample, mean = ctmax_c, sd = sd_lake)

    # Compare simulated river lenok CTmax values to observed lake lenok CTmax values
    test <- t.test(ctmax_sim, l_lake)

    # Extract p-value
    pval <- test$p.value

    # Record
    iter_df <- tibble(ctmax_c=ctmax_c,
                      nsample=nsample,
                      iter=j,
                      pval=pval)

    # Merge
    if(j==1){
      iter_df_combo <- iter_df}
    else{
      iter_df_combo <- bind_rows(iter_df_combo, iter_df)
    }

  }

  # Merge
  if(i==1){
    iter_out <- iter_df_combo}
  else{
    iter_out <- bind_rows(iter_out, iter_df_combo)
  }

}

# Summarize propotion of significant tests in each scenario
stats <- iter_out %>%
  group_by(ctmax_c, nsample) %>%
  summarize(n=n(),
            nsig=sum(pval<=0.05),
            psig=nsig/n) %>%
  ungroup()

# Plot results
g2 <- ggplot(stats, aes(x=nsample, y=psig, color=ctmax_c)) +
  geom_point() +
  # Ref line
  geom_hline(yintercept=0.05, linetype="dotted") +
  geom_vline(xintercept=12) +
  # Labels
  labs(x="Number of river lenok tested",
         y="Proportion of simulations\nwith sig. diff river vs. lake lenok",
         color="CTmax of river lenok") +
  # Legend
  scale_color_gradientn(colors=RColorBrewer::brewer.pal(9, "Spectral") %>% rev()) +
  # Theme
  theme_bw()
g2
g1

