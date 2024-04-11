
library(dplyr)
library(ggplot2)
library(gridExtra)

###########################################################

Pvaluecal <- function(modelsummary) {
  f_value <- modelsummary$fstatistic["value"]
  numdf <- modelsummary$fstatistic["numdf"]
  dendf <- modelsummary$fstatistic["dendf"]
  p <- pf(f_value, numdf, dendf, lower.tail = FALSE)
  return(p)
}

###############################################################

# Create linear models for each environmental factor against each diversity index
create_lm <- function(df, formula) {
  lm_results <- lapply(df, function(x) lm(formula, data = df))
  lm_summaries <- lapply(lm_results, summary)
  return(lm_summaries)
}

##############################################################

# Function to plot median diversity index against soil pH with regression line
plot_diversity_vs_soil_pH <- function(df, y, ylab) {
  ggplot(df, aes(x = X1.1.Soil.pH, y = !!rlang::sym(y))) +
    geom_point(size = 4) +
    geom_smooth(method = "lm", se = FALSE, fullrange = TRUE) +
    labs(x = "Soil pH", y = ylab) +
    theme(
      legend.text = element_text(color = "black", size = 30),
      legend.title = element_text(color = "black", size = 30),
      axis.title.x = element_text(size=30, face="bold"),
      axis.title.y = element_text(size=30, face="bold"),
      plot.title = element_text(size=30, face="bold"),
      axis.text.y = element_text(size=30, face="bold"),
      axis.text.x = element_text(size=30, face="bold"),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black")
    )
}

###########################################################



# Investigate the median alpha diversity in a location vs environmental factors
G2F_2019_median_alpha_diversity_by_location_data <- G2F_metadata_2019 %>%
  group_by(location) %>%
  summarize(
    median_shannon = median(Shannon, na.rm = TRUE),
    median_observed_features = median(Observed, na.rm = TRUE),
    median_simpson = median(Simpson, na.rm = TRUE)
  )

G2f_2019_soil_data <- read.csv("./g2f_2019_soil_data.csv", sep = ",")
G2f_2019_weather_data <- read.csv("./G2F_2019_weather_average_and_sum.tsv", sep = "\t")

# Joining soil and weather data with median alpha diversity data
G2F_2019_median_alpha_diversity_by_location_with_soil_weather_data <- G2F_2019_median_alpha_diversity_by_location_data %>%
  left_join(G2f_2019_soil_data, by = c("location" = "Location")) %>%
  left_join(G2f_2019_weather_data, by = "location") %>%
  filter(!is.na(X1.1.Soil.pH)) %>%
  as.data.frame()


shannon_result_summary <- create_lm(G2F_2019_median_alpha_diversity_by_location_with_soil_weather_data, median_shannon ~ .)
observed_features_result_summary <- create_lm(G2F_2019_median_alpha_diversity_by_location_with_soil_weather_data, median_observed_features ~ .)
simpson_result_summary <- create_lm(G2F_2019_median_alpha_diversity_by_location_with_soil_weather_data, median_simpson ~ .)



# Create plots for each diversity index
median_shannon_ph <- plot_diversity_vs_soil_pH(G2F_2019_median_alpha_diversity_by_location_with_soil_weather_data, "median_shannon", "Median Shannon Index")
median_observed_features_ph <- plot_diversity_vs_soil_pH(G2F_2019_median_alpha_diversity_by_location_with_soil_weather_data, "median_observed_features", "Median Observed Features")
median_simpson_ph <- plot_diversity_vs_soil_pH(G2F_2019_median_alpha_diversity_by_location_with_soil_weather_data, "median_simpson", "Median Simpson Index")

# Combine and save the plots
median_alpha_diversity_plot <- ggarrange(median_shannon_ph, median_simpson_ph, median_observed_features_ph, nrow = 1, ncol = 3)


median_alpha_diversity_plot

ggsave("./median_alpha_diversity_plot.png", plot = median_alpha_diversity_plot, height = 8, width = 24, device = "png")
