
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(Polychrome)
library(scales)
library(gridExtra)
library(patchwork)
library(reshape2)
library(grid)
library(vegan)


setwd("/Users/bross/Desktop/AIMS/Analysis/Measurements/SCF049")

data <- read_excel("/Users/bross/Desktop/AIMS/Analysis/Measurements/SCF049/SCF049_measures.xlsx", sheet = "Final")

set.seed(58584)
P50 <- createPalette(50,  c("#ff0000", "#00ff00", "#0000ff"))
names(P50)<-NULL

data <- data %>%
  mutate(Length = ifelse(Length >= Width, Length, Width),
         Width = ifelse(Length >= Width, Width, Length))

data$Length <- data$Length/1000
data$Width <- data$Width/1000

mean_points <- data %>%
  group_by(Species) %>%
  summarize(Mean_Length = mean(Length), Mean_Width = mean(Width), .groups = "drop")

data <- data %>%
  mutate(Species = case_when(
    Species == "SCF055" ~ "C. proliferum",
    Species == "SCF049" ~ "C. occultum",
    TRUE ~ Species  # Keep other values unchanged
  ))

mean_points <- mean_points %>%
  mutate(Species = case_when(
    Species == "SCF055" ~ "C. proliferum",
    Species == "SCF049" ~ "C. occultum",
    TRUE ~ Species  # Keep other values unchanged
  ))


p1 <- ggplot(data, aes(x = data$Length, y = data$Width, shape = data$Photo, color = data$Species)) +
  geom_point(size = 3, alpha = 0.5) +  # Plot points with triangles
  scale_shape_manual(values = c(17, 18, 19, 20, 21)) +  # Set shape to triangle
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_point(data = mean_points, 
             aes(x = Mean_Length, y = Mean_Width), 
             size = 5, 
             shape = 21, 
             fill = "black", 
             color = "black") +
  geom_text(data = mean_points, 
            aes(x = Mean_Length, y = Mean_Width, label = Species), 
            hjust = -0.5, # Adjust vertical position of labels
            color = "black", 
            fontface = "bold",
            size = 5) +
  guides(shape = guide_legend(title = NULL), 
         color = guide_legend(title = NULL)) +
  labs(x = expression("Length (" * mu * "m)"), y= expression("Width (" * mu * "m)"), color = "Photo") +# Add labels to axes and legend
  ggtitle("") + # Add title 
theme(
  axis.title = element_text(size = 14),# <-- Axis title font size
  axis.text = element_text(size = 12),        # Axis tick label size
  legend.text = element_text(size = 14),      # Legend label size
  legend.title = element_text(size = 13) )

p1

ggsave("measurements_scatter.pdf", plot =p1, device = "pdf")
ggsave("measurements_scatter.svg", plot =p1, device = "svg")

occ <- seq_along(data$Length)
data$Snap <- as.character(data$Snap)

# Identify the positions where Snap changes
sep <- which(data$Snap[-1] != data$Snap[-length(data$Snap)]) + 1

# Append the last index of the dataset to capture the final segment
change_points <- c(sep, length(data$Snap))

p2 <- ggplot (data, aes(x = occ, y=data$Length, color = data$Snap)) +
  geom_point() +
  geom_smooth() +
  scale_fill_manual(values = P50)+
  geom_vline(xintercept = sep, linetype = "dashed", color = "red")  

ggsave("measurements_photos_fittinglines.pdf", plot =p2, device = "pdf")

#########################################

mean_lenghh <- data %>%
  group_by(Species) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE))
sd_length <- data %>%
  group_by(Species) %>%
  summarize(sd_length = sd(Length, na.rm = TRUE))
mean_width <- data %>%
  group_by(Species) %>%
  summarize(mean_length = mean(Width, na.rm = TRUE))
sd_width <- data %>%
  group_by(Species) %>%
  summarize(sd_lenght = sd(Width, na.rm = TRUE))

#######################################################################


stati <- aggregate(cbind(Length, Width) ~ Snap + Species, data,  
                   FUN = function(x) c(mean = mean(x), sd = sd(x)))

p3 <- ggplot(stati, aes(x = stati$Length[,"mean"] , y = stati$Width[,"mean"], shape = stati$Snap, color = stati$Species)) +
  geom_point(shape = 17, size = 3) +  # Plot means with triangle shape
  geom_errorbar(aes(ymin = Width[,"mean"] - Width[,"sd"], ymax =Width[,"mean"] + Width[,"sd"]), 
                width = 0.1) +  # Error bars for width
  geom_errorbarh(aes(xmin = Length[,"mean"] - Length[,"sd"], xmax = Length[,"mean"] + Length[,"sd"]), 
                 height = 0.1) +  # Error bars for length
  geom_smooth(method = "lm", aes(group = Species), se = FALSE) +  # Add fitting lines
  labs(x = "Mean Length", y = "Mean Width") +  # Add labels to axes
  ggtitle("Scatter Plot of Mean Length vs Mean Width by Photo")  # Add title


ggsave("means_sd_species_trends_measures.pdf", plot =p3, device = "pdf")

########################################################################

# Define the bin width
bin_width <- 1

# Calculate the breaks for the histogram
breaks <- seq(min(data$Length), max(data$Length) + bin_width, by = bin_width)

# Create a histogram to count the frequency of length values in each bin
hist_data <- hist(data$Length, breaks = breaks, plot = FALSE)

# Create a data frame with the histogram data
hist_df <- data.frame(length = hist_data$breaks[-length(hist_data$breaks)], 
                      frequency = hist_data$counts)

# Plot the bar plot
ggplot(hist_df, aes(x = length, y = frequency)) +
  geom_bar(stat = "identity", width = bin_width - 0.1, fill = "skyblue") +
  labs(x = "Length", y = "Frequency") +
  ggtitle("Bar Plot of Length vs Frequency") +
  theme_minimal()

##############################################################################################################################################
##############################################################################################################################################
##############################################################################################################################################

# Calculate statistics grouped by Species
summary_stats <- data %>%
  group_by(Species) %>%
  summarize(
    avg_length = mean(Length, na.rm = TRUE),
    sd_length = sd(Length, na.rm = TRUE),
    sem_length = sd_length / sqrt(n()),
    ci_length = qt(0.975, df = n() - 1) * sem_length,
    lower_ci_length = avg_length - ci_length,
    upper_ci_length = avg_length + ci_length,
    avg_width = mean(Width, na.rm = TRUE),
    sd_width = sd(Width, na.rm = TRUE),
    sem_width = sd_width / sqrt(n()),
    ci_width = qt(0.975, df = n() - 1) * sem_width,
    lower_ci_width = avg_width - ci_width,
    upper_ci_width = avg_width + ci_width
  )

# Print the results
print(summary_stats)

# Optionally, print in a more detailed way
summary_stats %>%
  rowwise() %>%
  do({
    cat("Species:", .$Species, "\n")
    cat("Average length:", .$avg_length, "\n")
    cat("95% CI for length:", .$lower_ci_length, "-", .$upper_ci_length, "\n")
    cat("Average width:", .$avg_width, "\n")
    cat("95% CI for width:", .$lower_ci_width, "-", .$upper_ci_width, "\n")
    cat("Average length (CI):", .$avg_length, "+-", .$ci_length, "\n")
    cat("Average width (CI):", .$avg_width, "+-", .$ci_width, "\n\n")
  })

######################
z_value <- qnorm(0.975)

# Calculate summary statistics grouped by Species
summary_stats <- data %>%
  group_by(Species) %>%
  summarize(
    avg_length = mean(Length, na.rm = TRUE),
    sd_length = sd(Length, na.rm = TRUE),
    sem_length = sd_length / sqrt(n()),
    ci_length = z_value * sem_length,
    lower_ci_length = avg_length - ci_length,
    upper_ci_length = avg_length + ci_length,
    avg_width = mean(Width, na.rm = TRUE),
    sd_width = sd(Width, na.rm = TRUE),
    sem_width = sd_width / sqrt(n()),
    ci_width = z_value * sem_width,
    lower_ci_width = avg_width - ci_width,
    upper_ci_width = avg_width + ci_width
  )

# Print the results
print(summary_stats)

# Optionally, print in a more detailed way
summary_stats %>%
  rowwise() %>%
  do({
    cat("Species:", .$Species, "\n")
    cat("Average length:", .$avg_length, "\n")
    cat("95% CI for length:", .$lower_ci_length, "-", .$upper_ci_length, "\n")
    cat("Average width:", .$avg_width, "\n")
    cat("95% CI for width:", .$lower_ci_width, "-", .$upper_ci_width, "\n")
    cat("Average length (CI):", .$avg_length, "+-", .$ci_length, "\n")
    cat("Average width (CI):", .$avg_width, "+-", .$ci_width, "\n\n")
  })


SCF049 <- subset(data, Species == "SCF049")
SCF055 <- subset(data, Species == "SCF055")

t_test_length <- t.test(SCF049$Length, SCF055$Length)
t_test_length
t_test_width <- t.test(SCF049$Width, SCF055$Width)
t_test_width
