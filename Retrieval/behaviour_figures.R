#######################################################################################
# CFC behaviour figures for 2 and 5 hour
######################################################################################

library(tidyverse)
library(plotrix)

input_path <- "C:/Users/Natalie Meeson/Documents/PhD/CFC_2hr/Fear_data_analysis/Raw_data"

output_path <- "C:/Users/Natalie Meeson/Documents/PhD/Thesis/Results/CFC_chapter/Figures/"

samples <- c(1:36)
outliers <- c("6", "14", "21")
samples <- samples[!samples %in% outliers]

# logfile

log_file_name <- "C:/Users/Natalie Meeson/Documents/PhD/Thesis/Results/CFC_chapter/06.06.21_ANOVAs.txt"
log_description <- "log file for ANOVAs for 2 and 5 hour time points, 2hr exc outliers"

cat(str_glue("\n","Log file for : {log_description}",.sep = "\n"), 
    file = log_file_name, append = TRUE)

# read in data

cfc_data_2hr <- read_delim(paste0(input_path, "/conditioningdata.txt"), delim = "\t") %>%
  filter(RatIdentifier %in% samples)# only needed if have outliers

cfc_2hr_long <- gather(cfc_data_2hr, TimePoint, FreezingPercent, Pre_US:TenMin)

# filter to create 2 datasets for CFC and retrieval

cfc <- c("Pre_US", "Post_US")

retrieval <- c("TwoMin", "TenMin")

cfc_2hr <- filter(cfc_2hr_long, TimePoint %in% cfc)

retrieval_2hr <- filter(cfc_2hr_long, TimePoint %in% retrieval)

# Prep CFC data for plot

cfc_2hr$TimePoint <- ordered(cfc_2hr$TimePoint, levels = c("Pre_US", "Post_US"))

cfc_2hr$Group <- ordered(cfc_2hr$Group, levels = c("No Recall", "Recall", "Extinction"))

cfc_2hr_Sum <- group_by(cfc_2hr, Group, TimePoint) %>%
  summarise(MeanFreezing = mean(FreezingPercent), sem = std.error(FreezingPercent))

# prep retrieval data for plot

retrieval_2hr$TimePoint <- ordered(retrieval_2hr$TimePoint, levels = c("TwoMin", "TenMin"))

retrieval_2hr$Group <- ordered(retrieval_2hr$Group, levels = c("No Recall", "Recall", "Extinction"))

retrieval_2hr_Sum <- group_by(retrieval_2hr, Group, TimePoint) %>%
  summarise(MeanFreezing = mean(FreezingPercent), sem = std.error(FreezingPercent))

# plot

axis_labels <- c("Pre US", "Post US", "Recall", "Extinction")

cfc_2hr_plot <- ggplot() +
  geom_point(aes(x = TimePoint, y = MeanFreezing, col= Group, group = Group), cfc_2hr_Sum,
             size = 6, alpha = 0.8, shape = 16) +
  geom_errorbar(aes(x = TimePoint, ymin = MeanFreezing- sem, ymax = MeanFreezing + sem, col = Group,
                    group = Group), width = 0.1, size = 1.2,alpha = 0.8, cfc_2hr_Sum) +
  geom_line(aes(x = TimePoint, y = MeanFreezing, col= Group, group = Group), cfc_2hr_Sum, size = 2,
            alpha = 0.8) +
  geom_point(aes(x = TimePoint, y = MeanFreezing, col= Group, group = Group), retrieval_2hr_Sum,
             size = 6, alpha = 0.8, shape = 16) +
  geom_errorbar(aes(x = TimePoint, ymin = MeanFreezing- sem, ymax = MeanFreezing + sem, col = Group,
                    group = Group), width = 0.1, size = 1.2, alpha = 0.8, retrieval_2hr_Sum) +
  geom_line(aes(x = TimePoint, y = MeanFreezing, col= Group, group = Group), retrieval_2hr_Sum,
            size = 2, alpha = 0.8) +
  scale_x_discrete(limits = c("Pre_US", "Post_US", "TwoMin", 
                              "TenMin"), labels = function(x) str_wrap(axis_labels, width = 5),
                   name = "Time Point") +
  scale_colour_brewer(palette = "Dark2") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,100), breaks = c(0, 20, 40, 60, 80, 100), name = "Mean Freezing (%)")+ # forces axis at 0 and sets limits 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), # bg of the panel
        plot.background = element_blank(), #bg of plot
        axis.line = element_line(colour = "gray47", size = 1),
        axis.ticks = element_line(size = 1.2),
        axis.ticks.length = unit(0.25, "cm"),
        legend.key=element_blank(),# gets rid of background and grid lines
        legend.background = element_blank(), # get rid of legend bg
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(),
        axis.text.x = element_text(face = "bold", size = 12), # makes axis text bold
        axis.text.y = element_text(face = "bold", size= 12), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 15)))

#ggsave(paste0(output_path, "cfc_2hr.tiff"), plot = cfc_2hr_plot, dpi = 600,
       #width = 23, height = 15, units = "cm")

ggsave(paste0(output_path, "cfc_2hr_exc_outlier.tiff"), plot = cfc_2hr_plot, dpi = 600,
       width = 23, height = 15, units = "cm")

#################
# ANOVAs
#################
library(rstatix)

# Mixed ANOVA- between subjects= condition, within-subjects= timepoint

cfc_2hr$RatIdentifier <- as.factor(cfc_2hr$RatIdentifier)

res.aov <- anova_test(data = cfc_2hr, dv = FreezingPercent, wid = RatIdentifier,
                      between = Group, within = TimePoint)


cat("\n","Mixed ANOVA 2 hour exc outliers","\n" ,"between subjects= condition, within-subjects= timepoint (PreUS, PostUS)", "\n",.sep = "\n", 
    file = log_file_name, append = TRUE)

capture.output(get_anova_table(res.aov), file = log_file_name, append = TRUE)

# independent t-test between groups at 2 mins

cat("\n","Independent t-test at 2 min 2 hour exc outliers","\n", "\n",.sep = "\n", 
    file = log_file_name, append = TRUE)

# filter 2 min data

retrieval_2hr_2min <- filter(retrieval_2hr, TimePoint %in% "TwoMin" & 
                               Group %in% c("Recall", "Extinction")) 

# check homogeneity of variances

cat("\n","Homogeneity of variances","\n", "\n",.sep = "\n", 
    file = log_file_name, append = TRUE)

capture.output(var.test(FreezingPercent ~ Group, data = retrieval_2hr_2min), 
               file = log_file_name, append = TRUE)

# compute t-test

capture.output(t.test(formula = FreezingPercent ~ Group, data = retrieval_2hr_2min), 
               file = log_file_name, append = TRUE)

# within subject t-test 2 to 10 mins

cat("\n","Dependent t-test 2 hour exc outliers 2 to 10 mins","\n", "\n",.sep = "\n", 
    file = log_file_name, append = TRUE)

retrieval_2hr_ext <- filter(retrieval_2hr, TimePoint %in% c("TwoMin", "TenMin") & 
                               Group %in% "Extinction")

capture.output(t.test(formula = FreezingPercent ~ TimePoint, paired = TRUE, data = retrieval_2hr_ext), 
               file = log_file_name, append = TRUE)

#############################################################################################
# 5 hour data
#############################################################################################

input_path <- "C:/Users/Natalie Meeson/Documents/PhD/CFC_5hr/Behaviour/"

cfc_data_5hr <- read_delim(paste0(input_path, "/conditioningdata.txt"), delim = "\t")

cfc_5hr_long <- gather(cfc_data_5hr, TimePoint, FreezingPercent, Pre_US:TenMin)

# filter to create 2 datasets for CFC and retrieval

cfc <- c("Pre_US", "Post_US")

retrieval <- c("TwoMin", "TenMin")

cfc_5hr <- filter(cfc_5hr_long, TimePoint %in% cfc)

retrieval_5hr <- filter(cfc_5hr_long, TimePoint %in% retrieval)

# Prep CFC data for plot

cfc_5hr$TimePoint <- ordered(cfc_5hr$TimePoint, levels = c("Pre_US", "Post_US"))

cfc_5hr$Group <- ordered(cfc_5hr$Group, levels = c("No Recall", "Recall", "Extinction"))

cfc_5hr_Sum <- group_by(cfc_5hr, Group, TimePoint) %>%
  summarise(MeanFreezing = mean(FreezingPercent), sem = std.error(FreezingPercent))

# prep retrieval data for plot

retrieval_5hr$TimePoint <- ordered(retrieval_5hr$TimePoint, levels = c("TwoMin", "TenMin"))

retrieval_5hr$Group <- ordered(retrieval_5hr$Group, levels = c("No Recall", "Recall", "Extinction"))

retrieval_5hr_Sum <- group_by(retrieval_5hr, Group, TimePoint) %>%
  summarise(MeanFreezing = mean(FreezingPercent), sem = std.error(FreezingPercent))

# plot

axis_labels <- c("Pre US", "Post US", "Recall", "Extinction")

cfc_5hr_plot <- ggplot() +
  geom_point(aes(x = TimePoint, y = MeanFreezing, col= Group, group = Group), cfc_5hr_Sum,
             size = 6, alpha = 0.8, shape = 16) +
  geom_errorbar(aes(x = TimePoint, ymin = MeanFreezing- sem, ymax = MeanFreezing + sem, col = Group,
                    group = Group), width = 0.1, size = 1.2,alpha = 0.8, cfc_5hr_Sum) +
  geom_line(aes(x = TimePoint, y = MeanFreezing, col= Group, group = Group), cfc_5hr_Sum, size = 2,
            alpha = 0.8) +
  geom_point(aes(x = TimePoint, y = MeanFreezing, col= Group, group = Group), retrieval_5hr_Sum,
             size = 6, alpha = 0.8, shape = 16) +
  geom_errorbar(aes(x = TimePoint, ymin = MeanFreezing- sem, ymax = MeanFreezing + sem, col = Group,
                    group = Group), width = 0.1, size = 1.2, alpha = 0.8, retrieval_5hr_Sum) +
  geom_line(aes(x = TimePoint, y = MeanFreezing, col= Group, group = Group), retrieval_5hr_Sum,
            size = 2, alpha = 0.8) +
  scale_x_discrete(limits = c("Pre_US", "Post_US", "TwoMin", 
                              "TenMin"), labels = function(x) str_wrap(axis_labels, width = 5),
                   name = "Time Point") +
  scale_colour_brewer(palette = "Dark2") +
  scale_y_continuous(expand = c(0, 0), limits = c(0,100), breaks = c(0, 20, 40, 60, 80, 100), name = "Mean Freezing (%)")+ # forces axis at 0 and sets limits 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), # bg of the panel
        plot.background = element_blank(), #bg of plot
        axis.line = element_line(colour = "gray47", size = 1),
        axis.ticks = element_line(size = 1.2),
        axis.ticks.length = unit(0.25, "cm"),
        legend.key=element_blank(),# gets rid of background and grid lines
        legend.background = element_blank(), # get rid of legend bg
        legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(),
        axis.text.x = element_text(face = "bold", size = 12), # makes axis text bold
        axis.text.y = element_text(face = "bold", size= 12), 
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = "bold", size = 14, margin = margin(r = 15)))

ggsave(paste0(output_path, "cfc_5hr.tiff"), plot = cfc_5hr_plot, dpi = 600,
       width = 23, height = 15, units = "cm")

#################
# ANOVAs
#################

# Mixed ANOVA- between subjects= condition, within-subjects= timepoint

cfc_5hr$RatIdentifier <- as.factor(cfc_5hr$RatIdentifier)

res.aov <- anova_test(data = cfc_5hr, dv = FreezingPercent, wid = RatIdentifier,
                      between = Group, within = TimePoint)

cat("\n","Mixed ANOVA 5 hour outliers","\n" ,"between subjects= condition, within-subjects= timepoint (PreUS, PostUS)", "\n",.sep = "\n", 
    file = log_file_name, append = TRUE)

capture.output(get_anova_table(res.aov), file = log_file_name, append = TRUE)

# independent t-test between groups at 2 mins

cat("\n","Independent t-test at 2 min 5 hour","\n", "\n",.sep = "\n", 
    file = log_file_name, append = TRUE)

# filter 2 min data

retrieval_5hr_2min <- filter(retrieval_5hr, TimePoint %in% "TwoMin" & 
                               Group %in% c("Recall", "Extinction")) 

# compute t-test

capture.output(t.test(formula = FreezingPercent ~ Group, data = retrieval_5hr_2min), 
               file = log_file_name, append = TRUE)

# within subject t-test 2 to 10 mins

cat("\n","Dependent t-test 5 hour  2 to 10 mins","\n", "\n",.sep = "\n", 
    file = log_file_name, append = TRUE)

retrieval_5hr_ext <- filter(retrieval_5hr, TimePoint %in% c("TwoMin", "TenMin") & 
                              Group %in% "Extinction")

capture.output(t.test(formula = FreezingPercent ~ TimePoint, paired = TRUE, data = retrieval_5hr_ext), 
               file = log_file_name, append = TRUE)