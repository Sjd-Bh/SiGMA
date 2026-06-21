library(ggplot2)
library(dplyr)
library(readr)
library(ggpubr)

############################################################
# 1. Load data
############################################################

mda_sim  <- read_tsv("../1Esfand/MDA_summary.tsv", show_col_types = FALSE)
mda_real <- read_tsv("../1Esfand/real_MDA_summary.tsv", show_col_types = FALSE)
pta_sim  <- read_tsv("../1Esfand/PTA_initi_summary.tsv", show_col_types = FALSE)
pta_real <- read_tsv("../1Esfand/real_PTA_summary.tsv", show_col_types = FALSE)

# Filter real MDA and PTA for samples ending with "chr10.sort"
mda_real <- mda_real %>%
  filter(grepl("chr10\\.sort$", sample))

pta_real <- pta_real %>%
  filter(grepl("chr10\\.sort$", sample))


############################################################
# 2. Add metadata
############################################################

mda_sim  <- mda_sim  %>% mutate(Method = "MDA", Type = "Simulated")
mda_real <- mda_real %>% mutate(Method = "MDA", Type = "Real")

pta_sim  <- pta_sim  %>% mutate(Method = "PTA", Type = "Simulated")
pta_real <- pta_real %>% mutate(Method = "PTA", Type = "Real")

table(mda_real$Method, mda_real$Type)  # Should only show filtered samples
table(pta_real$Method, pta_real$Type)
############################################################
# 3. Merge all
############################################################

df <- bind_rows(mda_sim, mda_real, pta_sim, pta_real)

############################################################
# 4. Convert ADO to percentage
############################################################

df <- df %>%
  mutate(ADO_percent = ADO_rate * 100)

############################################################
# 5. Factor ordering (important for plot layout)
############################################################

df$Method <- factor(df$Method, levels = c("MDA", "PTA"))
df$Type   <- factor(df$Type, levels = c("Simulated", "Real"))

############################################################
# 6. Plot
############################################################
comparisons <- list(
  c("Simulated", "Real")
)

df_means <- df %>%
  group_by(Method, Type) %>%
  summarize(mean_ADO = mean(ADO_percent), .groups = "drop")

p <- ggplot(df, aes(x = Method, y = ADO_percent, fill = Type)) +
  
  # Boxplot
  geom_boxplot(
    width = 0.6,
    position = position_dodge(width = 0.7),
    alpha = 0.9,
    outlier.shape = NA   # hide default outliers (we show raw points)
  ) +
  
  # Jittered points
  geom_jitter(
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width  = 0.7
    ),
    color = "black",
    size = 2,
    alpha = 0.6
  ) +
  
  # Show mean as text on top of boxes
  geom_text(
    data = df_means,
    aes(x = Method, y = mean_ADO + 2, label = sprintf("%.1f", mean_ADO), group = Type),
    position = position_dodge(width = 0.7),
    color = "black",
    size = 4,
    fontface = "bold"
  ) +
  
  # Fill colors
  scale_fill_manual(values = c("Simulated" = "#4C72B0", "Real" = "#DD8452")) +
  
  # Labels
  labs(
    title = "Allelic Dropout (ADO) Comparison",
    subtitle = "MDA vs PTA in Simulated and Real Data (chr10.sort only)",
    x = "Amplification Method",
    y = "ADO (%)"
  ) +
  
  # Theme
  theme_classic(base_size = 14) +
  theme(
    plot.title      = element_text(face = "bold", size = 16),
    plot.subtitle   = element_text(size = 13),
    axis.title      = element_text(face = "bold"),
    axis.text       = element_text(color = "black"),
    legend.title    = element_blank(),
    legend.position = "top"
  ) 
#+
#  
#  stat_compare_means(
#    method = "wilcox.test",
#   aes(group = Type),
#    label = "p.signif"
#  )

############################################################
# 7. Save high-resolution figure
############################################################

ggsave(
  "ADO_boxplot_publication.png",
  p,
  width = 6,
  height = 5,
  dpi = 600
)

ggsave(
  "ADO_boxplot_publication.pdf",
  p,
  width = 6,
  height = 5
)

print(p)