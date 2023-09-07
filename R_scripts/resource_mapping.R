---
title: "Untitled"
author: "Christopher O'Connor"
date: "`r Sys.Date()`"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r}
# Create dataframe with processes used in pre-processing scRNAseq analysis
metrics <- data.frame(
  Process = c(
    "Load Module", "Import data", "Normalize", "Find Variables",
    "Scale", "PCA", "Neighbours", "Cluster/Leiden", "TSNE", "UMAP"
  ))

# Define elapsed time from timer_hook (R) and %%time function (Py)
seurat_time_sec <- c(4.9, 0.85, 0.86, 2.23, 12.83, 7.21, 2.8, 1.5, 32.47, 17.13)
scanpy_time_sec <- c(0.1, 0.8, 0.33, 2.0, 3.1, 9.5, 2.13, 0.3, 55.2, 14.7)

# Add times to dataframe
metrics$Seurat <- seurat_time_sec
metrics$Scanpy <- scanpy_time_sec
```

```{r}
library(ggplot2)
library(RColorBrewer)

metrics$Process <- factor(metrics$Process, levels = metrics$Process)
metrics$Process <- factor(metrics$Process, levels = rev(levels(metrics$Process)))

metrics_long <- tidyr::gather(metrics, Package, Time, -Process)
breaks_sequence <- seq(0, max(metrics_long$Time), by = 5)

ggplot(metrics_long, aes(x = Time, y = Process, fill = Package)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  labs(
    title = "Comparison of Process Times for Seurat and Scanpy Packages",
    x = "Time (Seconds)",
    y = "Process",
    fill = "Package"
  ) +
  scale_fill_brewer(palette = "Set2") +  # Colourblind friendly palette
  scale_x_continuous(breaks = breaks_sequence) +
  theme_minimal()
```

```{r}
# Using peakRAM module add metrics for RAM usage
seurat_Peak_RAM_Used_MiB = c(
                            49.5, 612.7, 517.1, 31.9,
                            4350.2, 634.4, 49.4, 10.8, 
                            90.3, 164.5)

# Using memory_profiler module add metrics for RAM usage
scanpy_Peak_RAM_Used_MiB = c(
                            168.2, 535.3, 536.3, 709.1, 709.1,
                            709.2, 707.9, 386.2, 1255.1, 1438.1,
                            1405.4, 1408.6, 1407.3, 1416.3) 

# Calculate average ram usage 
seu_average_peak_ram <- round(mean(seurat_Peak_RAM_Used_MiB), 2)
sc_average_peak_ram <- round(mean(scanpy_Peak_RAM_Used_MiB), 2)

RAM <- data.frame(
  Package = c("Seurat", "Scanpy"),
  "Average_Peak_RAM_Usage" = c(seu_average_peak_ram, sc_average_peak_ram )
)

# Plot
ggplot(RAM, aes(x = Package, y = Average_Peak_RAM_Usage, fill = Package)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(aes(label = Average_Peak_RAM_Usage), vjust = -0.5, size = 4) +
  labs(
    title = "Performance Comparison: Average Peak RAM Usage",
    x = "Package",
    y = "Average Peak RAM Usage",
    fill = "Package") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")
```

```{r}
seu_total_time <- 96 # Obtained with system.time() for R
sc_total_time <- 88 # Applied using %%time from the time module for Py.

Time <- data.frame(
  Package = c("Seurat", "Scanpy"),
  "Total_Elapsed_Time" = c(seu_total_time, sc_total_time ))

# Plot
ggplot(Time, aes(x = Package, y = Total_Elapsed_Time, fill = Package)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_text(aes(label = Total_Elapsed_Time), vjust = -0.5, size = 4) +
  labs(
    title = "Performance Comparison: Total Elapsed Time",
    x = "Package",
    y = "Total Elapsed Time (Seconds)",
    fill = "Package") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")
```
