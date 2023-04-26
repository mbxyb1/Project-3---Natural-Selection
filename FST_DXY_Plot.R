# Read in the first ENA Flatline file
LABvsNEN <- read.table("LABvsNEN15_BPM.txt", header = TRUE)

# Read in the second ENA Flatline file
LABvsODN <- read.table("LABvsODN15_BPM.txt", header = TRUE)

# Read in the third ENA Flatline file
NENvsODN <- read.table("NENvsODN15_BPM.txt", header = TRUE)

library(ggplot2)

# Plot the FST values for each file on one plot
ggplot(LABvsNEN, aes(x = start, y = dxy, color=dxy)) +
  geom_point() +
  labs(x = "start", y = "DXY", title = "Comparison of DXY values - LABvsNEN") +
  theme_minimal()

