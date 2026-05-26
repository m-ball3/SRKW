install.packages("vistime")
library(vistime)


timeline_data <- data.frame(
  event = c("Library Prep and Sequencing", "Chapter 1 + 2 Bioinformatics", "Chapter 1 + 2 Statistical Analysis", "Thesis Writing and Defense"),
  start = as.Date(c("2026-04-01", "2026-06-01", "2026-09-01", "2027-01-01")),
  end = as.Date(c("2026-05-31", "2026-08-31", "2027-03-31", "2027-12-31")),
  year = c("Year 1", "Year 1", "Year 1", "Year 2"),  # 2026=Year1, 2027=Year2
  group = c("Year 1", "Year 1", "Year 1", "Year 2")
)

# Interactive Plotly timeline
vistime(timeline_data)

# Static ggplot2 timeline
gg_vistime(timeline_data)



library(ggplot2)
library(dplyr)

timeline_data <- data.frame(
  event = c("Library Prep\n& Sequencing", 
            "Ch. 1+2\nBioinformatics", 
            "Ch. 1+2\nStats", 
            "Thesis\nWriting & Defense"),
  start = as.Date(c("2026-04-01", "2026-06-01", "2026-09-01", "2027-01-01")),
  end = as.Date(c("2026-05-31", "2026-08-31", "2027-03-31", "2027-12-31")),
  year = c("Year 1", "Year 1", "Year 1", "Year 2")
) %>%
  mutate(
    duration = as.numeric(end - start) / 30.44,  # months
    center = as.Date(start + (end - start)/2)    # middle of each bar
  )

ggplot(timeline_data, aes(x = duration, y = reorder(event, start))) +
  geom_col(fill = c("#2E86AB", "#A23B72", "#F18F01", "#C73E1D"), 
           alpha = 0.8, width = 0.7) +
  geom_vline(xintercept = c(3, 6, 9, 12), linetype = "dashed", alpha = 0.5) +
  labs(x = "Duration (months)", y = "", 
       title = "Master's Thesis Timeline") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 13),
                     breaks = c(3, 6, 9, 12)) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )


library(ggplot2)
library(dplyr)
library(scales)

timeline_data <- data.frame(
  event = c("Library Prep\n& Sequencing", 
            "Ch. 1+2\nBioinformatics", 
            "Ch. 1+2\nStats", 
            "Thesis\nWriting & Defense"),
  start = as.Date(c("2026-04-01", "2026-06-01", "2026-09-01", "2027-01-01")),
  end = as.Date(c("2026-05-31", "2026-08-31", "2027-03-31", "2027-12-31")),
  year = c("Year 1", "Year 1", "Year 1", "Year 2")
)

ggplot(timeline_data, aes(y = reorder(event, start), xmin = start, xmax = end)) +
  geom_rect(aes(xmin = start, xmax = end, ymin = -0.4, ymax = 0.4), 
            fill = c("#2E86AB", "#A23B72", "#F18F01", "#C73E1D"), 
            alpha = 0.8) +
  # Duration labels inside bars
  geom_text(aes(x = start + (end-start)/2, label = paste0(as.numeric(end-start)/30.44, "mo")), 
            hjust = 0.5, size = 3.5, fontface = "bold") +
  # Year labels
  geom_text(aes(x = end + 15, label = year), hjust = 0, size = 4, fontface = "bold") +
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y",
               limits = as.Date(c("2026-03-01", "2027-05-01"))) +
  labs(x = "", y = "", title = "Master's Thesis Timeline") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  guides(x = guide_axis(n.dodge = 2))
