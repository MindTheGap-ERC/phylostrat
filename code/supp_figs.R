load("converged_runs.RData")
library(scales)
p = df_converged |>
  group_by(analysis, nchars) |>
  summarise(n_converged = sum(converged), .groups = "drop") |>
  mutate(nchars = as.factor(nchars)) |>
  ggplot(aes(x = analysis, y = n_converged, color = nchars)) +
  geom_point(position = position_jitter(height = 0, width = 0.1)) +
  scale_y_continuous(
    breaks = breaks_width(1),
    labels = as.integer
  )
p
ggsave(filename = "figs/converged_runs.png",
       plot = p)
