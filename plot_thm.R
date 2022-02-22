
# Calendar year plots -----------------------------------------------------

thm1 <- theme_classic() +
  theme(
    # Format axes
    axis.title = element_text(family="Helvetica", size=16, color="black"),
    axis.text.y = element_text(family="Helvetica", size=14, color="black"),
    axis.text.x = element_text(family="Helvetica", size=14, color="black", angle=45, hjust=1),
    axis.line = element_line(size=0.75),
    axis.ticks = element_line(size=0.75),
    
    # Format legend
    legend.text = element_text(family="Helvetica", size=14, color="black", margin=margin(t=0.25,b=0.25, unit="lines")),
    legend.title = element_text(family="Helvetica", size=16, color="black"),
    legend.title.align = 0.5,
    legend.position = "bottom",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.direction = "horizontal",
    
    # Add space around plot
    plot.margin = unit(c(1, 1, 1, 1), "lines"),
    
    #Format facet title
    strip.text = element_text(family="Helvetica", size=14, color="black"),
    strip.background = element_blank(),
    panel.border = element_rect(size=0.75, fill=NA)
  )


# Other trend plots -------------------------------------------------------

thm2 <- theme_classic() +
  theme(
    # Format axes
    axis.title = element_text(family="Helvetica", size=16, color="black"),
    axis.text.y = element_text(family="Helvetica", size=14, color="black"),
    axis.line = element_line(size=0.75),
    axis.ticks = element_line(size=0.75),
    
    # Format legend
    legend.text = element_text(family="Helvetica", size=14, color="black", margin=margin(t=0.25,b=0.25, unit="lines")),
    legend.title = element_text(family="Helvetica", size=16, color="black"),
    legend.title.align = 0.5,
    legend.position = "bottom",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.direction="horizontal",
    
    # Add space around plot
    plot.margin = unit(c(1, 1, 1, 1), "lines")
  )


# Survival plots ----------------------------------------------------------

thm3 <- theme_classic() +
  theme(
    # Format axes
    axis.title = element_text(family="Helvetica", size=14, color="black"),
    axis.text = element_text(family="Helvetica", size=14, color="black"),
    axis.line = element_line(size=0.75),
    axis.ticks = element_line(size=0.75),
    
    # Format legend
    legend.text = element_text(family="Helvetica", size=14, color="black", margin=margin(t=0.25,b=0.25, unit="lines")),
    legend.title = element_text(family="Helvetica", size=14, color="black"),
    legend.title.align = 0.5,
    legend.position = "bottom",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.direction="horizontal",

    # Add space around plot
    plot.margin = unit(c(1, 1, 1, 1), "lines")
  )

