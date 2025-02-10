
# Make histogram

library(ggplot2)
library(extrafont)

make_hist <- function(x, symbol="sigma"){
  
  # Convert to data.frame
  df <- as.data.frame(x)
  
  # Construct the expression dynamically
  symbol_expr <- parse(text = paste0('"Estimates of " ~ hat(', symbol, ')^2'))
  # Plot
  ggplot(data=df, aes(x=x)) +
    geom_histogram(aes(y=..density..), color="black", fill="gray80", bins=15) + # lightblue
    geom_density(linewidth=0.5) +
    geom_vline(xintercept = mean(x), linewidth=0.7, linetype=2) + # dashed
    geom_vline(xintercept = median(x), linewidth=0.7, linetype=3) + # dotted
    labs(y="Probability Density",
         x=as.expression(symbol_expr)) +
    #scale_x_continuous(expand = c(0, 0)) + # Remove x-axis padding
    #scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +  # Remove y-axis padding
    theme_light() +
    theme(axis.title.x  = element_text(family = "Times New Roman"),
          axis.title.y  = element_text(family = "Times New Roman")
          ) -> fig
  # Output
  return(fig)
}

# ?plotmath
# x <- rnorm(n=10000)
# make_hist(x)
