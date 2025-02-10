

# Librerias
library(readxl)
library(ggplot2)
library(tidyr)   # to pivot_longer function

# Funcion
make_plot <- function(read_file, export_file, type_label="two.sample"){
  
  # Cargar datos
  sim <- read_excel(path=read_file)
  sim <- pivot_longer(sim, cols= -1, names_to="phi", values_to="power")
  #sim$phi <- factor(sim$phi)
  
  # X-axis label
  if(type_label=="one.sample"){xlabel = expression(mu[0])}
  if(type_label=="two.sample"){xlabel = expression(mu[2])}

  # Grafico
  fig <- ggplot(data=sim) +
    geom_line(aes(x=x, y=power, linetype=phi), linewidth = 0.6) +
    geom_vline(xintercept = 0,    linetype="dashed", linewidth=0.3) +
    geom_hline(yintercept = 0.05, linetype="dashed", linewidth=0.3) +
    annotate(geom="text", x=-1.7, y=0.08, label=expression(alpha==0.05)) +
    scale_linetype_manual(values=c("solid","longdash","dashed","dotdash","dotted"), 
                          labels=c("No correlation", expression(phi==10), expression(phi==30),expression(phi==60), expression(phi==90))) +
    labs(x = xlabel, y = expression(Pr~(~Reject ~ H[0])), title = "") +
    theme_light() +
    theme(
      # Legend
      legend.position = c(0.95, 0.5),
      legend.justification = c("right", "top"),
      legend.box.background = element_rect(color="black", linewidth=1),
      legend.title = element_blank(),
      legend.key.size =  unit(0.5, "cm"),
      legend.key.width = unit(1.1, "cm"), # largo de las lineas en la leyenda
      # Plot margin
      plot.margin = margin(t = 10,  # Top    margin
                           r = 20,  # Right  margin
                           b = 7,   # Bottom margin
                           l = 10), # Left   margin
      # Title
      plot.title = element_blank()
    )
  
  # Guardar grafico
  ggsave(filename=export_file, plot=fig, width=8, height=6, units="in")
  
  # Output
  return(fig)
  
  
} # fin de la funcion!

# Ejemplo
# read_file="C:/Users/driao/Documents/ARTICULOS/One and two sample test for geostatistical data/Cosas viejas/sim_exp_n50.xlsx"
# export_file="C:/Users/driao/Documents/ARTICULOS/One and two sample test for geostatistical data/Resultados/prueba2.pdf"
# 
# make_plot(read_file, export_file, type_label="one.sample")
