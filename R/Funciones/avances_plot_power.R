
library(tidyr)
library(ggplot2)

A <- output$sph_n100_tau0
Along <- pivot_longer(A, cols= -1, names_to="phi", values_to="power")
Along$phi <- factor(Along$phi, levels = paste0("phi",c(10,50,80,100)))

ggplot(data=Along) +
  geom_line(aes(x=x, y=power, linetype=factor(phi)), linewidth = 0.6) +
  geom_vline(xintercept = 0,    linetype="dashed", linewidth=0.3) +
  geom_hline(yintercept = 0.05, linetype="dashed", linewidth=0.3) +
  scale_linetype_manual(values=c("solid","longdash","dashed", "dotted"), # "dotdash"
                        labels=c(expression(phi==10), expression(phi==50),expression(phi==80), expression(phi==100))) +
  
  theme_light()
