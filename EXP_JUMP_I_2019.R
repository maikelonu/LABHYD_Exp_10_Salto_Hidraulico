# //////////////////////////////////////////////////////////////////////////////////
# INSTITUTO TECNOLOGICO DE COSTA RICA
# Escuela de Ingenieria en Construccion
# https://www.tec.ac.cr
# Session: FLUJO NO-UNIFORME/ SALTO HIDRAULICO

# M.Sc. Eng. Maikel Mendez M
# Water Resources + GIS + DataScience
# Instituto Tecnologico de Costa Rica
# https://www.tec.ac.cr
# https://orcid.org/0000-0003-1919-141X
# https://www.scopus.com/authid/detail.uri?authorId=51665581300
# https://scholar.google.com/citations?user=JnmSVFYAAAAJ&hl=en
# https://www.youtube.com/c/maikelmendez
# https://twitter.com/MaikelMendezM
# https://github.com/maikelonu
# //////////////////////////////////////////////////////////////////////////////////

# INFO:
# Analisis grafico avanzado
# ggplot2
# lattice
# Normalizacion y homogenizacion de variables
# Exportaci?n ASCII"
# //////////////////////////////////////////////////////////////////////////////////
# Workspace is cleared
rm(list = ls())

# Working directory is selected
setwd("C:/DATOS/R_ITC/R_LABHYD/EXP_JUMP")

# CRAN libraries are loaded
require(Agreement)
require(DescTools)
require(effects)
require(ggplot2)
require(MASS)
require(nls2)
require(nlstools)
require(pastecs)
require(reshape)
require(visreg)
require(gridExtra)
require(ggalt)

# /////////////////////////////////////////////////////////////
# BLOCK: Custom function, round data.frame to specif digits
# /////////////////////////////////////////////////////////////
round_df <- function(df, digits) {
  options(scipen = 0)
  options(scipen = -2)
  nums <- vapply(df, is.numeric, FUN.VALUE = logical(1))
  df[,nums] <- round(df[,nums], digits = digits)
  return(df)
}

# ////////////////////////////////////////////////////////
# BLOCK: Declarations
# ////////////////////////////////////////////////////////
base_m <- 0.086 # hydraulic flume base (m)
visco <- 1e-06 # water kynematic viscosity (m2/s)
     
# ////////////////////////////////////////////////////////
# BLOCK: Data input
# ////////////////////////////////////////////////////////

# Input data is loaded and a data.frame is created
df.base <- read.table("base.txt", header = TRUE)

# names {base} function is requested
names(df.base)

# An effective distance varible is created (m)
df.base$eff_cota <- df.base$Cota_m - 15.00 # OJO !!! this is a constant

# flow units are converted from m3/h to m3/s
df.base$q_m3_s <- (df.base$Q_m3h) / 3600

# An average water depth variable is created in m
df.base$y_m <- NA

# An averaging loop is created
counter <- length(df.base$Cota_m)

for(i in 1:counter) {
  
  df.base$y_m[i] <- (mean(c(df.base$Yi1_cm[i], df.base$Yi2_cm[i], df.base$Yi3_cm[i]))) / 100
  
}

# Delta Z is substracted from water depth
df.base$y_m <- df.base$y_m - ((df.base$DeltaZ_cm)/100)

# hydraulic area is calculated (m2)
df.base$area <- (df.base$y_m) * base_m

# hydraulic perimeter is calculated (m)
df.base$perimeter <- ((df.base$y_m) * 2) + base_m

# hydraulic radius is calculated
df.base$radius <- (df.base$area / df.base$perimeter)

# square root of hydraulic radius is calculated
df.base$radius.root <- (df.base$radius) ^ 0.5

# water velocity is calculated
df.base$vel <- (df.base$q_m3_s / df.base$area)

# Froude number is calculated
df.base$Froude <- df.base$vel / ((df.base$area * 9.81 / base_m) ^ 0.5)     

# Dynamic energy component is calculated (m)
df.base$dym <- ((df.base$vel)^2)/(2*9.81)

# Total energy is calculated (m)
df.base$energ_total <- df.base$y_m + df.base$dym

# A rule variable is created
df.base$rule <- NA

# If-statement & For-loop for Froude number
for(i in 1:counter) {
  if (df.base$Froude[i] > 1) {
    df.base$rule[i] = "SUPER"
  }
  else {
    df.base$rule[i] = "SUB"
  }
}

# A sequence variable is created
V.SEQ <- (1:counter)
df.base$SEQ <- V.SEQ

# A ggplot object is created
fg01 <- ggplot() +
 geom_point(aes(x = eff_cota,y = y_m),data=df.base,colour = '#0000ff', size = 3.5) +
 geom_point(aes(x = eff_cota,y = energ_total),data=df.base,shape = 8,colour = '#ff0000',size = 3.5) +
 geom_xspline(aes(x = eff_cota,y = y_m),data=df.base, spline_shape=-0.1, size=1.25, colour = '#0000ff') +
 geom_xspline(aes(x = eff_cota,y = energ_total),data=df.base, spline_shape=-0.1, size=1.25, colour = '#ff0000',linetype = 2) +
 scale_y_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
 scale_x_continuous(breaks = scales::pretty_breaks(n = 20.0,min.n = 20.0)) +
 geom_text(aes(x = eff_cota,y = energ_total,label = SEQ,vjust = -0.5),data=df.base,size = 8.0,parse = FALSE) +
 geom_text(aes(x = eff_cota,y = y_m,label = rule,vjust = 1.25),data=df.base,size = 5.5,parse = FALSE) +
 geom_vline(aes(xintercept = df.base$eff_cota[3]),data=df.base,colour = '#666600',linetype = 5, size = 0.75) +
 geom_vline(aes(xintercept = df.base$eff_cota[5]),data=df.base,colour = '#666600',linetype = 5, size = 0.75) +
 ggtitle("Evolucion del salto hidraulico con respecto de la distancia") +
 xlab("Cota-distancia (m)") +
 ylab("Tirante-energia (m)") +
 theme_bw(base_size = 14.0)

# A ggplot object is requested
fg01

# ////////////////////////////////////////////////////////
# BLOCK: Manual Calculations
# ////////////////////////////////////////////////////////

# Select Y1exp
Y1exp <- df.base$y_m[3]

# Select Y2exp
Y2exp <- df.base$y_m[5]

# Y2teor is calculated
Y2teor <- Y1exp*0.5*(-1 + sqrt(1+ (8*(df.base$Froude[3]^2))))

# Experimental energy loss (hl_exp) is calculated (m)
hl_exp <- ((Y2exp - Y1exp)^3) / (4*Y1exp*Y2exp)

# Theoretical energy loss (hl_teor) is calculated (m)
hl_teor <- ((Y2teor - Y1exp)^3) / (4*Y1exp*Y2teor)

# Experimental energy dissipiation is calculated (%)
E_perc_loss_exp <- ((df.base$energ_total[3] - df.base$energ_total[5]) / df.base$energ_total[3])*100

# ////////////////////////////////////////////////////////
# BLOCK: Manual Calculations for Y2teor
# ////////////////////////////////////////////////////////

# hydraulic area is calculated for Y2teor (m2)
area.Y2teor <- (Y2teor) * base_m

# hydraulic perimeter is calculated for Y2teor (m)
perimeter.Y2teor <- ((Y2teor) * 2) + base_m

# hydraulic radius is calculated for Y2teor
radius.Y2teor <- (area.Y2teor / perimeter.Y2teor)

# square root of hydraulic radius is calculated for Y2teor
radius.root.radius.Y2teor <- (radius.Y2teor) ^ 0.5

# water velocity is calculated for Y2teor
vel.Y2teor <- mean(df.base$q_m3_s) / area.Y2teor

# Froude number is calculated for Y2teor
Froude.Y2teor <- vel.Y2teor / ((area.Y2teor * 9.81 / base_m) ^ 0.5)     

# Dynamic energy component is calculated (m) for Y2teor
dym.Y2teor <- ((vel.Y2teor)^2)/(2*9.81)

# Total energy is calculated for for Y2teor ()
Energ.total.Y2teor <- dym.Y2teor + Y2teor

# Theoretical energy dissipiation is calculated (%)
E_perc_loss_teor <- ((df.base$energ_total[3] - Energ.total.Y2teor) / df.base$energ_total[3])*100

# ////////////////////////////////////////////////////////
# BLOCK: Y2/Y1 ratio
# ////////////////////////////////////////////////////////

# Experimental Y2/Y1 ratio is calculated
Y2_Y1_exp <- Y2exp / Y1exp

# Theoretical Y2/Y1 ratio is calculated
Y2_Y1_teor <- Y2teor / Y1exp

# ////////////////////////////////////////////////////////
# BLOCK: Hydraulic Jump Length
# ////////////////////////////////////////////////////////

# Effective Jump Lenght (m). You have to DO this manually!!
L_exp <- abs(df.base$eff_cota[3] - df.base$eff_cota[5])

# You have to compare this with:
# Kim, Y., Choi, G., Park, H., & Byeon, S. 2015. Hydraulic jump and energy dissipation with sluice gate. Water (Switzerland), 7(9), 5115-5133. 
# http://doi.org/10.3390/w7095115

# A compiled dta.frame is created
variable.total <- c("Fr01", "Fr02", "Y1exp", "Y2exp", "Y2teor", "hl_exp", "hl_teor", "E_perc_loss_exp", "E_perc_loss_teor", "Y2_Y1_exp", "Y2_Y1_teor", "L_exp")
valor.total <- c(df.base$Froude[3], df.base$Froude[5], Y1exp, Y2exp, Y2teor, hl_exp, hl_teor, E_perc_loss_exp, E_perc_loss_teor, Y2_Y1_exp, Y2_Y1_teor, L_exp)
df.compile <- data.frame(variable.total,valor.total)

# round_df function is applied to relevant data.frames
df.output <- round_df(df=df.base, digits=3)
df.output.02 <- round_df(df=df.compile, digits=3)

# Objects to export:
# df.output, df.output.02, fg01
write.csv(df.output, file = "df.output.csv")
write.csv(df.output.02, file = "df.output.02.csv")

# /////////////////////////////////////////////////////////////
# END OF SCRIPT
# /////////////////////////////////////////////////////////////
