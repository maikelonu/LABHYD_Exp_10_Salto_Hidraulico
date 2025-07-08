import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# Introduction
print("\n//////////////////////////////////////////////////////////////////////////////////")
print("INSTITUTO TECNOLÓGICO DE COSTA RICA")
print("Escuela de Ingeniería en Construcción")
print("https://www.tec.ac.cr")
print("Session: FLUJO NO-UNIFORME/ SALTO HIDRÁULICO")

print("\nM.Sc. Eng. Maikel Méndez M")
print("Water Resources + GIS + DataScience")
print("Instituto Tecnológico de Costa Rica")
print("https://www.tec.ac.cr")
print("https://orcid.org/0000-0003-1919-141X")
print("https://www.scopus.com/authid/detail.uri?authorId=51665581300")
print("https://scholar.google.com/citations?user=JnmSVFYAAAAJ&hl=en")
print("https://www.youtube.com/c/maikelmendez")
print("https://github.com/maikelonu")
print("//////////////////////////////////////////////////////////////////////////////////")

print("\n# INFO:")
print("Análisis gráfico avanzado")
print("Normalización y homogenización de variables")
print("Exportación ASCII")
print("//////////////////////////////////////////////////////////////////////////////////")

# /////////////////////////////////////////////////////////////
# BLOCK: Custom function, round DataFrame to specif digits
# /////////////////////////////////////////////////////////////
def round_df(df, digits):
    return df.round(decimals=digits)

print("\n////////////////////////////////////////////////////////")
print("BLOCK: Declarations")
print("////////////////////////////////////////////////////////\n")
base_m = 0.086  # Hydraulic flume base (m)
print("base_m <- 0.086")
visco = 1e-06  # Water kinematic viscosity (m2/s)
print("visco <- 1e-06")

print("\n////////////////////////////////////////////////////////")
print("BLOCK: Data input")
print("////////////////////////////////////////////////////////\n")

# Set the working directory
os.chdir("/home/shoe/Dropbox/Academics/LAB_Esencial/PYTHON_LABHYD_Exp/PYTHON_LABHYD_Exp_10")
print("\n# Working directory is selected:")
print(os.getcwd())

print("# Input data is loaded and a data.frame is created\n")
df_base = pd.read_csv("base.txt", sep="\t")

# Create an effective distance variable (m)
df_base['eff_cota'] = df_base['Cota_m'] - 15.00  # Constant value subtracted

# Convert flow units from m3/h to m3/s
df_base['q_m3_s'] = df_base['Q_m3h'] / 3600

# Calculate the average water depth (m)
df_base['y_m'] = df_base[['Yi1_cm', 'Yi2_cm', 'Yi3_cm']].mean(axis=1) / 100

# Subtract Delta Z from water depth
df_base['y_m'] = df_base['y_m'] - (df_base['DeltaZ_cm'] / 100)

# Calculate hydraulic area (m2)
df_base['area'] = df_base['y_m'] * base_m

# Calculate hydraulic perimeter (m)
df_base['perimeter'] = (df_base['y_m'] * 2) + base_m

# Calculate hydraulic radius
df_base['radius'] = df_base['area'] / df_base['perimeter']

# Calculate square root of hydraulic radius
df_base['radius_root'] = np.sqrt(df_base['radius'])

# Calculate water velocity (m/s)
df_base['vel'] = df_base['q_m3_s'] / df_base['area']

# Calculate Froude number
df_base['Froude'] = df_base['vel'] / np.sqrt(df_base['area'] * 9.81 / base_m)

# Calculate dynamic energy component (m)
df_base['dym'] = (df_base['vel'] ** 2) / (2 * 9.81)

# Calculate total energy (m)
df_base['energ_total'] = df_base['y_m'] + df_base['dym']

# Create a rule column for Froude number classification (SUPER/SUB)
df_base['rule'] = np.where(df_base['Froude'] > 1, 'SUPER', 'SUB')

# Create a sequence variable
df_base['SEQ'] = np.arange(1, len(df_base) + 1)

# ////////////////////////////////////////////////////////
# BLOCK: Plotting the data
# ////////////////////////////////////////////////////////

plt.figure(figsize=(10, 6))

# Plot y_m and energ_total
plt.scatter(df_base['eff_cota'], df_base['y_m'], color='#0000ff', s=50, label='Water Depth')
plt.scatter(df_base['eff_cota'], df_base['energ_total'], color='#ff0000', s=50, marker='x', label='Total Energy')

# Draw line connecting points
plt.plot(df_base['eff_cota'], df_base['y_m'], color='#0000ff', linestyle='-', linewidth=1.25)
plt.plot(df_base['eff_cota'], df_base['energ_total'], color='#ff0000', linestyle='--', linewidth=1.25)

# Add labels for SEQ and rule
for i in range(len(df_base)):
    plt.text(df_base['eff_cota'].iloc[i], df_base['energ_total'].iloc[i],
             df_base['SEQ'].iloc[i], fontsize=10, color='black', ha='center', va='bottom')
    plt.text(df_base['eff_cota'].iloc[i], df_base['y_m'].iloc[i],
             df_base['rule'].iloc[i], fontsize=10, color='black', ha='center', va='top')

# Add vertical lines
plt.axvline(df_base['eff_cota'].iloc[4], color='#666600', linestyle='--', linewidth=0.75)
plt.axvline(df_base['eff_cota'].iloc[4], color='#666600', linestyle='--', linewidth=0.75)

# Labels and title
plt.title("Evolución del salto hidráulico con respecto de la distancia")
plt.xlabel("Cota-distancia (m)")
plt.ylabel("Tirante-energia (m)")

plt.grid(True)
plt.legend()
plt.show()
print("\nPlot: 'Hydraulic Jump Evolution with Distance' has been generated.")

# ////////////////////////////////////////////////////////
# BLOCK: Manual Calculations
# ////////////////////////////////////////////////////////

# Select Y1exp
Y1exp = df_base['y_m'].iloc[2]

# Select Y2exp
Y2exp = df_base['y_m'].iloc[4]

# Y2teor is calculated
Y2teor = Y1exp * 0.5 * (-1 + np.sqrt(1 + (8 * (df_base['Froude'].iloc[2] ** 2))))

# Experimental energy loss (hl_exp) is calculated (m)
hl_exp = ((Y2exp - Y1exp) ** 3) / (4 * Y1exp * Y2exp)

# Theoretical energy loss (hl_teor) is calculated (m)
hl_teor = ((Y2teor - Y1exp) ** 3) / (4 * Y1exp * Y2teor)

# Experimental energy dissipiation is calculated (%)
E_perc_loss_exp = ((df_base['energ_total'].iloc[2] - df_base['energ_total'].iloc[4]) / df_base['energ_total'].iloc[2]) * 100

# ////////////////////////////////////////////////////////
# BLOCK: Manual Calculations for Y2teor
# ////////////////////////////////////////////////////////

# Hydraulic area is calculated for Y2teor (m2)
area_Y2teor = Y2teor * base_m

# hydraulic perimeter is calculated for Y2teor (m)
perimeter_Y2teor = (Y2teor * 2) + base_m

# Hydraulic radius is calculated for Y2teor
radius_Y2teor = area_Y2teor / perimeter_Y2teor

# Square root of hydraulic radius is calculated for Y2teor
radius_root_radius_Y2teor = (radius_Y2teor) ** 0.5

# Water velocity is calculated for Y2teor
vel_Y2teor = df_base['q_m3_s'].mean() / area_Y2teor

# Froude number is calculated for Y2teor
Froude_Y2teor = vel_Y2teor / np.sqrt(area_Y2teor * 9.81 / base_m)

# Dynamic energy component is calculated (m) for Y2teor
dym_Y2teor = (vel_Y2teor ** 2) / (2 * 9.81)

# Total energy is calculated for for Y2teor ()
Energ_total_Y2teor = dym_Y2teor + Y2teor

# Theoretical energy dissipiation is calculated (%)
E_perc_loss_teor = ((df_base['energ_total'].iloc[2] - Energ_total_Y2teor) / df_base['energ_total'].iloc[2]) * 100

# ////////////////////////////////////////////////////////
# BLOCK: Y2/Y1 ratio
# ////////////////////////////////////////////////////////

# Experimental Y2/Y1 ratio is calculated
Y2_Y1_exp = Y2exp / Y1exp

# Theoretical Y2/Y1 ratio is calculated
Y2_Y1_teor = Y2teor / Y1exp

# ////////////////////////////////////////////////////////
# BLOCK: Hydraulic Jump Length
# ////////////////////////////////////////////////////////

# Effective Jump Lenght (m). You have to DO this manually!!
L_exp = abs(df_base['eff_cota'].iloc[2] - df_base['eff_cota'].iloc[4])

# You have to compare this with:
# Kim, Y., Choi, G., Park, H., & Byeon, S. 2015. Hydraulic jump and energy dissipation with sluice gate. Water (Switzerland), 7(9), 5115-5133. 
# http://doi.org/10.3390/w7095115

# ////////////////////////////////////////////////////////
# BLOCK: Export results
# ////////////////////////////////////////////////////////

# Compile results into a DataFrame for export
variables_total = ['Fr01', 'Fr02', 'Y1exp', 'Y2exp', 'Y2teor', 'hl_exp', 'hl_teor', 'E_perc_loss_exp', 'E_perc_loss_teor', 'Y2_Y1_exp', 'Y2_Y1_teor', 'L_exp']
values_total = [df_base['Froude'].iloc[2], df_base['Froude'].iloc[4], Y1exp, Y2exp, Y2teor, hl_exp, hl_teor, E_perc_loss_exp, E_perc_loss_teor, Y2_Y1_exp, Y2_Y1_teor, L_exp]
df_compile = pd.DataFrame({'variable_total': variables_total, 'value_total': values_total})

# Round the DataFrames and export to CSV
df_output = round_df(df_base, digits=3)
df_output_02 = round_df(df_compile, digits=3)

df_output.to_csv("df.output.csv", index=False)
df_output_02.to_csv("df.output.02.csv", index=False)

print("\nData exported successfully.")

print("\n//////////////////////////////////////////////////////////////")
print("END OF SCRIPT")
print("/////////////////////////////////////////////////////////////\n")
