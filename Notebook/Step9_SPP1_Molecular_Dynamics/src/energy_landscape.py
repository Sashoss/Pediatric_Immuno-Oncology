import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Constants
k_B = 0.0083145  # kJ/mol·K
T = 300          # Kelvin, or use avg temp of your REMD sims

# Load data
df = pd.read_csv("./out/SPP1_replica_exchange/rmsd_rg_all.dat", delim_whitespace=True, header=None)
df.columns = ["Frame", "RMSD", "Rg", "Energy"]

df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["RMSD", "Rg"])
x = df["RMSD"]
y = df["Rg"]

hist, xedges, yedges = np.histogram2d(x, y, bins=100, density=True)

# Avoid log(0)
hist[hist == 0] = np.nan

# Compute Free Energy
F = -k_B * T * np.log(hist)
F = F - np.nanmin(F)  # Shift so minimum F = 0

# Plot
plt.figure(figsize=(6,5))
X, Y = np.meshgrid((xedges[:-1] + xedges[1:]) / 2,
                   (yedges[:-1] + yedges[1:]) / 2)

cp = plt.contourf(X, Y, F.T, levels=30, cmap='viridis')
plt.xlabel("RMSD (nm)", fontsize=16)
plt.ylabel("Radius of Gyration (nm)", fontsize=18)
plt.title("Free Energy Surface (kJ/mol)", fontsize=18)
plt.colorbar(cp, label="ΔG (kJ/mol)")
plt.subplots_adjust(left=0.15, bottom=0.15)

plt.xticks(ha='right', fontsize=12)
plt.yticks(rotation=0, fontsize=12)
plt.tight_layout()
plt.savefig("./out/SPP1_replica_exchange/Rg_RMSD.png", dpi=600)
plt.show()



# Get minima 1
target_rmsd_1 = 1.85
target_rg_1   = 4.09

# Get minima 2
target_rmsd_2 = 2.5
target_rg_2   = 4.09

tolerance   = 0.005 


near_minima_1 = df[(np.abs(df["RMSD"] - target_rmsd_1) < tolerance) &
                 (np.abs(df["Rg"] - target_rg_1) < tolerance)]

near_minima_2 = df[(np.abs(df["RMSD"] - target_rmsd_2) < tolerance) &
                 (np.abs(df["Rg"] - target_rg_2) < tolerance)]

#df.columns = ["Frame", "RMSD", "Rg", "Energy"]


print(near_minima_1[["Frame", "RMSD", "Rg", "Energy"]])

print(near_minima_2[["Frame", "RMSD", "Rg", "Energy"]])