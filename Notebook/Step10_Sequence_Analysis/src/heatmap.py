import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns



def sort_key(label):
    num = ''.join([c for c in label if c.isdigit()])
    return int(num) if num.isdigit() else 0

heatmap_data = heatmap_data.reindex(sorted(heatmap_data.index, key=sort_key))
annot = heatmap_data.applymap(lambda x: "*" if pd.notnull(x) and x > 0.8 else "")
plt.figure(figsize=(12, 8))

ax = sns.heatmap(heatmap_data, 
                 annot=annot,       
                 fmt='',          
                 cmap='viridis',    
                 cbar_kws={'label': 'Phosphorylation Score'},
                 linewidths=0.5,
                 linecolor='white',
                 square=True)

ax.set_xlabel("Phosphorylating Protein (Kinase)", fontsize=14)
ax.set_ylabel("Residue Position", fontsize=14)
ax.set_title("Phosphorylation Prediction Heatmap", fontsize=16, fontweight='bold')

plt.xticks(rotation=45, ha="right", fontsize=12)
plt.yticks(fontsize=12)

plt.tight_layout()

plt.savefig("phospho_heatmap.png", dpi=300)

plt.show()
