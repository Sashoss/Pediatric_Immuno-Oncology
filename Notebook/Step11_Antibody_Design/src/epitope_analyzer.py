import os
import requests
import pandas as pd


class EpitopeAnalyzer:
    def __init__(self, surface_res, epitope_dir="./in/epitope"):
        self.epitope_dir = epitope_dir
        os.makedirs(self.epitope_dir, exist_ok=True)
        self.freq_df = None
        self.df_subset = None
        self.df_ABC = None
        self.heavy_chain_hits = None
        self.heavy_chain_filtered = None
        self.store_surface_res = surface_res 

    def download_allele_frequency_file(self, url, filename):
        output_path = os.path.join(self.epitope_dir, filename)
        response = requests.get(url, allow_redirects=True)
        with open(output_path, "wb") as f:
            f.write(response.content)
        print(f"File downloaded to: {output_path}")
        return output_path

    def load_and_prepare_frequency_table(self, filepath, select_class=("HLA-A", "HLA-B", "HLA-C")):
        df = pd.read_excel(filepath)
        df = df.drop([0, 1])
        df.rename(columns={"Unnamed: 0": "Locus", "Unnamed: 1": "Allele"}, inplace=True)
        df["Median"] = df.iloc[:, 2:].median(axis=1)
        df["Standard_Allele"] = "HLA-" + df["Allele"].str.replace("g", "", regex=False)
        self.freq_df = df
        self.df_subset = df[df["Median"] > 0.1]  # adjustable threshold
        self.df_ret = self.df_subset[self.df_subset["Standard_Allele"].str.startswith(select_class)]
        return self.df_ret

    def load_heavy_chain_hits(self, filepath):
        self.heavy_chain_hits = pd.read_csv(filepath)
        common_alleles = self.df_ret["Standard_Allele"].unique()
        self.heavy_chain_filtered = self.heavy_chain_hits[self.heavy_chain_hits["Allele"].isin(common_alleles)]
        return self.heavy_chain_filtered

    def filter_epitope_hits(self, ic50_cutoff=150):
        df = self.heavy_chain_filtered[self.heavy_chain_filtered['IC50'] <= ic50_cutoff]
        df = df.drop_duplicates(subset="Peptide", keep="first").copy()
        self.heavy_chain_hits_filtered = df
        return df
    
    def read_classII_table(self, filepath):
        fdata = pd.read_csv(filepath)
        return fdata

    def calculate_surface_exposure(self, epitope_table, chain="H"):
        if not self.store_surface_res[chain]:
            raise ValueError(f"Surface residues for chain {chain} not set.")
        
        percentages = []
        for _, row in epitope_table.iterrows():
            start = int(row["Peptide start"])
            end = int(row["Peptide end"])
            peptide_positions = list(range(start, end + 1))
            exposed_count = sum([1 for pos in peptide_positions if pos in self.store_surface_res[chain]])
            exposed_percent = (exposed_count / len(peptide_positions)) * 100
            percentages.append(exposed_percent)

        epitope_table["% Surface Exposed"] = percentages
        sorted_df = epitope_table.sort_values(by="% Surface Exposed", ascending=False)
        return sorted_df
    
    