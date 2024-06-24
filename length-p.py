from datasets import data
from scipy.stats import spearmanr

import numpy as np
import plotly.express as px

def length():
    lengths = data.sv_df[["sv_id", "length"]]
    p_values = data.snp_sv_df[["sv_id", "P"]].groupby("sv_id").mean()
    df = lengths.merge(p_values, how="inner", left_on="sv_id", right_on="sv_id")

    pearson = np.corrcoef(df["length"], df["P"])
    spearman, p = spearmanr(df["length"], df["P"])

    print(f"Pearson: {pearson}")
    print(f"Spearman: {spearman}, P: {p}")

    fig = px.scatter(df, x="length", y="P", title=f"Correlation between length of an SV and<br>the average p of its interactions with SNPs.")

    fig.write_image("length-p.png")

    lengths = lengths.loc[lengths["length"] > 1]
    df = lengths.merge(p_values, how="inner", left_on="sv_id", right_on="sv_id")

    pearson = np.corrcoef(df["length"], df["P"])
    spearman, p = spearmanr(df["length"], df["P"])

    print(f"Pearson: {pearson}")
    print(f"Spearman: {spearman}, P: {p}")

    fig = px.scatter(df, x="length", y="P", title=f"Correlation between length of an SV and<br>the average p of its interactions with SNPs.")

    fig.write_image("length-p 2.png")

def distance():
    locations = data.sv_df[["sv_id", "start", "end"]]
    p_values = data.snp_sv_df[["sv_id", "POS", "-log10P"]]
    df = locations.merge(p_values, how="inner", left_on="sv_id", right_on="sv_id")

    def calculate_distance(row):
        if int(row["POS"]) < int(row["start"]):
            return int(row["start"]) - int(row["POS"])
        elif int(row["POS"]) > int(row["end"]):
            return int(row["POS"]) - int(row["end"])
        else:
            return 0
    
    df["distance"] = df.apply(calculate_distance, axis=1)

    pearson = np.corrcoef(list(df["distance"]), list(df["-log10P"]))
    spearman, p = spearmanr(df["distance"], df["-log10P"])

    print(f"Pearson: {pearson}")
    print(f"Spearman: {spearman}, P: {p}")

    fig = px.scatter(df, x="distance", y="-log10P", title=f"Correlation between distance between the SNP and SV and<br>the average -log(p) of their interaction.")

    fig.write_image("distance-p.png")

if __name__ == "__main__":
    distance()