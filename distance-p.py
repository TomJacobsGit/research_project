from datasets import data

import plotly.express as px

def main():
    locations = data.sv_df[["sv_id", "start", "end"]]
    p_values = data.snp_sv_df[["sv_id", "POS", "-log10P"]]
    df = locations.merge(p_values, how="inner", left_on="sv_id", right_on="sv_id")

    def calculate_distance(row):
        if row['position'] < row['start']:
            return row['start'] - row['position']
        elif row['position'] > row['end']:
            return row['position'] - row['end']
        else:
            return 0
    
    df["distance"] = df.apply(calculate_distance, axis=1)

    fig = px.scatter(df, x="distance", y="-log10P")

    fig.write_image("distance-p.png")

if __name__ == "main":
    main()