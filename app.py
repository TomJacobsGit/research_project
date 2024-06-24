# Libraries
from pathlib import Path
import matplotlib
from flask import (Flask, render_template, request)

from datasets import data
from miami_plot import Miami_plot

matplotlib.use('Agg')

# Set path to data
data_path = '/Data'
miami_plot = Miami_plot()

# Initialize the App
app = Flask(__name__)


@app.route('/', methods=["GET"])
def no_plot():
    return render_template("loading_data.html")


@app.route('/snp/', methods=["POST"])
def snp_specify():
    return render_template("viewtype_snp.html", snps=data.snps)


@app.route('/sv/', methods=["POST"])
def sv_specify():
    return render_template("viewtype_sv.html", sv_ids=[sv.sv_id for sv in data.svs])


@app.route('/trait/', methods=["POST"])
def trait_specify():
    return render_template("viewtype_trait.html", traits=data.traits)


@app.route('/region/', methods=["POST"])
def region_specify():
    return render_template("viewtype_region.html")


@app.route('/snp/plot/', methods=["POST"])
def snp_plot():
    window = int(request.form["window"])
    snps = request.form.getlist("selected_snps")
    color = request.form["color"]

    top_df = data.snp_sv_df.loc[data.snp_sv_df["POS"].astype("str").isin(snps)]
    mid_df = data.sv_df.loc[data.sv_df["sv_id"].astype("str").isin(top_df["sv_id"].astype("str").to_list())]
    bot_df = data.trait_df.loc[data.trait_df["POS"].astype("str").isin(snps)]

    plot = miami_plot.plot(top_df, mid_df, bot_df, color, "snp")

    return render_template("plot_snp.html", snps=data.snps, selected_snps=snps, plot=plot)


@app.route('/sv/plot/', methods=["POST"])
def sv_plot():
    window = int(request.form["window"])
    svs = request.form.getlist("selected_svs")
    color = request.form["color"]

    mid_df = data.sv_df.loc[data.sv_df["sv_id"].astype("str").isin(svs)]

    merged_df = mid_df.merge(data.snp_sv_df, left_on='sv_id', right_on='sv_id', how='inner')
    snp_sv_df = merged_df[(merged_df['POS'].astype(int) > merged_df['start'].astype(int) - window) & (merged_df['POS'].astype(int) < merged_df['end'].astype(int) + window)]

    top_df = snp_sv_df.loc[snp_sv_df["sv_id"].astype("str").isin(svs)]
    bot_df = data.trait_df.loc[data.trait_df["POS"].astype("str").isin(top_df["POS"].astype("str").to_list())]

    plot = miami_plot.plot(top_df, mid_df, bot_df, color, "snp")

    return render_template("plot_snp.html", svs=data.svs, selected_svs=svs, plot=plot)


@app.route('/trait/plot/', methods=["POST"])
def trait_plot():
    window = int(request.form["window"])
    traits = request.form.getlist("selected_traits")
    color = request.form["color"]

    bot_df = data.trait_df.loc[data.trait_df["RSID"].astype("str").isin(traits)]
    top_df = data.snp_sv_df.loc[data.snp_sv_df["POS"].astype("str").isin(bot_df["POS"].astype("str"))]
    mid_df = data.sv_df.loc[data.sv_df["sv_id"].astype("str").isin(top_df["sv_id"].astype("str").to_list())]

    plot = miami_plot.plot(top_df, mid_df, bot_df, color, "snp")

    return render_template("plot_trait.html", traits=data.traits, selected_traits=traits, plot=plot)


@app.route('/region/plot/', methods=["POST"])
def region_plot():
    window = int(request.form["window"])
    start = int(request.form["start"])
    end = int(request.form["end"])
    color = request.form["color"]
    snp = request.form.get("single_snp", None)

    mid_df = data.sv_df.loc[(data.sv_df["start"].astype("int") > start) & (data.sv_df["end"].astype("int") < end)]
    svs = mid_df["sv_id"]

    merged_df = mid_df.merge(data.snp_sv_df, left_on='sv_id', right_on='sv_id', how='inner')
    snp_sv_df = merged_df[(merged_df['POS'].astype(int) > merged_df['start'].astype(int) - window) & (merged_df['POS'].astype(int) < merged_df['end'].astype(int) + window)]

    top_df = snp_sv_df.loc[snp_sv_df["sv_id"].astype("str").isin(svs)]
    bot_df = data.trait_df.loc[data.trait_df["POS"].astype("str").isin(top_df["POS"].astype("str").to_list())]

    plot = miami_plot.plot(top_df, mid_df, bot_df, color, "snp")

    return render_template("plot_region.html", svs=list(mid_df["sv_id"]), start=start, end=end, plot=plot)


@app.after_request
def after_request(response):
    response.headers["Cache-Control"] = "no-store, max-age=0"
    return response

# Run the Flask application
if __name__ == '__main__':
    app.run(debug=True)
