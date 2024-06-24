from enum import StrEnum
import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.io as pio
import ipywidgets as widgets

from datasets import StructuralVariant
from plotly.subplots import make_subplots


class Miami_plot:
    def create_color_map(self, top_df: pd.DataFrame, mid_df: pd.DataFrame, bot_df: pd.DataFrame, color: str):
        if color == "snp":
            h_values = np.linspace(0, 360, len(top_df))
            np.random.shuffle(h_values)
            colors = ['hsl(' + str(h_values[i]) + ',50%' + ',50%)' for i in range(len(h_values))]

            top_df["color"] = colors
            bot_df["color"] = bot_df.merge(right=top_df, left_on="POS", right_on="POS", how="left")["color"]
            mid_df["color"] = mid_df.merge(right=top_df, left_on="sv_id", right_on="sv_id", how="left").drop_duplicates(subset=["sv_id"])["color"]  
        if color == "sv":
            h_values = np.linspace(0, 300, len(mid_df))
            np.random.shuffle(h_values)
            colors = ['hsl(' + str(h_values[i]) + ',50%' + ',50%)' for i in range(len(h_values))]

            mid_df["color"] = colors
            top_df["color"] = top_df.merge(right=mid_df, left_on="sv_id", right_on="sv_id", how="left")["color"]
            bot_df["color"] = bot_df.merge(right=top_df, left_on="POS", right_on="POS", how="left")["color"]
        if color == "trait":
            h_values = np.linspace(0, 360, len(top_df))
            np.random.shuffle(h_values)
            colors = ['hsl(' + str(h_values[i]) + ',50%' + ',50%)' for i in range(len(h_values))]

            top_df["color"] = colors
            bot_df["color"] = bot_df.merge(right=top_df, left_on="POS", right_on="POS", how="left")["color"]
            mid_df["color"] = mid_df.merge(right=top_df, left_on="sv_id", right_on="sv_id", how="left").drop_duplicates(subset=["sv_id"])["color"]
        if color == "snp_type":
            h_values = np.linspace(0, 360, len(top_df[["REF"]].drop_duplicates()))
            snp_types = list(top_df[["REF"]].drop_duplicates())
            color_dict = dict(zip(snp_types, h_values))
            np.random.shuffle(h_values)
            colors = ['hsl(' + str(color_dict[ref]) + ',50%' + ',50%)' for ref in top_df["REF"]]

            top_df["color"] = colors
            bot_df["color"] = bot_df.merge(right=top_df, left_on="POS", right_on="POS", how="left")["color"]
            mid_df["color"] = mid_df.merge(right=top_df, left_on="sv_id", right_on="sv_id", how="left").drop_duplicates(subset=["sv_id"])["color"]
        if color == "sv_type":
            mid_df["color"] = ["blue" if sv_type == "JOIN" else "red" for sv_type in mid_df[["variant"]]]
            top_df["color"] = top_df.merge(right=mid_df, left_on="sv_id", right_on="sv_id", how="left")["color"]
            bot_df["color"] = bot_df.merge(right=top_df, left_on="POS", right_on="POS", how="left")["color"]
        if color == "corr":
            top_df["color"] = ["blue" if corr >= 0 else "red" for corr in top_df[["T_STAT"]].astype(int)]
            bot_df["color"] = bot_df.merge(right=top_df, left_on="POS", right_on="POS", how="left")["color"]
            mid_df["color"] = mid_df.merge(right=top_df, left_on="sv_id", right_on="sv_id", how="left").drop_duplicates(subset=["sv_id"])["color"]

    def plot(self, top_df: pd.DataFrame, mid_df: pd.DataFrame, bot_df: pd.DataFrame, color: str, focus: str, line: None | str = None):
        top_df = top_df.reset_index(drop=True)
        mid_df = mid_df.reset_index(drop=True)
        bot_df = bot_df.reset_index(drop=True)
        
        self.create_color_map(top_df, mid_df, bot_df, color)
        
        mid_height = 0.04

        fig = make_subplots(
            rows=3,
            cols=1,
            shared_xaxes=True,
            row_heights=[(1 - mid_height) / 2, mid_height, (1 - mid_height) / 2],
            vertical_spacing=0.01,
        )

        top_scatter = go.Scatter(
            x=top_df["POS"],
            y=top_df["-log10P"],
            hovertext=top_df["POS"],
            hoverinfo="text",
            mode="markers",
            marker=dict(
                color=top_df["color"]
            ),
            showlegend=True,
        )

        mid_bar = go.Bar(
            y=[0 for _ in mid_df],
            base=mid_df["start"],
            x=mid_df["length"],
            hovertext=mid_df["sv_id"],
            hoverinfo="text",
            orientation="h",
            marker=dict(
                color=mid_df["color"]
            )
        )

        bot_scatter = go.Scatter(
            x=bot_df["POS"],
            y=bot_df["-log10P"],
            hovertext=bot_df["POS"],
            hoverinfo="text",
            mode="markers",
            marker=dict(
                color=bot_df["color"]
            ),
        )

        fig.add_trace(top_scatter, row=1, col=1)
        fig.add_trace(mid_bar, row=2, col=1)
        fig.add_trace(bot_scatter, row=3, col=1)

        fig.update_xaxes(showgrid=False)
        fig.update_yaxes(showgrid=False)

        if top_df.empty and bot_df.empty:
            top_max = bot_max = 60
        elif top_df.empty:
            top_max = bot_max = max(bot_df["-log10P"].to_list())
        elif bot_df.empty:
            top_max = bot_max = max(top_df["-log10P"].to_list())
        else:
            top_max = max(top_df["-log10P"].to_list())
            bot_max = max(bot_df["-log10P"].to_list())

        fig.update_layout(
            showlegend=False,
            plot_bgcolor = "rgba(0,0,0,0)",
            yaxis2=dict(showticklabels=False),
            yaxis3=dict(autorange="reversed"),
            hoversubplots="axis",
            hovermode="x",
            shapes=[
            dict(
                type="rect",
                xref="paper",
                yref="y",
                x0=0,
                y0=0,
                x1=1,
                y1=top_max,
                fillcolor="lightgray",
                opacity=0.5,
                layer="below",
                line_width=0,
            ),
            dict(
                type="rect",
                xref="paper",
                yref="y3",
                x0=0,
                y0=0,
                x1=1,
                y1=bot_max,
                fillcolor="lightgray",
                opacity=0.5,
                layer="below",
                line_width=0,
            )
        ]
        )

        if line:
            def scaley(y: float, loc: str) -> float:
                if loc == "top":
                    return 0.5 + (mid_height / 2) + ((y / top_max) * (0.5 - (mid_height / 2)))
                if loc == "bot":
                    return 0.5 - (mid_height / 2) - ((y / bot_max) * (0.5 + (mid_height / 2)))
                
            if line.startswith("chr"):
                top_x0 = top_df.loc[top_df["sv_id"] == str(line), "POS"]
                top_x1 = mid_df.loc[mid_df["sv_id"] ==  str(line), "center"]
                top_y0 = scaley(top_df.loc[top_df["sv_id"] == str(line), "-log10P"], "top")
                top_y1 = 0.5

                bot_x0 = bot_df.loc[bot_df["sv_id"] == str(line), "POS"]
                bot_x1 = mid_df.loc[mid_df["sv_id"] ==  str(line), "center"]
                bot_y0 = scaley(bot_df.loc[bot_df["POS"] == str(bot_x0), "-log10P"], "bot")
                bot_y1 = 0.5

                color = top_df.loc[top_df["POS"] == str(top_x0), "color"]
            else:
                sv_id = top_df.loc[top_df["POS"] == str(line), "sv_id"]

                top_x0 = int(line)
                top_x1 = mid_df.loc[mid_df["sv_id"] == sv_id, "center"]
                top_y0 = scaley(top_df.loc[top_df["POS"] == str(line), "-log10P"], "top")
                top_y1 = 0.5

                bot_x0 = int(line)
                bot_x1 = mid_df.loc[mid_df["sv_id"] == sv_id, "center"]
                bot_y0 = scaley(bot_df.loc[bot_df["POS"] == str(line), "-log10P"], "bot")
                bot_y1 = 0.5

                color = top_df.loc[top_df["POS"] == str(line), "color"]
                
            fig.add_shape(
                type="line",
                xref="paper", yref="paper",
                x0=int(top_x0), y0=int(top_y0),
                x1=int(top_x1), y1=int(top_y1),
                line=dict(
                    color=color,
                    width=3,
                ),
            )

            fig.add_shape(
                type="line",
                xref="paper", yref="paper",
                x0=int(bot_x0), y0=int(bot_y0),
                x1=int(bot_x1), y1=int(bot_y1),
                line=dict(
                    color=color,
                    width=3,
                ),
            )

        if False:
            top_middle = (0.5 * (0.5 - (mid_height / 2))) / top_max
            
            if focus == "snp":
                for i, row in top_df.iterrows():
                    sv = StructuralVariant.parse(row["sv_id"])
                    fig.add_annotation(
                        ax=row["POS"], ay=row["-log10P"],
                        axref="x1", ayref="y1",
                        x=(int(sv.start) + int(sv.end)) / 2, y=-row["-log10P"],
                        xref="x1", yref="y1",
                        xclick=row["POS"], yclick=row["-log10P"],
                        showarrow=True,
                        # visible=False,
                        arrowhead=0,
                        arrowsize=1,
                        arrowwidth=2,
                        arrowcolor=row["color"],
                        clicktoshow="onoff",
                    )

                    traits = bot_df[bot_df["POS"] == row["POS"]]

                    if not traits.empty:
                        trait = traits.iloc[0]
                        fig.add_annotation(
                            ax=trait["POS"], ay=trait["-log10P"],
                            axref="x3", ayref="y3",
                            x=(int(sv.start) + int(sv.end)) / 2, y=-trait["-log10P"],
                            xref="x3", yref="y3",
                            xclick=trait["POS"], yclick=trait["-log10P"],
                            showarrow=True,
                            # visible=False,
                            arrowhead=0,
                            arrowsize=1,
                            arrowwidth=2,
                            arrowcolor=row["color"],
                            clicktoshow="onoff",
                        )

        plot_json = pio.to_json(fig)
        return plot_json