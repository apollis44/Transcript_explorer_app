import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import math


def create_localization_plot(localization_data, all_transcripts):
    print("Creating localization plot...")
    if len(all_transcripts) != 0:
        df_exploded = localization_data.copy().reset_index()
        df_exploded["Protein_ID"] = df_exploded["Protein_ID"].str.split("<br>")
        df_exploded = df_exploded.explode("Protein_ID")
        df_exploded.set_index("Protein_ID", inplace=True)
        localization_data_plot = df_exploded
    else:
        localization_data_plot = localization_data.copy()
        localization_data_plot.index = [
            index.split("<br>")[0] for index in localization_data_plot.index
        ]

    # Sort by alphabetical order
    localization_data_plot.sort_index(inplace=True)

    fig = go.Figure(
        data=go.Heatmap(
            z=localization_data_plot.values,
            x=localization_data_plot.columns.tolist(),
            y=localization_data_plot.index.tolist(),
            colorscale=[[0.0, "#2563eb"], [0.5, "#f1f5f9"], [1.0, "#d32f2f"]],
            zmin=0,
            zmax=1,
            xgap=1.5,
            ygap=1.5,
            colorbar=dict(
                title=dict(
                    text="Probability",
                    font=dict(color="#eceef4", size=14, family="Inter"),
                    side="top",
                ),
                tickfont=dict(color="#8f9bb3", size=12, family="Inter"),
                thickness=30,
                lenmode="pixels",
                len=300,
                yanchor="top",
                y=1,
                xanchor="left",
                x=1.02,
            ),
        )
    )

    fig.update_layout(
        height=len(localization_data_plot) * 50 + 100,
        yaxis=dict(
            autorange="reversed",
            tickfont=dict(
                color="#8f9bb3", size=11, family="Inter"
            ),  # Soft blue-gray text color
            showgrid=False,
            ticklabelstandoff=15,
        ),
        xaxis=dict(
            tickfont=dict(
                color="#8f9bb3", size=11, family="Inter"
            ),  # Soft blue-gray text color
            showgrid=False,
            side="top",
        ),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        margin=dict(l=150, r=150, t=50, b=50),
    )

    return fig


def create_topology_plot(
    mapping, sequences_data, unique_transcripts, x_label, all_transcripts
):
    print("Creating topology plot...")

    # High contrast vibrant color map
    color_map = {
        "S": "#ff7675",  # Salmon / Coral Red
        "O": "#fbbf24",  # Amber (Extracellular)
        "I": "#60a5fa",  # Blue (Intracellular)
        "-": "#64748b",  # Slate (Alignment gap)
        "M": "#10b981",  # Emerald (Transmembrane)
    }

    letter_to_label = {
        "S": "Signal peptide",
        "O": "Extracellular",
        "I": "Intracellular",
        "-": "Alignment gap",
        "M": "Transmembrane",
    }

    fig = go.Figure()
    isoforms = list(dict.fromkeys(mapping.values()))
    y_labels = [
        transcript_id for isoform in isoforms for transcript_id in isoform.split("<br>")
    ]
    y_labels, sequences_data = zip(
        *sorted(zip(y_labels, sequences_data), key=lambda x: x[0])
    )

    y_label_available = [y_label in unique_transcripts for y_label in y_labels]

    # Group data by feature to minimize Plotly traces
    feature_data = {
        feat: {"x": [], "base": [], "y": [], "hovertext": []}
        for feat in color_map.keys()
    }

    # Loop through sequences (isoforms)
    for i, seq_data in enumerate(sequences_data):
        if len(all_transcripts) == 0 and not y_label_available[i]:
            continue

        y_val = y_labels[i]

        for feature, ranges in seq_data.items():
            if feature not in feature_data:
                feature_data[feature] = {"x": [], "base": [], "y": [], "hovertext": []}
            for start, width in ranges:
                feature_data[feature]["x"].append(width)
                feature_data[feature]["base"].append(start)
                feature_data[feature]["y"].append(y_val)
                feat_label = letter_to_label.get(feature, feature)
                feature_data[feature]["hovertext"].append(
                    f"<b>{feat_label}</b><br>"
                    f"Range: {start} - {start + width}<br>"
                    f"Length: {width}"
                )

    for feature, data in feature_data.items():
        if not data["x"]:
            continue
        fig.add_trace(
            go.Bar(
                name=letter_to_label.get(feature, feature),
                x=data["x"],
                base=data["base"],
                y=data["y"],
                orientation="h",
                marker_color=color_map.get(feature, "#000000"),
                hovertext=data["hovertext"],
                hovertemplate="%{hovertext}<extra></extra>",
                legendgroup=feature,
                showlegend=True,
            )
        )

    # Calculate height to match original logic
    if len(all_transcripts) == 0:
        num_isoforms = len(unique_transcripts)
    else:
        num_isoforms = len(y_labels)
    calculated_height = num_isoforms * 50 + 100

    fig.update_layout(
        xaxis_title=x_label,
        barmode="stack",
        bargap=0.15,
        height=calculated_height,
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="rgba(0,0,0,0)",
        legend=dict(
            title=dict(
                text="Features", font=dict(color="#eceef4", size=11, family="Inter")
            ),
            font=dict(color="#8f9bb3", size=10, family="Inter"),
            orientation="v",
            yanchor="top",
            y=1,
            xanchor="left",
            x=1.02,
            bgcolor="rgba(50, 55, 78, 0.8)",
            bordercolor="rgba(255,255,255,0.15)",
            borderwidth=1,
        ),
        xaxis=dict(
            title=dict(font=dict(color="#8f9bb3", size=12, family="Inter")),
            tickfont=dict(color="#8f9bb3", size=10, family="Inter"),
            gridcolor="rgba(255,255,255,0.15)",
            showgrid=True,
            zeroline=False,
        ),
        yaxis=dict(
            autorange="reversed",
            tickfont=dict(color="#8f9bb3", size=11, family="Inter"),
            showgrid=False,
            zeroline=False,
            ticklabelstandoff=15,
        ),
        margin=dict(l=150, r=150, t=50, b=50),
    )

    return fig


def plot_expression_data(expression_df):
    print("Creating expression plot...")

    tissue_types = expression_df.loc[:, "tissue_type"].unique().tolist()

    if len(tissue_types) == 0:
        return []

    # Calculate global y-axis range across all tissues to share the same y-axis scale
    global_min = None
    global_max = None
    for _, row in expression_df.iterrows():
        for lf in row["lowerfence"]:
            if lf is not None and not pd.isna(lf):
                if global_min is None or lf < global_min:
                    global_min = lf
        for uf in row["upperfence"]:
            if uf is not None and not pd.isna(uf):
                if global_max is None or uf > global_max:
                    global_max = uf
        for outliers in row["y"]:
            for val in outliers:
                if val is not None and not pd.isna(val):
                    if global_min is None or val < global_min:
                        global_min = val
                    if global_max is None or val > global_max:
                        global_max = val

    if global_min is not None and global_max is not None:
        y_padding = (global_max - global_min) * 0.05 if global_max != global_min else 1.0
        y_range = [global_min - y_padding, global_max + y_padding]
    else:
        y_range = None

    # Coordinated theme-aligned high contrast colors from the reference image
    vibrant_colors = [
        "#ff7675",  # Salmon / Coral Red
        "#2ec4b6",  # Mint / Turquoise
        "#ffa62b",  # Amber / Orange
        "#9b5de5",  # Violet / Purple
        "#00b4d8",  # Cyan / Light Blue
        "#82c91e",  # Lime / Yellow-Green
    ]

    figures = []

    for i, tissue_type in enumerate(tissue_types):
        data_for_each_tissue_type = expression_df.loc[
            (expression_df.loc[:, "tissue_type"] == tissue_type), :
        ]
        current_color = vibrant_colors[i % len(vibrant_colors)]

        fig = go.Figure()

        # Add Box plot
        fig.add_trace(
            go.Box(
                x=data_for_each_tissue_type["protein"].iloc[0],
                q1=data_for_each_tissue_type["q1"].iloc[0],
                q3=data_for_each_tissue_type["q3"].iloc[0],
                median=data_for_each_tissue_type["median"].iloc[0],
                lowerfence=data_for_each_tissue_type["lowerfence"].iloc[0],
                upperfence=data_for_each_tissue_type["upperfence"].iloc[0],
                showlegend=False,
                marker_color=current_color,
                line=dict(width=1.5, color=current_color),
                fillcolor="rgba(0, 0, 0, 0)",
            )
        )

        proteins = data_for_each_tissue_type["protein"].iloc[0]
        outliers_lists = data_for_each_tissue_type["y"].iloc[0]

        # Combine coordinates for Scatter
        all_x = []
        all_y = []

        for protein, values in zip(proteins, outliers_lists):
            all_x.extend([protein] * len(values))
            all_y.extend(values)

        fig.add_trace(
            go.Scatter(
                x=all_x,
                y=all_y,
                mode="markers",
                marker=dict(
                    size=4.5,
                    symbol="circle-open",
                    color=current_color,
                    opacity=0.8,
                ),
                showlegend=False,
            )
        )

        # Style individual plot layout
        fig.update_layout(
            title=dict(
                text=tissue_type,
                font=dict(size=13, color="#eceef4", family="Inter"),
                x=0.5,
                xanchor="center",
            ),
            height=380,
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)",
            margin=dict(l=50, r=20, t=50, b=50),
            xaxis=dict(
                tickfont=dict(color="#8f9bb3", size=10, family="Inter"),
                showgrid=False
            ),
            yaxis=dict(
                title=dict(
                    text="log2(TPM+1)", font=dict(color="#8f9bb3", size=11, family="Inter")
                ),
                tickfont=dict(color="#8f9bb3", size=10, family="Inter"),
                gridcolor="rgba(255,255,255,0.15)",
                showgrid=True,
                zeroline=False,
                range=y_range,
            )
        )

        figures.append((tissue_type, fig))

    return figures
