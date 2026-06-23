import dash
import dash_bootstrap_components as dbc
from dash import Input, Output, dcc, html, State, ctx, ALL, no_update
from scripts.plots_generation import (
    create_localization_plot,
    plot_expression_data,
    create_topology_plot,
)
import pandas as pd
import os
from functools import lru_cache
import shelve
import numpy as np
import urllib.parse
import plotly.graph_objects as go
import bisect

# Initial values for the dropdowns
tissue_types_inital_value = []

app = dash.Dash(
    external_stylesheets=[dbc.themes.BOOTSTRAP], suppress_callback_exceptions=True
)
server = app.server
base_dir = os.path.dirname(os.path.abspath(__file__))

# Open shelve databases globally at startup in read-only mode for performance
db_names = shelve.open(f"{base_dir}/files_for_plots/genes", flag="r")
db_deeploc = shelve.open(f"{base_dir}/files_for_plots/deeploc2_output", flag="r")
db_tcga = shelve.open(f"{base_dir}/files_for_plots/TCGA_GTEx_plotting_data", flag="r")
db_mapping = shelve.open(
    f"{base_dir}/files_for_plots/transcripts_to_isoforms_mapping", flag="r"
)
db_topology = shelve.open(
    f"{base_dir}/files_for_plots/membrane_topology_objects", flag="r"
)

# Sort protein options once at startup for fast binary search
protein_options = sorted([p.upper() for p in db_names.keys()])


def getting_gene_names(gene_name):
    if gene_name not in db_names:
        return None
    gene_names = db_names[gene_name]
    return gene_names


header = html.Div(
    [
        html.H3("Transcript Explorer", className="header-title me-auto"),
        html.Div(
            [
                dbc.Input(
                    id="protein-input",
                    placeholder="Search protein (e.g. TP53)...",
                    autoComplete="off",
                    className="custom-input me-2",
                    style={"width": "280px"},
                ),
                dbc.ListGroup(
                    id="protein-options",
                    className="custom-autocomplete",
                    style={
                        "position": "absolute",
                        "top": "100%",
                        "left": 0,
                        "width": "280px",
                        "zIndex": 1000,
                        "maxHeight": "200px",
                        "overflowY": "auto",
                        "display": "none",
                    },
                ),
            ],
            style={"position": "relative"},
            className="d-flex align-items-center me-2",
        ),
        dbc.Button("Search", id="protein-submit", n_clicks=0, className="custom-btn"),
    ],
    className="custom-header d-flex align-items-center mb-4",
)

sidebar = html.Div(
    [
        dbc.Nav(
            [
                dbc.NavLink(
                    "Description", id="link-description", href="/", active="exact"
                ),
                dbc.NavLink(
                    "Localization",
                    id="link-localization",
                    href="/Localization",
                    active="exact",
                ),
                dbc.NavLink(
                    "Topology", id="link-topology", href="/Topology", active="exact"
                ),
                dbc.NavLink(
                    "Expression",
                    id="link-expression",
                    href="/Expression",
                    active="exact",
                ),
            ],
            vertical=True,
            pills=True,
        ),
    ],
    className="custom-sidebar",
)

content = html.Div(id="page-content")

app.layout = dbc.Container(
    [
        dcc.Location(id="url", refresh=False),
        dcc.Store(
            id="checklist-store",
            data={"localization": ["all"], "topology": ["all"], "protein": None},
        ),
        header,
        dbc.Row(
            [
                # Left Column (Sidebar) - Width 2/12
                dbc.Col(sidebar, width=2),
                # Right Column (Content) - Width 10/12
                dbc.Col(
                    [
                        content,
                    ],
                    width=10,
                ),
            ],
            className="main-row g-4",
        ),
    ],
    fluid=True,
    className="p-0",
)


@app.callback(
    Output("protein-options", "children"),
    Output("protein-options", "style"),
    Input("protein-input", "value"),
    State({"type": "result-item", "index": ALL}, "n_clicks"),
)
def update_protein_options(value, n_clicks):
    if value is None:
        return [], {"display": "none"}

    if n_clicks != [] and (np.array(n_clicks) != None).any():
        return no_update, no_update

    options = []

    prefix = value.upper()
    start_idx = bisect.bisect_left(protein_options, prefix)
    end_idx = bisect.bisect_left(protein_options, prefix + "\xff")
    valid_protein_options = protein_options[start_idx:end_idx]

    for protein in valid_protein_options[:10]:
        options.append(
            dbc.ListGroupItem(
                protein, id={"type": "result-item", "index": protein}, action=True
            )
        )

    if len(valid_protein_options) > 10:
        options.append(
            dbc.ListGroupItem(
                f"{len(valid_protein_options) - 10} more options",
                disabled=True,
                color="light",
            )
        )

    style = {
        "position": "absolute",
        "top": "100%",
        "left": 0,
        "width": "100%",
        "zIndex": 1000,
        "display": "block",
    }
    return options, style


@app.callback(
    Output("protein-input", "value"),
    Output("protein-options", "style", allow_duplicate=True),
    Output("protein-submit", "n_clicks", allow_duplicate=True),
    Input({"type": "result-item", "index": ALL}, "n_clicks"),
    State("protein-submit", "n_clicks"),
    prevent_initial_call=True,
)
def select_item(n_clicks, protein_submit_n_clicks):
    if not ctx.triggered:
        return no_update, no_update, no_update

    if (np.array(n_clicks) == None).all():
        return no_update, no_update, no_update

    triggered_id = ctx.triggered_id
    selected_value = triggered_id["index"]

    return selected_value, {"display": "none"}, protein_submit_n_clicks + 1


@lru_cache(maxsize=10)
def get_localization_data(protein):
    df = db_deeploc.get(protein)
    if df is None:
        return None
    df = df.iloc[
        :, :-4
    ]  # We exclude the last 4 columns that contains data we don't use
    return df


@lru_cache(maxsize=10)
def get_expression_data(protein):
    df = db_tcga.get(protein)
    if df is None:
        return None
    df = df.groupby(["study", "tissue_type"], sort=False).agg(list).reset_index()
    return df


@lru_cache(maxsize=10)
def get_topology_data(protein):
    mapping = db_mapping.get(protein)
    sequences_data = db_topology.get(protein)

    if mapping is None or sequences_data is None:
        return None, None, None

    unique_transcripts = []
    for unique_transcript in set(mapping.values()):
        unique_transcripts.append(unique_transcript.split("<br>")[0])

    return mapping, sequences_data, unique_transcripts


@lru_cache(maxsize=10)
def get_query_data(search):
    search_str = search.lstrip("?") if "?" in search else ""
    parsed_search = urllib.parse.parse_qs(search_str)
    if "protein" not in parsed_search:
        return None
    return parsed_search["protein"][0]


@app.callback(
    Output("url", "search"),
    Input("protein-input", "n_submit"),
    Input("protein-submit", "n_clicks"),
    State("url", "search"),
    State("protein-input", "value"),
    prevent_initial_call=True,
)
def update_query(n_submit, n_clicks, search, protein):
    if not ctx.triggered:
        return no_update

    search_str = search.lstrip("?") if "?" in search else ""
    parsed_search = urllib.parse.parse_qs(search_str)
    parsed_search["protein"] = [protein]
    search_str = urllib.parse.urlencode(parsed_search, doseq=True)
    return "?" + search_str


@app.callback(
    [
        Output("link-description", "href"),
        Output("link-localization", "href"),
        Output("link-topology", "href"),
        Output("link-expression", "href"),
    ],
    Input("url", "search"),
)
def update_nav_links(search):
    # If search is None or empty, just return the base paths
    query = search if search else ""

    # We return the base path + the current query string for each link
    return (
        f"/{query}",
        f"/Localization{query}",
        f"/Topology{query}",
        f"/Expression{query}",
    )


@app.callback(
    Output("page-content", "children"),
    Output("protein-options", "style", allow_duplicate=True),
    Input("url", "search"),
    Input("url", "pathname"),
    State("checklist-store", "data"),
    prevent_initial_call=True,
)
def render_page_content(query, pathname, checklist_data):

    protein = get_query_data(query)

    if protein is None:
        return html.Div(
            [
                html.Div(
                    [
                        html.H4(
                            "Welcome to Transcript Explorer", className="card-title"
                        ),
                        html.P(
                            "Please search and select a protein in the search bar above to begin analysis.",
                            className="desc-text",
                        ),
                    ],
                    className="custom-card text-center py-5",
                )
            ]
        ), {"display": "none"}

    protein = getting_gene_names(protein.upper())

    if protein is None:
        return html.Div(
            [
                html.Div(
                    [
                        html.H4(
                            "Protein Not Found", className="card-title text-danger"
                        ),
                        html.P(
                            "The requested gene/protein is not in the database. Please try another search.",
                            className="desc-text",
                        ),
                    ],
                    className="custom-card text-center py-5",
                )
            ]
        ), {"display": "none"}

    # Check if the gene encodes any valid protein
    has_valid_protein = db_mapping.get(protein) is not None

    if not has_valid_protein:
        return html.Div(
            [
                html.Div(
                    [
                        html.H4(f"{protein} Profile Overview", className="card-title"),
                        html.P(
                            f"This gene ({protein}) does not encode for any known protein coding variant in this database.",
                            className="desc-text text-warning",
                        ),
                    ],
                    className="custom-card",
                )
            ]
        ), {"display": "none"}

    loc_value = ["all"]
    top_value = ["all"]
    if checklist_data and checklist_data.get("protein") == protein:
        loc_value = checklist_data.get("localization", ["all"])
        top_value = checklist_data.get("topology", ["all"])

    if pathname == "/":
        return html.Div(
            [
                html.Div(
                    [
                        html.H4(f"{protein} Profile Overview", className="card-title"),
                        html.P(
                            f"Welcome to the interactive Transcript Explorer for the gene {protein}.",
                            className="desc-text mb-4",
                        ),
                        html.Div(
                            [
                                html.Span(
                                    "Status: Active Entry", className="badge-info me-2"
                                ),
                                html.Span(
                                    "Available: Localization, Topology, Expression",
                                    className="badge-info",
                                ),
                            ],
                            className="mb-4",
                        ),
                        html.Hr(style={"borderColor": "rgba(255,255,255,0.05)"}),
                        html.H5(
                            "Explore Transcript Variants",
                            className="mt-4 mb-3",
                            style={"fontWeight": "600", "fontSize": "1.1rem"},
                        ),
                        html.P(
                            "Select from the options in the sidebar navigation to view cellular localization predictions, "
                            "membrane topology alignments, and transcript expression comparison levels.",
                            className="desc-text",
                        ),
                    ],
                    className="custom-card",
                )
            ]
        ), {"display": "none"}

    elif pathname == "/Localization":
        return html.Div(
            [
                html.Div(
                    [
                        html.H4("Subcellular Localization", className="card-title"),
                        html.P(
                            "This heatmap displays the predicted probability scores of different subcellular localization destinations "
                            "for the transcript isoforms of the protein.",
                            className="desc-text mb-4",
                        ),
                        dbc.Checklist(
                            options=[
                                {
                                    "label": "Show all transcripts, including duplicates mapping to the same protein sequence",
                                    "value": "all",
                                },
                            ],
                            value=loc_value,
                            switch=True,
                            id="localization-checklist",
                            className="mb-4",
                        ),
                        html.Div(
                            dbc.Spinner(
                                children=dcc.Graph(id="localization-plot"),
                                size="lg",
                                color="primary",
                                type="border",
                                fullscreen=False,
                                id="localization-spinner",
                                spinner_style={
                                    "position": "absolute",
                                    "top": "50px",
                                    "left": "50%",
                                    "transform": "translate(-50%, -50%)",
                                },
                            ),
                            style={"position": "relative", "minHeight": "200px"},
                        ),
                    ],
                    className="custom-card",
                )
            ]
        ), {"display": "none"}

    elif pathname == "/Topology":
        return html.Div(
            [
                html.Div(
                    [
                        html.H4("Membrane Topology Alignment", className="card-title"),
                        html.P(
                            "Visualizes the predicted structural features (transmembrane regions, intracellular/extracellular segments) "
                            "along the amino acid alignments for each transcript isoform.",
                            className="desc-text mb-4",
                        ),
                        dbc.Checklist(
                            options=[
                                {
                                    "label": "Show all transcripts, including duplicates mapping to the same protein sequence",
                                    "value": "all",
                                },
                            ],
                            value=top_value,
                            switch=True,
                            id="topology-checklist",
                            className="mb-4",
                        ),
                        html.Div(
                            dbc.Spinner(
                                children=dcc.Graph(id="topology-plot"),
                                size="lg",
                                color="primary",
                                type="border",
                                fullscreen=False,
                                id="topology-spinner",
                                spinner_style={
                                    "position": "absolute",
                                    "top": "50px",
                                    "left": "50%",
                                    "transform": "translate(-50%, -50%)",
                                },
                            ),
                            style={"position": "relative", "minHeight": "200px"},
                        ),
                    ],
                    className="custom-card",
                )
            ]
        ), {"display": "none"}

    elif pathname == "/Expression":
        return html.Div(id="expression-container"), {"display": "none"}


@app.callback(
    Output("expression-container", "children"),
    Input("expression-container", "id"),
    State("url", "search"),
    optional=True,
)
def manage_expression_page(expression_container_id, search):
    protein = get_query_data(search)
    protein = getting_gene_names(protein.upper())
    expression_df = get_expression_data(protein)

    if expression_df is None:
        return html.Div(
            [
                html.Div(
                    [
                        html.H4(
                            "Transcript Expression Profile", className="card-title"
                        ),
                        html.P(
                            "Expression data is not available for this protein's transcripts.",
                            className="desc-text",
                        ),
                    ],
                    className="custom-card",
                )
            ]
        )

    return html.Div(
        [
            html.Div(
                [
                    html.H4("Transcript Expression Profile", className="card-title"),
                    html.P(
                        "Analyze and compare expression distributions across various healthy tissues and cancer types. "
                        "Select target cancer/tissue types below or leave empty to plot all datasets.",
                        className="desc-text mb-4",
                    ),
                    html.Div(
                        [
                            dcc.Dropdown(
                                options=expression_df.loc[:, "tissue_type"].unique(),
                                value=tissue_types_inital_value,
                                multi=True,
                                id="expression-cancer-type-dropdown",
                                placeholder="Select tissue/cancer types...",
                                className="dash-dropdown mb-3",
                                style={"color": "#1e293b"},
                            ),
                            # Load button
                            dbc.Button(
                                "Load Expression Plot",
                                id="expression-load-button",
                                className="custom-btn w-100",
                            ),
                        ],
                        id="expression-parameters-container",
                    ),
                    html.Div(
                        [
                            html.Div(
                                [
                                    dbc.Button(
                                        "Filter Parameters",
                                        id="expression-reset-button",
                                        className="custom-btn btn-secondary mb-4",
                                        style={
                                            "background": "#1e293b",
                                            "border": "1px solid rgba(255,255,255,0.08)",
                                            "color": "#f8fafc",
                                        },
                                    ),
                                ],
                                className="d-flex justify-content-end",
                            ),
                            html.Div(
                                dbc.Spinner(
                                    children=dcc.Graph(id="expression-plot"),
                                    size="lg",
                                    color="primary",
                                    type="border",
                                    id="expression-spinner",
                                    # Supprime ou vide spinner_style car Flexbox va gérer le centrage automatiquement
                                    spinner_style={},
                                ),
                                style={"position": "relative", "minHeight": "150px"},
                            ),
                        ],
                        id="expression-container",
                        style={"display": "none"},
                    ),
                    dcc.Store(id="expression-parameters"),
                ],
                className="custom-card",
            )
        ]
    )


@app.callback(
    Output("expression-container", "style"),
    Output("expression-parameters-container", "style"),
    Output("expression-parameters", "data"),
    Output("expression-cancer-type-dropdown", "value"),
    Input("expression-load-button", "n_clicks"),
    Input("expression-reset-button", "n_clicks"),
    State("expression-cancer-type-dropdown", "value"),
    prevent_initial_call=True,
    optional=True,
)
def expression_container_style(_1, _2, tissue_types):
    ctx = dash.callback_context
    trigger_id = ctx.triggered[0]["prop_id"].split(".")[0]
    if trigger_id == "expression-load-button":
        return {"display": "block"}, {"display": "none"}, tissue_types, tissue_types
    elif trigger_id == "expression-reset-button":
        return (
            {"display": "none"},
            {"display": "block"},
            tissue_types_inital_value,
            tissue_types_inital_value,
        )


@app.callback(
    Output("localization-plot", "figure"),
    Input("url", "search"),
    Input("localization-checklist", "value"),
    optional=True,
)
def localization_plot(search, all_transcripts):
    protein = get_query_data(search)
    if not protein:
        return go.Figure()
    protein = getting_gene_names(protein.upper())
    if not protein:
        return go.Figure()
    localization_data = get_localization_data(protein)
    if localization_data is None:
        return go.Figure()
    fig = create_localization_plot(localization_data, all_transcripts)
    return fig


@app.callback(
    Output("topology-plot", "figure"),
    Input("url", "search"),
    Input("topology-checklist", "value"),
    optional=True,
)
def topology_plot(search, all_transcripts):
    protein = get_query_data(search)
    if not protein:
        return go.Figure()
    protein = getting_gene_names(protein.upper())
    if not protein:
        return go.Figure()
    x_label = "Amino acid position in MSA"
    mapping, sequences_data, unique_transcripts = get_topology_data(protein)

    fig = create_topology_plot(
        mapping, sequences_data, unique_transcripts, x_label, all_transcripts
    )
    return fig


@app.callback(
    Output("expression-plot", "figure"),
    Input("expression-parameters", "data"),
    State("url", "search"),
    Input("expression-reset-button", "n_clicks"),
    prevent_initial_call=True,
    optional=True,
)
def update_expression_plot(parameters, search, is_clicked):
    ctx = dash.callback_context
    trigger_id = ctx.triggered[0]["prop_id"].split(".")[0]
    if trigger_id == "expression-reset-button":
        return

    protein = get_query_data(search)
    protein = getting_gene_names(protein.upper())
    expression_df = get_expression_data(protein)
    tissue_types = parameters

    # Filter cancer type
    if len(tissue_types) > 0:
        expression_df = expression_df.loc[
            expression_df.loc[:, "tissue_type"].isin(tissue_types), :
        ]

    # Generate the plot
    fig = plot_expression_data(expression_df)

    # Reset the buttons and close the popover
    return fig


@app.callback(
    Output("checklist-store", "data", allow_duplicate=True),
    Input("localization-checklist", "value"),
    State("checklist-store", "data"),
    State("url", "search"),
    prevent_initial_call=True,
)
def update_localization_store(value, current_data, search):
    protein = get_query_data(search)
    if protein:
        protein = getting_gene_names(protein.upper())

    if current_data is None or current_data.get("protein") != protein:
        return {
            "localization": value if value is not None else [],
            "topology": ["all"],
            "protein": protein,
        }

    current_data["localization"] = value if value is not None else []
    return current_data


@app.callback(
    Output("checklist-store", "data", allow_duplicate=True),
    Input("topology-checklist", "value"),
    State("checklist-store", "data"),
    State("url", "search"),
    prevent_initial_call=True,
)
def update_topology_store(value, current_data, search):
    protein = get_query_data(search)
    if protein:
        protein = getting_gene_names(protein.upper())

    if current_data is None or current_data.get("protein") != protein:
        return {
            "localization": ["all"],
            "topology": value if value is not None else [],
            "protein": protein,
        }

    current_data["topology"] = value if value is not None else []
    return current_data


if __name__ == "__main__":
    app.run(debug=True)
