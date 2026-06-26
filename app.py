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
        dcc.Store(id="last-selected-protein", data=None),
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
    State("last-selected-protein", "data"),
)
def update_protein_options(value, last_selected):
    if value is None or value == "":
        return [], {"display": "none"}

    if value == last_selected:
        return no_update, {"display": "none"}

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
    Output("last-selected-protein", "data", allow_duplicate=True),
    Output("protein-submit", "n_clicks", allow_duplicate=True),
    Input({"type": "result-item", "index": ALL}, "n_clicks"),
    State("protein-submit", "n_clicks"),
    prevent_initial_call=True,
)
def select_item(n_clicks, protein_submit_n_clicks):
    if not ctx.triggered:
        return no_update, no_update, no_update, no_update

    if (np.array(n_clicks) == None).all():
        return no_update, no_update, no_update, no_update

    triggered_id = ctx.triggered_id
    selected_value = triggered_id["index"]

    return selected_value, {"display": "none"}, selected_value, protein_submit_n_clicks + 1


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
                                maxHeight=400,
                            ),
                            html.Br(),
                            html.P(
                                "Select Study Datasets:",
                                className="mb-2 text-white",
                                style={"fontWeight": "600"},
                            ),
                            dbc.Checklist(
                                options=[
                                    {"label": "GTEx (healthy)", "value": "GTEX"},
                                    {"label": "TCGA (cancer)", "value": "TCGA"},
                                ],
                                value=["GTEX", "TCGA"],
                                id="expression-study-filter",
                                inline=True,
                                switch=True,
                                className="mb-4 text-white",
                            ),
                            html.P("Sort order:", className="mb-2 text-white"),
                            dcc.Dropdown(
                                options=[
                                    "Alphabetical sort",
                                    "Median expression (descending)",
                                ],
                                value="Alphabetical sort",
                                id="expression-sorting-dropdown",
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
                                        "Reset",
                                        id="expression-reset-to-main-button",
                                        className="custom-btn btn-secondary mb-4 me-2",
                                        style={
                                            "background": "#dc2626",
                                            "border": "1px solid rgba(255,255,255,0.08)",
                                            "color": "#f8fafc",
                                        },
                                    ),
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
                                    children=html.Div(id="expression-plot"),
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
                        id="expression-plot-container",
                        style={"display": "none"},
                    ),
                    dbc.Modal(
                        [
                            dbc.ModalHeader(
                                dbc.ModalTitle("Filter & Sorting Parameters"),
                                close_button=True,
                            ),
                            dbc.ModalBody(
                                [
                                    # 1. Select tissue/cancer types
                                    html.P(
                                        "Select tissue/cancer types:",
                                        className="mb-2 text-white",
                                        style={"fontWeight": "600"},
                                    ),
                                    dcc.Dropdown(
                                        options=expression_df.loc[
                                            :, "tissue_type"
                                        ].unique(),
                                        value=tissue_types_inital_value,
                                        multi=True,
                                        id="expression-cancer-type-dropdown-modal",
                                        placeholder="Select tissue/cancer types...",
                                        className="dash-dropdown mb-4",
                                        style={"color": "#1e293b"},
                                        maxHeight=400,
                                    ),
                                    # 2. Select healthy vs cancer
                                    html.P(
                                        "Select Study Datasets:",
                                        className="mb-2 text-white",
                                        style={"fontWeight": "600"},
                                    ),
                                    dbc.Checklist(
                                        options=[
                                            {
                                                "label": "GTEx (healthy)",
                                                "value": "GTEX",
                                            },
                                            {"label": "TCGA (cancer)", "value": "TCGA"},
                                        ],
                                        value=["GTEX", "TCGA"],
                                        id="expression-study-filter-modal",
                                        inline=True,
                                        switch=True,
                                        className="mb-4 text-white",
                                    ),
                                    # 3. Sort order dropdown
                                    html.P(
                                        "Sort order:",
                                        className="mb-2 text-white",
                                        style={"fontWeight": "600"},
                                    ),
                                    dcc.Dropdown(
                                        options=[
                                            "Alphabetical sort",
                                            "Median expression (descending)",
                                        ],
                                        value="Alphabetical sort",
                                        id="expression-sorting-dropdown-modal",
                                        className="dash-dropdown mb-3",
                                        style={"color": "#1e293b"},
                                    ),
                                ]
                            ),
                            dbc.ModalFooter(
                                dbc.Button(
                                    "Apply",
                                    id="expression-modal-apply",
                                    className="custom-btn ms-auto",
                                )
                            ),
                        ],
                        id="expression-parameters-modal",
                        is_open=False,
                        className="custom-modal",
                    ),
                    dcc.Store(id="expression-parameters"),
                ],
                className="custom-card",
            )
        ]
    )


@app.callback(
    Output("expression-plot-container", "style"),
    Output("expression-parameters-container", "style"),
    Input("expression-load-button", "n_clicks"),
    Input("expression-reset-to-main-button", "n_clicks"),
    prevent_initial_call=True,
)
def toggle_expression_containers(load_clicks, reset_clicks):
    ctx = dash.callback_context
    if not ctx.triggered:
        return no_update, no_update

    trigger_id = ctx.triggered[0]["prop_id"].split(".")[0]
    if trigger_id == "expression-load-button":
        return {"display": "block"}, {"display": "none"}
    elif trigger_id == "expression-reset-to-main-button":
        return {"display": "none"}, {"display": "block"}
    return no_update, no_update


@app.callback(
    Output("expression-parameters", "data"),
    Input("expression-load-button", "n_clicks"),
    Input("expression-modal-apply", "n_clicks"),
    Input("expression-reset-to-main-button", "n_clicks"),
    State("expression-cancer-type-dropdown", "value"),
    State("expression-sorting-dropdown", "value"),
    State("expression-study-filter", "value"),
    State("expression-cancer-type-dropdown-modal", "value"),
    State("expression-sorting-dropdown-modal", "value"),
    State("expression-study-filter-modal", "value"),
    prevent_initial_call=True,
)
def update_parameters(
    load_clicks,
    apply_clicks,
    reset_clicks,
    initial_tissues,
    initial_sorting,
    initial_studies,
    modal_tissues,
    modal_sorting,
    modal_studies,
):
    ctx = dash.callback_context
    if not ctx.triggered:
        return no_update

    trigger_id = ctx.triggered[0]["prop_id"].split(".")[0]

    if trigger_id == "expression-load-button":
        return {
            "tissue_types": initial_tissues if initial_tissues is not None else [],
            "sorting": initial_sorting
            if initial_sorting is not None
            else "Alphabetical sort",
            "studies": initial_studies
            if initial_studies is not None
            else ["GTEX", "TCGA"],
        }
    elif trigger_id == "expression-modal-apply":
        return {
            "tissue_types": modal_tissues if modal_tissues is not None else [],
            "sorting": modal_sorting
            if modal_sorting is not None
            else "Alphabetical sort",
            "studies": modal_studies if modal_studies is not None else ["GTEX", "TCGA"],
        }
    elif trigger_id == "expression-reset-to-main-button":
        return None
    return no_update


@app.callback(
    Output("expression-parameters-modal", "is_open"),
    Output("expression-cancer-type-dropdown-modal", "value", allow_duplicate=True),
    Output("expression-sorting-dropdown-modal", "value"),
    Output("expression-study-filter-modal", "value"),
    Output("expression-cancer-type-dropdown-modal", "options", allow_duplicate=True),
    Input("expression-reset-button", "n_clicks"),
    Input("expression-modal-apply", "n_clicks"),
    State("expression-parameters-modal", "is_open"),
    State("expression-parameters", "data"),
    State("url", "search"),
    prevent_initial_call=True,
)
def toggle_modal(reset_clicks, apply_clicks, is_open, current_params, search):
    ctx = dash.callback_context
    if not ctx.triggered:
        return no_update, no_update, no_update, no_update, no_update

    trigger_id = ctx.triggered[0]["prop_id"].split(".")[0]

    if trigger_id == "expression-reset-button":
        tissues = current_params.get("tissue_types", []) if current_params else []
        sorting = (
            current_params.get("sorting", "Alphabetical sort")
            if current_params
            else "Alphabetical sort"
        )
        studies = (
            current_params.get("studies", ["GTEX", "TCGA"])
            if current_params
            else ["GTEX", "TCGA"]
        )

        protein = get_query_data(search)
        protein = getting_gene_names(protein.upper())
        expression_df = get_expression_data(protein)
        if expression_df is not None and studies:
            upper_studies = [s.upper() for s in studies]
            filtered_df = expression_df.loc[
                expression_df.loc[:, "study"].str.upper().isin(upper_studies), :
            ]
            options = list(filtered_df.loc[:, "tissue_type"].unique())
        else:
            options = []

        return True, tissues, sorting, studies, options
    elif trigger_id == "expression-modal-apply":
        return False, no_update, no_update, no_update, no_update

    return is_open, no_update, no_update, no_update, no_update


def update_tissue_dropdown_options(studies, search, current_values):
    protein = get_query_data(search)
    protein = getting_gene_names(protein.upper())
    expression_df = get_expression_data(protein)
    if expression_df is None or not studies:
        return [], []

    # Filter by study (GTEx/TCGA)
    upper_studies = [s.upper() for s in studies]
    filtered_df = expression_df.loc[
        expression_df.loc[:, "study"].str.upper().isin(upper_studies), :
    ]

    options = list(filtered_df.loc[:, "tissue_type"].unique())

    # Filter selected value to only retain valid options
    valid_values = [v for v in current_values if v in options] if current_values else []

    return options, valid_values


@app.callback(
    Output("expression-cancer-type-dropdown", "options"),
    Output("expression-cancer-type-dropdown", "value", allow_duplicate=True),
    Input("expression-study-filter", "value"),
    State("url", "search"),
    State("expression-cancer-type-dropdown", "value"),
    prevent_initial_call=True,
)
def update_main_tissue_dropdown(studies, search, current_values):
    return update_tissue_dropdown_options(studies, search, current_values)


@app.callback(
    Output("expression-cancer-type-dropdown-modal", "options"),
    Output("expression-cancer-type-dropdown-modal", "value", allow_duplicate=True),
    Input("expression-study-filter-modal", "value"),
    State("url", "search"),
    State("expression-cancer-type-dropdown-modal", "value"),
    prevent_initial_call=True,
)
def update_modal_tissue_dropdown(studies, search, current_values):
    return update_tissue_dropdown_options(studies, search, current_values)


@app.callback(
    Output("expression-cancer-type-dropdown", "value"),
    Output("expression-sorting-dropdown", "value"),
    Output("expression-study-filter", "value"),
    Input("expression-reset-to-main-button", "n_clicks"),
    prevent_initial_call=True,
)
def reset_main_page_inputs(n_clicks):
    if n_clicks:
        return tissue_types_inital_value, "Alphabetical sort", ["GTEX", "TCGA"]
    return no_update, no_update, no_update


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
    Output("expression-plot", "children"),
    Input("expression-parameters", "data"),
    State("url", "search"),
    prevent_initial_call=True,
    optional=True,
)
def update_expression_plot(parameters, search):
    if parameters is None:
        return html.Div()

    protein = get_query_data(search)
    protein = getting_gene_names(protein.upper())
    expression_df = get_expression_data(protein)
    if expression_df is None:
        return html.Div()

    tissue_types = parameters.get("tissue_types", [])
    sorting = parameters.get("sorting", "Alphabetical sort")
    studies = parameters.get("studies", ["GTEX", "TCGA"])

    # Filter cancer type
    if len(tissue_types) > 0:
        expression_df = expression_df.loc[
            expression_df.loc[:, "tissue_type"].isin(tissue_types), :
        ]

    # Filter study (GTEx/TCGA)
    if len(studies) > 0:
        upper_studies = [s.upper() for s in studies]
        expression_df = expression_df.loc[
            expression_df.loc[:, "study"].str.upper().isin(upper_studies), :
        ]
    else:
        expression_df = expression_df.iloc[0:0]

    # Write a message if GTEX and TCGA are unselected or if no matching data is found
    if len(studies) == 0:
        fig = go.Figure()
        fig.add_annotation(
            text="Please select at least one study dataset (GTEx or TCGA) to view expression data.",
            xref="paper",
            yref="paper",
            x=0.5,
            y=0.5,
            showarrow=False,
            font=dict(size=14, color="#cbd5e1", family="Inter"),
        )
        fig.update_layout(
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)",
            xaxis=dict(showgrid=False, showticklabels=False, zeroline=False),
            yaxis=dict(showgrid=False, showticklabels=False, zeroline=False),
        )
        return dcc.Graph(figure=fig)

    if expression_df.empty:
        fig = go.Figure()
        fig.add_annotation(
            text="No expression data matches the current selection.",
            xref="paper",
            yref="paper",
            x=0.5,
            y=0.5,
            showarrow=False,
            font=dict(size=14, color="#cbd5e1", family="Inter"),
        )
        fig.update_layout(
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)",
            xaxis=dict(showgrid=False, showticklabels=False, zeroline=False),
            yaxis=dict(showgrid=False, showticklabels=False, zeroline=False),
        )
        return dcc.Graph(figure=fig)

    if sorting == "Alphabetical sort":
        expression_df = expression_df.sort_values(by="tissue_type", ascending=True)
    elif sorting == "Median expression (descending)":
        expression_df["average_median"] = [
            np.mean(medians) for medians in expression_df.loc[:, "median"]
        ]
        expression_df = expression_df.sort_values(by="average_median", ascending=False)
        expression_df = expression_df.drop(columns=["average_median"])

    # Generate the plots
    tissue_figures = plot_expression_data(expression_df)

    graphs = []
    for tissue_type, fig in tissue_figures:
        # Determine column width based on the number of transcripts in the plot
        num_transcripts = 0
        if fig.data and hasattr(fig.data[0], "x") and fig.data[0].x is not None:
            num_transcripts = len(set(fig.data[0].x))

        if num_transcripts <= 2:
            md_val = 6  # 2 plots per line
            lg_val = 3  # 4 plots per line
            xl_val = 2  # 6 plots per line
        elif num_transcripts <= 5:
            md_val = 6  # 2 plots per line
            lg_val = 4  # 3 plots per line
            xl_val = 3  # 4 plots per line
        elif num_transcripts <= 8:
            md_val = 12  # 1 plot per line
            lg_val = 6  # 2 plots per line
            xl_val = 4  # 3 plots per line
        elif num_transcripts <= 12:
            md_val = 12  # 1 plot per line
            lg_val = 12  # 1 plot per line
            xl_val = 6  # 2 plots per line
        else:
            md_val = 12  # 1 plot per line
            lg_val = 12  # 1 plot per line
            xl_val = 12  # 1 plot per line

        graphs.append(
            dbc.Col(
                dcc.Graph(
                    figure=fig,
                    config={"displayModeBar": False},
                ),
                xs=12,  # 1 plot per line on mobile/narrow screens
                md=md_val,  # 1 plot per line on medium screens (not large)
                lg=lg_val,
                xl=xl_val,
            )
        )

    return dbc.Row(graphs, className="g-4")


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
