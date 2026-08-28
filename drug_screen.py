#!/usr/bin/env python

import os
import io
import dash
from dash import Dash, html, dcc, callback, Output, Input, no_update, State, ALL, MATCH, get_asset_url
import dash_ag_grid as dag
import pandas as pd
import plotly.express as px
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import dash_auth
from dash_auth import BasicAuth
import math
from scipy import stats
import numpy as np
from sqlalchemy import create_engine, text
import dash_bio as dashbio
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, leaves_list
import plotly.figure_factory as ff
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Uncommnent this line out for multiple pages web app
dash.register_page(__name__)

# drug name, cell line and mysql table dictionary
DATASETS2 = {
    'MCF7': ['CDK2i_ARTNET-100nM_VK061625'],
    'T47D-WT': ['Atrimo-50nM_ARTNET-100nM_VK0001', 'Palbo-50nM_ARTNET-100nM_VK0001'],
    'T47D-CDK6-OE': ['Atrimo-50nM_ARTNET-100nM_VK0002', 'Palbo-50nM_ARTNET-100nM_VK0002'],
    'MDA-MB-468': ['shCCNE2_ARTNET-250nM_KS0001', 'shSKP2_ARTNET-200nM_VK060426', 'shCCNE1_ARTNET-250nM_KS0005'],
    'MDA-MB-436': ['shCCNE2_ARTNET-250nM_KS0002'],
    'HCC-1937': ['shCCNE2_ARTNET-250nM_KS0004', 'shCCNE2_ARTNET-500nM_KS0003', 'shCCNE1_ARTNET-250nM_KS0005'],
    'HCC-1806': ['shSKP2_ARTNET-100nM_VK18061911']
}

# Cell size datasets
DATASETS3 = {
    'MCF7': ['CDK2i_ARTNET-100nM_VK061625'],
}

#  DATABASE CONNECTION SETUP
USER = "vectra"
PASSWORD = "ritho9Ng"
HOST = "10.126.0.41"
PORT = "3306"
DATABASE = "drug_screen"

# Create a reusable SQLAlchemy connection engine
engine = create_engine(f"mysql+pymysql://{USER}:{PASSWORD}@{HOST}:{PORT}/{DATABASE}")

# Get available drug options
with engine.connect() as conn:
    try:
        query = text(f"SELECT distinct Drugs FROM `ARTNET` order by Drugs")

        df_table = pd.read_sql(query, conn)

        if not df_table.empty:
            drugs = df_table['Drugs'].to_list()

    except Exception as e:
        print(f"Could not read table. Error: {e}")

# Get available drug class options
with engine.connect() as conn:
    try:
        query = text(f"SELECT distinct Class FROM `drug_annotation` order by Class")

        df_table = pd.read_sql(query, conn)

        if not df_table.empty:
            drug_class = df_table['Class'].to_list()

    except Exception as e:
        print(f"Could not read table. Error: {e}")
        
# Get available screens
with engine.connect() as conn:
    try:
        query = text("""
                SELECT table_name
                FROM information_schema.tables
                WHERE table_schema = DATABASE()
                  AND table_name LIKE '%_ARTNET_%'
                  AND table_name NOT LIKE '%Cell-size%'
                """)

        table_names = pd.read_sql(query, conn)
        

        if not table_names.empty:
            screen_names = table_names.iloc[:, 0].tolist()

    except Exception as e:
        print(f"Could not read table. Error: {e}")
        

# footer_markdown=dcc.Markdown("""
#     This is the footer text...
# """)

layout = dbc.Container([
    dbc.Row(html.H1("Explore Drug Screen Data Sets"), style = {'textAlign' : 'center'}),
    dbc.Row([
        # Left Panel (Input Placeholder)
        dbc.Col([
            html.H4("Controls"),
            html.Div(id="left-panel-inputs")  
        ], width=2, style={"background-color": "#f8f9fa", "padding": "20px"}),
        
        # Right Panel (Tabs & Graph Placeholder)
        dbc.Col([
            dbc.Tabs(id="right-panel-tabs", active_tab="growth_tab", children=[
                dbc.Tab(label="Cell Growth Curves", tab_id="growth_tab"),
                dbc.Tab(label="Drug Effect Similarity", tab_id="drug_similarity_tab"),
                dbc.Tab(label="Cell Size Change", tab_id="cell_size_tab"),
                #dbc.Tab(label="Drug Correlation", tab_id="drug_correlation_tab"),
            ]),
            html.Div(id="right-panel-content", style={"padding": "20px"})
        ], width=10),
        dcc.Store(id='mean-cell-size-data'),
    ]),
    # dbc.Row(footer_markdown),
], fluid=True)


# 3. Callback 1: Dynamically update the Left Panel based on selected tab
@callback(
    Output("left-panel-inputs", "children"),
    Input("right-panel-tabs", "active_tab")
)
def render_left_panel(active_tab):
    if active_tab == "growth_tab":
        return [
            html.P(
                "1). Select a Cell Line Model",
            ),
            dcc.Dropdown(
                id={"type": "cell-line-dropdown2", "index": "cellline"},
                options=[{'label': cell, 'value': cell} for cell in DATASETS2.keys()],
                placeholder="Select a cell line model...",
                clearable=True,
            ),
            html.Br(),
            html.P(
                "2). Select Treatment(s)/Condition(s)",
            ),
            dcc.Dropdown(
                id={"type": "condition-dropdown", "index": "condition"},
                placeholder="Select a treatment/condition...",
                clearable=True,
            ),
            html.Br(),
            html.P(
                "3). Select Drug(s)",
            ),
            dcc.Dropdown(
                id={"type": "drug-name-dropdown", "index": "drug"},
                options=[{'label': drug, 'value': drug} for drug in drugs],
                placeholder="Select by drug(s)/treatment...",
                clearable=True,
                multi=True,
            ),
            html.Br(),
            html.P(
                "or Drug Class(es)",
            ),
            html.Br(),
            dcc.Dropdown(
                id={"type": "drug-class-dropdown", "index": "drug-class"},
                options=[{'label': dc, 'value': dc} for dc in drug_class],
                placeholder="Select by drug class(es)...",
                clearable=True,
                multi=True,
            ),
            html.Br(),            
            html.Button("Reset", id="reset-button", n_clicks=0, 
                style={
                    "backgroundColor": "#f3ac0c",
                    "color": "white",
                    "border": "none",
                    "padding": "5px 15px",
                    "borderRadius": "8px",
                    "cursor": "pointer",
                    "display": "flex",
                    "justifyContent": "flex-end"

                    }
            ),
            html.Br(),
            dcc.Markdown(
                """
                Click the Reset button above to make new drug selection choice.
                """
            ),
            html.H4("Filters"),
            html.Label("Filters results by end-point fold-change difference", style={
                    "textAlign": "center", 
                    "display": "block",
                }),
            dcc.Dropdown(
                id={"type": "delta-threshold", "index": "delta-th"},
                options=[{'label': delta, 'value': delta} for delta in [round(x * 0.5, 1) for x in range(30)]],
                placeholder="Select delta threshold...",
                clearable=True,
            ),
        ]
        
    elif active_tab == "drug_similarity_tab":
        return [
            html.Label("Select Cell Line/Condition:"),
            dcc.Dropdown(
                id={"type": "drug-screen-dropdown", "index": "drug-screen"},
                options=[{'label': screen, 'value': screen} for screen in screen_names],
                placeholder="Select a drug screen...",
                clearable=True,
                multi=True
            ),
            html.Br(),
            html.Label("Select Drug Class(es):"),
            dcc.Dropdown(
                id={"type": "drug-class-dropdown2", "index": "drug-class2"},
                options=[{'label': dc, 'value': dc} for dc in drug_class],
                placeholder="Select by drug class(es)...",
                clearable=True,
                multi=True,
            ),
            html.Br(),
            html.H4("Refine Selections"),
            html.Label("Filter data by viability no greater than", style={
                    "textAlign": "center", 
                    "display": "block",
                    #"fontWeight": "bold"  # Optional: makes it stand out as a clean separator
                }),
            dcc.Dropdown(
                id={"type": "fc-threshold", "index": "fc-th"},
                options=[{'label': ft, 'value': ft} for ft in [round(x * 0.1, 1) for x in range(11)]],
                placeholder="Select viability threshold...",
                clearable=True,
            ),
            html.Label("and in at at least", style={
                    "textAlign": "center", 
                    "display": "block",
                    #"fontWeight": "bold"  # Optional: makes it stand out as a clean separator
                }),
            dcc.Dropdown(
                id={"type": "screen-count-threshold", "index": "screen-count-th"},
                options=[{'label': sc, 'value': sc} for sc in range(3, len(screen_names)+1)],
                placeholder="Select by count threshold...",
                clearable=True,
            ),
            html.Label("screens"),
            html.Br(),
            
        ]
    elif active_tab == "cell_size_tab":
        return [
            html.P(
                "1). Select a Cell Line Model",
            ),
            dcc.Dropdown(
                id="cell-line-dropdown3",
                options=[{'label': cell, 'value': cell} for cell in DATASETS3.keys()],
                placeholder="Select a cell line model...",
                clearable=True,
            ),
            html.Br(),
            html.P(
                "2). Select Treatment(s)/Condition(s)",
            ),
            dcc.Dropdown(
                id="condition-dropdown3",
                placeholder="Select a treatment/condition...",
                clearable=True,
            ),
            html.Br(),
            dcc.Dropdown(
                id={"type": "drug-class-dropdown"},
                options=[{'label': dc, 'value': dc} for dc in drug_class],
                style={"display": "none"}
            ),
        ]
        
    return "No tab selected"



# Drop down options for conditions, tab1 and tab2
@callback(
    # FIX: This must match the dictionary format from your render_left_panel layout
    Output({"type": "condition-dropdown", "index": "condition"}, "options"),
    Input({"type": "cell-line-dropdown2", "index": "cellline"}, "value")
)
def update_condition_dropdown(selected_cell):
    
    if not selected_cell:
        return dash.no_update
        #return []
    
    # Perform the dictionary lookup
    val = DATASETS2.get(selected_cell, [])
    
    # Generate options safely
    options = [{'label': cell, 'value': cell} for cell in ([val] if isinstance(val, str) else val)]
    
    return options

# Dropdown for tab3, top pane
@callback(
    Output("order-by-dropdown", "options"),
    Input("cell-line-dropdown3", "value"),
    Input("condition-dropdown3", "value"),
)
def update_order_by_dropdown(selected_cell, selected_condition):
    
    if not selected_cell or not selected_condition:
        return []
    
    table_name = f"Cell-size_{selected_cell}_{selected_condition}"
    
    with engine.connect() as conn:
        try:
            query = text(f"""
                SELECT DISTINCT `Condition`
                FROM `{table_name}`
            """)

            result = conn.execute(query)
            conditions = [row[0] for row in result]

        except Exception as e:
            print(f"Could not read table. Error: {e}")
            conditions = []
    
    # Generate options safely 
    options = [{'label': val, 'value': val} for val in conditions]
    
    # add 'Difference' to the options
    options.append('Difference')
    
    return options

# Dropdown for tab3, bottom pane
@callback(
    # FIX: This must match the dictionary format from your render_left_panel layout
    Output("condition-dropdown3", "options"),
    Input("cell-line-dropdown3", "value"),
)
def update_condition_dropdown3(selected_cell):
    
    if not selected_cell:
        return dash.no_update
    
    # Perform the dictionary lookup
    val = DATASETS3.get(selected_cell, [])
    
    # Generate options safely
    options = [{'label': cell, 'value': cell} for cell in ([val] if isinstance(val, str) else val)]
    
    return options


# Dropdown for tab3 drug selection, make it inactive if parents dropdown has not be selected
@callback(
    Output("drug-name-dropdown3", "options"),
    Output("drug-name-dropdown3", "disabled"),
    Input("condition-dropdown3", "value")
)
def enable_drug_name_dropdown(condition_dropdown3):

    if not condition_dropdown3:

        return dash.no_update, True
    
    options=[{'label': drug, 'value': drug} for drug in drugs]
    
    return options, False


@callback(
    # FIX: Updated to match the dictionary IDs from your layout
    Output({"type": "drug-name-dropdown", "index": "drug"}, "disabled"),
    Output({"type": "drug-class-dropdown", "index": "drug-class"}, "disabled"),
    
    # Inputs (These look good as they match the layout IDs perfectly)
    Input({"type": "condition-dropdown", "index": "condition"}, "value"),
    Input({"type": "drug-name-dropdown", "index": "drug"}, "value"),
    Input({"type": "drug-class-dropdown", "index": "drug-class"}, "value"),
)
def toggle_dropdowns(val0, val1, val2):
 
    if not val0:
        disabled_1 = True
        disabled_2 = True
    else:
        # Note on handling multi-select lists:
        # Since drug and drug-class are multi=True dropdowns, 
        # Dash returns an empty list [] instead of None when cleared.
        has_val1 = val1 is not None and len(val1) > 0
        has_val2 = val2 is not None and len(val2) > 0
        
        # If dropdown-2 (drug class) has a value → disable dropdown-1 (drug name)
        # If dropdown-1 (drug name) has a value → disable dropdown-2 (drug class)
        disabled_1 = has_val2
        disabled_2 = has_val1

    return disabled_1, disabled_2


@callback(
    Output({"type": "drug-name-dropdown", "index": "drug"}, "value"),
    Output({"type": "drug-class-dropdown", "index": "drug-class"}, "value"),
    Input("reset-button", "n_clicks"),
    prevent_initial_call=True
)
def reset_dropdowns(n_clicks):
    return None, None

@callback(
    Output("right-panel-content", "children"),
    Input("right-panel-tabs", "active_tab"),
    # CHANGE: Swapped out specific indices for the ALL wildcard selector
    Input({"type": "cell-line-dropdown2", "index": ALL}, "value"),
    Input({"type": "condition-dropdown", "index": ALL}, "value"),
    Input({"type": "drug-name-dropdown", "index": ALL}, "value"),
    Input({"type": "drug-class-dropdown", "index": ALL}, "value"),
)
def render_right_panel(active_tab, cell_line_list, condition_list, drug_list, class_list):
    # CRITICAL: Because we used ALL, Dash bundles inputs into lists. 
    # Let's extract the actual string values safely. If the tab hasn't rendered them yet, 
    # the list will be empty, so we default to None.
    cell_line = cell_line_list[0] if cell_line_list else None
    condition = condition_list[0] if condition_list else None
    selected_drugs = drug_list[0] if drug_list else None
    selected_classes = class_list[0] if class_list else None
    # 1. If we are on the growth tab, render the card layout structure
    if active_tab == "growth_tab":
        
        card_growth_rate = dbc.Card([
            dbc.CardHeader(html.H4("Cell Growth Line Plots")),
            dbc.CardBody([
                dbc.Row([
                    dbc.Col([
                        dcc.Loading(
                            dcc.Graph(id='cell-growth-lineplot'),
                        ),
                    ], width=12),
                ]),
                dbc.Row([
                    dbc.Col([
                        dcc.Loading(
                            dcc.Graph(id='endpoint-heatmap'),
                        ),
                    ], width=12),
                ]),
            ]),
            dbc.CardFooter([
                dcc.Markdown(
                    """
                    Difference is the fold-change value calculated at the same time-point between experimental conditions.
                    """
                ),
            ]),
        ], className="shadow")
        
        return card_growth_rate

    # 2. Render content for the other tab
    elif active_tab == "drug_similarity_tab":

        card_drug_similarity = html.Div([

            # -------------------------------------------------------------
            # ROW 1: Top Container Card
            # -------------------------------------------------------------
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardHeader(html.H4("Hierarchical Clustering")),
                        dbc.CardBody([
                            # Inside the card body: Exactly one row and one column
                            dbc.Row([
                                dbc.Col([
                                    # html.H5("Top Card Title"),
                                    dcc.Loading(dcc.Graph(id='endpoint-viability-heatmap')),
                                ], width=12) # Spans full internal card width
                            ])
                        ])
                    ], className="shadow-sm")
                ], width=12) # Spans full horizontal screen width
            ], className="mb-3"), # mb-3 adds margin-bottom spacing to separate it from Row 2

            # -------------------------------------------------------------
            # ROW 2: Middle Container Card
            # -------------------------------------------------------------
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardHeader(html.H4("PCA Analysis")),
                        dbc.CardBody([
                            # Inside the card body: Exactly one row and one column
                            dbc.Row([
                                dbc.Col([
                                    #html.H5("Middle Card Title"),
                                    dcc.Loading(                      
                                        dcc.Graph(id='pca-plot'),
                                    ),
                                ], width=12)
                            ])
                        ])
                    ], className="shadow")
                ], width=12)
            ], className="mb-3"), # mb-3 separates Row 2 from Row 3

            # -------------------------------------------------------------
            # ROW 3: Bottom Container Card
            # -------------------------------------------------------------
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardHeader(html.H4("Condition Correlation Scatterplots")),
                        dbc.CardBody([
                            dbc.Row([
                                dbc.Col([
                                    #html.H5("Bottom Card Title"),
                                    dcc.Loading(                      
                                        dcc.Graph(id='master_scatterplot'),
                                    ),
                                ], width=12)
                            ])
                        ]),
                        dbc.CardFooter([
                            dcc.Markdown(
                                """
                                Values on x- and y-axis are fold-change (cell count ratio at end-time point vs. time zero).
                                Diagonal line is x = y. Any drug resides on this line has no effect on cell growth under the represive condition used.
                                """
                            ),
                        ]),
                    ], className="shadow-sm"),
                    
                    
                ], width=12)
            ], className="mb-3")

        ])

        return card_drug_similarity
    
    elif active_tab == "cell_size_tab":
        card_cell_size = html.Div([
            # -------------------------------------------------------------
            # ROW 1: Top Container Card
            # -------------------------------------------------------------
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardHeader(html.H4("Heatmap of Mean Cell Size")),
                        dbc.CardBody([
                            dbc.Row([
                                dbc.Col([
                                    dcc.Loading(
                                        
                                        children=[
                                            dcc.Store(id="mean-cell-size-data"),
                                            dcc.Graph(id="cell-size-heatmap")
                                        ]
                                    ),                                  
                                ], width=10) # Spans full internal card width
                            ]),
                            dbc.Row([
                                dbc.Col([
                                    html.P(
                                        "3). Re-order heatmap column by",
                                    ),
                                    dcc.Dropdown(
                                        id="order-by-dropdown",
                                        options=[],
                                        clearable=True,
                                        multi=False,
                                        disabled=True,
                                    ),       
                                    html.Br(),
                                ], width=3),
                            ]),
                        ]),
                    ], className="shadow-sm")
                ], width=12) # Spans full horizontal screen width
            ], className="mb-3"), # mb-3 adds margin-bottom spacing to separate it from Row 2

            # -------------------------------------------------------------
            # ROW 2: Second Container Card
            # -------------------------------------------------------------
            dbc.Row([
                dbc.Col([
                    dbc.Card([
                        dbc.CardHeader(html.H4("Violinplots and Segmentation Images")),
                        dbc.CardBody([
                            dbc.Row([
                                dbc.Col([
                                    html.P(
                                        "3). Select a Drug",
                                    ),
                                    dcc.Dropdown(
                                        id="drug-name-dropdown3",
                                        placeholder="Select a drug",
                                        clearable=True,
                                        multi=False,
                                    ),

                                    html.Br(),
                                ], width=3),
                            ]),
                            # Inside the card body: Exactly one row and one column
                            dbc.Row([
                                dbc.Col([
                                    #html.H5("Middle Card Title"),
                                    dcc.Loading(                      
                                        dcc.Graph(id='cell-size-violinplot'),
                                    
                                    ),
                                ], width=5)
                            ]),

                            dbc.Row([
                                dbc.Col([
                                    html.H5("Control", style={"textAlign": "center", "marginBottom": "5px"}),
                                    dcc.Loading(
                                        html.Img(
                                            id="cell-image-control",
                                            style={"width": "100%"}
                                        ),
                                    ),
                                ], width=12),
                            ], style={"marginTop": "0px", "marginBottom": "0px"}),

                            # Row 3: Treated (className="mt-0 g-0" removes top spacing and gutters)
                            dbc.Row([
                                dbc.Col([
                                    html.H5(
                                        id="image-title",
                                        style={"textAlign": "center", "marginBottom": "5px"}
                                    ),
                                    dcc.Loading(
                                        html.Img(
                                            id="cell-image-treated",
                                            style={"width": "100%"}
                                        ),
                                    ),
                                ], width=12),
                            ], style={"marginTop": "0px", "marginBottom": "0px"}),
                        ]),
                    ], className="shadow")
                ], width=12)
            ], className="mb-3"),
        ])

        
        return card_cell_size
    
    return "Select a tab to view content."



@callback(
    Output('cell-growth-lineplot', 'figure'),
    Output('endpoint-heatmap', 'figure'),
    Input({"type": "cell-line-dropdown2", "index": "cellline"}, "value"),
    Input({"type": "condition-dropdown", "index": "condition"}, "value"),
    Input({"type": "drug-name-dropdown", "index": "drug"}, "value"),
    Input({"type": "drug-class-dropdown", "index": "drug-class"}, "value"),
    Input({"type": "delta-threshold", "index": "delta-th"}, "value"),
)
def update_tab1_graphs(selected_cell, selected_condition, selected_drugs, selected_drug_class, delta_th):
   
    if not selected_cell or not selected_condition or not (selected_drugs or selected_drug_class):
        return dash.no_update, dash.no_update
        
    # Initialize a list to pool data frames fetched from individual tables
    collected_data = {}
    
    # Establish a clean context connection to the database
    with engine.connect() as conn:
        
        # 1. get control data table name
        parts = selected_condition.split("_")
        
        if selected_condition.startswith('sh'):
            parts[0] = 'Ctrl'
        else:
            parts[0] = 'DMSO'

        table_name1 = selected_cell + '_' + "_".join(parts)
        table_name2 = selected_cell + '_' + selected_condition
        
        if selected_drugs:
            # 1. Dynamically create placeholders: ":d0, :d1, :d2"
            placeholders = ", ".join(f":d{i}" for i in range(len(selected_drugs)))

            # 3. Construct the parameter dictionary: {"d0": "aspirin", "d1": "ibuprofen"}
            params_dict = {f"d{i}": drug for i, drug in enumerate(selected_drugs)}
            
        elif selected_drug_class:
            # 1. Dynamically create placeholders: ":d0, :d1, :d2"
            placeholders = ", ".join(f":d{i}" for i in range(len(selected_drug_class)))

            # 3. Construct the parameter dictionary: {"d0": "aspirin", "d1": "ibuprofen"}
            params_dict = {f"d{i}": dc for i, dc in enumerate(selected_drug_class)}
            query_str = f"SELECT * FROM drug_annotation WHERE Class IN ({placeholders})"
            query = text(query_str)
            dc_table = pd.read_sql(query, conn, params=params_dict)
            
            placeholders = ", ".join(f":d{i}" for i in range(len(dc_table['Drugs'])))
            params_dict = {f"d{i}": drug for i, drug in enumerate(dc_table['Drugs'])}
                                 
        for idx, table_name in enumerate([table_name1, table_name2]):
            try:
                # 2. Build the query string using the generated placeholders
                query_str = f"SELECT * FROM `{table_name}` WHERE Drugs IN ({placeholders})"
                query = text(query_str)

                # 4. Pass the dictionary explicitly to the params keyword
                df_table = pd.read_sql(query, conn, params=params_dict)

                if not df_table.empty:
                    # Inject label tracking metadata back into the frame
                    if idx == 0:
                        df_table['Condition'] = parts[0]       
                    else:
                        df_table['Condition'] = selected_condition   
                        
                    collected_data[idx] = df_table

            except Exception as e:
                fig = go.Figure()
                fig.update_layout(title=f"Database error {e}", template="plotly_white")
                return fig, fig
            
    # Pull a high-contrast qualitative color palette array
    palette = px.colors.qualitative.Plotly

    # Define a custom hex-to-rgba converter for trace shading
    def hex_to_rgba(hex_str, opacity):
        hex_str = hex_str.lstrip('#')
        rgb = tuple(int(hex_str[i:i+2], 16) for i in (0, 2, 4))
        
        return f"rgba({rgb[0]}, {rgb[1]}, {rgb[2]}, {opacity})"
    
     # Configuration for subplots, we will make three subplots: 1. DMSO, 2. Treated, 3. difference
    cols = 3 
    rows = 1
    
    # Initialize subplots
    fig = make_subplots(
        rows=rows, 
        cols=cols, 
        subplot_titles=['DMSO', selected_condition.split("_")[0], 'Difference'],
        vertical_spacing=0.2, 
        horizontal_spacing=0.05
    )
    
    # Make growth curves for control and drug treated separately first
    
    # collecte summary data for later use
    df_summary_saved = []
    df_joined = None
    
    for idx, df_table in collected_data.items():
                                 
        df_summary = (
            df_table.groupby(["Time", "Drugs"])["fold_change"]
            .agg(Mean="mean", StdErr=lambda x: stats.sem(x, ddof=1) if len(x) > 1 else 0)
            .reset_index()
        )
        df_summary["StdErr"] = df_summary["StdErr"].fillna(0)
        
        if len(df_summary_saved) == 0:
            df_summary_saved = df_summary
        else:
            df_joined = pd.merge(df_summary_saved, df_summary, on=["Time", "Drugs"], how="inner")
            
            # Calculate the fold change differece between treated and control
            df_joined['Delta'] = df_joined['Mean_x'] - df_joined['Mean_y']
            
            # Filter data if a delta threshold is chosen
            if delta_th:
                to_keep = df_joined[df_joined['Delta'] > delta_th]
                df_joined = df_joined[df_joined['Drugs'].isin(to_keep['Drugs'])]
          
    
    legend_identifier = 'legend'
        
    # Repeat to make line plots for each condition
    for idx, df_table in collected_data.items():
                                 
        df_summary = (
            df_table.groupby(["Time", "Drugs"])["fold_change"]
            .agg(Mean="mean", StdErr=lambda x: stats.sem(x, ddof=1) if len(x) > 1 else 0)
            .reset_index()
        )
        df_summary["StdErr"] = df_summary["StdErr"].fillna(0)      
            
        if delta_th:
            df_summary = df_summary[df_summary['Drugs'].isin(to_keep['Drugs'])]
            
        if idx == 0:
            current_col = 1
        else:
            current_col = 2
        
        for i, drug in enumerate(df_summary['Drugs'].unique()):

            group_data = df_summary[df_summary['Drugs'] == drug].sort_values('Time')

            # Assign specific matching palette colors to this category group
            base_color_hex = palette[i % len(palette)]
            fill_color_rgba = hex_to_rgba(base_color_hex, opacity=0.15)

            fig.add_trace(go.Scatter(
                x=group_data['Time'],
                y=group_data['Mean'],
                mode='lines+markers',
                name=drug,
                showlegend=False,
                line=dict(color=base_color_hex, width=2.5),
                marker=dict(size=6, symbol='circle'),

                # --- CRITICAL LEGEND ROUTING CONFIGURATION ---
                legend=legend_identifier,       # Binds this trace to a custom legend box
                legendgroup=drug,          # Groups related trace properties together
                error_y=dict(
                    type='data',                  # Tells Plotly to read values from an array
                    array=group_data['StdErr'],   # Sets the error bar heights above/below the mean
                    visible=True,                 # Renders the error bars visible
                    thickness=1.5,                # Line width of the error bars
                    width=4,                      # Width of the horizontal crossbar cap (T-bar)
                    color=base_color_hex               # Automatically matches the line color
                )
            ), row=1, col = current_col)
            

    current_col = 3

    # Make the third plot to show the treatment difference
    for i, drug in enumerate(df_joined['Drugs'].unique()):
        group_data = df_joined[df_joined['Drugs'] == drug].sort_values('Time')

        # Assign specific matching palette colors to this category group
        base_color_hex = palette[i % len(palette)]
        fill_color_rgba = hex_to_rgba(base_color_hex, opacity=0.15)

        fig.add_trace(go.Scatter(
            x=group_data['Time'],
            y=group_data['Delta'],
            mode='lines+markers',
            name=drug,
            showlegend=True,        
            
            # --- Single Line Hover Configuration ---
            text=[drug] * len(group_data),
            hovertemplate="%{text}<extra></extra>", 
           
            line=dict(color=base_color_hex, width=2.5),
            marker=dict(size=6, symbol='circle'),

            # --- CRITICAL LEGEND ROUTING CONFIGURATION ---
            legend=legend_identifier,       # Binds this trace to a custom legend box
            legendgroup=drug,          # Groups related trace properties together

        ), row=1, col = current_col)
        
    fig.update_layout(hovermode='closest')
    
    
    fig.update_layout(
        height=400,
        width=400* 3,       
        template="plotly_white",
        #hovermode="x unified",
    )
    fig.update_xaxes(title_text="Time Points")
    fig.update_yaxes(title_text="Fold Change Difference")
    df_joined
 
    # Heatmap
    data_for_heatmap = df_joined[df_joined['Time'] == df_joined.Time.max()]
    data_for_heatmap = data_for_heatmap[['Drugs', 'Mean_x', 'Mean_y', 'Delta']]
    data_for_heatmap.rename(columns = {'Mean_x': 'DMSO', 'Mean_y': selected_condition}, inplace = True)
    data_for_heatmap = data_for_heatmap.set_index('Drugs').T

    # Order by Delta
    data_for_heatmap = data_for_heatmap.sort_values(by="Delta", ascending=False, axis=1)

    heatmap = px.imshow(data_for_heatmap, height = 300)
 
    heatmap.update_layout(
        title='Drug Effect Heatmap',
        title_x=0.5,  # Centers the title
        template="plotly_white"
    )

    return fig, heatmap
   

@callback(
    Output('endpoint-viability-heatmap', 'figure'),
    Output('pca-plot', 'figure'),
    Output('master_scatterplot', 'figure'),
    Input({"type": "drug-screen-dropdown", "index": "drug-screen"}, "value"),
    Input({"type": "drug-class-dropdown2", "index": "drug-class2"}, "value"),
    State({"type": "fc-threshold", "index": "fc-th"}, "value"),
    Input({"type": "screen-count-threshold", "index": "screen-count-th"}, "value"),
)
def update_tab2_heatmap(selected_drug_screen, selected_drug_class, fold_change, screen_count):    
    if not selected_drug_screen or not selected_drug_class:
        return dash.no_update
    if len(selected_drug_screen) < 3:
        return dash.no_update
        
    # retrieve end-point viability data from database
    result = pd.DataFrame()

    for table_name in selected_drug_screen:
        with engine.connect() as conn:
            try:
                # read db table to get last time point values
                query = text(f"select MAX(Time) from `{table_name}`")
                df = pd.read_sql(query, conn)
                max_time = df.iloc[0,0]

                # retrived data from db based on last time point
                query = text(f"select Wells,fold_change from `{table_name}` where Time = {max_time}")

                df_table = pd.read_sql(query, conn)

            except Exception as e:
                print(f"Could not read table. Error: {e}")

            # rename column
            df_table.rename(columns={"fold_change": table_name}, inplace = True)

            if result.empty:
                result = df_table
            else:
                result = pd.merge(result, df_table, on = 'Wells', how='inner')

    
    with engine.connect() as conn:
        # Get drug library data
        try:
            query = text(f"select * from ARTNET")
            drug_library = pd.read_sql(query, conn)
        except Exception as e:
            print(f"Could not read table. Error: {e}")
        # get drug annotation data
        try:
            query = text(f"select * from drug_annotation")
            drug_annotation = pd.read_sql(query, conn)
        except Exception as e:
            print(f"Could not read table. Error: {e}")

    # merge with drug screen data
    result = pd.merge(drug_library, result, on='Wells', how='inner')
    result = result.drop(columns=['Wells'])

    # Some drugs have replicates, so take mean values for each drug
    result = result.groupby('Drugs', as_index=False).mean()
    result = result.set_index('Drugs')

    # get reference cell counts for calculating fold change or viability
    control_row = result.loc['DMSO']
    normalized = result.div(control_row, axis=1)
    normalized = pd.merge(drug_annotation, normalized, on='Drugs', how='right').set_index('Drugs')
    
    # Fix NaN values
    normalized['Class'] = normalized['Class'].fillna('Other').astype(str)
    normalized['Mechanism'] = normalized['Mechanism'].fillna('Other').astype(str)
    
    normalized_jsonified = normalized.to_json(date_format='iso', orient='split')
    normalized = normalized[normalized.Class.isin(selected_drug_class)]
    
    if normalized.shape[0] < 3:
         return  px.scatter(title="Minimum number of drugs is 3, please select more drugs"),px.scatter(title="")

    if screen_count and fold_change:
        sample_cols = normalized.select_dtypes(include='number').columns.tolist()
        tmp = normalized[sample_cols]
        to_keep = (tmp <= fold_change).sum(axis = 1)
        to_keep = to_keep >= screen_count
        normalized = normalized[to_keep]

    sample_cols = normalized.select_dtypes(include='number').columns.tolist()

    # Isolate strictly numeric data
    matrix_df = normalized[sample_cols]
    
    #=======================================
    # Heatmap
    #=======================================
    # --------------------------------------------------------------------------
    # A. COMPUTE MATHEMATICAL HIERARCHICAL CLUSTERING
    # --------------------------------------------------------------------------
    # Row Linkage (using complete linkage & Euclidean distance)
    row_dist = pdist(matrix_df.values, metric='euclidean')
    row_linkage = linkage(row_dist, method='complete')
    row_order = leaves_list(row_linkage)
    
    # Column Linkage
    col_dist = pdist(matrix_df.values.T, metric='euclidean')
    col_linkage = linkage(col_dist, method='complete')
    col_order = leaves_list(col_linkage)
    
    # Reorder our source labels and datasets based on leaf sorting
    ordered_rows = [matrix_df.index[i] for i in row_order]
    ordered_cols = [matrix_df.columns[i] for i in col_order]
    
    clustered_matrix = matrix_df.values[row_order, :][:, col_order]
    
    # Reorder Metadata matching row tree changes
    ordered_classes = normalized['Class'].iloc[row_order]
    
    # --------------------------------------------------------------------------
    # B. GENERATE TREE STRUCTURE DENDROGRAM FIGURES
    # --------------------------------------------------------------------------
    fig_col_dendro = ff.create_dendrogram(matrix_df.values.T, orientation='bottom', linkagefun=lambda x: col_linkage)
    fig_row_dendro = ff.create_dendrogram(matrix_df.values, orientation='left', linkagefun=lambda x: row_linkage)
    
    # --------------------------------------------------------------------------
    # C. MAP COVARIATE CATEGORIES TO COLOR INTEGRATION
    # --------------------------------------------------------------------------
    unique_classes = sorted(list(normalized['Class'].unique()))
    class_to_int = {cls: idx for idx, cls in enumerate(unique_classes)}
    int_classes = ordered_classes.map(class_to_int).values.reshape(-1, 1)
    
    color_palette = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6','#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080','#cccccc']
    plotly_colorscale = []
    num_classes = len(unique_classes)
    for i, cls in enumerate(unique_classes):
        plotly_colorscale.append([i / num_classes, color_palette[i % len(color_palette)]])
        plotly_colorscale.append([(i + 1) / num_classes, color_palette[i % len(color_palette)]])

    # --------------------------------------------------------------------------
    # D. ASSEMBLE SUBPLOT GRID SYSTEM
    # --------------------------------------------------------------------------
    # dynamically specify column with ratios
    width=matrix_df.shape[1]*25 
    annot_width=20 # set to fixed pixel width
    proportion_col1=annot_width/width
    proportion_col2=(width - annot_width)/width
    
    fig = make_subplots(
        #rows=2, cols=3,
        rows=1, cols=2,
        shared_xaxes=False,
        shared_yaxes=True,
        #row_heights=[0.12, 0.88],
        #column_widths=[0.04, 0.82, 0.14],
        # column_widths=[0.04, 0.96],
        column_widths=[proportion_col1, proportion_col2],
        #vertical_spacing=0.001,
        horizontal_spacing=0.002,
    )

#     # 1. Add Column Dendrogram Traces (Top Plot, Center Right)
#     for trace in fig_col_dendro['data']:
#         fig.add_trace(go.Scatter(trace, mode='lines', marker=dict(color='#7f7f7f'), showlegend=False), row=1, col=2)

    # 2. Add Row Covariate Track (Bottom Plot, Left Column)
    fig.add_trace(
        go.Heatmap(
            z=int_classes,
            x=["Class"],
            y=ordered_rows,
            colorscale=plotly_colorscale,
            showscale=False,
            hoverongaps=False,
            hovertemplate="Drug: %{y}<br>Class: %{customdata}<extra></extra>",
            customdata=ordered_classes.values.reshape(-1, 1)
        ),
        row=1, col=1
    )

    # 3. Add Primary Expression Matrix Heatmap (Bottom Plot, Middle Column)
    fig.add_trace(
        go.Heatmap(
            z=clustered_matrix,
            zmin=0,
            zmax=1,
            x=ordered_cols,
            y=ordered_rows,
            colorscale="Viridis",
            showscale=True,
            colorbar=dict(title="Fold_Change<br>(Viability)", thickness=15, len=0.4, y=0.3)
        ),
        row=1, col=2
    )

#     # 4. Add Row Dendrogram Traces (Bottom Plot, Right Column)
#     for trace in fig_row_dendro['data']:
#         fig.add_trace(go.Scatter(trace, mode='lines', marker=dict(color='#7f7f7f'), showlegend=False), row=2, col=3)

    # --------------------------------------------------------------------------
    # E. ADD DISCRETE EXPLANATORY LEGENDS
    # --------------------------------------------------------------------------
    for i, cls in enumerate(unique_classes):
        fig.add_trace(
            go.Scatter(
                x=[None], y=[None], mode='markers',
                marker=dict(color=color_palette[i % len(color_palette)], symbol='square', size=14),
                name=cls, legendgroup='Drug Class Metadata', showlegend=True
            ), row=1, col=1
        )

    # --------------------------------------------------------------------------
    # F. STYLING OVERRIDES & CLEANUP
    # --------------------------------------------------------------------------
    fig.update_layout(
        height=matrix_df.shape[0]*20 + 150 if matrix_df.shape[0] < 100 else matrix_df.shape[0]*5 + 150,
        width=matrix_df.shape[1]*25 + 600,
        template="plotly_white",
        legend=dict(
            title="Drug Class",
            yanchor="top", y=0.95,
            xanchor="left", x=1.0 #1.12
        ),
        # Clean up tree background visibility
        yaxis=dict(visible=True),   # Column tree vertical axis 
        xaxis5=dict(visible=True),  # Row tree horizontal axis
        xaxis6=dict(visible=False),
    )
    
    # Clear visual chart grid lines for the dendrogram areas
    fig.update_xaxes(showgrid=False, row=1, col=2)
    fig.update_yaxes(showgrid=False, row=1, col=2)

    
    fig.update_xaxes(
        showticklabels=True,      # Forces the column names to render
        tickangle=45,             # Angles the text at 45 degrees so long names don't overlap
        tickfont=dict(size=10),   # Adjust text font size to fit your sample names cleanly
        row=1, col=2
    )


    # Optional: Clean up layout padding to ensure long column labels don't get cut off
    fig.update_layout(
        margin=dict(b=150) # Adds 150px of extra whitespace buffer at the very bottom
    )

    fig.update_yaxes(tickson="boundaries", row=1, col=1)

    
    # Force the bottom heatmap row subplot (Row 1, Col 2) to use your exact text strings
    fig.update_xaxes(
        showticklabels=True,
        type='category',            # Forces Plotly to treat columns as discrete categories
        tickmode='array',            # Tells Plotly we are passing manual arrays
        tickvals=list(range(len(ordered_cols))), # Maps a numerical coordinate slot to each column
        ticktext=ordered_cols,       # Inject your actual sample string names list here
        tickangle=45,               # Rotates the long labels diagonally
        tickfont=dict(size=10),      # Sets a clean text size
        row=1, col=2                 # Targets the Main Heatmap subplot grid explicitly
    )
    
    # 3. Add ample bottom margin buffer so the rotated labels don't clip off the screen
    fig.update_layout(
        title_text="Drug Similarity Heatmap",
        margin=dict(b=180, l=80, r=80, t=50) 
    )


    #=======================================
    # PCA plot
    #=======================================
    scaled_data = StandardScaler().fit_transform(matrix_df.values)
    pca = PCA(n_components=2)
    pca_features = pca.fit_transform(scaled_data)
    
    # Calculate the exact variance explained by PC1 and PC2 for plot labeling
    exp_var_cum = pca.explained_variance_ratio_ * 100
    
    # Merge the PCA coordinates back with your qualitative text columns for visualization mapping
    pca_df = pd.DataFrame(
        pca_features, 
        columns=['PC1', 'PC2'], 
        index=matrix_df.index
    )
    pca_df['Class'] = normalized['Class']
    pca_df['Mechanism'] = normalized['Mechanism']
    pca_df['Drug_Name'] = pca_df.index  # Keep index labels for hover boxes
    fig_pca = px.scatter(
        pca_df, 
        x='PC1', 
        y='PC2', 
        color='Class',             # Automatically assigns categorical color tags
        symbol='Mechanism',        # Assigns distinct shapes to different drug mechanisms
        #text='Drug_Name',          # Optional: Displays the row label right next to the node
        hover_data={
            'Drug_Name': True, 
            'Class': True, 
            'Mechanism': True,
            'PC1': ':.2f',         # Format coordinates to 2 decimal places
            'PC2': ':.2f'
        },
        labels={
            'PC1': f'PC1 ({exp_var_cum[0]:.1f}% Variance Explained)',
            'PC2': f'PC2 ({exp_var_cum[1]:.1f}% Variance Explained)'
        },
        title="PCA plot"
    )
    
    # Clean up node visibility and layout scaling
    fig_pca.update_traces(
        textposition='top center', 
        marker=dict(size=8, line=dict(width=1, color='DarkSlateGrey'))
    )
    fig_pca.update_layout(
        height=650, 
        width=1200, 
        template="plotly_white",
        legend=dict(title="Metadata Mappings")
    )

    #==========================================
    # Scatterplot to show pairwise correlation
    #==========================================
    # group by cell line

    cell_lines = [s.split('_')[0] for s in matrix_df.columns.tolist()]
    cell_lines = list(set(cell_lines))

    plot_registry = []
    seen = dict() # to track plot list

    for cl in cell_lines:
        df_subset = matrix_df.filter(regex=f"^{cl}")
        
        # put a limit of values to make plot working properly
        df_subset = df_subset.clip(upper=2.5)

        # get sub-groups based on library concentration
        drug_lib_conc = [s.split('_')[2] for s in df_subset.columns.tolist()]
        drug_lib_conc = list(set(drug_lib_conc))

        for conc in drug_lib_conc:
            data_to_plot = df_subset.filter(regex=f"{conc}")

            for col1 in data_to_plot.columns.tolist():
                for col2 in data_to_plot.columns.tolist():
                    key = col1 + col2
                    if col1!= col2 and not key in seen.keys():
                        key2 = col2 + col1
                        seen[key2] = ''
                        
                        drug_names = data_to_plot.index
                        subfig = px.scatter(data_to_plot, 
                                            x=f"{col1}", 
                                            y=f"{col2}", 
                                            hover_name=drug_names,
                                            #trendline="ols"
                                           )
    
                        # add a x=y diagonal line
                        max_x = max(data_to_plot[col1])
                        max_y = max(data_to_plot[col2])
                        max_z = max([max_x, max_y])
                    
                        subfig.add_scatter(
                            x=[0, 1.2],
                            y=[0, 1.2],
                            mode="lines",  
                            line=dict(
                                color="red",
                                width=2,
                                dash="dash"
                            ),
                        )

                        # Save the trace data AND a descriptive title marker as a tuple block
                        x = col1.split('_')[1]
                        y = col2.split('_')[1]
                        title = ' '.join([col1.split('_')[0], col1.split('_')[2]])
                        max_y = max(data_to_plot[col2])
                        plot_registry.append((subfig.data, f"{title}", f"{x}", f"{y}", max_y))


    # 1. Configuration
    cols = 4  # Set your fixed number of columns
    total_plots = len(plot_registry)
    rows = math.ceil(total_plots / cols)

    # Dynamically build a grid layout: 1 Column, and as many rows as saved plots
    master_scatterplot = make_subplots(
        rows=rows, cols=4,
        vertical_spacing=0.08,
        horizontal_spacing=0.05
    )

    for idx, (traces, panel_title, x_title, y_title, max_y) in enumerate(plot_registry):

        curr_row = (idx // cols) + 1
        curr_col = (idx % cols) + 1
        
        # Map the visual elements
        for trace in traces:
            #trace.line.color = "black"
            #trace.showlegend = True if idx == 0 else False 
            trace.showlegend = False
            master_scatterplot.add_trace(trace, row=curr_row, col=curr_col)
            
        # 3. Dynamic Axis Labels Injection
        master_scatterplot.update_xaxes(title_text=x_title, row=curr_row, col=curr_col)
        master_scatterplot.update_yaxes(title_text=y_title, row=curr_row, col=curr_col)
   
        master_scatterplot.add_annotation(
            text=f"<b>{panel_title}</b>",  # HTML formatting tags work perfectly here
            xref="paper", yref="paper",
            x=0.5,                        # Centers the title horizontally over the subplot panel
            y=max_y + 0.1,                       # Positions the title slightly above the cell frame border
            showarrow=False,
            font=dict(size=14, color="black"),
            row=curr_row, col=curr_col            # Bind this annotation explicitly to the target grid cell
        )

        # 5. Finalize Master Canvas Dimensions
        master_scatterplot.update_layout(
            height=380 * rows, 
            width=350 * cols,
            template="plotly_white"
        )
    
        #master_scatterplot.update_xaxes(range=[0, 1.5])
        
    return fig, fig_pca, master_scatterplot


@callback(
    Output('mean-cell-size-data', 'data'),
    Input("cell-line-dropdown3", "value"),
    Input("condition-dropdown3", "value"),
)
def prepare_tab3_heatmap_data(selected_cell, selected_condition):
   
    if not selected_cell or not selected_condition:
        return dash.no_update
        
    table_name = f"Cell-size_{selected_cell}_{selected_condition}"
    
    with engine.connect() as conn:
        try:
            query = text(f"select * from `{table_name}` where Area < 4000")
            df = pd.read_sql(query, conn)

        except Exception as e:
            print(f"Could not read table. Error: {e}")
            
    # Based on the image with scale bar, 1pixel = 400/320 = 1.25 micron
    # One square pixel is equavalent to 1.25^2 = 1.5625 micron
    df['Area'] = df['Area'] * 1.5625

    # Aggregate Area by Condition and Well
    df_pivot = df.pivot_table(
        values='Area',
        index='Condition',
        columns='Wells',
        aggfunc='mean'
    )

    # Add a difference row
    treatment = selected_condition.split('_')[0]
    df_pivot.loc['Difference'] = df_pivot.loc[treatment] - df_pivot.loc['Ctrl']

    # convert to z-score
    df_zscore = df_pivot.sub(df_pivot.mean(axis=1), axis=0).div(
        df_pivot.std(axis=1), axis=0
    )

    return df_zscore.to_json(date_format='iso', orient='split')


@callback(
    Output('cell-size-heatmap', 'figure'),
    Output("order-by-dropdown", "disabled"),
    Input("mean-cell-size-data", "data"),
    Input("order-by-dropdown", "value"),
)
def update_tab3_heatmap(jsonified_aggregated_data, order_by):
   
    if jsonified_aggregated_data is None:
        return dash.no_update, True
    
    df_zscore = pd.read_json(jsonified_aggregated_data, orient='split')  
            
    with engine.connect() as conn:
        try:
            query = text(f"select * from ARTNET")
            drug_library = pd.read_sql(query, conn)
        except Exception as e:
            print(f"Could not read table. Error: {e}")
            
        # get drug annotation data
        try:
            query = text(f"select * from drug_annotation")
            drug_annotation = pd.read_sql(query, conn)
        except Exception as e:
            print(f"Could not read table. Error: {e}")


    #=======================================
    # Heatmap
    #=======================================
    # data for main heatmap
    # Sort the DataFrame columns based on 'row_label'
    #order_by = selected_condition.split("_")[0]
    
    if not order_by:
        order_by = 'Difference'
        
    df_sorted = df_zscore.sort_values(by=order_by, axis=1, ascending=False)

    # Get the ordered list of column names
    sorted_column_names = df_sorted.columns.tolist()

    #---------------------------------------------------------------------------
    # data for top covariate bar
    #---------------------------------------------------------------------------
    # Merge with drug and drug class data
    drug_df = pd.merge(drug_library, drug_annotation, on = 'Drugs', how="left")

    # Convert the target column to a Categorical type custom order
    drug_df["Wells"] = pd.Categorical(
        drug_df["Wells"], categories=sorted_column_names, ordered=True
    )

    # 3. Sort the DataFrame by that column
    drug_df_sorted = drug_df.sort_values(by="Wells")

    # --------------------------------------------------------------------------
    # MAP COVARIATE CATEGORIES TO COLOR INTEGRATION
    # --------------------------------------------------------------------------
    unique_classes = sorted(list(drug_df_sorted['Class'].unique()))
    class_to_int = {cls: idx for idx, cls in enumerate(unique_classes)}
    int_classes = drug_df_sorted['Class'].map(class_to_int).values.reshape(1, -1)

    color_palette = ['#e6194b', '#3cb44b', '#ffe119', '#4363d8', '#f58231', '#911eb4', '#46f0f0', '#f032e6','#bcf60c', '#fabebe', '#008080', '#e6beff', '#9a6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#808080','#cccccc']

    plotly_colorscale = []
    num_classes = len(unique_classes)
    for i, cls in enumerate(unique_classes):
        plotly_colorscale.append([i / num_classes, color_palette[i % len(color_palette)]])
        plotly_colorscale.append([(i + 1) / num_classes, color_palette[i % len(color_palette)]])

    # --------------------------------------------------------------------------
    # ASSEMBLE SUBPLOT GRID SYSTEM
    # --------------------------------------------------------------------------
    fig = make_subplots(
        rows=2, cols=1,
        shared_xaxes=True,
        shared_yaxes=False,
        vertical_spacing=0.005,
        row_heights=[0.3, 0.7]
    )

    # Add Row Covariate Track (top plot)
    z = df_sorted.values
    x = x =  [str1 + '_' + str2 for str1, str2 in zip(drug_df_sorted.Drugs, drug_df_sorted.Wells)]
    y = df_sorted.index.to_list()

    fig.add_trace(
        go.Heatmap(
            z=int_classes,
            x=x,
            y=["Drug Class"],
            colorscale=plotly_colorscale,
            showscale=False,
            hoverongaps=False,
            hovertemplate="Drug: %{y}<br>Class: %{customdata}<extra></extra>",
            customdata=drug_df_sorted.values.reshape(-1, 1)
        ),
        row=1, col=1
    )
    fig.update_yaxes(autorange="reversed")
    
    # Add Primary Expression Matrix Heatmap (bottom plot)
    fig.add_trace(
        go.Heatmap(
            z=z,
            x=x,
            y=y,
            colorscale=[
                [0.0, "blue"],   # Deepest blue for your lowest negative value
                [0.5, "white"],  # Midpoint color
                [1.0, "red"]     # Deepest red for your highest positive value
            ],
            zmid=0,  
            showscale=True,
            colorbar=dict(title="Zscore of<br>Mean<br>cell size", thickness=15, len=0.8, y=0.55)
        ),
        row=2, col=1
    )

    # --------------------------------------------------------------------------
    # ADD DISCRETE EXPLANATORY LEGENDS
    # --------------------------------------------------------------------------
    for i, cls in enumerate(unique_classes):
        fig.add_trace(
            go.Scatter(
                x=[None], y=[None], mode='markers',
                marker=dict(color=color_palette[i % len(color_palette)], symbol='square', size=14),
                name=cls, legendgroup='Drug Class', showlegend=True
            ), row=2, col=1
        )


    # --------------------------------------------------------------------------
    # E. CONFIGURE BOTTOM LEGEND WITH TEXT WRAPPING COLUMNS
    # --------------------------------------------------------------------------
    fig.update_layout(
        legend=dict(
            title="Drug Class",
            orientation="h",         
            x=0.5,                   
            xanchor="center",       

            # --- ADJUST THIS TO INCREASE THE GAP ---
            # Changing from -0.15 to -0.30 pushes the legend lower down the canvas
            y=-0.7,                 
            yanchor="top",

            entrywidth=380,          
            entrywidthmode="pixels",
        ),

        # This prevents the x-axis text tick labels from overlapping your new gap
        xaxis2=dict(
            tickangle=45,
            title=dict(standoff=30) # Pushes the axis title down if you have one
        ),

        margin=dict(
            b=450,                   
            t=50,                    
            l=50,
            r=50
        ),
        height=1000,                  
        width=1200
    )  
    
    return fig, False

@callback(
    Output('cell-size-violinplot', 'figure'),
    Output('cell-image-control', 'src'),
    Output('cell-image-treated', 'src'),
    Output("image-title", "children"),
    Input("cell-line-dropdown3", "value"),
    Input("condition-dropdown3", "value"),
    Input("drug-name-dropdown3", "value"),
)
def update_tab3_graphs(selected_cell, selected_condition, selected_drug):
   
    if not selected_cell or not selected_condition or not selected_drug:
        return dash.no_update, dash.no_update, dash.no_update, dash.no_update
        
    # Check if data exsits in database
    table_name = f"Cell-size_{selected_cell}_{selected_condition}"
    
    from sqlalchemy import inspect

    inspector = inspect(engine)

    if not inspector.has_table(table_name):
        return  px.scatter(title="Data not availble for this cell line and screen condition yet.")     
            
    # get drug name to plate well mapping
    with engine.connect() as conn:
        try:
            query = text(f"select * from ARTNET where Drugs = '{selected_drug}'")
            drug_to_wells = pd.read_sql(query, conn)
            
            # Some wells are used for a single drug, so only use the first occurence to reduce redundancy
            well = drug_to_wells["Wells"].iloc[0]
            
           
            if well[1] == '0':
                well2 = well[0] + well[2]
            else:
                well2 = well
                
            parts = selected_condition.split("_")
        
            parts[0] = 'DMSO'
            control_condition = "_".join(parts)

            image_file_control = f"images/{selected_cell}/{control_condition}/VID2303_{well2}_1_04d20h00m.png"
            image_file_treated = f"images/{selected_cell}/{selected_condition}/VID2305_{well2}_1_04d20h00m.png"
     
        except Exception as e:
            print(f"Could not read table. Error: {e}")
    
    
    
    with engine.connect() as conn:
        try:
            query = text(f"select * from `{table_name}` where Wells = '{well}'")
            df = pd.read_sql(query, conn)
            
            # Based on the image with scale bar, 1pixel = 400/320 = 1.25 micron
            # One square pixel is equavalent to 1.25^2 = 1.5625 micron
            df['Area2'] = df['Area'] * 1.5625
            
            # fitler to remove exreamly large cells, this is likely artfacts
            df = df[df.Area2 < 6000]
            
            violinplots = px.violin(
                df, 
                y="Area2", 
                x="Condition", 
                color="Condition",         # Splits violins side-by-side by gender
                box=True,            # Draws a box plot inside the violin
                points="all",        # Shows individual data points next to violin
                hover_data=df.columns
            )
    
            # Adjust layout formatting
            violinplots.update_layout(
                # title=f"Cell Area Violinplots",
                yaxis_title="Cell Area, square micron",
                violinmode='group'  # Groups colored violins next to each other
            )
            
        except Exception as e:
            print(f"Could not read table. Error: {e}")
            
        image_title = selected_condition.split("_")[0]
         
    return violinplots, get_asset_url(image_file_control), get_asset_url(image_file_treated), image_title




