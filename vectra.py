#!/usr/bin/env python
# coding: utf-8

# In[1]:


# %load /mnt/Jason/venvs/scimap/labdb/pages/vectra.py
import os
import dash
from dash import Dash, html, dcc, callback, Output, Input, no_update, State
import dash_ag_grid as dag
import pandas as pd
import plotly.express as px
import dash_bootstrap_components as dbc
import dash_bio as dashbio
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import dash_auth
from dash_auth import BasicAuth
import scimap as sm
import anndata as ad
import math
from sqlalchemy import create_engine, text

dash.register_page(__name__)

DATAPATH = '/mnt/Jason/venvs/scimap/data'

# Datasets
pheno_brca81 = pd.read_csv(os.path.join(DATAPATH, 'RPCI_BrCa81_pheno_percent.csv'))
pheno_brca85 = pd.read_csv(os.path.join(DATAPATH, 'RPCI_BrCa85c_22_CC2_pheno_percent.csv'))
pheno_bcm_tnbc = pd.read_csv(os.path.join(DATAPATH, 'BCM_TNBC_CC2_pheno_percent.csv'))

# Format: mysql+pymysql://username:password@host:port/database_name
DATABASE_URL = "mysql+pymysql://vectra:ritho9Ng@10.126.0.41:3306/Vectra"
engine = create_engine(DATABASE_URL)

# Color pallett for consistant cluster colors
distinct20 = {
    '0' : '#e6194b', 
    '1' : '#3cb44b', 
    '2' : '#ffe119', 
    '3' : '#4363d8', 
    '4' : '#f58231', 
    '5' : '#911eb4', 
    '6' : '#46f0f0', 
    '7' : '#f032e6', 
    '8' : '#bcf60c', 
    '9' : '#fabebe', 
    '10': '#008080', 
    '11': '#e6beff', 
    '12': '#9a6324', 
    '13': '#fffac8', 
    '14': '#800000', 
    '15': '#aaffc3', 
    '16': '#808000', 
    '17': '#ffd8b1', 
    '18': '#000075', 
    '19': '#808080',
    'Others': '#cccccc'
}

DATASETS = {
#     'BRCA81': {
#         'Pheno': pheno_brca81,
#         'Consolid': "rpci_brca81_expression"
#     },
    'BRCA85': {
        'Pheno': pheno_brca85,
        'Consolid': "rpci_brca85_cc2_expression",
        'SciMap': "brca85_scimap",
        'Trajectory': "brca85_trajectory",
        'Profile': "brca85_profile"
    },
    'BCM_TNBC': {
        'Pheno': pheno_bcm_tnbc,
        'Consolid': "bcm_tnbc_cc2_expression",
        'SciMap': "bcm_tnbc_scimap",
        'Trajectory': "bcm_tnbc_trajectory",
        'Profile': "bcm_tnbc_profile"
    }
}

# Data for different sources
data_sources = {dataset : DATASETS[dataset]['Pheno'].columns.to_list()[2:]  for dataset in DATASETS.keys() }


# In[11]:




left_pane = [
    dbc.CardBody(
        [
            
            html.P(
                "Select a Data Source",
            ),
            dcc.Dropdown(
                id='data-source-dropdown',
                options=[{'label': source, 'value': source} for source in data_sources.keys()],
                placeholder="Select the data source...",
                clearable=True,
            ),
            html.Br(),
            dcc.Markdown(
                """
                Navigate by clicking the tabs on the right side panel. Currently, only the 
                aggregated pheno type data and marker intensity data at single cell level 
                are availble. This resource is still under active development...
                """
            ),
        ]
    ),
]

card_pheno_correlation = dbc.Card([
    dbc.CardHeader(html.H5("Visualize Correlations of Markers")),
    dbc.CardBody([
        dbc.Row([
            dbc.Col([                  
                dbc.Row([
                    dbc.Label("1. Select X Axis:"),
                    dcc.Dropdown(
                        id='x-dropdown',
                        placeholder="Select a marker...",
                        clearable=True,
                    ),           
                ]),
                dbc.Row([
                    dbc.Label("2. Select Y Axis:"),
                    dcc.Dropdown(
                        id='y-dropdown',
                        placeholder="Select a marker...",
                        clearable=True,
                    ),
                                   
                ]),
            ], width=3, align="center"),    
            dbc.Col(
                dcc.Graph(id='scatter-plot'), 
            width=7),
        ]),
    ]),
    dbc.CardFooter("Some notes here..."),
])
card_heatmaps = dbc.Card([
    dbc.CardHeader(html.H4("Cluster Markers with Heatmap")),
    dbc.CardBody([
        dbc.Row([                 
            dbc.Col([
                dbc.Label("Select/Deselect Markers:"),
                
                # This container starts empty and is filled by the callback            
                dcc.Dropdown(
                    id="markers-dropdown-heatmap",
                    options=[],
                    value=[],
                    placeholder="Select markers...",
                    clearable=True,
                    multi=True,
                                        
                ),
            ], width=10),
        ]),
        dbc.Row([
            dbc.Col([
                #dbc.Label("Tumor"),
                html.H6(
                        html.B("Tumor"),
                        className="text-center mt-2 mb-1",
                        style={"color": "Purple", "text-decoration": "None",},
                    ),
                dcc.Loading(
                dcc.Graph(id="heatmap1"),
                ),
            ], width=5),
            dbc.Col([
                #dbc.Label("Stroma"),
                html.H6(
                        html.B("Stroma"),
                        className="text-center mt-2 mb-1",
                        style={"color": "Blue", "text-decoration": "None",},
                    ),
                dcc.Loading(
                dcc.Graph(id="heatmap2"),
                ),
                              
            ], width=5),
        ], align="center"),          
    ]),
    dbc.CardFooter("Some notes here..."),
])
card_pheno_table = dbc.Card([
    dbc.CardHeader(html.H4('Table View')),
    dbc.CardBody([
        html.H5("Table: Percent of cells positive for each marker"),
        dag.AgGrid(
            id="dynamic-grid",
            dashGridOptions={"pagination": True},
            defaultColDef={"sortable": True, "filter": True, "resizable": True},                 
        ),
    ]),
    dbc.CardFooter([dbc.Button("Download as CSV", id="btn-download"),
    
        # The Download component (invisible)
        dcc.Download(id="download-dataframe-csv")]
    )
], className="shadow")
card_consolid_table = dbc.Card([
    dbc.CardHeader(html.H4("Single Cell Level Measurement")),
    dbc.CardBody([
        html.H5("Table: consolidated data"),
        dag.AgGrid(
            id="dynamic-grid2",
            dashGridOptions={"pagination": True},
            defaultColDef={"sortable": True, "filter": True, "resizable": True},   
        ),
    ]),
    dbc.CardFooter([dbc.Button("Download as CSV", id="btn-download2"),
    
        # The Download component (invisible)
        dcc.Download(id="download-dataframe-csv2")])
], className="shadow")

card_scimap = dbc.Card([
    dbc.CardHeader(html.H4("Sci-Map Analysis")),
    dbc.CardBody([
        dbc.Row([
            dbc.Col([
                dcc.Loading(
                    dcc.Graph(id='leiden-heatmap'),
                ),
            ], width=6),
            dbc.Col([
                dcc.Loading(
                    dcc.Graph(id='leiden-umap')
                ),
            ], width=6),
        ]),
        dbc.Row([
            dcc.Loading(
            dcc.Graph(id='markers-umap'),
            ),
        ]),
    ]),
    dbc.CardFooter("In UMAPs, data is down-sampled to include 5,000 randomly selected cells to speed up the display")
])

card_trajectory_scatter = dbc.Card([
    dbc.CardHeader(html.H4("Trajectory Analysis: LMDS plot")),
    dbc.CardBody([
        dbc.Row([
            dbc.Col([
                dbc.Label("Select a Sample/TMA:"),
                
            # This container starts empty and is filled by the callback            
            dcc.Dropdown(
                id="tma-dropdown",
                options=[],
                value=[],
                placeholder="Select a Sample/TMA...",
                
                clearable=True,
                multi=True,                          
            ),], width=2)
            
        ]),
        dbc.Row([
            dcc.Loading(
            dcc.Graph(id='trajectory-scatter'),
            ),
        ]),
    ]),
    
])

card_trajectory_profile = dbc.Card([
    dbc.CardHeader(html.H4("Trajectory Analysis: pseudo-time profile plot")),
    dbc.CardBody([
        dbc.Row([
            dbc.Col([
                dbc.Label("Select a Sample/TMA:"),
                
            # This container starts empty and is filled by the callback            
            dcc.Dropdown(
                id="tma-dropdown-profile",
                options=[],
                value=[],
                placeholder="Select a Sample/TMA...",
                clearable=True,
                multi=False,                          
            ),], width=2),
            
            dbc.Col([
                dbc.Label("Show/Hide Markers:"),
                
            # This container starts empty and is filled by the callback            
            dcc.Dropdown(
                id="show-hide-markers",
                options=[],
                value=[],
                placeholder="Select Markers to Show/Hide.",
                clearable=True,
                multi=True,                          
            ),], width=2)     
        ]),
        dbc.Row([
            dcc.Loading(
            dcc.Graph(id='trajectory-profile'),
            ),
        ]),
    ]),
    
])
footer_markdown=dcc.Markdown("""
    This is the footer text...
""")

layout = dbc.Container([
    dbc.Row(html.H1("Explore Vectra Polaris Cancer Sample Data Sets"), style = {'textAlign' : 'center'}),
     dbc.Row([
        dbc.Col(dbc.Row(dbc.Card(left_pane, color="light")), width=2),
        dbc.Col([
            dbc.Tabs([
                dbc.Tab(label='Aggregated Pheno Data', tab_id="pheno_tab", children=[
                    card_pheno_correlation,
                    card_heatmaps, 
                    card_pheno_table,                   
                ]),
                dbc.Tab(label='Sci-Map analysis', tab_id="scimap_tab", children=[
                    card_consolid_table,
                    card_scimap,
                ]),
                dbc.Tab(label='Trajectory Analysis', tab_id="trajectory_tab", children=[
                    card_trajectory_scatter,
                    card_trajectory_profile,
                ]),
            ],
            id="card-tabs",
            active_tab="pheno_tab", 
            ),
            
        ], width=10),
    ]),
    dbc.Row(footer_markdown),
], fluid=True)                          


# Callback 1: Populate X and Y dropdowns based on selected data source
@callback(
    [Output('x-dropdown', 'options'),
     Output('x-dropdown', 'value'),
     Output('y-dropdown', 'options'),
     Output('y-dropdown', 'value')],
    Input('data-source-dropdown', 'value')
)
def update_dropdowns(data_source):
    if not data_source:
        return [], None, [], None
    
    df = DATASETS[data_source]['Pheno']
    
    numeric_cols = df.select_dtypes(include=['number']).columns.tolist()
    
    options = [{'label': col, 'value': col} for col in numeric_cols]
    
    # Default to first two numeric columns if available
    default_x = numeric_cols[0] if numeric_cols else None
    default_y = numeric_cols[1] if len(numeric_cols) > 1 else default_x
    
    return options, default_x, options, default_y


# Callback 2: Update scatter plot based on selections
@callback(
    Output('scatter-plot', 'figure'),
    [Input('data-source-dropdown', 'value'),
     Input('x-dropdown', 'value'),
     Input('y-dropdown', 'value')]
)
def update_scatter_plot(data_source, x_col, y_col):
    if not data_source or not x_col or not y_col:
        # Return empty figure if inputs are missing
        return px.scatter(title="")
    
    df = DATASETS[data_source]['Pheno']
    
    # Create scatter plot
    fig = px.scatter(
        df,
        x=x_col,
        y=y_col,
        #title=f"{data_source.capitalize()} Dataset: {y_col} vs {x_col}",
        hover_data=df.columns.tolist() if len(df.columns) < 10 else None,
        color=df.columns[1] if df.columns[0] not in [x_col, y_col] else None
    )
    
    fig.update_layout(
        transition_duration=300,
        xaxis_title=x_col,
        yaxis_title=y_col
    )

    fig.update_layout(legend=dict(
        yanchor="top",
        y=0.99,
        xanchor="left",
        x=0.01
        ),
        template='plotly_white'
    )

    return fig

@callback(
    Output("dynamic-grid", "rowData"),
    Output("dynamic-grid", "columnDefs"),
    Input("data-source-dropdown", "value"),
    prevent_initial_call=True,
)
def update_grid(selected_dataset):
    if not selected_dataset:
        return [], []
    else: 
        df = DATASETS[selected_dataset]['Pheno']

    # Generate Column Definitions dynamically
    # This maps each column name to the format AG Grid expects
    column_defs = [{"field": i} for i in df.columns]

    # Convert DataFrame to records (List of Dicts)
    row_data = df.to_dict('records')
    
    # Return the data using the dcc.send_data_frame helper
    return row_data, column_defs
   
@callback(
    Output("download-dataframe-csv", "data"),
    State("data-source-dropdown", "value"),
    Input("btn-download", "n_clicks"),
    prevent_initial_call=True,
)
def download_data(selected_dataset, n_clicks):
    if not selected_dataset:
        return [], []
    else: 
        df = DATASETS[selected_dataset]['Pheno']

    # Return the data using the dcc.send_data_frame helper
    return dcc.send_data_frame(df.to_csv, f"data_{selected_dataset}.csv", index=False)
   
    
@callback(
    Output("dynamic-grid2", "rowData"),
    Output("dynamic-grid2", "columnDefs"),
    Input("data-source-dropdown", "value"),
    prevent_initial_call=True,
)
def update_grid2(selected_dataset):
    if not selected_dataset:
        return [], []    
        
    table_name = DATASETS[selected_dataset]['Consolid']
    
    try:
        with engine.connect() as conn:
            
            query = text(f"SELECT * FROM `{table_name}` LIMIT 1000")
            df = pd.read_sql(query, conn)
        
            # Generate Column Definitions dynamically
            # This maps each column name to the format AG Grid expects
            column_defs = [{"field": i} for i in df.columns]

            # Convert DataFrame to records (List of Dicts)
            row_data = df.to_dict('records')
    
            # Return the data using the dcc.send_data_frame helper
            return row_data, column_defs
         
    except Exception as e:
        print(f"Error reading table '{table_name}': {e}")
        
        return [], []
    
    
   
@callback(
    Output("download-dataframe-csv2", "data"),
    State("data-source-dropdown", "value"),
    Input("btn-download2", "n_clicks"),
    prevent_initial_call=True,
)
def download_data2(selected_dataset, n_clicks2):
    if not selected_dataset:
        return [], []
    
    table_name = DATASETS[selected_dataset]['Consolid']
    
    try:
        with engine.connect() as conn:
            query = text(f"SELECT * FROM `{table_name}` LIMIT 1000")
            df = pd.read_sql(query, conn)
    
            # Return the data using the dcc.send_data_frame helper
            return dcc.send_data_frame(df.to_csv, f"Consolidated_data_{selected_dataset}.csv", index=False)
        
    except Exception as e:
        print(f"Error reading table '{table_name}': {e}")
        
        return [], []
    
    
   
    
@callback(
    Output('markers-dropdown-heatmap', 'options'),
    Output('markers-dropdown-heatmap', 'value'),
    Input('data-source-dropdown', 'value')
)
def create_dropdown(data_source):
    if not data_source:
        return [], []
    
    # Define options based on selection   
    df = DATASETS[data_source]['Pheno']
    
    numeric_cols = df.select_dtypes(include=['number']).columns.tolist()
    
    options = [{'label': col, 'value': col} for col in numeric_cols]
    
    #default_values = [options[0]['value'], options[1]['value']] if len(options) >= 2 else []
    default_values = [i['value'] for i in options] if len(options) >= 2 else []
    
    return options, default_values

@callback(
    Output('heatmap1', 'figure'),
    Output('heatmap2', 'figure'),
    Input('data-source-dropdown', 'value'),
    Input('markers-dropdown-heatmap', 'value')
)
def update_heatmap(data_source, selected_markers):
    if not data_source:
        # Return empty figure if inputs are missing
        #return go.Figure() # emptyb figure
        return  px.scatter(title=""), px.scatter(title="")
    else:
        df = DATASETS[data_source]['Pheno']
    
    if not selected_markers:
        return no_update
    
    if len(selected_markers) < 2:
        return "Please select at least two markers to display."
     
    # split df by tissue category
    tumor = df[df["Tissue Category"]=='Tumor']
    stroma = df[df["Tissue Category"]=='Stroma']
    
    tumor = tumor.select_dtypes(include=['number'])
    tumor = tumor.loc[:, selected_markers]
     
    clustergram1 = dashbio.Clustergram(
        data=tumor.values,
        column_labels=list(tumor.columns),
        row_labels=list(tumor.index),
        height=700,
        width=600,
        color_threshold={
            'row': 2,
            'col': 2
        },
        hidden_labels=['row'], # Hides individual row names if the list is too long
        display_ratio=0.05,     # Proportion of the plot used for dendrograms
        color_map=[
            [0.0, 'blue'],
            [0.5, 'white'],
            [1.0, 'red']
        ]
    )
    stroma = stroma.select_dtypes(include=['number'])
    stroma = stroma.loc[:, selected_markers]
     
    clustergram2 = dashbio.Clustergram(
        data=stroma.values,
        column_labels=list(stroma.columns),
        row_labels=list(stroma.index),
        height=700,
        width=600,
        color_threshold={
            'row': 2,
            'col': 2
        },
        hidden_labels=['row'], # Hides individual row names if the list is too long
        display_ratio=0.05,     # Proportion of the plot used for dendrograms
        color_map=[
            [0.0, 'blue'],
            [0.5, 'white'],
            [1.0, 'red']
        ]
    )
    
    return clustergram1, clustergram2

@callback(
    Output('leiden-heatmap', 'figure'),
    Output('leiden-umap', 'figure'),
    Output('markers-umap', 'figure'),
    Input('data-source-dropdown', 'value'),
    prevent_initial_call=True,
)
def update_scimap(data_source):
    if not data_source:      
        return  px.scatter(title="Please select data source"), px.scatter(title=""), px.scatter(title="")
    
    if not 'SciMap' in DATASETS[data_source]:
        return  px.scatter(title="No Data"), px.scatter(title="No Data"), px.scatter(title="No Data")
    
    table_name = DATASETS[data_source]['SciMap']
    
    try:
        with engine.connect() as conn:
            query = text(f"SELECT * FROM `{table_name}` ORDER BY RAND() limit 10000")
            df = pd.read_sql(query, conn)    
    except Exception as e:
        print(f"Error reading table '{table_name}': {e}")
        return  px.scatter(title=f"database error: {table_name}"), px.scatter(title=""), px.scatter(title="")
    
    df2 = df.iloc[:, :df.columns.get_loc('X_centroid')]
    markers = df2.columns.to_list()
    markers.append('leiden')
    
    # add leiden column
    df2 = df[markers]
    df_heatmap = df2.groupby(by = 'leiden').mean()
    
    # get the markers name for later use
    markers.remove('leiden')
    del df2
    
    # Leiden cluster heatmap
    leiden_heatmap = dashbio.Clustergram(
        data=df_heatmap.values,
        column_labels=list(df_heatmap.columns),
        row_labels=list(df_heatmap.index),
        height=350,
        width=600,
        color_threshold={
            'row': 1,
            'col': 1
        },
        display_ratio=0.05,     # Proportion of the plot used for dendrograms
        color_map=[
            [0.0, 'blue'],
            [0.5, 'white'],
            [1.0, 'red']
        ]
    )
    leiden_heatmap.update_layout(
        title='Leiden Cluster Heatmap',
        title_x=0.5  # Centers the title
    )
  
    # Leiden Cluster UMAP colored by clusters
    
    df_umap = df[['UMAP_1', 'UMAP_2', 'leiden']]
    
    leiden_umap = px.scatter(
        df_umap,
        x='UMAP_1',
        y='UMAP_2',
        title="Leiden Clusters UMAP",
        # hover_data=df_umap_subset.columns.tolist() if len(df_umap_subset.columns) < 10 else None,
        color=df_umap.leiden,
        color_discrete_map=distinct20,
        width=410,
        height=380
    )
    leiden_umap.update_traces(marker_size=2)
    leiden_umap.update_layout(
        legend_title_text='Cluster', 
        title_x=0.5, 
        template='plotly_white',
        legend=dict(
            font=dict(size=16),
            itemsizing='constant'  # ensures legend markers don’t shrink
        )
    )
    
    # 1. Configuration
    cols = 4  # Set your fixed number of columns
    total_plots = len(markers)
    rows = math.ceil(total_plots / cols)

    # 2. Initialize subplots
    marker_umap = make_subplots(
        rows=rows, 
        cols=cols, 
        subplot_titles=list(markers),
        vertical_spacing=0.1,  # Add space between rows
        horizontal_spacing=0.05
    )
      
    # 4. Loop with Row/Col Calculation
    
    for i, marker in enumerate(markers):
        # Calculate current row and column (1-based)
        curr_row = (i // cols) + 1
        curr_col = (i % cols) + 1

        marker_umap.add_trace(
            go.Scatter(
                x=df_umap['UMAP_1'], 
                y=df_umap['UMAP_2'],
                name=marker,
                mode='markers',
                
                marker=dict(
                    size=3,
                    color=df[marker].values,
                    colorscale=[[0, 'blue'], [0.5, 'white'], [1, 'red']],  # blue-white-red
                    cmin=df[marker].min(),
                    cmax=df[marker].max(),
                    showscale=True if i == total_plots - 1 else False
                )

#                 marker=dict(
#                     size=3,
#                     color=df[marker].values,
#                     colorscale='Viridis',
#                     # Show scale only for the very last plot
#                     showscale=True if i == total_plots - 1 else False 
#                 )
            ),
            row=curr_row, col=curr_col
        )

    # 5. Adjust layout height dynamically based on rows
    marker_umap.update_layout(
        height=300 * rows, 
        width=1200, 
        showlegend=False,
        title_text="Marker Intensity Distribution",
        title_x=0.5,
        template='plotly_white',
        
    )

    return leiden_heatmap, leiden_umap, marker_umap

# Trajectory
@callback(
    Output('tma-dropdown', 'options'),
    Input('data-source-dropdown', 'value')
)
def update_tma_dropdowns(data_source):
    if not data_source:
        return []
    
    if not 'Trajectory' in DATASETS[data_source]:
        return  px.scatter(title="No data")
    
    table_name_lmds = DATASETS[data_source]['Trajectory']
    
    try:
        with engine.connect() as conn:
            query = text(f"SELECT distinct(SampleID) FROM `{table_name_lmds}`")
            df_lmds = pd.read_sql(query, conn)    
    except Exception as e:
        print(f"Error reading table '{table_name_lmds}': {e}")
        return  px.scatter(title=f"database error: {table_name_lmds}")
    
    options = [{'label': s, 'value': s} for s in df_lmds.SampleID.to_list()]
    
    return options

# Trajectory
@callback(
    Output('tma-dropdown-profile', 'options'),
    Input('data-source-dropdown', 'value')
)
def update_tma_dropdowns2(data_source):
    if not data_source:
        return []
    
    if not 'Profile' in DATASETS[data_source]:
        return  px.scatter(title="No data")
    
    table_name = DATASETS[data_source]['Profile']
    
    try:
        with engine.connect() as conn:
            query = text(f"SELECT distinct(source) FROM `{table_name}`")
            df = pd.read_sql(query, conn)    
    except Exception as e:
        print(f"Error reading table '{table_name}': {e}")
        return  px.scatter(title=f"database error: {table_name}")
    
    options = [{'label': s, 'value': s} for s in df.source.to_list()]
    
    return options
@callback(
    Output('trajectory-scatter', 'figure'),
    Input('data-source-dropdown', 'value'),
    Input('tma-dropdown', 'value'),
    prevent_initial_call=True,
)
def update_trajectory_scatter(data_source, tma_selected):
    if not data_source:      
        return  px.scatter(title="Please select data source")
    
    if not 'Trajectory' in DATASETS[data_source]:
        return  px.scatter(title="Analysis not present")
    
    if not tma_selected:
        return  px.scatter(title="No TMA selected")
    
    table_name_scatter = DATASETS[data_source]['Trajectory']
    
    try:
        with engine.connect() as conn:
            # Use a colon (:variable) to define a parameter placeholder
            query = text(f"SELECT * FROM `{table_name_scatter}` WHERE SampleID IN :sample_id ORDER BY RAND() LIMIT 5000")
        
            # Pass the parameter value safely inside pd.read_sql using params
            df_lmds = pd.read_sql(query, conn, params={"sample_id": tma_selected})   
            
    except Exception as e:
        print(f"Error reading table '{table_name_scatter}': {e}")
        return  px.scatter(title=f"database error: {table_name_scatter}")
    
    # Get column names for markers for coloring the scatterplot
    try:
        with engine.connect() as conn:
            query = text(f"desc `{table_name_scatter}`")
            df_col = pd.read_sql(query, conn)    
    except Exception as e:
        print(f"Error reading table '{table_name_scatter}': {e}")
        
    markers = df_col.Field.to_list()
    markers = markers[6:len(markers)]
    
    # 1. Configuration
    cols = 4 # Set your fixed number of columns
    total_plots = len(markers)
    rows = math.ceil(total_plots/cols)

    # 2. Initialize subplots
    marker_scatter = make_subplots(
        rows=rows, 
        cols=cols, 
        subplot_titles=markers,
        vertical_spacing=0.1,  # Add space between rows
        horizontal_spacing=0.05
    )

    # 4. Loop with Row/Col Calculation
    for i, col_name in enumerate(markers):
        # Calculate current column (1-based)
        curr_row = (i // cols) + 1
        curr_col = (i % cols) + 1

        marker_scatter.add_trace(
            go.Scatter(
                x=df_lmds['comp_1'], 
                y=df_lmds['comp_2'],
                name=col_name,
                mode='markers',               
                marker=dict(
                        size=3,
                        color=df_lmds[col_name].values,
                        colorscale='Turbo',
                        colorbar=dict(
                            title="Scaled<br>Intensity"
                        ),
                        showscale=True if i == total_plots - 1 else False 
                    )

            ),
            row=curr_row, col=curr_col
        )

    # 5. Adjust layout height dynamically based on rows
    marker_scatter.update_layout(
        height=320*rows, 
        width=310*cols, 
        showlegend=False,
        title_text="Trajectory Plots",
        title_x=0.5,
        template='plotly_white'
    )
    
    return marker_scatter

@callback(
    Output('show-hide-markers', 'options'),
    Input('data-source-dropdown', 'value'),
    Input('tma-dropdown-profile', 'value'),
)
def update_markers_dropdowns(data_source, tma_selected):
    if not data_source or not tma_selected:
        return []
    
    if not 'Profile' in DATASETS[data_source]:
        return  px.scatter(title="No data")
    
    table_name = DATASETS[data_source]['Profile']
    
    try:
        with engine.connect() as conn:
            query = text(f"SELECT distinct(Marker) FROM `{table_name}`")
            df = pd.read_sql(query, conn)    
    except Exception as e:
        print(f"Error reading table '{table_name}': {e}")
        return  px.scatter(title=f"database error: {table_name}")
    
    options = [{'label': m, 'value': m} for m in df.Marker.to_list()]
    
    return options

@callback(
    Output('trajectory-profile', 'figure'),
    Input('data-source-dropdown', 'value'),
    Input('tma-dropdown-profile', 'value'),
    Input('show-hide-markers', 'value'),
    prevent_initial_call=True,
)
def update_trajectory_profile(data_source, tma_selected, markers_selected):
    if not data_source:      
        return  px.scatter(title="Please select data source")
    
    if not 'Profile' in DATASETS[data_source]:
        return  px.scatter(title="Analysis not present")
    
    if not tma_selected:
        return  px.scatter(title="No TMA selected")
    
    table_name_profile = DATASETS[data_source]['Profile']
    
    try:
        with engine.connect() as conn:
            # Use a colon (:variable) to define a parameter placeholder
            query = text(f"SELECT * FROM `{table_name_profile}` WHERE source = :sample_id")
        
            # Pass the parameter value safely inside pd.read_sql using params
            df_profile = pd.read_sql(query, conn, params={"sample_id": tma_selected})
            
            if markers_selected:
                df_profile = df_profile.loc[~df_profile['Marker'].isin(markers_selected)]
            
    except Exception as e:
        print(f"Error reading table '{table_name_profile}': {e}")
        return  px.scatter(title=f"database error: {table_name_profile}")
    
    fig = px.line(df_profile, x='Pseudotime', y='Expression', color='Marker', template="plotly_white")
    
    return fig
    





