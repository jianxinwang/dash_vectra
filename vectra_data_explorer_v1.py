#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
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


# In[2]:


VALID_USERNAME_PASSWORD_PAIRS = {
    'dbuser': 'xxxxxxxx'
}


# In[3]:


os.path.join(os.path.dirname(os.path.abspath('data')), 'data')


# In[4]:


DATAPATH = '/mnt/Jason/venvs/scimap/data'

# Datasets
pheno_brca81 = pd.read_csv(os.path.join(DATAPATH, 'RPCI_BrCa81_pheno_percent.csv'))
pheno_brca85 = pd.read_csv(os.path.join(DATAPATH, 'RPCI_BrCa85c_22_CC2_pheno_percent.csv'))
pheno_bcm_tnbc = pd.read_csv(os.path.join(DATAPATH, 'BCM_TNBC_CC2_pheno_percent.csv'))

consolidated_brca81 = pd.read_csv(os.path.join(DATAPATH, 'RPCI_BrCa81_expression.csv'))
consolidated_brca85 = pd.read_csv(os.path.join(DATAPATH, 'RPCI_BrCa85c_22_CC2_expression.csv'))
consolidated_bcm_tnbc = pd.read_csv(os.path.join(DATAPATH, 'BCM_TNBC_CC2_expression.csv'))


# In[5]:


adata_brca85 = ad.read_h5ad(os.path.join(DATAPATH,'BrCa85c_res_02_data.h5ad'))
adata_bcm_tnbc = ad.read_h5ad(os.path.join(DATAPATH,'BCM_res_02_data.h5ad'))                         


# In[6]:


DATASETS = {
    'BRCA81': {
        'Pheno': pheno_brca81,
        'Consolid': consolidated_brca81[1:1000]
    },
    'BRCA85': {
        'Pheno': pheno_brca85,
        'Consolid': consolidated_brca85.iloc[1:1000],
        'Adata': adata_brca85
    },
    'BCM_TNBC': {
        'Pheno': pheno_bcm_tnbc,
        'Consolid': consolidated_bcm_tnbc.iloc[1:1000],
        'Adata': adata_bcm_tnbc
    }
}


# In[7]:


# Data for different sources
data_sources = {dataset : DATASETS[dataset]['Pheno'].columns.to_list()[2:]  for dataset in DATASETS.keys() }

colors = {
    'background': '#111111',
    'text': '#7FDBFF'
}

app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])


# In[8]:


# auth = dash_auth.BasicAuth(
#     app,
#     VALID_USERNAME_PASSWORD_PAIRS,
#     #secret_key="WLJWJCPVrb3n1VJ6FuKCmTRU-BH3CVlBzRP91CdrWOg"
# )


# In[9]:


data_source_dropdown = dcc.Dropdown(
    id='data-source-dropdown',
    options=[{'label': source, 'value': source} for source in data_sources.keys()],
    placeholder="Select the data source...",
    clearable=True,
)

about_markdown = dcc.Markdown(
"""
Description of this webpage: This resource is about...
"""
)
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
                dcc.Graph(id="heatmap1"),
            ], width=5),
            dbc.Col([
                #dbc.Label("Stroma"),
                html.H6(
                        html.B("Stroma"),
                        className="text-center mt-2 mb-1",
                        style={"color": "Blue", "text-decoration": "None",},
                    ),
                dcc.Graph(id="heatmap2"),
                              
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
                
                dcc.Graph(id='leiden-heatmap'),], width=6),
            dbc.Col([
                
                dcc.Graph(id='leiden-umap')], width=6),
        ]),
        dbc.Row([dcc.Graph(id='markers-umap'),]),
        
    ]),
    
])

footer_markdown=dcc.Markdown("""
    This is the footer text...
""")


# In[10]:


app.layout = dbc.Container([
    dbc.Row(html.H1("Explore Vectra Polaris Cancer Sample Data Sets"), style = {'textAlign' : 'center'}),
     dbc.Row([
        dbc.Col([data_source_dropdown,
                about_markdown,
                ], width=2),
        dbc.Col([
            dbc.Tabs([
                dbc.Tab(label='Aggregated Pheno Data', tab_id="pheno_tab", children=[
                    card_pheno_correlation,
                    card_heatmaps, 
                    card_pheno_table,                   
                ]),
                dbc.Tab(label='Single Cell Level Data', tab_id="single_tab", children=[
                    card_consolid_table,
                    card_scimap,
                ]),
                
            ],
            id="card-tabs",
            active_tab="pheno_tab", 
            ),
            
        ], width=10),
    ]),
    dbc.Row(footer_markdown),
], fluid=True)                          


# In[11]:


# Callback 1: Populate X and Y dropdowns based on selected data source
@app.callback(
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
@app.callback(
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
    ))

    return fig

@app.callback(
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
   
@app.callback(
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
   
    
@app.callback(
    Output("dynamic-grid2", "rowData"),
    Output("dynamic-grid2", "columnDefs"),
    Input("data-source-dropdown", "value"),
    prevent_initial_call=True,
)
def update_grid2(selected_dataset):
    if not selected_dataset:
        return [], []
    else: 
        df = DATASETS[selected_dataset]['Consolid']

    # Generate Column Definitions dynamically
    # This maps each column name to the format AG Grid expects
    column_defs = [{"field": i} for i in df.columns]

    # Convert DataFrame to records (List of Dicts)
    row_data = df.to_dict('records')
    
    # Return the data using the dcc.send_data_frame helper
    return row_data, column_defs
   
@app.callback(
    Output("download-dataframe-csv2", "data"),
    State("data-source-dropdown", "value"),
    Input("btn-download2", "n_clicks"),
    prevent_initial_call=True,
)
def download_data2(selected_dataset, n_clicks2):
    if not selected_dataset:
        return [], []
    else: 
        df = DATASETS[selected_dataset]['Consolid']

    # Return the data using the dcc.send_data_frame helper
    return dcc.send_data_frame(df.to_csv, f"Consolidated_data_{selected_dataset}.csv", index=False)
   
    
@app.callback(
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

@app.callback(
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

@app.callback(
    Output('leiden-heatmap', 'figure'),
    Output('leiden-umap', 'figure'),
    Output('markers-umap', 'figure'),
    Input('data-source-dropdown', 'value') 
)
def update_scimap(data_source):
    if not data_source:      
        return  px.scatter(title=""), px.scatter(title=""), px.scatter(title="")
    
    if not 'Adata' in DATASETS[data_source]:
        return  px.scatter(title="No Data"), px.scatter(title="No Data"), px.scatter(title="No Data")
        
    adata = DATASETS[data_source]['Adata']
    df = adata.to_df()
    df['Leiden'] = adata.obs.leiden
    df2 = df.groupby(by = 'Leiden').mean()
    
    # Leiden cluster heatmap
    leiden_heatmap = dashbio.Clustergram(
        data=df2.values,
        column_labels=list(df2.columns),
        row_labels=list(df2.index),
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
    
    # Leiden cluster UMAP
    df2 = pd.DataFrame(adata.obsm['X_umap'], columns=['UMAP_1', 'UMAP_2'])
    df2_subset = df2.sample(n=5000)
    Cluster=adata.obs.leiden.iloc[df2_subset.index]
    
    
    leiden_umap = px.scatter(
        df2_subset,
        x='UMAP_1',
        y='UMAP_2',
        title="Leiden Clusters UMAP",
        hover_data=df2_subset.columns.tolist() if len(df2_subset.columns) < 10 else None,
        color=Cluster
    )
    leiden_umap.update_traces(marker_size=5)
    leiden_umap.update_layout(legend_title_text='Cluster', title_x=0.5)
    
    # 1. Ensure subplots are defined correctly
    df = df.drop('Leiden', axis=1)
    
    # 1. Configuration
    cols = 4  # Set your fixed number of columns
    total_plots = len(df.columns)
    rows = math.ceil(total_plots / cols)

    # 2. Initialize subplots
    marker_umap = make_subplots(
        rows=rows, 
        cols=cols, 
        subplot_titles=list(df.columns),
        vertical_spacing=0.1,  # Add space between rows
        horizontal_spacing=0.05
    )
    
    
    # 2. Ensure data alignment
    df_subset = df.iloc[df2_subset.index] 

    
    # 3. Loop with Row/Col Calculation
    for i, col_name in enumerate(df.columns):
        # Calculate current row and column (1-based)
        curr_row = (i // cols) + 1
        curr_col = (i % cols) + 1

        marker_umap.add_trace(
            go.Scatter(
                x=df2_subset['UMAP_1'], 
                y=df2_subset['UMAP_2'],
                name=col_name,
                mode='markers',
                marker=dict(
                    size=3,
                    color=df_subset[col_name].values,
                    colorscale='Viridis',
                    # Show scale only for the very last plot
                    showscale=True if i == total_plots - 1 else False 
                )
            ),
            row=curr_row, col=curr_col
        )

    # 4. Adjust layout height dynamically based on rows
    marker_umap.update_layout(
        height=300 * rows, 
        width=1200, 
        showlegend=False,
        title_text="Marker Intensity Distribution",
        title_x=0.5
    )

    
    #marker_umap.update_layout(height=300, width=1800, title_text="Marker Intensity Distribution")
    
    return leiden_heatmap, leiden_umap, marker_umap

if __name__ == '__main__':
    app.run(debug=True, host="10.126.0.41", port=8787)


# In[ ]:





# In[ ]:




