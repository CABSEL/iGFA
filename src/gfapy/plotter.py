import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


class res_plotter_generic():
    def __init__(self, results, multi_result=False, met_meta=None, rxn_meta=None, name=None, time_col_name='Time (WD)'):
        self.name = name
        self.time_col_name = time_col_name
        if multi_result is False:
            self.results = {'Main': results.copy()}
            self.curr_ind = 'Main'
        else:
            self.results = results.copy()
            self.curr_ind = list(results.keys())[0]

        if rxn_meta is not None:
            self.rxn_meta = rxn_meta

        if met_meta is not None:
            self.met_meta = met_meta

    def summarize_runs(self):
        color_cols = sorted(set(self.results[self.curr_ind]['objectives'].keys()) - {'message'})
        self.model_summary = (pd.DataFrame({k: v['objectives'] for k, v in self.results.items()}).
                              T.reset_index(names='Runs').
                              sort_values(by=['obj'] + sorted(set(color_cols) - {'obj'})).
                              head(100))
        multistart_report = (self.model_summary.
                             style.set_caption(f'<b>Objective in {len(self.results)} multistart runs sorted by objective (showing top {np.min([100, len(self.results)])})</b>').
                             hide(axis="index").
                             set_table_styles([{'selector': 'th:not(.index_name)',
                                                'props': 'background-color: #000000; color: white;border-left: 1px solid white;border-top: 1px solid white;'},
                                               {'selector': 'td', 'props': 'border: 1px solid black;'}]).
                             highlight_quantile(q_right=0.05, axis='index', color='yellow', subset=color_cols,
                                                props='border-left: 2px dashed red;background-color:yellow;').
                             highlight_min(axis='index', props='border: 2px dashed red; background-color: green;', subset=color_cols).
                             highlight_null())

        return multistart_report

    def sel_result(self, res_index='Main'):
        if res_index in list(self.results.keys()):
            self.curr_ind = res_index
            print(self.results[self.curr_ind]['status']['Solver'])
        else:
            raise ValueError('Key not found in results!')

        self.curr_result = self.results[self.curr_ind]

    def choose_top_n(self, obj_col, n=10, ascending=True):
        return self.model_summary.sort_values(by=obj_col, ascending=ascending).head(n)['Runs'].tolist()

    def get_params_n(self, n_res, param_col):
        return pd.concat({e: self.results[e][param_col] for i, e in enumerate(n_res)})

    def sel_compartment(self, comp=1):
        self.curr_result['beta'] = self.curr_result['beta_allcomps'][comp]
        self.curr_result['gamma'] = self.curr_result['gamma_allcomps'][comp]
        self.curr_result['v_ref'] = self.curr_result['v_ref_allcomps'][comp]
        self.curr_result['secreted_flux'] = self.curr_result['secreted_flux_allcomps'][comp]
        self.curr_result['internal_flux'] = self.curr_result['internal_flux_allcomps'][comp]
        self.curr_result['secreted_flux_ratio'] = self.curr_result['secreted_flux_ratio_allcomps'][comp]

    def plot_alphas(self, time_col, meas_cols=None, ncols=5, figsize=None, res_ids=None, **plot_kwargs):
        if res_ids is None:
            res_ids = [self.curr_ind]

        orig_data = self.get_params_n(res_ids, 'alpha')
        if meas_cols is None:
            meas_cols = orig_data.columns.tolist()

        ncols = min([len(meas_cols), ncols])
        nrows = np.ceil(len(meas_cols)/ncols).astype(int)

        if figsize is None:
            figsize = (ncols*4, nrows*4.3)

        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, layout='tight', sharey=False, sharex=False)
        if (nrows > 1) | (ncols > 1):
            axes = axes.ravel()
        else:
            axes = [axes]

        default_plot_opts = {'c': 'k', 'label': "Alpha",
                             'marker': 'o', 'markersize': 7, 'markerfacecolor': 'r',
                             'errorbar': 'sd', 'estimator': 'mean', 'err_style': 'bars',
                             'err_kws': {'capsize': 2.5}}
        default_plot_opts.update(plot_kwargs)
        print(plot_kwargs)

        for i, ax in enumerate(axes):
            if i < len(meas_cols):
                sns.lineplot(data=orig_data, x=time_col, y=meas_cols[i], ax=axes[i],
                             **default_plot_opts)
                axes[i].get_legend().remove()
                axes[i].set_xlabel(axes[i].get_xlabel(), fontsize=13)
                axes[i].set_title(meas_cols[i], fontsize=13, fontweight='bold')
                axes[i].set_xticks(axes[i].get_xticks())
                axes[i].set_xticklabels(axes[i].get_xticklabels(), fontsize=12)
                axes[i].set_yticks(axes[i].get_yticks())
                axes[i].set_yticklabels(axes[i].get_yticklabels(), fontsize=12)
                axes[i].set_ylabel('Alpha', fontsize=13)
            else:
                axes[i].set_axis_off()

        return fig

    def plot_gamma(self, time_col='Time (WD)', meas_cols=None, ncols=5, figsize=None, res_ids=None, **plot_kwargs):
        if res_ids is None:
            res_ids = [self.curr_ind]

        orig_data = self.get_params_n(res_ids, 'gamma')
        if meas_cols is None:
            meas_cols = orig_data.columns.tolist()

        ncols = min([len(meas_cols), ncols])
        nrows = np.ceil(len(meas_cols)/ncols).astype(int)

        if figsize is None:
            figsize = (ncols*4, nrows*4.3)

        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, layout='tight', sharey=False, sharex=False)
        if (nrows > 1) | (ncols > 1):
            axes = axes.ravel()
        else:
            axes = [axes]

        default_plot_opts = {'c': 'k', 'label': "Gamma",
                             'marker': 'o', 'markersize': 7, 'markerfacecolor': 'r',
                             'errorbar': 'sd', 'estimator': 'mean', 'err_style': 'bars',
                             'err_kws': {'capsize': 2.5}}
        default_plot_opts.update(plot_kwargs)
        print(plot_kwargs)

        for i, ax in enumerate(axes):
            if i < len(meas_cols):
                sns.lineplot(data=orig_data, x=time_col, y=meas_cols[i], ax=axes[i], **default_plot_opts)
                axes[i].get_legend().remove()
                axes[i].set_xlabel(axes[i].get_xlabel(), fontsize=13)
                axes[i].set_xticks(axes[i].get_xticks())
                axes[i].set_xticklabels(axes[i].get_xticklabels(), fontsize=12)
                axes[i].set_yticks(axes[i].get_yticks())
                axes[i].set_yticklabels(axes[i].get_yticklabels(), fontsize=12)
                axes[i].set_ylabel('Gamma', fontsize=13)
            else:
                axes[i].set_axis_off()

        return fig

    def plot_entry_flux(self, time_col='Time (WD)', meas_cols=None, ncols=5, figsize=None, res_ids=None, **plot_kwargs):
        if res_ids is None:
            res_ids = [self.curr_ind]

        orig_data = self.get_params_n(res_ids, 'entry_flux')
        if meas_cols is None:
            meas_cols = orig_data.columns.tolist()

        ncols = min([len(meas_cols), ncols])
        nrows = np.ceil(len(meas_cols)/ncols).astype(int)

        if figsize is None:
            figsize = (ncols*4, nrows*4.3)

        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, layout='tight', sharey=False, sharex=False)
        if (nrows > 1) | (ncols > 1):
            axes = axes.ravel()
        else:
            axes = [axes]

        default_plot_opts = {'c': 'k', 'label': "Entry Flux",
                             'marker': 'o', 'markersize': 7, 'markerfacecolor': 'r',
                             'errorbar': 'sd', 'estimator': 'mean', 'err_style': 'bars',
                             'err_kws': {'capsize': 2.5}}
        default_plot_opts.update(plot_kwargs)
        print(plot_kwargs)

        for i, ax in enumerate(axes):
            if i < len(meas_cols):
                sns.lineplot(data=orig_data, x=time_col, y=meas_cols[i], ax=axes[i], **default_plot_opts)
                axes[i].get_legend().remove()
                axes[i].set_xlabel(axes[i].get_xlabel(), fontsize=13)
                axes[i].set_title(meas_cols[i], fontsize=13, fontweight='bold')
                axes[i].set_xticks(axes[i].get_xticks())
                axes[i].set_xticklabels(axes[i].get_xticklabels(), fontsize=12)
                axes[i].set_yticks(axes[i].get_yticks())
                axes[i].set_yticklabels(axes[i].get_yticklabels(), fontsize=12)
                axes[i].set_ylabel('Entry Flux', fontsize=13)
            else:
                axes[i].set_axis_off()

        return fig

    def plot_vref(self, x_col='Reactions', meas_cols=None, ncols=5, figsize=None, res_ids=None, **plot_kwargs):
        if res_ids is None:
            res_ids = [self.curr_ind]

        orig_data = self.get_params_n(res_ids, 'v_ref')
        if orig_data.index.names[1] is None:
            print(orig_data.index.names[1])
            orig_data = orig_data.rename_axis(index=['trial', x_col])  # The rename step needs to be removed in future

        if meas_cols is None:
            meas_cols = orig_data.columns.tolist()

        ncols = min([len(meas_cols), ncols])
        nrows = np.ceil(len(meas_cols)/ncols).astype(int)

        if figsize is None:
            figsize = (ncols*4, nrows*4.3)

        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, layout='tight', sharey=False, sharex=False)
        if (nrows > 1) | (ncols > 1):
            axes = axes.ravel()
        else:
            axes = [axes]

        default_plot_opts = {'capsize': 0.15, 'label': 'Reference Flux', 'fill': True, 'facecolor': '#F7F6F6', 'edgecolor': 'k', 'linewidth': 1, 'errorbar': 'sd', 'err_kws': {'linewidth': 1}, 'hatch': '\\'}
        default_plot_opts.update(plot_kwargs)
        print(plot_kwargs)

        for i, ax in enumerate(axes):
            if i < len(meas_cols):
                sns.barplot(data=orig_data, x=x_col, y=meas_cols[i], ax=axes[i], **default_plot_opts)
                axes[i].get_legend().remove()
                axes[i].set_xlabel(axes[i].get_xlabel(), fontsize=13)
                axes[i].set_title(meas_cols[i], fontsize=13, fontweight='bold')
                axes[i].set_xticks(axes[i].get_xticks())
                axes[i].set_xticklabels(axes[i].get_xticklabels(), fontsize=12)
                axes[i].set_yticks(axes[i].get_yticks())
                axes[i].set_yticklabels(axes[i].get_yticklabels(), fontsize=12)
                axes[i].set_ylabel('Reference Flux', fontsize=13)
                axes[i].set_ylim((0, axes[i].get_ylim()[1]))
            else:
                axes[i].set_axis_off()
        
        return fig

    def plot_secretedflux(self, time_col, meas_cols=None, ncols=5, figsize=None, res_ids=None, **plot_kwargs):
        if res_ids is None:
            res_ids = [self.curr_ind]

        orig_data = self.get_params_n(res_ids, 'secreted_flux')
        if meas_cols is None:
            meas_cols = orig_data.columns.tolist()

        ncols = min([len(meas_cols), ncols])
        nrows = np.ceil(len(meas_cols)/ncols).astype(int)

        if figsize is None:
            figsize = (ncols*4, nrows*4.3)

        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, layout='tight', sharey=False, sharex=False)
        if (nrows > 1) | (ncols > 1):
            axes = axes.ravel()
        else:
            axes = [axes]

        default_plot_opts = {'c': 'k', 'label': "Secreted Flux",
                             'marker': 'o', 'markersize': 7, 'markerfacecolor': 'r',
                             'errorbar': 'sd', 'estimator': 'mean', 'err_style': 'bars',
                             'err_kws': {'capsize': 2.5}}
        default_plot_opts.update(plot_kwargs)
        print(plot_kwargs)

        for i, ax in enumerate(axes):
            if i < len(meas_cols):
                sns.lineplot(data=orig_data, x=time_col, y=meas_cols[i], ax=axes[i],
                             **default_plot_opts)
                axes[i].get_legend().remove()
                axes[i].set_xlabel(axes[i].get_xlabel(), fontsize=13)
                axes[i].set_title(meas_cols[i], fontsize=13, fontweight='bold')
                axes[i].set_xticks(axes[i].get_xticks())
                axes[i].set_xticklabels(axes[i].get_xticklabels(), fontsize=12)
                axes[i].set_yticks(axes[i].get_yticks())
                axes[i].set_yticklabels(axes[i].get_yticklabels(), fontsize=12)
                axes[i].set_ylabel('Secreted Flux', fontsize=13)
            else:
                axes[i].set_axis_off()

        return fig

    def plot_internalflux(self, time_col, meas_cols=None, ncols=5, figsize=None, res_ids=None, **plot_kwargs):
        if res_ids is None:
            res_ids = [self.curr_ind]

        orig_data = (self.get_params_n(res_ids, 'internal_flux').reset_index(names=['trial', time_col]).
                     melt(id_vars=[time_col, 'trial'], var_name='Reaction ID', value_name='Internal Flux').
                     merge(self.rxn_meta.reset_index(), on='Reaction ID', how='left').
                     groupby(['trial', time_col, 'Enzymes'])['Internal Flux'].sum().
                     reset_index(level=2).pivot(columns='Enzymes', values='Internal Flux'))

        if meas_cols is None:
            meas_cols = orig_data.columns.tolist()

        ncols = min([len(meas_cols), ncols])
        nrows = np.ceil(len(meas_cols)/ncols).astype(int)

        if figsize is None:
            figsize = (ncols*4, nrows*4.3)

        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, layout='tight', sharey=False, sharex=False)
        if (nrows > 1) | (ncols > 1):
            axes = axes.ravel()
        else:
            axes = [axes]

        default_plot_opts = {'c': 'k', 'label': "Enzymes",
                             'marker': 'o', 'markersize': 7, 'markerfacecolor': 'r',
                             'errorbar': 'sd', 'estimator': 'mean', 'err_style': 'bars',
                             'err_kws': {'capsize': 2.5}}
        default_plot_opts.update(plot_kwargs)
        print(plot_kwargs)

        for i, ax in enumerate(axes):
            if i < len(meas_cols):
                sns.lineplot(data=orig_data, x=time_col, y=meas_cols[i], ax=axes[i],
                             **default_plot_opts)
                axes[i].get_legend().remove()
                axes[i].set_xlabel(axes[i].get_xlabel(), fontsize=13)
                axes[i].set_title(meas_cols[i], fontsize=13, fontweight='bold')
                axes[i].set_xticks(axes[i].get_xticks())
                axes[i].set_xticklabels(axes[i].get_xticklabels(), fontsize=12)
                axes[i].set_yticks(axes[i].get_yticks())
                axes[i].set_yticklabels(axes[i].get_yticklabels(), fontsize=12)
                axes[i].set_ylabel('Internal Flux', fontsize=13)
            else:
                axes[i].set_axis_off()

        return fig

    def plot_betas(self, time_col, meas_cols=None, ncols=5, figsize=None, res_ids=None, orig_kwargs=None, smooth_kwargs=None):
        if res_ids is None:
            res_ids = [self.curr_ind]

        orig_data = self.get_params_n(res_ids, 'secreted_flux_ratio')
        smooth_data = self.get_params_n(res_ids, 'beta')
        if meas_cols is None:
            meas_cols = orig_data.columns.tolist()

        ncols = min([len(meas_cols), ncols])
        nrows = np.ceil(len(meas_cols)/ncols).astype(int)

        if figsize is None:
            figsize = (ncols*4, nrows*4.3)

        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, layout='tight', sharey=False, sharex=False)
        if (nrows > 1) | (ncols > 1):
            axes = axes.ravel()
        else:
            axes = [axes]

        orig_plot_opts = {'c': 'r', 'label': "Secreted Flux Ratio",
                          'marker': 'o', 'markersize': 7, 'markerfacecolor': 'r',
                          'errorbar': 'sd', 'estimator': 'mean', 'err_style': 'bars',
                          'err_kws': {'capsize': 2.5}}
        smooth_plot_opts = {'c': 'k', 'label': "Beta",
                            'marker': 'o', 'markersize': 7, 'markerfacecolor': 'k',
                            'errorbar': 'sd', 'estimator': 'mean', 'err_style': 'bars',
                            'err_kws': {'capsize': 2.5}}

        if orig_kwargs is not None:
            orig_plot_opts.update(orig_kwargs)
            print(orig_kwargs)

        if smooth_kwargs is not None:
            smooth_plot_opts.update(smooth_kwargs)
            print(smooth_kwargs)

        for i, ax in enumerate(axes):
            if i < len(meas_cols):
                sns.lineplot(data=orig_data, x=time_col, y=meas_cols[i], ax=axes[i], **orig_plot_opts)
                sns.lineplot(data=smooth_data, x=time_col, y=meas_cols[i], ax=axes[i], **smooth_plot_opts)
                axes[i].get_legend().remove()
                axes[i].set_xlabel(axes[i].get_xlabel(), fontsize=13)
                axes[i].set_title(meas_cols[i], fontsize=13, fontweight='bold')
                axes[i].set_xticks(axes[i].get_xticks())
                axes[i].set_xticklabels(axes[i].get_xticklabels(), fontsize=12)
                axes[i].set_yticks(axes[i].get_yticks())
                axes[i].set_yticklabels(axes[i].get_yticklabels(), fontsize=12)
                axes[i].set_ylabel('Beta', fontsize=13)
            else:
                axes[i].set_axis_off()

        handles, labels = axes[len(meas_cols)-1].get_legend_handles_labels()
        fig.legend(handles=handles, loc='lower center', ncol=2, bbox_to_anchor=(0.5, 1),
                   frameon=True, fontsize=14, title_fontproperties={'weight': 'bold', 'size': 14})

        return fig

    def plot_fracs(self, time_col, meas_cols=None, ncols=5, figsize=None, res_ids=None, orig_kwargs=None, smooth_kwargs=None):
        if res_ids is None:
            res_ids = [self.curr_ind]

        orig_data = self.get_params_n([self.curr_ind], 'orig_glycofracs')
        smooth_data = self.get_params_n(res_ids, 'pred_glycofracs')
        if meas_cols is None:
            meas_cols = orig_data.columns.tolist()

        ncols = min([len(meas_cols), ncols])
        nrows = np.ceil(len(meas_cols)/ncols).astype(int)

        if figsize is None:
            figsize = (ncols*4, nrows*4.3)

        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, layout='tight', sharey=False, sharex=False)
        if (nrows > 1) | (ncols > 1):
            axes = axes.ravel()
        else:
            axes = [axes]

        orig_plot_opts = {'c': 'r', 'label': "Measured",
                          'marker': 'o', 's': 7**2}
        smooth_plot_opts = {'c': 'k', 'label': "Predicted (GFA)",
                            'marker': 'o', 'markersize': 7, 'markerfacecolor': 'k',
                            'errorbar': 'sd', 'estimator': 'mean', 'err_style': 'bars',
                            'err_kws': {'capsize': 2.5}}

        if orig_kwargs is not None:
            orig_plot_opts.update(orig_kwargs)
            print(orig_kwargs)

        if smooth_kwargs is not None:
            smooth_plot_opts.update(smooth_kwargs)
            print(smooth_kwargs)

        for i, ax in enumerate(axes):
            if i < len(meas_cols):
                sns.scatterplot(data=orig_data, x=time_col, y=meas_cols[i], ax=axes[i], **orig_plot_opts)
                sns.lineplot(data=smooth_data, x=time_col, y=meas_cols[i], ax=axes[i], **smooth_plot_opts)
                axes[i].get_legend().remove()
                axes[i].set_xlabel(axes[i].get_xlabel(), fontsize=13)
                axes[i].set_title(meas_cols[i], fontsize=13, fontweight='bold')
                axes[i].set_xticks(axes[i].get_xticks())
                axes[i].set_xticklabels(axes[i].get_xticklabels(), fontsize=12)
                axes[i].set_yticks(axes[i].get_yticks())
                axes[i].set_yticklabels(axes[i].get_yticklabels(), fontsize=12)
                axes[i].set_ylabel('Fractions', fontsize=13)
            else:
                axes[i].set_axis_off()

        handles, labels = axes[len(meas_cols)-1].get_legend_handles_labels()
        fig.legend(handles=handles, loc='lower center', ncol=2, bbox_to_anchor=(0.5, 1),
                   frameon=True, fontsize=14, title_fontproperties={'weight': 'bold', 'size': 14})

        return fig

    def plot_interactive(self, init_model, port=5000):

        from dash import Dash, html, Input, Output, State, dcc, ctx, no_update  # , callback
        import dash_cytoscape as cyto
        import dash_bootstrap_components as dbc
        import plotly.graph_objects as go
        from glycowork.motif.draw import GlycoDraw
        import urllib.parse

        results = self.results[self.curr_ind]

        palette = sns.color_palette().as_hex()
        edge_meta = init_model.feat_meta.loc[init_model.feat_meta['internal'], :].copy()
        edge_meta.loc[:, 'enz_colors'] = edge_meta.loc[:, 'Enzymes'].map({e: palette[i] for i, e in enumerate(edge_meta['Enzymes'].unique())})

        node_meta = init_model.met_meta.assign(diag_locations=lambda df: df.apply(lambda row: (GlycoDraw(row['Structure'], show_linkage=False, vertical=True).
                                                                                               as_svg(header='<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE svg>')),
                                                                                  axis=1))

        diag = {idx: f"data:image/svg+xml;utf8,{urllib.parse.quote(row['diag_locations'])}" for idx, row in node_meta.iterrows()}

        app = Dash(external_stylesheets=[dbc.themes.BOOTSTRAP])
        app.title = "Glycosylation Flux Analysis"

        ##########################################################################################################
        banner_url = 'url("https://raw.githubusercontent.com/srira25/GlycoCARTA/master/glycocarta_banner.jpg")'
        app_banner = dbc.Col(children=html.H1(children=['GFA Wizard'],
                                              style={'background-image': banner_url,
                                                     'background-fit': 'contain',
                                                     'color': 'white'}),
                             width='100%')

        ##########################################################################################################
        # Create nodes and edges
        nodes_compartment = [{'classes': 'Compartments',
                              'data': {'id': name,
                                       'label': name}} for name in node_meta['Compartments'].unique()]
        nodes = [{'classes': 'glycoforms',
                  'data': {'id': name,
                           'label': name,
                           'parent': node_meta.loc[name, 'Compartments'],
                           'url': diag[name]}} for name in node_meta.index]
        edges = [{'data': {'id': rxn,
                           'source': vals.index[vals < 0][0],
                           'target': vals.index[vals > 0][0],
                           'label': rxn+'-'+edge_meta.loc[rxn, 'Enzymes'],
                           'colors': edge_meta.loc[rxn, 'enz_colors']}}
                 for rxn, vals in init_model.stoich_internal.items()]

        stylesheet = [{'selector': 'node', 'style': {'content': 'data(label)'}},
                      {'selector': '.glycoforms', 'style': {'background-fit': 'contain',
                                                            'background-color': 'white',
                                                            'border-width': 2,
                                                            'border-color': 'black',
                                                            'width': 100,
                                                            'height': 100,
                                                            'background-image': 'data(url)'}},
                      {'selector': 'edge',
                       'style': {'curve-style': 'bezier',  # The default curve style does not work with certain arrows
                                 'width': 2,
                                 'target-arrow-color': 'data(colors)',
                                 'target-arrow-shape': 'vee',
                                 'line-color': 'data(colors)',
                                 'content': 'data(label)',
                                 'text-margin-x': 10,
                                 'text-margin-y': 10,
                                 'edge-text-rotation': 'autorotate'}}]

        main_graph = html.Div(cyto.Cytoscape(id='GFA',
                                             layout={'name': 'breadthfirst', 'roots': ['Man7']},
                                             style={'height': '75vh', 'margin-bottom': '10px', 'padding': '0px'},
                                             stylesheet=stylesheet,
                                             elements=nodes_compartment + nodes + edges,
                                             minZoom=0.01, responsive=True, pan={'x': 0, 'y': 0},
                                             zoomingEnabled=True,
                                             userZoomingEnabled=True,
                                             userPanningEnabled=True,
                                             wheelSensitivity=0.2),
                              style={'margin-left': '10px', 'margin-top': '10px', 'border': '2px black solid', 'width': '45%'},
                              className="p-0 rounded-3")

        ##########################################################################################################
        # Create Dropdowns and Buttons
        def create_dropdown(label_text, dropdown_id, options, value):
            return dbc.Row(children=[html.B(children=[label_text],
                                            style={'width': '20%', 'margin-top': '10px'}),
                                     dcc.Dropdown(id=dropdown_id, value=value, clearable=False,
                                                  options=[{'label': name.capitalize(), 'value': name} for name in options],
                                                  style={'width': '50%', 'margin-top': '2px', 'margin-left': '0px'})])

        layout_dropdown = create_dropdown(label_text='Choose Network Layout:',
                                          dropdown_id='dropdown-update-layout',
                                          options=['breadthfirst', 'grid', 'random', 'circle', 'cose', 'concentric'],
                                          value='breadthfirst')
        layout_root_dropdown = create_dropdown(label_text='Breadthfirst Layout Root:',
                                               dropdown_id='dropdown-update-root',
                                               options=node_meta.index.tolist(),
                                               value='Man7')

        def create_button(label_text, button_id):
            return html.Button(children=[label_text],
                               id=button_id,
                               style={'width': 'auto', 'margin-top': '10px', 'margin-left': '10px'})

        button_save = create_button(label_text='Save as png', button_id='btn-get-png')
        button_reset = create_button(label_text='Reset layout', button_id='btn-reset')
        buttons = dbc.Row(children=[button_save, button_reset])

        ctrl_panel = dbc.Col(children=[layout_dropdown,
                                       layout_root_dropdown,
                                       buttons],
                             style={'margin-left': '10px', 'margin-bottom': '10px'})

        ##########################################################################################################
        jumbotron = html.Div(children=dbc.Container(fluid=True,
                                                    className="py-3",
                                                    id='jumbo-cont'),
                             className="p-0 bg-body-secondary rounded-3",
                             style={'width': '45%', 'margin-left': '10px'})

        ##########################################################################################################

        def create_figure(xaxis_title, yaxis_title):
            return go.Figure(layout={'margin': {"t": 0, "b": 0, "l": 0, "r": 0},
                                     'xaxis': {'showgrid': False},
                                     'xaxis_title': xaxis_title,
                                     'yaxis': {'showgrid': True, 'gridcolor': 'black', 'gridwidth': 1, 'griddash': 'dot'},
                                     'yaxis_title': yaxis_title,
                                     'plot_bgcolor': 'white',
                                     'paper_bgcolor': 'white'})

        def create_tab(figure, idx, label, tab_id):
            return dbc.Tab(children=dbc.Card(children=dbc.CardBody(children=[dcc.Graph(figure=figure,
                                                                                       id=idx,
                                                                                       style={'width': '100%', 'height': '25vh'})],
                                                                   style={'height': '100%'}),
                                             className="mt-2",
                                             style={'width': '100%', 'height': '100%'}),
                           label=label,
                           tab_id=tab_id,
                           style={'height': '100%'})

        fig1 = create_figure('Time (h)', 'Gamma')
        fig1.add_trace(go.Scatter(x=results['gamma'].index,
                                  y=results['gamma']['gamma'],
                                  mode='lines+markers',
                                  name='Gamma'))
        tab1_layout = create_tab(figure=fig1, idx='gamma', label='Gamma', tab_id='tab-gamma')

        fig2 = create_figure('Reactions', 'Internal Flux at t_ref')
        fig2.add_trace(go.Bar(x=results['v_ref'].index,
                              y=results['v_ref']['v_ref']))
        tab2_layout = create_tab(figure=fig2, idx='v_ref', label='v_ref', tab_id='tab-vref')

        results_tabs_general = html.Div(children=dbc.Tabs([tab1_layout, tab2_layout],
                                                          id='tab-layout-gen',
                                                          style={'width': '100%'}),
                                        style={'margin-left': '0px', 'width': '50%', 'height': '100%'})

        static_plot_sections = dbc.Row(children=[jumbotron, results_tabs_general],
                                       justify='left',
                                       style={'height': '48%'})

        ##########################################################################################################
        results_tabs = html.Div(dbc.Tabs(id='tab-layout'),
                                style={'margin-top': '10px', 'margin-left': '0px', 'width': '90%', 'height': '48%'})

        ##########################################################################################################
        sidebar_content = dbc.Col(children=[static_plot_sections,
                                            results_tabs],
                                  style={'margin-left': '10px', 'margin-top': '10px', 'height': '75vh', 'padding': '0px', 'width': '90%'})

        ##########################################################################################################
        # Declare app layout
        app.layout = html.Div(children=dbc.Stack(children=[app_banner,
                                                           ctrl_panel,
                                                           dbc.Row(children=[main_graph,
                                                                             sidebar_content],
                                                                   style={'margin-left': '0px'})]))

        ##########################################################################################################
        # Save image as png/jpg/svg
        @app.callback(Output("GFA", "generateImage"),
                      Input("btn-get-png", "n_clicks"))
        def get_image(get_png_clicks):
            action = 'store'
            if ctx.triggered:
                action = 'download'

            return {'type': 'png', 'action': action}

        # Alter layouts
        @app.callback(Output('GFA', 'layout'),
                      Output('GFA', 'zoom'),
                      Input('dropdown-update-layout', 'value'),
                      Input('dropdown-update-root', 'value'),
                      Input('btn-reset', 'n_clicks'))
        def update_layout(layout, root, n_clicks):
            if layout == 'breadthfirst':
                return {'name': layout, 'roots': [root], 'animate': True}, 1

            return {'name': layout, 'roots': [root], 'animate': True}, 1  # Don't know why but reset button only works if roots is provided

        # Change color/width of selected component on tap of node/edge
        @app.callback(Output('GFA', 'stylesheet'),
                      Input('GFA', 'tapNode'),
                      Input('GFA', 'tapEdge'))
        def changeColorData(node, edge):
            triggered_id = [*ctx.triggered_prop_ids][0].split(".")[-1] if ctx.triggered else ''

            if triggered_id == 'tapNode':
                node_id = node['data']['id']
                new_styles = [{'selector': f'node[id = "{node_id}"]',
                               'style': {'border-color': 'blue'}}]
                return stylesheet+new_styles

            elif triggered_id == 'tapEdge':
                edge_id = edge['data']['id']
                new_styles = [{'selector': f'edge[id = "{edge_id}"]',
                               'style': {'width': 5}}]
                return stylesheet+new_styles

            return stylesheet

        @app.callback(Output('jumbo-cont', 'children'),
                      Input('GFA', 'tapNode'),
                      Input('GFA', 'tapEdge'))
        def displayJumbotronData(node, edge):
            triggered_id = [*ctx.triggered_prop_ids][0].split(".")[-1] if ctx.triggered else ''
            # empty_imgstyle = {'height':'0vh', 'width':'auto'}
            notempty_imgstyle = {'height': '20vh', 'width': 'auto'}
            comps = [html.B(children='Instructions', id='glycoform-name', className="display-6"),
                     html.Hr(className="my-2"),
                     html.P(children=['Click on a node or edge to know more'], id='glycoform-desc')]

            if triggered_id == 'tapNode' and not node['isParent']:
                comps = [html.B(children=node['data']['id'], id='glycoform-name', className="display-6"),
                         html.Hr(className="my-2"),
                         dbc.Stack(children=[html.Img(src=diag[node['data']['id']],
                                                      id='tapNode-image-src',
                                                      style=notempty_imgstyle)],
                                   style={'width': '100%'},
                                   direction='horizontal'),  # , 'transform': 'rotate(90deg)', 'margin-left': '1em'
                         html.P(children=[''], id='glycoform-desc')]

            if triggered_id == 'tapEdge':
                edge_desc = f"Enzyme: {edge_meta.loc[edge['data']['id'], 'Enzymes']}\nGenes: {edge_meta.loc[edge['data']['id'], 'Genes']}"
                comps = [html.B(children=edge['data']['id'], id='glycoform-name', className="display-6"),
                         html.Hr(className="my-2"),
                         dbc.Stack(children=[html.Img(src=diag[edge['data']['source']],
                                                      id='tapNode-image-src',
                                                      style=notempty_imgstyle),
                                             html.B(children=[" \u2192"],
                                                    id='glycoform-edge',
                                                    className="display-6",
                                                    style={'width': 'auto'}),
                                             html.Img(src=diag[edge['data']['target']],
                                                      id='tapNode-image-tgt',
                                                      style=notempty_imgstyle)],
                                   style={'width': '100%'},
                                   direction='horizontal'),  # , 'transform': 'rotate(90deg)', 'margin-left': '1em'
                         html.P(children=[edge_desc], id='glycoform-desc')]

            return comps

        # Show results (from results dictionary) corresponding to selected component on tap of node/edge
        @app.callback(Output('tab-layout', 'children'),
                      Output('tab-layout', 'active_tab'),
                      Input('GFA', 'tapNode'),
                      Input('GFA', 'tapEdge'),
                      State('tab-layout', 'active_tab'))
        def show_results(node, edge, active):
            triggered_id = [*ctx.triggered_prop_ids][0].split(".")[-1] if ctx.triggered else ''
            tabs = []
            if triggered_id == 'tapNode' and not node['isParent']:
                fig1 = create_figure('Time (h)', 'Fractions')
                fig1.add_trace(go.Scatter(x=results['orig_glycofracs'].index,
                                          y=results['orig_glycofracs'].loc[:, node['data']['id']],
                                          mode='lines+markers',
                                          name='Original',
                                          line={'dash': 'dash'}))
                fig1.add_trace(go.Scatter(x=results['pred_glycofracs'].index,
                                          y=results['pred_glycofracs'].loc[:, node['data']['id']],
                                          mode='lines+markers',
                                          name='GFA Predicted'))
                tab1_content = create_tab(figure=fig1, idx='fracs', label='Fractions', tab_id='tab-frac')
                tabs.append(tab1_content)

                fig2 = create_figure('Time (h)', 'Secreted Flux')
                fig2.add_trace(go.Scatter(x=results['secreted_flux'].index,
                                          y=results['secreted_flux'].loc[:, node['data']['id']],
                                          mode='lines+markers',
                                          name='Flux'))
                tab2_content = create_tab(figure=fig2, idx='secflux', label='Secreted Fluxes', tab_id='tab-sec')
                tabs.append(tab2_content)

                if node_meta.loc[node['data']['id'], 'beta_mets'] is True:
                    fig3 = create_figure('Time (h)', 'Beta')
                    fig3.add_trace(go.Scatter(x=results['beta'].index,
                                              y=results['beta'].loc[:, node['data']['id']],
                                              mode='lines+markers',
                                              name='Beta'))
                    tab3_content = create_tab(figure=fig3, idx='beta', label='Beta', tab_id='tab-beta')
                    tabs.append(tab3_content)

                active = active if active in [tab.tab_id for tab in tabs] else 'tab-frac'

            if triggered_id == 'tapEdge':
                fig1 = create_figure('Time (h)', 'Internal Flux')
                fig1.add_trace(go.Scatter(x=results['internal_flux'].index,
                                          y=results['internal_flux'].loc[:, edge['data']['id']],
                                          mode='lines+markers',
                                          name='Internal Flux',
                                          line={'dash': 'dash'}))
                tab1_content = create_tab(figure=fig1, idx='int_flux', label='Internal Flux', tab_id='tab-int')
                tabs.append(tab1_content)

                if edge_meta.loc[edge['data']['id'], 'nonlin_rxns'] is True:
                    fig2 = create_figure('Time (h)', 'Alpha')
                    fig2.add_trace(go.Scatter(x=results['alpha'].index,
                                              y=results['alpha'].loc[:, edge_meta.loc[edge['data']['id'], 'Enzymes']],
                                              mode='lines+markers',
                                              name='Alpha',
                                              line={'dash': 'dash'}))
                    tab2_content = create_tab(figure=fig2, idx='alpha', label='Alpha', tab_id='tab-alpha')
                    tabs.append(tab2_content)

                active = active if active in [tab.tab_id for tab in tabs] else 'tab-int'

            if len(tabs) == 0:
                active = no_update

            return tabs, active

        app.run(debug=True, port=port)

        return nodes + nodes_compartment, edges

    def __repr__(self):
        printer = f'Results class (name={self.name}) with {len(self.results)} experiments (Active={self.curr_ind}):\n' + ', '.join(map(str, self.results.keys()))

        return printer


def plot_meas(time_col, orig_data, smooth_data=None, meas_cols=None, ncols=3, figsize=(5, 5)):
    assert (time_col in orig_data.columns) or (time_col in orig_data.index.names)
    if smooth_data is not None:
        assert (time_col in smooth_data.columns) or (time_col in smooth_data.index.names)

    if meas_cols is None:
        meas_cols = sorted(set(orig_data.columns) - set([time_col]))
    else:
        assert all(col in orig_data.columns for col in meas_cols)

    if smooth_data is not None:
        assert all(col in smooth_data.columns for col in meas_cols)

    ncols = min([len(meas_cols), ncols])
    nrows = np.ceil(len(meas_cols)/ncols).astype(int)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, layout='tight')
    if (nrows > 1) | (ncols > 1):
        axes = axes.ravel()
    else:
        axes = [axes]

    for i, ax in enumerate(axes):
        if i < len(meas_cols):
            if smooth_data is not None:
                sns.scatterplot(data=orig_data, x=time_col, y=meas_cols[i], c='k', label="Original Noised Data", ax=axes[i])
                sns.lineplot(data=smooth_data, x=time_col, y=meas_cols[i], c='r', label="Fitted Curve", ax=axes[i],
                             errorbar='sd', estimator='mean',
                             err_style='bars', err_kws={'capsize': 2.5})
            else:
                sns.lineplot(data=orig_data, x=time_col, y=meas_cols[i], c='r', marker='o', markerfacecolor='k', label="Original Noised Data", ax=axes[i],
                             errorbar='sd', estimator='mean',
                             err_style='bars', err_kws={'capsize': 2.5})

            axes[i].get_legend().remove()
        else:
            axes[i].set_axis_off()

    if smooth_data is not None:
        handles, labels = axes[len(meas_cols)-1].get_legend_handles_labels()
        fig.legend(handles, labels, loc='lower center', ncol=2, bbox_to_anchor=(0.5, 1), frameon=True)

    return fig


def plot_perturb(time_col, x_col, orig_data, meas_cols=None, ncols=3, figsize=(5, 5)):
    assert (time_col in orig_data.columns) or (time_col in orig_data.index.names)
    assert (x_col in orig_data.columns) or (x_col in orig_data.index.names)
    if meas_cols is None:
        meas_cols = sorted(set(orig_data.columns) - set([time_col]))
    else:
        assert all(col in orig_data.columns for col in meas_cols)

    ncols = min([len(meas_cols), ncols])
    nrows = np.ceil(len(meas_cols)/ncols).astype(int)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharey=False, sharex=False, layout='tight')
    if (nrows > 1) | (ncols > 1):
        axes = axes.ravel()
    else:
        axes = [axes]

    for i, ax in enumerate(axes):
        if i < len(meas_cols):
            sns.barplot(orig_data, x=x_col, y=meas_cols[i], hue=time_col, ax=axes[i])
            axes[i].get_legend().remove()
            axes[i].set_xticks(axes[i].get_xticks())
            axes[i].set_xticklabels(axes[i].get_xticklabels(), rotation=90)
        else:
            axes[i].set_axis_off()

    handles, labels = axes[len(meas_cols)-1].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', ncol=2, bbox_to_anchor=(0.5, 1), frameon=True)

    fig.subplots_adjust(top=0.92)  # adjust the Figure in rp

    return fig


def plot_fig(time_col, data, meas_col, smooth_col=None, ax=None):
    assert (time_col in data.columns) or (time_col in data.index.names)
    assert meas_col in data.columns
    if smooth_col is not None:
        assert smooth_col in data.columns

    if ax is None:
        ax = sns.scatterplot(data=data, x=time_col, y=meas_col, marker='o', c='k', label="Original Noised Data")
    else:
        sns.scatterplot(data=data, x=time_col, y=meas_col, marker='o', c='k', label="Original Noised Data", ax=ax)

    if smooth_col is not None:
        sns.lineplot(data=data, x=time_col, y=smooth_col, c='r', label="Fitted Curve", ax=ax,
                     errorbar='sd', estimator='mean',
                     err_style='bars', err_kws={'capsize': 2.5})
    else:
        sns.lineplot(data=data, x=time_col, y=meas_col, c='r', ax=ax,
                     errorbar='sd', estimator='mean',
                     err_style='bars', err_kws={'capsize': 2.5})
        ax.get_legend().remove()
    ax.set_xlabel(time_col)
    fig = ax.get_figure()

    return fig