import warnings
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from sympy import symbols, lambdify, parse_expr, diff
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn import metrics

# warnings.filterwarnings("ignore", "use_inf_as_na")  # Generated by seaborn


class curve_fitter:
    """
    A class for fitting data to various mathematical models such as Polynomial, Exponential, Power, Logarithmic, 
    Fourier, Gaussian, Weibull, Hill-type, and Sigmoidal models. It also allows custom model expressions.
    """
    def __init__(self):
        """
        Initializes the CurveFitter object by defining a set of mathematical models and their variations. The user
        can select from these models or provide a custom function.
        """
        self.funcs = {"Polynomial": {'model_expr': {i: parse_expr(' + '.join([f'a{i-j}*x**{i-j}' for j in np.arange(i+1)])) for i in np.arange(1, 10)},
                                     'label': 'Polynomial',
                                     'subopt_name': 'Degree'},
                      "Exponential": {'model_expr': {i: parse_expr(' + '.join([f'a{j} * exp(b{j}*x)' for j in np.arange(1, i+1)])) for i in np.arange(1, 10)},
                                      'label': 'Exponential',
                                      'subopt_name': 'Number of Terms'},
                      "Power": {'model_expr': {'with deviation': parse_expr('a * x**b + c'),
                                               'without deviation': parse_expr('a * x**b')},
                                'label': 'Power Law',
                                'subopt_name': 'Deviation'},
                      "Logarithmic": {'model_expr': {'2': parse_expr('a * log2(b*x) + c'),
                                                     'e': parse_expr('a * log(b*x) + c'),
                                                     '10': parse_expr('10 * log10(b*x) + c')},
                                      'label': 'Logarithmic',
                                      'subopt_name': 'Base'},
                      "Fourier": {'model_expr': {i: parse_expr(' + '.join([f'a{j}*sin({j}*b*x) + c{j}*cos({j}*b*x)' for j in np.arange(1, i+1)]) + f' + a{0}') for i in np.arange(1, 10)},
                                  'label': 'Fourier Series',
                                  'subopt_name': 'Number of Terms'},
                      "Gaussian": {'model_expr': {i: parse_expr(' + '.join([f'a{j} * exp(-(x-b{j})**2 / (c{j}**2))' for j in np.arange(i)])) for i in np.arange(1, 10)},
                                   'label': 'Gaussian',
                                   'subopt_name': 'Number of Terms'},
                      "Weibull": {'model_expr': parse_expr('a * b * x**(b-1) * exp(-a*(x**b))'),
                                  'label': 'Weibull'},
                      "Hill-type": {'model_expr': parse_expr('A / (D + (x / B)**(-n))'),
                                    'label': 'Hill-type'},
                      "Sigmoidal": {'model_expr': {"Logistic": parse_expr('a / (1 + exp(-b*(x - c)))'),
                                                   "4PL": parse_expr('d + (a - d) / (1 + (x / c)**b)'),
                                                   "Alternate Logistic": parse_expr('A / (exp(B * x) + C * exp(-D * x))'),
                                                   "Gompertz": parse_expr('d + (a - d) * exp(exp(-b * (x - c)))')},
                                    'label': 'Sigmoidal',
                                    'subopt_name': 'model'}}

        print('Choose a model from below keys or "Custom:"')
        print(list(self.funcs.keys()))
        #: Name of model to be chosen by user using choose_model.
        self.curr_model_name = None
        #: Model to be chosen by user using choose_model.
        self.model = None
        #: Dictionary of details about fit populated . You can retrieve the fitted and derivative function by using the methods get_fitmodel and get_diffmodel
        self.current_stats = None

    def ingest_data(self, data, x_col):
        """
        Ingests the dataset to be used for curve fitting.

        :param data: The dataset containing the independent and dependent variables.
        :type data: pd.DataFrame
        :param x_col: The column name representing the independent variable (x).
        :type x_col: str
        """
        self.data = data
        self.x_col = x_col

    def choose_model(self, model_name):
        """
        Selects a predefined model or allows the user to choose a custom model. If the model has suboptions, choose_submodel needs to be run.

        :param model_name: The name of the model to choose from the predefined set of models or 'Custom'.
        :type model_name: str 
        """
        self.model = None
        assert model_name in list(self.funcs.keys())+['Custom']
        self.curr_model_name = model_name
        if model_name in list(self.funcs.keys()):
            if isinstance(self.funcs[self.curr_model_name]['model_expr'], dict):
                print(f'{self.curr_model_name} selected. Choose submodel using choose_submodel function from below keys.')
                self.suboptions = list(self.funcs[self.curr_model_name]['model_expr'].keys())
                print(f"{self.funcs[self.curr_model_name]['subopt_name']}:{self.suboptions}")
            else:
                print(f'{self.curr_model_name} selected')
                self.suboptions = None
                self.model = self.funcs[self.curr_model_name]['model_expr']
        elif model_name == 'Custom':
            self.suboptions = None
            print(f'{self.curr_model_name} selected. Enter custom function in Python readable string using custom_func function.')

    def choose_submodel(self, submodel_name):
        """
        Selects a submodel from the chosen model if applicable.

        :param submodel_name: The name of the submodel to select.
        :type submodel_name: str
        """
        assert submodel_name in list(self.funcs[self.curr_model_name]['model_expr'].keys())
        self.model = self.funcs[self.curr_model_name]['model_expr'][submodel_name]

    def custom_func(self, model_string):
        """
        Defines a custom model expression when the 'Custom' model option is selected.

        :param model_string: A string representing the custom mathematical model in Python syntax.
        :type model_string: str
        """
        if self.curr_model_name == 'Custom':
            self.model = parse_expr(model_string)
        else:
            raise ValueError('Run choose_model first')

    def calc_diff(self, fitted_params):
        """
        Calculates the derivative of the fitted model with respect to the independent variable.

        :param fitted_params: The parameters of the fitted model.
        :type fitted_params: dict
        :return: The derivative values of the fitted model for the data's x-values.
        :rtype: np.ndarray
        """
        diff_model = self.get_diffmodel(fitted_params)
        return diff_model(self.data[self.x_col].values)

    def get_diffmodel(self, fitted_params):
        """
        Retrieves the symbolic derivative of the model.

        :param fitted_params: The parameters of the fitted model.
        :type fitted_params: dict
        :return: A function that computes the derivative values.
        :rtype: callable
        """
        assert self.model is not None
        x = symbols('x')
        diff_model = diff(self.model, x)
        temp = diff_model.subs(fitted_params)
        print(f'Derivative of function: {temp}')

        return lambdify([x], temp, modules=['numpy'])

    def get_fitmodel(self, fitted_params):
        """
        Retrieves the symbolic model after substituting the fitted parameters.

        :param fitted_params: The parameters of the fitted model.
        :type fitted_params: dict
        :return: A function representing the fitted model.
        :rtype: callable
        """
        assert self.model is not None
        x = symbols('x')
        model = self.model.subs(fitted_params)

        return lambdify([x], model, modules=['numpy'])

    def fit(self, y_col, start_params=None, **kwargs):
        """
        Fits the chosen model to the dataset.

        :param y_col: The column name representing the dependent variable (y).
        :type y_col: str
        :param start_params: (optional) Initial guesses for the model parameters. If None, default guesses are used.
        :type start_params: list or None
        :param kwargs: Additional keyword arguments to pass to the curve fitting function.
        :return: A tuple containing the fitted parameters, fitted y-values, the statistical results, and the plot axis.
        :rtype: tuple
        """
        xvals = self.data[self.x_col].values
        yvals = self.data[y_col].values
        assert self.model is not None
        x = symbols('x')
        params = sorted(self.model.free_symbols - {x}, key=lambda symbol: symbol.name)
        func = lambdify([x] + params, self.model, modules=['numpy'])
        if start_params is None:
            start_params = [1] * len(params)
        popt, _ = curve_fit(func, xvals, yvals, p0=start_params, method='trf', **kwargs)

        fitted_y = func(xvals, *popt)
        self.equation = str(self.model.subs({e: round(popt[i], 3) for i, e in enumerate(params)}))

        # Calculate R-squared
        fitted_params = {e.name: popt[i] for i, e in enumerate(params)}
        stats_df = pd.DataFrame(fitted_params, index=['Value'])
        stats_df.loc[:, 'R-squared'] = metrics.r2_score(yvals, fitted_y)
        stats_df.loc[:, 'RMSE'] = metrics.root_mean_squared_error(yvals, fitted_y)
        stats_df.loc[:, 'Equation'] = str(self.model)

        stats_df = stats_df.T.reset_index(names='Parameter')

        ax = sns.scatterplot(x=xvals, y=yvals, label='Original', c='r')
        ax.set_xlabel(self.x_col)
        ax.set_ylabel(y_col)
        fit_x_plot = np.linspace(xvals.min(), xvals.max(), num=1000)
        fit_y_plot = func(fit_x_plot, *popt)
        sns.lineplot(x=fit_x_plot, y=fit_y_plot, ax=ax, label='Fitted')

        return fitted_params, fitted_y, stats_df, ax

    def fit_jupyter(self, y_col):
        """
        Fits the chosen model to the dataset and displays the fitted results interactively in a Jupyter notebook environment.

        :param y_col: The column name representing the dependent variable (y).
        :type y_col: str
        """
        # Create widgets
        import ipywidgets as widgets
        from IPython.display import display, clear_output

        model_dropdown = widgets.Dropdown(options=list(self.funcs.keys()) + ['Custom'], value=None, description='Model:')
        submodel_dropdown = widgets.Dropdown(description='Submodel:')
        custom_func_text = widgets.Text(description='Custom Func:')
        fit_button = widgets.Button(description='Fit Curve', layout=widgets.Layout(margin='10px 0px 0px 90px'))

        submodel_dropdown.layout.display = 'none'
        custom_func_text.layout.display = 'none'
        output = widgets.Output()

        widgets_row = widgets.VBox([widgets.HBox([model_dropdown, submodel_dropdown, custom_func_text]), fit_button],
                                   layout=widgets.Layout(align_items='flex-start'))

        def on_model_change(change):
            self.choose_model(change['new'])
            if self.curr_model_name in list(self.funcs.keys()):
                if isinstance(self.funcs[self.curr_model_name]['model_expr'], dict):
                    submodel_dropdown.options = list(self.funcs[self.curr_model_name]['model_expr'].keys())
                    submodel_dropdown.layout.display = 'block'
                    custom_func_text.layout.display = 'none'
                else:
                    submodel_dropdown.layout.display = 'none'
                    custom_func_text.layout.display = 'none'
            else:
                submodel_dropdown.layout.display = 'none'
                custom_func_text.layout.display = 'block'

        def on_submodel_change(change):
            self.choose_submodel(change['new'])

        def on_custom_func_submit(change):
            self.custom_func(change['new'])

        model_dropdown.observe(on_model_change, names='value')
        submodel_dropdown.observe(on_submodel_change, names='value')
        custom_func_text.observe(on_custom_func_submit, names='value')

        display(widgets_row)

        @output.capture()
        def fit_curve(b):
            clear_output()
            stats, fitted_y, df, ax = self.fit(y_col, max_nfev=5000)
            diff_y = self.calc_diff(stats)
            display(df)
            plt.show()
            self.current_stats = {'params': stats, 'fitted': fitted_y, 'deriv': diff_y, 'detailed': df}

        fit_button.on_click(fit_curve)
        display(output)

    # def interactive_fitter(self, y_col):
    #     import plotly.graph_objs as go
    #     import dash
    #     from dash import dcc, html, Input, Output, State, dash_table
    #     import dash_bootstrap_components as dbc

    #     xvals = self.data[self.x_col].values
    #     yvals = self.data[y_col].values

    #     model_options = [{'label': v['label'], 'value': k}for k, v in self.funcs.items()] + [{'label': 'Custom', 'value': 'Custom'}]
    #     x = symbols('x')

    #     app = dash.Dash(__name__, suppress_callback_exceptions=True)
    #     option_selector = dbc.Row([html.B(children=['Model: '],
    #                                       style={'width': '20%', 'margin-top': '10px'}),
    #                                dcc.Dropdown(id='function-select',
    #                                             options=model_options,
    #                                             value=model_options[0]['label'],
    #                                             clearable=False,
    #                                             style={'width': '50%', 'margin-top': '2px', 'margin-left': '0px'})],
    #                               justify='left')

    #     fit_button = html.Button('Fit', id='fit-button', n_clicks=0, style={'marginTop': 20})

    #     def create_figure(xaxis_title, yaxis_title):
    #         return go.Figure(layout={'margin': {"t": 0, "b": 0, "l": 0, "r": 0},
    #                                  'xaxis': {'showgrid': False},
    #                                  'xaxis_title': xaxis_title,
    #                                  'yaxis': {'showgrid': True, 'gridcolor': 'black', 'gridwidth': 1, 'griddash': 'dot'},
    #                                  'yaxis_title': yaxis_title,
    #                                  'plot_bgcolor': 'white',
    #                                  'paper_bgcolor': 'white',
    #                                  'legend': {'x': 0, 'y': 1},
    #                                  'hovermode': 'closest'})
    #     fig = create_figure('Time', 'Data')
    #     fig.add_trace(go.Scatter(x=xvals, y=yvals, mode='markers', name='Original Data'))
    #     main_graph = dbc.Row([dcc.Graph(figure=fig, id='curve-fitting-plot', style={'display': 'inline-block', 'width': '49%'}),
    #                           html.Div(id='fit-stats', style={'display': 'inline-block', 'width': '49%', 'verticalAlign': 'top'})])
    #     fitted_equation = html.Div(id='fitted-equation',
    #                                style={'marginTop': 20})

    #     app.layout = html.Div(dbc.Stack([option_selector,
    #                                      dbc.Row([], id='suboptions', style={'marginTop': 20}),
    #                                      dbc.Row([fit_button, fitted_equation]),
    #                                      main_graph]))

    #     @app.callback(Output('suboptions', 'children'),
    #                   Input('function-select', 'value'))
    #     def display_subfunc(option):
    #         if option not in ['Weibull', 'Custom']:
    #             suboption_label = html.B(children=[self.funcs[option]['subopt_name']+':'],
    #                                      style={'width': '20%', 'margin-top': '10px'})
    #             suboption_selector = dcc.Dropdown(id='subfunction-select',
    #                                               options=list(self.funcs[option]['model_expr'].keys()),
    #                                               clearable=False,
    #                                               style={'width': '50%', 'margin-top': '2px', 'margin-left': '0px'})
    #             return [suboption_label, suboption_selector]
    #         elif option == 'Custom':
    #             custom_func_label = html.B(children=['Enter custom function:'],
    #                                        style={'width': '20%', 'marginTop': '10px'})
    #             custom_func_input = dcc.Input(id='subfunction-select',
    #                                           type='text',
    #                                           placeholder='e.g., a*x**3 + b*x + c',
    #                                           style={'marginTop': 20, 'width': '15%'})
    #             return [custom_func_label, custom_func_input]
    #         else:
    #             return dcc.Input(id='subfunction-select', style={'display': 'none'})

    #     @app.callback([Output('curve-fitting-plot', 'figure'),
    #                    Output('fitted-equation', 'children'),
    #                    Output('fit-stats', 'children')],
    #                   [Input('fit-button', 'n_clicks')],
    #                   [State('function-select', 'value'),
    #                    State('subfunction-select', 'value')], prevent_initial_call=True)
    #     def update_figure(n_clicks, selected_func_name, selected_subfunc):
    #         fig = create_figure('Time', 'Data')
    #         fig.add_trace(go.Scatter(x=xvals, y=yvals, mode='markers', name='Original Data'))

    #         if selected_func_name not in ['Weibull', 'Custom'] and selected_subfunc:
    #             selected_func = self.funcs[selected_func_name]['model_expr'][selected_subfunc]
    #         elif selected_func_name in ['Weibull']:
    #             selected_func = self.funcs[selected_func_name]['model_expr']
    #         elif selected_func_name == "Custom" and selected_subfunc:
    #             selected_func = parse_expr(selected_subfunc)
    #         else:
    #             return fig, "", ""

    #         params = sorted(selected_func.free_symbols - {x}, key=lambda symbol: symbol.name)
    #         func = lambdify([x] + params, selected_func, modules=['numpy'])
    #         try:
    #             popt, _ = curve_fit(func, xvals, yvals, p0=[1] * len(params), method='trf')
    #         except Exception as e:
    #             equation = f"Fit Error: {str(e)}"
    #             return fig, equation, ""

    #         fitted_y = func(xvals, *popt)
    #         equation = str(selected_func.subs({e: round(popt[i], 3) for i, e in enumerate(params)}))

    #         # Calculate R-squared
    #         ss_res = np.sum((yvals - fitted_y) ** 2)
    #         ss_tot = np.sum((yvals - np.mean(yvals)) ** 2)
    #         r_squared = 1 - (ss_res / ss_tot)

    #         # Create the plot
    #         fig.add_trace(go.Scatter(x=xvals, y=fitted_y, mode='lines', name='Fitted Curve'))

    #         stats_df = pd.DataFrame({e.name: popt[i] for i, e in enumerate(params)}, index=['Value'])
    #         stats_df.loc[:, 'R-squared'] = r_squared
    #         stats_df.loc[:, 'Equation'] = str(selected_func)

    #         stats_df = stats_df.T.reset_index(names='Parameter')
    #         stats_table = dash_table.DataTable(id="table",
    #                                            columns=[{"name": i, "id": i} for i in stats_df.columns],
    #                                            data=stats_df.to_dict("records"),
    #                                            export_format="csv",
    #                                            export_columns='all')

    #         return fig, equation, stats_table

    #     app.run_server(debug=True, port=8265)