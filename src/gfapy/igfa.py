import numpy as np
import pandas as pd
from tqdm.autonotebook import tqdm
import pyomo.environ as pyo
import pyomo.dae as pydae
import seaborn as sns
import matplotlib.pyplot as plt
from pyomo.opt import SolverStatus, TerminationCondition
from .model import parse_model_generic
from .plotter import res_plotter_generic


class parse_model(parse_model_generic):
    def __init__(self, stoich, feat_meta, met_meta, n_comp, lin_rxns, name=None, author=None):
        super().__init__(stoich, feat_meta, met_meta, n_comp, lin_rxns, name, author)

    def read_measurements(self, data, time_col_name='Time (WD)'):
        '''
        met_meta and frac should have same glycoform IDs
        All keys in data should have timepoints (float/int) as index
        '''
        self.time_col_name = time_col_name
        self.spec_prod = data['q_prod'].sort_index(axis='index')
        self.glyco_fractions = data['Fractions'].sort_index(axis='index')
        if 'fit_Fractions' in data.keys():
            self.glyco_fractions_smooth = data['fit_Fractions'].sort_index(axis='index')

        if 'diff_Fractions' in data.keys():
            self.glyco_fractions_dot = data['diff_Fractions'].sort_index(axis='index')

        self._read_core(data, time_data=self.glyco_fractions.index.tolist())

    def fetch_data(self, vcd_col='VCD (1E6 VC/mL)', titer_col='Titer (g/L)', spec_prod_col='Spec Prod (pg/cells/day)', get_diff=False):
        titer = self.titer.loc[self.time_col, :]
        vcd = self.vcd.loc[self.time_col, :]
        spec_prod = self.spec_prod.loc[self.time_col, :]
        glyco_fractions = self.glyco_fractions.loc[self.time_col, :]
        org_data = {None: {'timepoints': {None: self.time_col},
                           'timepoints_data': {None: self.time_col},
                           'titer': titer[titer_col].to_dict(),
                           'vcd': vcd[vcd_col].to_dict(),
                           'spec_prod': spec_prod[spec_prod_col].to_dict(),
                           'glyco_frac_data': (glyco_fractions.reset_index(names=self.time_col_name).
                                               melt(id_vars=self.time_col_name).
                                               set_index([self.time_col_name, 'variable']).squeeze().to_dict())}}

        if get_diff is True:
            glyco_fractions_dot = self.glyco_fractions_dot.loc[self.time_col, :]
            org_data[None]['glyco_frac_dot'] = (glyco_fractions_dot.reset_index(names=self.time_col_name).
                                                melt(id_vars=self.time_col_name).
                                                set_index([self.time_col_name, 'variable']).squeeze().to_dict())
        return org_data

    def _init_pyomomodel(self):
        model = pyo.AbstractModel()
        model.name = 'iGFA'
        model.n_comp = self.n_comp
        model.substrate = pyo.Set(initialize=self.substrate, dimen=1)
        model.mets = pyo.Set(initialize=self.met_meta.index.tolist(), dimen=1)
        model.enzymes = pyo.Set(initialize=self.enzymes.cat.categories.tolist(), dimen=1)
        model.beta_mets = pyo.Set(initialize=self.met_meta[self.met_meta.beta_mets].index.tolist(), dimen=1)

        model.timepoints_data = pyo.Set(dimen=1)
        model.timepoints = pydae.ContinuousSet()

        model.nonlin_rxns = pyo.Set(initialize=self.feat_meta[(self.feat_meta.nonlin_rxns)].index.tolist(), dimen=1)
        model.internal_rxns = pyo.Set(initialize=self.feat_meta[self.feat_meta.internal].index.tolist(), dimen=1)
        model.lin_rxns = pyo.Set(initialize=model.internal_rxns-model.nonlin_rxns, dimen=1)

        # Parameters
        @model.Param(model.mets, model.internal_rxns)
        def stoich(m, i, j):
            return self.stoich_internal.loc[i, j]

        # Variables
        model.entry_flux = pyo.Var(model.timepoints, model.substrate, domain=pyo.NonNegativeReals, bounds=(0, 100))
        model.lin_flux = pyo.Var(model.timepoints, model.lin_rxns, domain=pyo.NonNegativeReals)
        model.v_ref = pyo.Var(model.nonlin_rxns, domain=pyo.NonNegativeReals, bounds=(0, 100))

        def alpha_bounds(m, i, j):
            return (1, 1) if i == m.timepoints.at(1) else (0, 10)
        model.alpha = pyo.Var(model.timepoints, model.enzymes, domain=pyo.NonNegativeReals, bounds=alpha_bounds)

        def beta_bounds(m, i, j):
            return (1, 1) if i == m.timepoints.at(1) else (0, 15)
        model.beta = pyo.Var(model.timepoints, model.beta_mets, domain=pyo.NonNegativeReals, bounds=beta_bounds)

        def gamma_bounds(m, i):
            return (1, 1) if i == m.timepoints.at(1) else (0, 10)
        model.gamma = pyo.Var(model.timepoints, domain=pyo.NonNegativeReals, bounds=gamma_bounds)

        # Expressions
        @model.Expression(model.timepoints, model.internal_rxns)
        def internal_flux(m, i, j):
            if j in m.nonlin_rxns:
                return m.alpha[i, self.enzymes[j]]*m.beta[i, self.producers[j]]*m.gamma[i]*m.v_ref[j]
            else:
                return m.lin_flux[i, j]

        @model.Expression(model.timepoints, model.mets)
        def secreted_flux(m, i, j):
            entry_flux = m.entry_flux[i, j] if j in m.substrate else 0
            return sum(m.stoich[j, k] * m.internal_flux[i, k] for k in m.internal_rxns if m.stoich[j, k] != 0) + entry_flux

        return model

    def _formulate_measured_flux(self, model, formulation='ode'):
        TITER_UNIT_CONV = 1000  # Convert titer units so that Q prod is in pg/cell/day
        if formulation == 'ode':
            model.glyco_frac_dot = pydae.DerivativeVar(model.glyco_frac, wrt=model.timepoints)  # Derivative

            # Constraints
            @model.Constraint(model.timepoints_data, model.mets)  # ODE Equation
            def Xdot(m, i, j):
                return m.secreted_flux[i, j] == (m.titer[i]*TITER_UNIT_CONV/m.vcd[i])*m.glyco_frac_dot[i, j] + m.glyco_frac[i, j]*m.spec_prod[i]
        elif formulation == 'integral':
            @model.Expression(model.timepoints, model.mets)
            def glyco_frac_int(m, i, j):
                return (m.secreted_flux[i, j] - m.glyco_frac[i, j]*m.spec_prod[i])*(m.vcd[i]/(m.titer[i]*TITER_UNIT_CONV))

            @model.Constraint(model.timepoints_data, model.mets)  # ODE Equation
            def Xint(m, i, j):
                if i != m.timepoints_data.at(1):
                    k = m.timepoints_data.prev(i)
                    return m.glyco_frac[i, j] == m.glyco_frac[k, j] + ((i - k)*(m.glyco_frac_int[i, j] + m.glyco_frac_int[k, j])/2)
                else:
                    return pyo.Constraint.Skip
        elif formulation == 'precalc_diff':
            model.glyco_frac_dot = pyo.Param(model.timepoints_data, model.mets)

            @model.Constraint(model.timepoints_data, model.mets)  # ODE Equation
            def Xdot(m, i, j):
                return m.secreted_flux[i, j] == (m.titer[i]*TITER_UNIT_CONV/m.vcd[i])*m.glyco_frac_dot[i, j] + m.glyco_frac[i, j]*m.spec_prod[i]

        return model

    def create_pyomomodel(self, formulation='ode', fit_beta=True, regularize_params=True):
        model = self._init_pyomomodel()

        # Data
        model.titer = pyo.Param(model.timepoints_data)
        model.vcd = pyo.Param(model.timepoints_data)
        model.spec_prod = pyo.Param(model.timepoints_data)
        model.glyco_frac_data = pyo.Param(model.timepoints_data, model.mets)

        model.glyco_frac = pyo.Var(model.timepoints, model.mets, domain=pyo.NonNegativeReals)

        model = self._formulate_measured_flux(model, formulation=formulation)

        @model.Constraint(model.mets)  # Initialize fractions at timepoint t
        def frac_init(m, i):
            return m.glyco_frac[m.timepoints.at(1), i] == m.glyco_frac_data[m.timepoints.at(1), i]

        @model.Constraint(model.timepoints, model.mets)
        def pos_secretions(m, i, j):
            return m.secreted_flux[i, j] >= 0

        # Regularization
        if regularize_params is True:
            model.reg_param = pyo.Param(mutable=True)

            @model.Constraint()
            def min_params(m):
                return (sum((m.alpha[i, j] - 1)**2 for i in m.timepoints for j in m.enzymes) +
                        # sum(abs(m.beta[i, j] - 1) for i in m.timepoints for j in m.beta_mets) +
                        sum((m.gamma[i] - 1)**2 for i in m.timepoints)) <= m.reg_param

            # @model.Expression()
            # def min_params(m):
            #     return (sum((m.internal_flux[i, j]**2)/m.spec_prod[i] for i in m.timepoints for j in m.internal_rxns))/len(m.internal_flux)

        # Expressions needed for objective
        @model.Expression()
        def fit_fracs(m):
            return sum(((m.glyco_frac[i, j]-m.glyco_frac_data[i, j])**2)/abs(m.glyco_frac_data[i, j]) for i in m.timepoints_data for j in m.mets)

        @model.Expression()
        def fit_spec_prod(m):
            return sum(((sum(m.secreted_flux[i, j] for j in m.mets) - m.spec_prod[i])**2)/m.spec_prod[i] for i in m.timepoints_data)

        @model.Expression()
        def fit_betas(m):
            return sum(((m.secreted_flux[i, j] - (m.beta[i, j]*m.secreted_flux[m.timepoints.at(1), j]))**2)/m.spec_prod[i]
                       for i in m.timepoints for j in m.beta_mets)

        # @model.Constraint(model.timepoints, model.beta_mets)
        # def fit_betas(m, i, j):
        #     return m.secreted_flux[i, j] == m.beta[i, j]*m.secreted_flux[m.timepoints.at(1), j]

        # @model.Constraint()
        # model.bilevel_weight = pyo.Param(mutable=True, initialize=1e-4)
        # def bilevel_objective(m):
        #     return (m.fit_fracs + m.fit_spec_prod) <= m.bilevel_weight

        # Objective
        @model.Objective(sense=pyo.minimize)
        def obj(m):
            added_objective = m.fit_betas if fit_beta is True else 0
            return (m.fit_fracs + m.fit_spec_prod + added_objective)

        # @model.Objective(sense=pyo.minimize)
        # def obj(m):
        #     return (m.fit_fracs + m.fit_spec_prod)

        return model

    def get_results(self, instance):
        results = {}
        results['objectives'] = {e: pyo.value(instance.find_component(e)) for e in ['fit_betas', 'fit_fracs', 'fit_spec_prod', 'obj', 'min_params']}  #
        results['pred_glycofracs'] = (pd.DataFrame({j: {i: pyo.value(instance.glyco_frac[i, j]) for i in instance.timepoints_data} for j in instance.mets}).
                                      rename_axis(index=self.time_col_name))
        results['orig_glycofracs'] = (pd.DataFrame({j: {i: pyo.value(instance.glyco_frac_data[i, j]) for i in instance.timepoints_data} for j in instance.mets}).
                                      rename_axis(index=self.time_col_name))

        glyco_frac_dot = instance.find_component('glyco_frac_dot')
        if (glyco_frac_dot is not None) and (glyco_frac_dot.active is True):
            results['diff_glycofracs'] = (pd.DataFrame({j: {i: pyo.value(glyco_frac_dot[i, j]) for i in instance.timepoints_data} for j in instance.mets}).
                                          rename_axis(index=self.time_col_name))

        results['alpha'] = (pd.DataFrame({j: {i: pyo.value(instance.alpha[i, j]) for i in instance.timepoints_data} for j in instance.enzymes}).
                            rename_axis(index=self.time_col_name))
        results['beta'] = (pd.DataFrame({j: {i: pyo.value(instance.beta[i, j]) for i in instance.timepoints_data} for j in instance.beta_mets}).
                           rename_axis(index=self.time_col_name))
        results['gamma'] = pd.Series({i: pyo.value(instance.gamma[i]) for i in instance.timepoints_data}).to_frame(name='gamma').rename_axis(index=self.time_col_name)
        results['v_ref'] = (pd.Series({i: pyo.value(instance.v_ref[i]) for i in instance.nonlin_rxns}).to_frame(name='v_ref').
                            rename_axis(index='Reaction ID'))
        results['secreted_flux'] = (pd.DataFrame({j: {i: pyo.value(instance.secreted_flux[i, j]) for i in instance.timepoints_data} for j in instance.mets}).
                                    rename_axis(index=self.time_col_name))
        results['internal_flux'] = (pd.DataFrame({j: {i: pyo.value(instance.internal_flux[i, j]) for i in instance.timepoints_data} for j in instance.internal_rxns}).
                                    rename_axis(index=self.time_col_name))
        results['entry_flux'] = (pd.DataFrame({j: {i: pyo.value(instance.entry_flux[i, j]) for i in instance.timepoints_data} for j in instance.substrate}).
                                 rename_axis(index=self.time_col_name))
        results['secreted_flux_ratio'] = (pd.DataFrame({j: {i: pyo.value(instance.secreted_flux[i, j]/instance.secreted_flux[instance.timepoints.at(1), j]) for i in instance.timepoints_data} for j in instance.beta_mets}).
                                          rename_axis(index=self.time_col_name))

        return results

    def create_perturb_model(self):
        model = pyo.AbstractModel()
        model.name = 'iGFA_perturb'
        model.n_comp = self.n_comp
        model.substrate = pyo.Set(initialize=self.substrate, dimen=1)
        model.mets = pyo.Set(initialize=self.met_meta.index.tolist(), dimen=1)
        model.enzymes = pyo.Set(initialize=self.enzymes.cat.categories.tolist(), dimen=1)
        model.beta_mets = pyo.Set(initialize=self.met_meta[self.met_meta.beta_mets].index.tolist(), dimen=1)

        model.nonlin_rxns = pyo.Set(initialize=self.feat_meta[(self.feat_meta.nonlin_rxns)].index.tolist(), dimen=1)
        model.internal_rxns = pyo.Set(initialize=self.feat_meta[self.feat_meta.internal].index.tolist(), dimen=1)
        model.lin_rxns = pyo.Set(initialize=model.internal_rxns-model.nonlin_rxns, dimen=1)

        # Parameters
        @model.Param(model.mets, model.internal_rxns)
        def stoich(m, i, j):
            return self.stoich_internal.loc[i, j]

        model.timepoints = pydae.ContinuousSet()
        model.alpha_deltapert = pyo.Param(model.timepoints, model.enzymes, domain=pyo.NonNegativeReals, initialize=1)
        model.vref_deltapert = pyo.Param(model.internal_rxns, domain=pyo.NonNegativeReals, initialize=1)

        # Data about unperturbed conditions
        model.internal_flux_nominal = pyo.Param(model.timepoints, model.internal_rxns, domain=pyo.NonNegativeReals)
        model.secreted_flux_nominal = pyo.Param(model.timepoints, model.mets, domain=pyo.NonNegativeReals)
        model.entry_flux_nominal = pyo.Param(model.timepoints, model.substrate, domain=pyo.NonNegativeReals)

        model.beta_deltapert = pyo.Var(model.timepoints, model.beta_mets, domain=pyo.NonNegativeReals, bounds=(0, 15))

        @model.Expression(model.timepoints, model.internal_rxns)
        def internal_flux_perturb(m, i, j):
            return m.alpha_deltapert[i, self.enzymes[j]]*m.beta_deltapert[i, self.producers[j]]*m.vref_deltapert[j]*m.internal_flux_nominal[i, j]

        @model.Expression(model.timepoints, model.mets)
        def secreted_flux_perturb(m, i, j):
            entry_flux = m.entry_flux_nominal[i, j] if j in m.substrate else 0
            return sum(m.stoich[j, k] * m.internal_flux_perturb[i, k] for k in m.internal_rxns if m.stoich[j, k] != 0) + entry_flux

        @model.Constraint(model.timepoints, model.beta_mets)
        def fit_betas(m, i, j):
            return m.secreted_flux_perturb[i, j] == m.beta_deltapert[i, j]*m.secreted_flux_nominal[i, j]

        @model.Constraint(model.timepoints, model.mets)
        def pos_secretions(m, i, j):
            return m.secreted_flux_perturb[i, j] >= 0

        @model.Expression()
        def min_betas(m):
            return sum(((m.beta_deltapert[i, j] - 1)**2) for i in m.timepoints for j in m.beta_mets)

        @model.Objective(sense=pyo.minimize)
        def obj(m):
            return m.min_betas

        return model


class res_plotter(res_plotter_generic):
    def __init__(self, results, multi_result=False, met_meta=None, rxn_meta=None, name=None, time_col_name='Time (WD)'):
        super().__init__(results, multi_result, met_meta, rxn_meta, name, time_col_name)

    def fetch_nominal_vars(self, result):
        secreted_flux = result['secreted_flux']
        self.time_col = secreted_flux.index.tolist()
        internal_flux = result['internal_flux'].loc[self.time_col, :]
        entry_flux = result['entry_flux'].loc[self.time_col, :]

        org_data = {None: {'timepoints': {None: self.time_col},
                           'secreted_flux_nominal': (secreted_flux.reset_index(names=self.time_col_name).
                                                     melt(id_vars=self.time_col_name).
                                                     set_index([self.time_col_name, 'variable']).squeeze().to_dict()),
                           'internal_flux_nominal': (internal_flux.reset_index(names=self.time_col_name).
                                                     melt(id_vars=self.time_col_name).
                                                     set_index([self.time_col_name, 'variable']).squeeze().to_dict()),
                           'entry_flux_nominal': (entry_flux.reset_index(names=self.time_col_name).
                                                  melt(id_vars=self.time_col_name).
                                                  set_index([self.time_col_name, 'variable']).squeeze().to_dict())}}

        return org_data

    @staticmethod
    def get_pertresults(instance, time_col_name):
        perturb_results = {}
        perturb_results['beta_delta'] = (pd.DataFrame({j: {i: pyo.value(instance.beta_deltapert[i, j]) for i in instance.timepoints} for j in instance.beta_mets}).
                                         rename_axis(index=time_col_name))
        perturb_results['secreted_flux_ratio'] = (pd.DataFrame({j: {i: pyo.value(instance.secreted_flux_perturb[i, j]/instance.secreted_flux_nominal[i, j]) for i in instance.timepoints} for j in instance.mets}).
                                                  rename_axis(index=time_col_name))
        perturb_results['secreted_flux'] = (pd.DataFrame({j: {i: pyo.value(instance.secreted_flux_perturb[i, j]) for i in instance.timepoints} for j in instance.mets}).
                                            rename_axis(index=time_col_name))
        perturb_results['internal_flux'] = (pd.DataFrame({j: {i: pyo.value(instance.internal_flux_perturb[i, j]) for i in instance.timepoints} for j in instance.internal_rxns}).
                                            rename_axis(index=time_col_name))
        perturb_results['secreted_flux_diff'] = (pd.DataFrame({j: {i: pyo.value(instance.secreted_flux_perturb[i, j] - instance.secreted_flux_nominal[i, j]) for i in instance.timepoints} for j in instance.mets}).
                                                 rename_axis(index=time_col_name))
        perturb_results['internal_flux_diff'] = (pd.DataFrame({j: {i: pyo.value(instance.internal_flux_perturb[i, j] - instance.internal_flux_nominal[i, j]) for i in instance.timepoints} for j in instance.internal_rxns}).
                                                 rename_axis(index=time_col_name))
        perturb_results['entry_flux'] = (pd.DataFrame({j: {i: pyo.value(instance.entry_flux_nominal[i, j]) for i in instance.timepoints} for j in instance.substrate}).
                                         rename_axis(index=time_col_name))

        return perturb_results

    def perturb_alpha(self, model, perc_keys, perc_val=0.1, nominal_result=None, solver_opts=None, get_fullresults=False):
        if nominal_result is None:
            nominal_result = self.curr_result

        data = self.fetch_nominal_vars(nominal_result)
        results = {}
        status_run = {'Optimal': [], 'MaxIterations': [], 'Error': []}
        for enz in (pbar := tqdm(perc_keys, total=len(perc_keys), desc='Progress',
                                 bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]')):

            data[None]['alpha_deltapert'] = {(t, enz): (1+perc_val) for t in self.time_col}
            instance = model.create_instance(data=data)
            instance.alpha_deltapert.pprint()

            solver = pyo.SolverFactory('ipopt')
            default_solver_options = {'max_iter': 100000,
                                      'linear_solver': 'mumps',
                                      'warm_start_init_point': 'yes',
                                      'print_level': 4}  # ,
            if solver_opts is not None:
                default_solver_options.update(solver_opts)

            solver.options = default_solver_options
            status = solver.solve(instance, tee=True)

            results_enz = self.get_pertresults(instance, self.time_col_name)
            results_enz['dsecreted_flux'] = results_enz['secreted_flux_diff'].div(perc_val, axis='index')
            results_enz['dinternal_flux'] = results_enz['internal_flux_diff'].div(perc_val, axis='index')
            results_enz['norm_dsecreted_flux'] = results_enz['secreted_flux_diff'].div(perc_val*nominal_result['secreted_flux'])
            results_enz['norm_dinternal_flux'] = results_enz['internal_flux_diff'].div(perc_val*nominal_result['internal_flux'])
            results_enz['status'] = status
            results[enz] = results_enz
            if (status.solver.status == SolverStatus.ok) and (status.solver.termination_condition == TerminationCondition.optimal):
                status_run['Optimal'].append(enz)
            elif (status.solver.status == SolverStatus.warning) and (status.solver.termination_condition == TerminationCondition.maxIterations):
                status_run['MaxIterations'].append(enz)
            else:
                status_run['Error'].append(enz)

            pbar.set_description(f"Progress ({', '.join(f'{k}: {len(v)}' for k, v in status_run.items())})", refresh=True)

        if get_fullresults is True:
            return results
        else:
            perturb_values = {e: (pd.concat({k: v[e] for k, v in results.items() if e in v}, names=['Enzymes'])) for i, e in enumerate(['norm_dsecreted_flux', 'norm_dinternal_flux'])}

            return perturb_values['norm_dsecreted_flux'], perturb_values['norm_dinternal_flux']

    def perturb_vref(self, model, perc_keys, perc_val=0.1, nominal_result=None, solver_opts=None, get_fullresults=False):
        if nominal_result is None:
            nominal_result = self.curr_result

        data = self.fetch_nominal_vars(nominal_result)
        results = {}
        status_run = {'Optimal': [], 'MaxIterations': [], 'Error': []}
        for rxn in (pbar := tqdm(perc_keys, total=len(perc_keys), desc='Progress',
                                 bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]')):
            data[None]['vref_deltapert'] = {rxn: (1+perc_val)}
            instance = model.create_instance(data=data)
            instance.vref_deltapert.pprint()

            solver = pyo.SolverFactory('ipopt')
            default_solver_options = {'max_iter': 100000,
                                      'linear_solver': 'mumps',
                                      'warm_start_init_point': 'yes',
                                      'print_level': 4}  # ,
            if solver_opts is not None:
                default_solver_options.update(solver_opts)

            solver.options = default_solver_options
            status = solver.solve(instance, tee=True)

            results_enz = self.get_pertresults(instance, self.time_col_name)
            results_enz['dsecreted_flux'] = results_enz['secreted_flux_diff'].div(perc_val, axis='index')
            results_enz['dinternal_flux'] = results_enz['internal_flux_diff'].div(perc_val, axis='index')
            results_enz['norm_dsecreted_flux'] = results_enz['secreted_flux_diff'].div(perc_val*nominal_result['secreted_flux'])
            results_enz['norm_dinternal_flux'] = results_enz['internal_flux_diff'].div(perc_val*nominal_result['internal_flux'])
            results_enz['status'] = status
            results[rxn] = results_enz
            if (status.solver.status == SolverStatus.ok) and (status.solver.termination_condition == TerminationCondition.optimal):
                status_run['Optimal'].append(rxn)
            elif (status.solver.status == SolverStatus.warning) and (status.solver.termination_condition == TerminationCondition.maxIterations):
                status_run['MaxIterations'].append(rxn)
            else:
                status_run['Error'].append(rxn)

            pbar.set_description(f"Progress ({', '.join(f'{k}: {len(v)}' for k, v in status_run.items())})", refresh=True)

        if get_fullresults is True:
            return results
        else:
            perturb_values = {e: (pd.concat({k: v[e] for k, v in results.items() if e in v}, names=['Reaction ID'])) for i, e in enumerate(['norm_dsecreted_flux', 'norm_dinternal_flux'])}

            return perturb_values['norm_dsecreted_flux'], perturb_values['norm_dinternal_flux']

    def plot_perturb(self, perturb_data, x_col, time_col, meas_cols=None, ncols=5, figsize=None, **plot_kwargs):
        if meas_cols is None:
            meas_cols = perturb_data.columns.tolist()

        ncols = min([len(meas_cols), ncols])
        nrows = np.ceil(len(meas_cols)/ncols).astype(int)

        if figsize is None:
            figsize = (ncols*4, nrows*4.3)

        fig, axes = plt.subplots(nrows, ncols, figsize=figsize, layout='tight', sharey=False, sharex=False)
        if (nrows > 1) | (ncols > 1):
            axes = axes.ravel()
        else:
            axes = [axes]

        default_plot_opts = {'capsize': 0.15, 'fill': True,
                             'edgecolor': 'k', 'linewidth': 1, 'errorbar': 'sd', 'err_kws': {'linewidth': 1}}
        default_plot_opts.update(plot_kwargs)
        print(plot_kwargs)

        for i, ax in enumerate(axes):
            if i < len(meas_cols):
                sns.barplot(data=perturb_data, x=x_col, y=meas_cols[i], hue=time_col, ax=axes[i], **default_plot_opts)
                axes[i].get_legend().remove()
                axes[i].set_xlabel(axes[i].get_xlabel(), fontsize=13)
                axes[i].set_ylabel('Perturbation', fontsize=13)
                axes[i].set_title(meas_cols[i], fontsize=13, fontweight='bold')
                axes[i].set_xticks(axes[i].get_xticks())
                axes[i].set_xticklabels(axes[i].get_xticklabels(), fontsize=12, rotation=90)
                axes[i].set_yticks(axes[i].get_yticks())
                axes[i].set_yticklabels(axes[i].get_yticklabels(), fontsize=12)
            else:
                axes[i].set_axis_off()

        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles=handles, loc='lower center', ncol=4, bbox_to_anchor=(0.5, 1), title=time_col,
                   frameon=True, fontsize=14, title_fontproperties={'weight': 'bold', 'size': 14})

        return fig