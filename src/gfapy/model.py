import numpy as np
import pandas as pd
from tqdm.autonotebook import tqdm
from collections import namedtuple
import pyomo.environ as pyo
from pyomo.opt import SolverStatus, TerminationCondition
from pyomo.contrib.multistart.reinit import reinitialize_variables


class parse_model_generic():
    def __init__(self, stoich, feat_meta, met_meta, n_comp, lin_rxns, name=None, author=None):
        """
        Initializes the parse_model_generic class.

        Necessary columns in excel file:
        - stoich (index - glycoform id): Should be all numbers.
        - feat_meta (index - reactions): Enzymes, Interactive plotting (Genes).
        - met_meta (index - glycoform id): Interactive plotting (Structure, Compartment).

        :param stoich: Stoichiometry matrix of reactions.
        :type stoich: pd.DataFrame
        :param feat_meta: Feature metadata including enzyme information.
        :type feat_meta: pd.DataFrame
        :param met_meta: Metabolite metadata for glycoform id.
        :type met_meta: pd.DataFrame
        :param n_comp: Number of compartments.
        :type n_comp: int
        :param lin_rxns: List of linear reactions.
        :type lin_rxns: list[str]
        :param name: Name of the model (optional).
        :type name: str, optional
        :param author: Author of the model (optional).
        :type author: str, optional
        """

        self.name = name if name else None
        self.author = author if author else None
        self.n_comp = n_comp
        self.met_meta = met_meta.copy()
        self.feat_meta = feat_meta.copy()
        self.stoich = stoich.values
        self.feat_meta['internal'] = np.abs(self.stoich).sum(axis=0) > 1
        self.feat_meta['nonlin_rxns'] = (self.feat_meta['internal']) & (~np.in1d(self.feat_meta.index.tolist(), lin_rxns))
        self.n_linear = sum((~self.feat_meta['nonlin_rxns']) & self.feat_meta['internal'])

        self.feat_meta['Enzymes'] = self.feat_meta['Enzymes'].astype('category')
        nonlin_stoich = self.stoich[:, self.feat_meta['nonlin_rxns']]
        stoich_prod = (np.abs(nonlin_stoich)-nonlin_stoich)/2
        self.met_meta['beta_mets'] = np.sum(stoich_prod, axis=1) > 0
        self.met_meta['substrate_glycoform'] = np.sum(self.stoich[:, ~self.feat_meta['internal']] > 0, axis=1) > 0
        self.stoich_prod = pd.DataFrame(stoich_prod[self.met_meta['beta_mets'], :],
                                        index=self.met_meta.index[self.met_meta['beta_mets']],
                                        columns=self.feat_meta.index[self.feat_meta['nonlin_rxns']])

        self.enzymes = self.feat_meta.Enzymes[self.feat_meta.nonlin_rxns].cat.remove_unused_categories()
        self.producers = pd.Series({j: self.stoich_prod.index[self.stoich_prod[j] > 0].values[0] for j in self.feat_meta.index[self.feat_meta['nonlin_rxns']]})
        self.substrate = self.met_meta.index[self.met_meta['substrate_glycoform']].tolist()
        self.stoich_internal = pd.DataFrame(self.stoich[:, self.feat_meta['internal']],
                                            index=self.met_meta.index,
                                            columns=self.feat_meta.index[self.feat_meta['internal']])

    def _read_core(self, data, time_data):
        """
        Reads core data including titer and VCD.

        :param data: Data containing titer and VCD values.
        :type data: pd.DataFrame
        :param time_data: Time data for the experiments.
        :type time_data: pd.Series
        """
        self.time_col = time_data
        self.titer = data['Titer'].sort_index(axis='index')
        self.vcd = data['VCD'].sort_index(axis='index')

    def run_multistart(self, instance, strategy='rand', iterations=100, suppress_warning=True, solver_options=None):
        """
        Runs a multi-start optimization for the given instance.

        :param instance: The Pyomo model instance to be solved.
        :type instance: pyo.ConcreteModel
        :param strategy: The strategy for multi-start optimization (default is 'rand').
        :type strategy: str, optional
        :param iterations: The number of iterations to run (default is 100).
        :type iterations: int, optional
        :param suppress_warning: Flag to suppress warnings (default is True).
        :type suppress_warning: bool, optional
        :param solver_options: Additional solver options to configure.
        :type solver_options: dict[str, any], optional
        :return: A tuple containing the results dictionary and the status run dictionary.
        :rtype: tuple[dict[int, any], dict[str, list[int]]]
        """
        # Should be moved to generic parse_model?
        config = namedtuple("CONFIG", ["strategy", "iterations", "suppress_unbounded_warning"])
        CONFIG = config(strategy=strategy, iterations=iterations, suppress_unbounded_warning=suppress_warning)

        default_solver_options = {'max_iter': 100000,
                                  'linear_solver': 'mumps',
                                  'warm_start_init_point': 'yes',
                                  'print_level': 4}

        if solver_options is not None:
            default_solver_options.update(solver_options)

        solver = pyo.SolverFactory('ipopt')
        solver.options = default_solver_options
        results_dict = {}
        status_run = {'Optimal': [], 'MaxIterations': [], 'Error': []}
        for i in (pbar := tqdm(np.arange(CONFIG.iterations), total=CONFIG.iterations, desc='Progress',
                               bar_format='{desc}: {percentage:3.0f}%|{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]')):
            print('------------------------------------------------------------------------------')
            if i > 0:
                reinitialize_variables(instance, CONFIG)
                initial_vars = pd.DataFrame({v.name: v.value for v in instance.component_data_objects(pyo.Var, descend_into=True)},
                                            index=['variable name']).T
            try:
                status = solver.solve(instance, tee=True)
                if (status.solver.status == SolverStatus.ok) and (status.solver.termination_condition == TerminationCondition.optimal):
                    status_run['Optimal'].append(i)
                elif (status.solver.status == SolverStatus.warning) and (status.solver.termination_condition == TerminationCondition.maxIterations):
                    status_run['MaxIterations'].append(i)
                print(status)
                result = self.get_results(instance)
                if i > 0:
                    result['initial_vars'] = initial_vars.copy()

                result['objectives']['message'] = status['Solver'].message
                result['status'] = status.copy()
                results_dict[i] = result.copy()
                print(i, status.solver.status, status.solver.termination_condition)
                pbar.set_description(f"Progress ({', '.join(f'{k}: {len(v)}' for k, v in status_run.items())})", refresh=True)
            except:
                status_run['Error'].append(i)
                pbar.set_description(f"Progress ({', '.join(f'{k}: {len(v)}' for k, v in status_run.items())})", refresh=True)
                continue

        # solver = pyo.SolverFactory('multistart')
        # status = solver.solve(instance, solver='ipopt',
        #                       solver_args={'tee': True,
        #                                    'options': default_solver_options},
        #                       strategy=CONFIG.strategy, iterations=CONFIG.iterations, suppress_unbounded_warning=CONFIG.suppress_unbounded_warning)
        # result = self.get_results(instance)
        # result['status'] = status.copy()

        return results_dict, status_run

    def run_singlestart(self, instance, solver_options=None):
        """
        Runs a single-start optimization for the given instance.

        :param instance: The Pyomo model instance to be solved.
        :type instance: pyo.ConcreteModel
        :param solver_options: Additional solver options to configure.
        :type solver_options: dict[str, any], optional
        :return: A dictionary containing the results of the optimization.
        :rtype: dict[str, any]
        """
        solver = pyo.SolverFactory('ipopt')
        default_solver_options = {'max_iter': 100000,
                                  'linear_solver': 'mumps',
                                  'warm_start_init_point': 'yes',
                                  'print_level': 4}  # ,
        # 'tol': 1e-4,
        # 'constr_viol_tol': 1e-4}
        if solver_options is not None:
            default_solver_options.update(solver_options)

        solver.options = default_solver_options
        status = solver.solve(instance, tee=True)
        result = self.get_results(instance)
        result['status'] = status.copy()

        return result


def summarize_instance(instance):
    """
    Summarizes the instance by printing the statistics, variable names, constraint names, and objective names.

    :param instance: The Pyomo model instance to summarize.
    :type instance: pyo.ConcreteModel
    """
    instance.compute_statistics()
    print(instance.statistics)
    print('\n')
    print('Names of variables in model:')
    for v in instance.component_objects(pyo.Var, active=True):
        print(f'{v.name} : {len(list(v.keys()))} variables')

    print('\n')
    print('Names of constraints in model:')
    for v in instance.component_objects(pyo.Constraint, active=True):
        print(f'{v.name} : {len(list(v.keys()))} constraints')

    print('\n')
    print('Names of objectives in model:')
    for v in instance.component_objects(pyo.Objective, active=True):
        print(f'{v.name} : {len(list(v.keys()))} objective')
    print('\n')

    print(instance.timepoints.pprint())
    print(instance.mets.pprint())
    print(instance.enzymes.pprint())
    print(instance.internal_rxns.pprint())