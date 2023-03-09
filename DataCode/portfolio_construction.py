'''
This program creates the testing assets
'''


#%%
from os.path import join, expanduser
import numpy as np
import pandas as pd
from pandas.tseries.offsets import *
from functools import reduce


main_path = expanduser("~/Dropbox/Projects/Implied Interest Rates/Estimation")
raw_path = join(main_path, "Raw Data")
process_path = join(main_path, "Processed Data")
out_path = join(main_path, "Input")



def generate_portfolio(segment, drop_20):
    """generate portfolio according to Fama French 1993

    Args:
        segment (list): A list of number between 0 and 1 to indicate the quantiles used to divide stocks into groups.
        drop_20 (boolean): Whether to drop the smallest 20% ME stocks in the data.

    Returns:
        DataFrame: Two DataFrames
    """

    ######################################################
    # Generation
    ######################################################

    betas = pd.read_csv(
        join(raw_path, "beta_estimated.csv"), engine="pyarrow", index_col=0
    ).assign(date=lambda x: pd.to_datetime(x["date"]) + MonthEnd(0))

    crsp = (
        pd.read_csv(join(raw_path, "CCM", "crsp_m.csv"), index_col=0, engine='pyarrow')
        .assign(date = lambda x: pd.to_datetime(x['date']) + MonthEnd(0))
    )

    mseall = (
        pd.read_csv(join(raw_path, "CCM", "mseall.csv"), index_col=0, engine='pyarrow')
        .assign(date = lambda x: pd.to_datetime(x['date']) + MonthEnd(0))
        .drop_duplicates(subset=["date", "permno", "permco"])
    )

    merged = reduce(
        lambda left, right: pd.merge(
            left, right, on=["permno", "permco", "date"], how="inner"
        ),
        [crsp, betas, mseall[["permno", "permco", "date", "exchcd", "shrcd"]]],
    )

    merged = merged.query(
        "((exchcd == 1 and shrcd == 10) or (shrcd == 11)) and (date.dt.month == 6)"
    )

    beta_breakpoints = merged.groupby("date")["beta"].describe(
        percentiles=np.arange(5, 105, 5) / 100
    )
    beta_breakpoints = beta_breakpoints.iloc[:, 4:-1]
    column_names = []
    for i in range(20):
        column_names.append("beta_p{0}".format((i + 1) * 5))
    beta_breakpoints.columns = column_names
    beta_breakpoints = beta_breakpoints.reset_index()





    # %%
    '''
    ######################################################
    # CRSP Compustat Merged
    ######################################################
    '''
    ccm = (
        pd.read_csv(join(process_path, "ccm.csv"), engine='pyarrow')
        .assign(
            ddate=lambda x: pd.to_datetime(x["ddate"]) + MonthEnd(0),
            me=lambda x: x["me"] / 1000
        )
        .rename(columns={"ddate": "date"})
    )



    if drop_20 is True:
        # Get rid of the 20% smallest stocks
        size_cutoff = ccm.groupby(['date'])['me'].quantile(0.2).reset_index().rename(columns={"me":"size_cutoff"})
        ccm = pd.merge(ccm, size_cutoff, on='date', how='inner')
        ccm = ccm.query("me >= size_cutoff")


    # %%

    ccm_june = ccm[ccm["date"].dt.month == 6]

    ccm_june2 = ccm_june.copy()
    for var in ['bm', 'me', 'op', 'inv']:
        breaks = (
            ccm_june
            .query("exchcd == 1 and bm > 0 and me > 0 and (shrcd == 10 or shrcd == 11)")
            .replace([np.inf, -np.inf], np.nan)
            .dropna(how='all')
            .groupby('date')[var]
            .quantile(segment)
            .reset_index()
            .rename(columns={'level_1': 'group'})
            )
        group_num = 1
        breaks['group2'] = 0
        for i in segment:
            breaks['group2'] = np.where(breaks['group']==i, group_num, breaks['group2'])
            group_num = group_num + 1
        
        breaks['group'] = breaks['group2'].astype(int)
        breaks = breaks.drop('group2', axis=1)
        breaks = breaks.pivot(index='date', columns='group', values=var).reset_index()
        names = ['date']
        for ix, i in enumerate(segment):
            names.append("{0}_p{1}".format(var, int(segment[ix]*100) ))
        breaks = breaks.set_axis(names, axis = 1)
        
        ccm_june2 = pd.merge(ccm_june2, breaks, how = 'inner', on = ['date'])

    ccm_june2 = pd.merge(ccm_june2, betas, on=["permno", "permco", "date"], how="inner")
    ccm_june2 = pd.merge(ccm_june2, beta_breakpoints, on="date", how="inner")
            

    for var in ["bm", "me", "beta", "op", 'inv']:
        conds = []
        for cut in [30, 70, 100]:
            new_var = var + "_p" + str(cut)
            conds.append(ccm_june2[var] <= ccm_june2[new_var])

        vals = [1, 2, 3]
        port_var = var + "_port"
        ccm_june2[port_var] = np.select(conds, vals)
            
            


    # %%
    ccm_june2 = (
        ccm_june2.pipe(pd.DataFrame.drop_duplicates, subset=["permco", "permno", "date"])
        .assign(year=lambda x: x["date"].dt.year + 1)
        .loc[:, ["permco", "permno", "year", "bm_port", "me_port", "beta_port", "op_port", "inv_port"]]
    )


    ccm["year"] = ccm["date"].dt.year
    ccm2 = pd.merge(ccm, ccm_june2, on=["permno", "permco", "year"], how="inner")


    # %%
    def wavg(group, avg_name, weight_name):
        """
        Compute weighted average by group.
        
        INPUT
        group - group instance produced by the pandas groupby method.
        avg_name - target variable to compute the statistics
        weight_name - variable to be used as weights
        
        OUTPUT
        weighted_mean - weighted mean by group.
        
        by Qitong Wang - 01.28.2022
        """
        d = group[avg_name]
        w = group[weight_name]
        try:
            return (d * w).sum() / w.sum()
        except ZeroDivisionError:
            return np.nan


    # %%
    ccm3 = ccm2.query("bm_port != 0 and me_port != 0 and beta_port != 0 and op_port != 0 and inv_port != 0")
    ccm4 = ccm2.query("bm_port != 0 and me_port != 0 and beta_port != 0")


    def generate_portfolio(dataframe, ports):
        group_data = (
            dataframe.groupby(["date"] +ports)
            .apply(wavg, "retadj", "wt")
            .reset_index()
            .rename(columns={0: "mret"})
            )
        group_data["mret"] = group_data["mret"] * 100
        
        group_data = group_data.pivot_table(index=ports, columns="date")
        new_index = []
        port_tags = [sub.replace("_port", "") for sub in ports]
        for i in group_data.index:
            new_var_name = ""
            for ix, tag in enumerate(port_tags):
                new_var_name = new_var_name + tag + "_" + str(i[ix])
            
            new_index.append(new_var_name)

        group_data = group_data.set_index(pd.Index(new_index))
        group_data.columns = group_data.columns.get_level_values(1)
        group_data.columns.name = None
        group_data = group_data.T.reset_index()
        group_data.rename(columns={"index": "date"}, inplace=True)
        
        return group_data


    beta_size_op = generate_portfolio(ccm3, ['beta_port', 'me_port', 'op_port'])
    beta_size_inv = generate_portfolio(ccm3, ['beta_port', 'me_port', 'inv_port'])
    bm_size_beta = generate_portfolio(ccm3, ['bm_port', 'me_port', 'beta_port'])

    # only 27 portfolio from FF3
    beta_size_bm = generate_portfolio(ccm4, ['beta_port', 'me_port', 'bm_port'])


    # %%
    industry_portfolios = pd.read_csv(
        join(raw_path, "49_Industry_Portfolios.CSV"), skiprows=11, nrows=1150
    )
    industry_portfolios.rename(
        columns={industry_portfolios.columns[0]: "date"}, inplace=True
    )
    industry_portfolios["date"] = pd.to_datetime(
        industry_portfolios["date"].astype(str) + "01"
    ) + MonthEnd(0)

    # %%

    output_portfolio = reduce(
        lambda left, right: pd.merge(left, right, how = 'inner', on = ['date']),
        [bm_size_beta, beta_size_op, beta_size_inv, industry_portfolios]
    )

    output_portfolio_27 = reduce(
        lambda left, right: pd.merge(left, right, how = 'inner', on = ['date']),
        [beta_size_bm, industry_portfolios]
    )
    
    return output_portfolio, output_portfolio_27

    
    

# %%
