'''
This program creates the testing assets
'''


#%%
from os.path import join
import numpy as np
import pandas as pd
from pandas.tseries.offsets import *
from functools import reduce


######################################################
# Hyperparameters
######################################################



def generate_groups(main_path, segment, beta_segment, drop_20, share_code_restriction):
    raw_path = join(main_path, "Raw Data")
    process_path = join(main_path, "Processed Data")



    ######################################################
    # Generation
    ######################################################

    betas = pd.read_csv(
        join(raw_path, "beta_estimated.csv"), engine="pyarrow", index_col=0
    ).assign(date=lambda x: pd.to_datetime(x["date"]) + MonthEnd(0))





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


    # %%

    ccm_june = ccm[ccm["date"].dt.month == 6]
    ccm_june = ccm_june.merge(betas, on=["permno", "permco", "date"], how="inner")
    
    if drop_20 is True:
        # Get rid of the 20% smallest stocks
        size_cutoff = ccm_june.groupby(['date'])['me'].quantile(0.2).reset_index().rename(columns={"me":"size_cutoff"})
        ccm_june = pd.merge(ccm_june, size_cutoff, on='date', how='inner')
        ccm_june = ccm_june.query("me >= size_cutoff")

    ccm_june2 = ccm_june.copy()
    for var in ['bm', 'me', 'op', 'inv', 'beta']:
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

    
    # 
    for var in ["bm", "me", "beta", "op", 'inv']:
        conds = []
        for cut in [30, 70]:
            new_var = var + "_p" + str(cut)
            conds.append(ccm_june2[var] <= ccm_june2[new_var])

        vals = [1, 2]
        port_var = var + "_port"
        ccm_june2[port_var] = np.select(conds, vals)
        
        ccm_june2.loc[ccm_june2[var] > ccm_june2[new_var], port_var] = 3
    
    # A separate section to compute the beta groups by size groups.
    # This is because the small stocks may have inaccurate betas due to
    # liquidity reasons.
    for group_var in [
        ['me_port', 'op_port'], ['me_port', 'inv_port'], ['me_port', 'bm_port']
    ]:  
        two_var = "_".join([str(item).replace("_port", "") for item in group_var])
        
        if share_code_restriction is True:
            beta_breaks = (
                ccm_june2
                .query("bm > 0 and me > 0 and (shrcd == 10 or shrcd == 11)")
                .replace([np.inf, -np.inf], np.nan)
                .dropna(how='all')
                .groupby(['date']+group_var)['beta']
                .quantile(beta_segment)
                .reset_index()
                .rename(columns={'level_3': 'group'})
                )
        else:
            beta_breaks = (
                ccm_june2
                .query("bm > 0 and me > 0")
                .replace([np.inf, -np.inf], np.nan)
                .dropna(how='all')
                .groupby(['date']+group_var)['beta']
                .quantile(beta_segment)
                .reset_index()
                .rename(columns={'level_3': 'group'})
                )

        group_num = 1
        beta_breaks['group2'] = 0
        for i in beta_segment:
            beta_breaks['group2'] = np.where(beta_breaks['group']==i, group_num, beta_breaks['group2'])
            group_num = group_num + 1

        beta_breaks['group'] = beta_breaks['group2'].astype(int)
        beta_breaks = beta_breaks.drop('group2', axis=1)
        beta_breaks = beta_breaks.pivot(index=['date']+group_var, columns='group', values='beta').reset_index()

        names = ['date']+group_var
        for ix, i in enumerate(beta_segment):
            names.append("{0}_p{1}".format('beta_' + two_var, int(beta_segment[ix]*100) ))
        beta_breaks = beta_breaks.set_axis(names, axis = 1)

        ccm_june2 = pd.merge(ccm_june2, beta_breaks, how = 'left', on = ['date']+group_var)
        conds = []
        for cut in [x * 100 for x in beta_segment[1:-1]]:
            new_var = 'beta_' + two_var + "_p" + str(int(cut))
            conds.append(ccm_june2['beta'] <= ccm_june2[new_var])

        vals = list(range(1, len(segment)-1))
        port_var = 'beta' + "_port_" + two_var
        ccm_june2[port_var] = np.select(conds, vals)      
        
        ccm_june2.loc[ccm_june2['beta'] > ccm_june2[new_var], port_var] = len(segment)-1
        
        # An useful object for debugging.
        # beta_counts = (
        #     ccm_june2
        #     .groupby(['date']+ group_var + [port_var])['beta']
        #     .count()
        # )


            
            


    # %%
    
    port_list = []
    for col in ccm_june2.columns:
        if "_port" in col:
            port_list.append(col)
            
                
    ccm_june2 = (
        ccm_june2.pipe(pd.DataFrame.drop_duplicates, subset=["permco", "permno", "date"])
        .assign(year=lambda x: x["date"].dt.year + 1)
        .loc[:, ["permco", "permno", "year"]+port_list]
    )


    ccm["year"] = ccm["date"].dt.year
    ccm2 = pd.merge(ccm, ccm_june2, on=["permno", "permco", "year"], how="inner")
    
    return ccm2


    
def generate_portfolio(main_path, segment, beta_segment, drop_20, share_code_restriction):
    """generate portfolio according to Fama French 1993

    Args:
        segment (list): A list of number between 0 and 1 to indicate the quantiles used to divide stocks into groups.
        drop_20 (boolean): Whether to drop the smallest 20% ME stocks in the data.
        
    Returns:
        DataFrame: Two DataFrames
    """    
    
    raw_path = join(main_path, "Raw Data")
    
    ccm2 = generate_groups(main_path, segment, beta_segment, drop_20, share_code_restriction)
    
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


    def generate_portfolio(dataframe, ports):
        group_data = (
            dataframe[(dataframe[ports]!=0).all(1)]
            .groupby(["date"] +ports)
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


    # FF5 portfolios with size-dependent beta groups
    beta_size_op = generate_portfolio(ccm2, ['beta_port_me_op', 'me_port', 'op_port'])
    beta_size_inv = generate_portfolio(ccm2, ['beta_port_me_inv', 'me_port', 'inv_port'])
    bm_size_beta = generate_portfolio(ccm2, ['bm_port', 'me_port', 'beta_port_me_bm'])

    # only 27 portfolio from FF3 with size-dependent beta groups
    beta_size_bm = generate_portfolio(ccm2, ['beta_port_me_bm', 'me_port', 'bm_port'])


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

if __name__ == "__main__":
    with open('path_variables.md', 'r') as f:
        main_path = f.readline().split('"')[1]
    generate_portfolio(main_path, [0, 0.3, 0.7, 1], [0, 0.33, 0.66, 1], True)