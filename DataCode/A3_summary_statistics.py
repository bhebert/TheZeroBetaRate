######################################################
# Summary Statistics
# Produce summary statistics for 
# 1. Instruments
# 2. Equity portfolio constructed
# 
# Produce table in the appendix and useful for debugging.
######################################################



# %%
from os.path import join
import os
import numpy as np
import pandas as pd
from pandas.tseries.offsets import *
from functools import reduce
from portfolio_construction import generate_groups



with open('path_variables.md', 'r') as f:
    main_path = f.readline().split('"')[1]
    
    
raw_path = join(main_path, "Raw Data")
process_path = join(main_path, "Processed Data")
input_path = join(main_path, "Input")
output_path = join(main_path, "Output", "Data Output")
if not os.path.exists(output_path):
    os.mkdir(output_path)


# %%
######################################################
# Other Instruments
######################################################
instruments = pd.read_csv(join(input_path, "Instruments.csv")).query("Date >= '1973-02-01' and Date <= '2020-11-30'")
summary_vars = ['RF', 'UMP', 'SAHM', 'EBP', 'CAPE', 'TSP', 'shadow_spread', 'BAAS']

instruments_summary = (
    instruments[summary_vars]
    .describe()
    .set_axis([r'T-bill Yield (\%, Monthly)', r'Unemployment (\%)', r'Sahm Rate (\%)', r'Excess Bond Premium (\%, Annual)', r'CAPE', r'Term Spread (\%, Annual)', r'Shadow Spread (\%, Monthly)', r'Corporate Bond Spread (\%, Annual)'], axis=1)
    .set_axis(['Count', 'Mean', 'SD', 'Min', r'25\%', r'50\%', r'75\%', 'Max'], axis=0)
    )

######################################################
# Rolling Inflation 
######################################################
inflation_summary = (
    instruments.query("Date >= '1973-01-01' and Date <= '2020-10-31'")
    .loc[:, ['CPI_rolling']]
    .describe()
    .set_axis(['Count', 'Mean', 'SD', 'Min', r'25\%', r'50\%', r'75\%', 'Max'], axis=0)
)



######################################################
# Consumption
######################################################
jw_m = pd.read_csv(join(input_path, "Consumption", "jw_m.csv"))
jw_m = (
    jw_m
    .assign(log_cg = lambda x: 100*np.log(1 + x['cg_nd_sv_r_pc']/100) )
    .query("date >= '1973-03-31' and date <= '2020-12-31'")
)

cons_summary = (
    jw_m['log_cg']
    .describe()
    .set_axis(['Count', 'Mean', 'SD', 'Min', r'25\%', r'50\%', r'75\%', 'Max'], axis=0)
    )


summary_table = (
    pd.concat([instruments_summary, inflation_summary, cons_summary], axis=1)
    .rename(columns={'log_cg':r'Consumption Growth (\%, Monthly)', 'CPI_rolling':r'Rolling Inflation (\%, Monthly)'})
    .iloc[1:].map('{:.3f}'.format)
    .T
    )


summary_table.style.to_latex(join(output_path, 'summary_table.tex'), hrules=True)






######################################################
# Other Instruments
######################################################
instruments = pd.read_csv(join(input_path, "Instruments.csv")).query("Date >= '1973-02-01' and Date <= '2019-11-30'")
summary_vars = ['RF', 'UMP', 'SAHM', 'EBP', 'CAPE', 'TSP', 'shadow_spread', 'BAAS']

instruments_summary = (
    instruments[summary_vars]
    .describe()
    .set_axis([r'T-bill Yield (\%, Monthly)', r'Unemployment (\%)', r'Sahm Rate (\%)', r'Excess Bond Premium (\%, Annual)', r'CAPE', r'Term Spread (\%, Annual)', r'Shadow Spread (\%, Monthly)', r'Corporate Bond Spread (\%, Annual)'], axis=1)
    .set_axis(['Count', 'Mean', 'SD', 'Min', r'25\%', r'50\%', r'75\%', 'Max'], axis=0)
    )

######################################################
# Rolling Inflation 
######################################################
inflation_summary = (
    instruments.query("Date >= '1973-01-01' and Date <= '2019-10-31'")
    .loc[:, ['CPI_rolling']]
    .describe()
    .set_axis(['Count', 'Mean', 'SD', 'Min', r'25\%', r'50\%', r'75\%', 'Max'], axis=0)
)



######################################################
# Consumption
######################################################
jw_m = pd.read_csv(join(input_path, "Consumption", "jw_m.csv"))
jw_m = (
    jw_m
    .assign(log_cg = lambda x: 100*np.log(1 + x['cg_nd_sv_r_pc']/100) )
    .query("date >= '1973-03-31' and date <= '2019-12-31'")
)

cons_summary = (
    jw_m['log_cg']
    .describe()
    .set_axis(['Count', 'Mean', 'SD', 'Min', r'25\%', r'50\%', r'75\%', 'Max'], axis=0)
    )


summary_table = (
    pd.concat([instruments_summary, inflation_summary, cons_summary], axis=1)
    .rename(columns={'log_cg':r'Consumption Growth (\%, Monthly)', 'CPI_rolling':r'Rolling Inflation (\%, Monthly)'})
    .iloc[1:].map('{:.3f}'.format)
    .T
    )


summary_table.style.to_latex(join(output_path, 'summary_table_precovid.tex'), hrules=True)






# %%
######################################################
# Test Assets
######################################################
# Hyperparameters (should mirror the settings in 00_create_data)

segment = [0, 0.3, 0.7, 1]
beta_segment = [0, 0.33, 0.66, 1]
drop_20 = True
share_code_restriction = False 

# Generation
ccm2 = generate_groups(main_path, segment, beta_segment, drop_20, share_code_restriction)
ccm2 = ccm2.query("date >= '1973-01-01'")


# %%
######################################################
# Summary Statistics
###################################################### 

def format_result(group_data, ports):
    
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



def portfolio_tables(ports):

    count = (
        ccm2[(ccm2[ports]!=0).all(1)]
        .groupby(["date"] +ports)
        .size()
        .reset_index()
        .rename(columns={0:'N'})
        )
    
    beta_equal = (
        ccm2[(ccm2[ports]!=0).all(1)]
        .groupby(["date"] +ports)
        .agg({'beta':'mean'})
        .reset_index()
        .rename(columns={'beta':'Beta'})
        )
    
    beta_weighted = (
        ccm2[(ccm2[ports]!=0).all(1)]
        .groupby(["date"] +ports)
        .apply(wavg, "beta", "wt", include_groups=False)
        .reset_index()
        .rename(columns={0:'Beta'})
        )
    
    return_equal = (
        ccm2[(ccm2[ports]!=0).all(1)]
        .groupby(["date"] +ports)
        .agg({'retadj':'mean'})
        .reset_index()
        .rename(columns={'retadj':'Return'})
        .assign(Return = lambda x: x['Return'] * 100)
        )
    
    return format_result(count, ports), format_result(beta_equal,ports), format_result(beta_weighted,ports), format_result(return_equal,ports)
        
        
    



# 3x3x3 Beta x Size x Book-Market portfolio
ports = ['beta_port_me_bm', 'me_port', 'bm_port']
count_27, beta_equal_27, beta_weighted_27, return_equal_27 = portfolio_tables(ports)

count_27.to_csv(join(output_path, "count_27.csv"), index=False)
beta_equal_27.to_csv(join(output_path, "beta_equal_27.csv"), index=False)
beta_weighted_27.to_csv(join(output_path, "beta_weighted_27.csv"), index=False)
return_equal_27.to_csv(join(output_path, "return_equal_27.csv"), index=False)

# 3x3x3 x 3 FF5 portfolio 
ports_list = [['beta_port_me_op', 'me_port', 'op_port'], ['beta_port_me_inv', 'me_port', 'inv_port'], ['bm_port', 'me_port', 'beta_port_me_bm']]
counts = []
beta_equals = []
beta_weighteds = []
return_equals = []
for ports in ports_list:
    counts.append(portfolio_tables(ports)[0])
    beta_equals.append(portfolio_tables(ports)[1])
    beta_weighteds.append(portfolio_tables(ports)[2])
    return_equals.append(portfolio_tables(ports)[3])


count_ff5 = reduce(
    lambda left, right: pd.merge(left, right, how = 'inner', on = ['date']),
    counts
)

beta_equal_ff5 = reduce(
    lambda left, right: pd.merge(left, right, how = 'inner', on = ['date']),
    beta_equals
)

beta_weighted_ff5 = reduce(
    lambda left, right: pd.merge(left, right, how = 'inner', on = ['date']),
    beta_weighteds
)

return_equal_ff5 = reduce(
    lambda left, right: pd.merge(left, right, how = 'inner', on = ['date']),
    return_equals
)


count_ff5.to_csv(join(output_path, "count_ff5.csv"), index=False)
beta_equal_ff5.to_csv(join(output_path, "beta_equal_ff5.csv"), index=False)
beta_weighted_ff5.to_csv(join(output_path, "beta_weighted_ff5.csv"), index=False)
return_equal_ff5.to_csv(join(output_path, "return_equal_ff5.csv"), index=False)

# %%
