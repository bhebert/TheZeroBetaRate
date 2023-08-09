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
summary_vars = ['RF', 'UMP', 'EBP', 'CAPE', 'TSP', 'shadow_spread', 'BAAS']

instruments_summary = (
    instruments[summary_vars]
    .describe()
    .set_axis(['T-bill Yield (\%, Monthly)', 'Unemployment (\%)', 'Excess Bond Premium (\%, Monthly)', 'CAPE', 'Term Spread (\%, Annual)', 'Shadow Spread (\%, Monthly)', 'Corporate Bond Spread (\%, Annual)'], axis=1)
    .set_axis(['Count', 'Mean', 'SD', 'Min', '25\%', '50\%', '75\%', 'Max'], axis=0)
    )

######################################################
# Rolling Inflation 
######################################################
inflation_summary = (
    instruments.query("Date >= '1973-01-01' and Date <= '2020-10-31'")
    .loc[:, ['CPI_rolling']]
    .describe()
    .set_axis(['Count', 'Mean', 'SD', 'Min', '25\%', '50\%', '75\%', 'Max'], axis=0)
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
    .set_axis(['Count', 'Mean', 'SD', 'Min', '25\%', '50\%', '75\%', 'Max'], axis=0)
    )


summary_table = (
    pd.concat([instruments_summary, inflation_summary, cons_summary], axis=1)
    .rename(columns={'log_cg':'Consumption Growth (\%, Monthly)', 'CPI_rolling':'Rolling Inflation (\%, Monthly)'})
    .iloc[1:].applymap('{:.3f}'.format)
    .T
    )


summary_table.style.to_latex(join(output_path, 'summary_table.tex'))








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


# 3x3x3 Beta x Size x Book-Market portfolio
ports = ['beta_port_me_bm', 'me_port', 'bm_port']
count_27 = (
    ccm2[(ccm2[ports]!=0).all(1)]
    .groupby(["date"] +ports)
    .size()
    .reset_index()
    .rename(columns={0:'N'})
    )

count_27 = format_result(count_27, ports)




# 3x3x3 x 3 FF5 portfolio 
ports_list = [['beta_port_me_op', 'me_port', 'op_port'], ['beta_port_me_inv', 'me_port', 'inv_port'], ['bm_port', 'me_port', 'beta_port_me_bm']]
counts = []
for ports in ports_list:
    count = (
    ccm2[(ccm2[ports]!=0).all(1)]
    .groupby(["date"] +ports)
    .size()
    .reset_index()
    .rename(columns={0:'N'})
    )
    counts.append(format_result(count, ports))
    
count_ff5 = reduce(
    lambda left, right: pd.merge(left, right, how = 'inner', on = ['date']),
    counts
)


# %%
