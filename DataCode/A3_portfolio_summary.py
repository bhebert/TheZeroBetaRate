######################################################
# Equity Portfolio Summary
# Produce summary statistics for equity portfolio constructed
# in portfolio_construction.py
# Produce table in the appendix and useful for debugging.
######################################################



# %%
from os.path import join
import numpy as np
import pandas as pd
from pandas.tseries.offsets import *
from functools import reduce
from portfolio_construction import generate_groups



with open('path_variables.md', 'r') as f:
    main_path = f.readline().split('"')[1]
    
    
raw_path = join(main_path, "Raw Data")
process_path = join(main_path, "Processed Data")

######################################################
# Hyperparameters (should mirror the settings in 00_create_data)
######################################################

segment = [0, 0.3, 0.7, 1]
beta_segment = [0, 0.33, 0.66, 1]
drop_20 = True
share_code_restriction = False 


######################################################
# Generation
######################################################

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
