'''
Created on Thu Sep 8

@author: Qitong Wang

Load in monetary shock datasets.

'''


# %% 
import pandas as pd
from os.path import join
from pandas.tseries.offsets import *

def generate_MS(main_path):

    data_path = join(main_path, "Raw Data", "Monetary_shocks")
    output_path = join(main_path, "Input")


    # %% Nakamura and Steinsson Monetary Policy Shock
    ns = (
        pd.read_excel(join(data_path, "NS.xlsx"), sheet_name='shocks')
        .assign(Date = lambda x: pd.to_datetime(x['fomc']) + MonthEnd(0))
        .drop(columns = ['fomc'])
        .set_index(['Date'])
        .asfreq(freq='1M', fill_value=0)
        .reset_index()
    )

    # %% Romer and Romer Shock
    rr = (
        pd.read_stata(join(data_path,"RR_monetary_shock_monthly.dta"))
        .assign(Date = lambda x: pd.to_datetime(x['date'] + MonthEnd(0))) 
        .drop(columns = ['date'])
        .set_index(['Date'])
        .asfreq(freq='1M', fill_value=0)
        .reset_index()
    )

    # %%

    ns.to_csv(join(output_path, "NS_shocks.csv"), index=False)
    rr.to_csv(join(output_path, "RR_shocks.csv"), index=False)