#!/usr/bin/env python
# coding: utf-8

'''
This file estimated the betas on individual stocks using time series data on stock returns 
and market returns.

CAPM model:
return_{it} = alpha + b Market Return_{it} + e_{it}

i - stock
t - monthly data

This procedure replicates Fama French procedure of beta estimation. 
β for June of year t is estimated using the preceding five years (two minimum) of past monthly returns.

by Qitong Wang, 03.18.2022
'''

# # Beta Estimation

# This file estimates the betas for firms in the data.



# %%
import pandas as pd
import numpy as np
import datetime as dt
import os

# to handle dates
from dateutil.relativedelta import *
from pandas.tseries.offsets import *

# for estimation
import statsmodels.formula.api as smf

# better loop
from joblib import Parallel, delayed

# progress tracking
import contextlib
import joblib
from tqdm import tqdm

@contextlib.contextmanager
def tqdm_joblib(tqdm_object):
    """Context manager to patch joblib to report into tqdm progress bar given as argument"""
    class TqdmBatchCompletionCallback(joblib.parallel.BatchCompletionCallBack):
        def __call__(self, *args, **kwargs):
            tqdm_object.update(n=self.batch_size)
            return super().__call__(*args, **kwargs)

    old_batch_callback = joblib.parallel.BatchCompletionCallBack
    joblib.parallel.BatchCompletionCallBack = TqdmBatchCompletionCallback
    try:
        yield tqdm_object
    finally:
        joblib.parallel.BatchCompletionCallBack = old_batch_callback
        tqdm_object.close()


with open('path_variables.md', 'r') as f:
    main_path = f.readline().split('"')[1]

raw_path = os.path.join(main_path, "Raw Data")

# keep track of the time
print("The program starts at {0}".format(dt.datetime.now()))


# ## 2. Merge with CRSP dataset
# Load the monthly CRSP dataset
crsp_m = (
    pd.read_csv(os.path.join(raw_path, "CCM", "crsp_m.csv"), index_col=0, engine='pyarrow')
    .assign(date = lambda x: pd.to_datetime(x['date']) + MonthEnd(0))
    .drop_duplicates(subset=['permno', 'permco', 'date'])
    )

ff = (
    pd.read_csv(os.path.join(raw_path, "Factors", "F-F_Research_Data_5_Factors_2x3.csv"), skiprows=3)
    .set_axis(['date', 'mktrf', 'smb', 'hml', 'rmw', 'cma', 'rf'], axis=1)
    .assign(date = lambda x: pd.to_datetime(x['date'].astype(str)+'01') + MonthEnd(0),
            mktrf = lambda x: x['mktrf']/100,
            rf = lambda x:x['rf']/100)
)
# merge with the monthly index data
crsp = pd.merge(crsp_m, ff, how='inner', on=['date'])

# exclude extreme values
# crsp = crsp[(crsp['ret'] > -0.5) & (crsp['ret'] < 1.5)]

# compute excess return
crsp['eret'] = crsp['ret'] - crsp['rf']

# %%
# ## 3. Estimation

# Estimate monthly beta for individual stock using monthly returns. OLS regress stock return on a constant and the contemporaneous return on the market portfolio return. 
# - The estimation uses one to three years of return data. 
# - Skip the firm missing values, but require the included value to be in a four-year window.
# - Censor the stock return to the range (-50%, 100%) to limit the influence of extreme firm-specific outliers.



# betas = pd.DataFrame(columns=['permno', 'date', 'beta'])

# select list of permno code to experiment
codes = crsp['permno'].unique()

def beta_estimation(permno):
    
    firm_returns = crsp.loc[crsp['permno'] == permno, ['permno', 'permco', 'date', 'eret', 'mktrf']].dropna()
    
    firm_returns = firm_returns.sort_values(by=['date'])
    
    outputs = pd.DataFrame(columns=['permno', 'permco', 'date', 'beta'])

    for dt in firm_returns['date']:

        # this is the date currently looking at
        current_dt = dt
        
        # this is the firm code
        firm_code = firm_returns.loc[firm_returns['date']==dt].iat[0,1]

        # set the earliest date to include data (5Y)
        min_dt = current_dt + MonthEnd(-60)

        est_df = firm_returns.loc[(firm_returns['date'] <= current_dt) & (firm_returns['date']>=min_dt)].sort_values(by=['date']).reset_index()

        # need at least 24 months of obs to do estimation
        if len(est_df) < 24:
            #betas = betas.append({'permno':permno, 'date':current_dt, 'beta':np.nan}, ignore_index=True)
            
            if outputs.empty:
                outputs = pd.DataFrame.from_records([{'permno':permno, 'permco':firm_code, 'date':current_dt, 'beta':np.nan}])
            else:
                outputs = pd.concat([outputs, pd.DataFrame.from_records([{'permno':permno, 'permco':firm_code, 'date':current_dt, 'beta':np.nan}])], ignore_index=True)
            continue
        
        # One factor model        
        
        result = smf.ols(formula='eret ~ mktrf', data = est_df).fit()
        beta = result.params['mktrf']
    
        if outputs.empty:
            outputs = pd.DataFrame.from_records([{'permno':permno, 'permco':firm_code, 'date':current_dt, 'beta':beta}])
        else:
            outputs = pd.concat([outputs, pd.DataFrame.from_records([{'permno':permno, 'permco':firm_code, 'date':current_dt, 'beta':beta}])], ignore_index=True)
        
        
    return outputs

        
        
# %%


# run results in parallel
with tqdm_joblib(tqdm(desc="Beta Estimation", total=len(codes))) as progress_bar:
    results = Parallel(n_jobs=8, verbose=0)([delayed(beta_estimation)(code) for code in codes])

# unpack the results
results = [result for result in results if not result.empty]
output = pd.concat(results, ignore_index=True)
    
# output = output.rename(columns={0:'permno', 1:'date', 2:'beta'})

# %%
output.convert_dtypes().to_csv(os.path.join(raw_path, "beta_estimated.csv"))


# keep track of the time
print("The program ends at {0}".format(dt.datetime.now()))