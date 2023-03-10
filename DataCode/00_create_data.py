#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 09:35:46 2022

@author: Qitong Wang


Edit: Load all data at once, them output separate files for testing assets, factors, instruments, 
with the same length (time span).

Instruments: 
1. RF (t-bill)
2. CPI (Inflation)
3. UMP (Unemployment)
4. (CAPE) (Shiller PE ratio)
5. TSP (Term Spread Yield)
6. BAAS (Corporate Bond Premium)
7. EBP (Excess Bond Premium)
8. CPI Rolling (12 month rolling average)

Factors:
1. Market Return
2. SMB
3. HML
4. Bond Return Term Spread

Portfolios:
1. 27 Beta x Size x Book-to-Market sorted
2. 49 Industries
3. FF (2015)

"""

################################################################
#                          Main Path
################################################################

main_path = "YOUR DIRECTORY"
username = "YOUR WRDS USERNAME" # to access the WRDS database

################################################################



# %%
import pandas as pd
from pandas.tseries.offsets import MonthEnd

import os
os.chdir(main_path)
import sys
sys.path.insert(1, 'Data Scripts')
input_path = os.path.join(main_path, "Input")
if not os.path.exists(input_path):
    os.mkdir(input_path)

# %%
'''
Hyper-parameters
'''

convert_to_real = False
if convert_to_real is True:
    print("Constructing portfolios with real returns.")
elif convert_to_real is False:
    print("Constructing portfolios with nominal returns.")



# %% 
################################################################
#                          Instruments
################################################################

CAPE = (
    pd.read_excel("Raw Data/Instruments/CAPE_Shiller.xls", sheet_name="Data", skiprows=7)
    .query("Date == Date")
    .assign(Year = lambda x: x['Date'].astype('int'), 
            Month = lambda x: ((x['Date'] - x['Date'].astype(int))*100).round().astype(int),
            Date = lambda x: pd.to_datetime(x['Year'].round().astype(str) +'-' + x['Month'].astype(str)) + MonthEnd()
            )
    .loc[:, ['Date', 'CAPE']]
            
    )


BAA = (
    pd.read_csv("Raw Data/Instruments/BAA_FRED.csv")
    .rename(columns={'DATE':'Date'})
    .assign(Date = lambda x: pd.to_datetime(x['Date']) + MonthEnd(0))
    )

AAA = (
    pd.read_csv("Raw Data/Instruments/AAA_FRED.csv")
    .rename(columns={'DATE':'Date'})
    .assign(Date = lambda x: pd.to_datetime(x['Date']) + MonthEnd(0))
)


EBP = (
    pd.read_csv("Raw Data/Instruments/ExcessBondPremium_FRED.csv")
    .loc[:, ['date', 'ebp']]
    .assign(Date = lambda x: pd.to_datetime(x['date']) + MonthEnd(0))
    .drop(columns=['date'])
    .rename(columns={'ebp':'EBP'})
    )

TY10 = (
    pd.read_csv("Raw Data/Instruments/GS10_FRED.csv")
    .set_axis(['Date', 'TY10'], axis = 1)
    .assign(Date = lambda x: pd.to_datetime(x['Date']) + MonthEnd(0))
    )

TY3 = (
    pd.read_csv("Raw Data/Instruments/TB3MS_FRED.csv")
    .set_axis(['Date', 'TY3'], axis = 1)
    .assign(Date = lambda x: pd.to_datetime(x['Date']) + MonthEnd(0))
    )


CPI = (
    pd.read_csv("Raw Data/Instruments/CPILFESL_monthly_change.csv")
    .set_axis(['Date', 'CPI'], axis=1)
    .assign(Date = lambda x: pd.to_datetime(x['Date']) + MonthEnd(0),
            CPI_rolling = lambda x: x['CPI'].rolling(12).mean())
    )


UMP = (
    pd.read_csv("Raw Data/Instruments/UNRATE.csv")
    .set_axis(['Date', 'UMP'], axis = 1)
    .assign(Date = lambda x: pd.to_datetime(x['Date']) + MonthEnd(0))
    )


from shadow_spread import generate_shadow_spread
generate_shadow_spread(main_path)

shadow_spread = (
    pd.read_csv("Processed Data/shadow_spread.csv", index_col=0)
    .set_axis(['Date', 'shadow_spread'], axis=1)
    .assign(Date = lambda x: pd.to_datetime(x['Date']) + MonthEnd(0))
)




RF = (
    pd.read_csv("Raw Data/Instruments/betaff3_1963.csv")
    .loc[:, ['ym', 'RF']]
    .assign(Date = lambda x: pd.to_datetime(x['ym'],  format='%Ym%m') + MonthEnd(0),
            RF = lambda x: x['RF'].shift(-1))
    .drop(columns=['ym'])
    )


# notice that TY10 is yearly yield, and RF is monthly yield.
TSP = (
    pd.merge(TY10, RF, on='Date', how='inner')
    .assign(TSP = lambda x: x['TY10'] - x['RF']*12)
    .loc[:, ['Date', 'TSP']]
    )   

to_merge = [CPI, UMP, EBP, CAPE, TSP, BAA, AAA, shadow_spread]

instruments = RF

for df in to_merge:
    instruments = pd.merge(instruments, df, on='Date', how='left')

instruments = (
    instruments
    .loc[:, ['Date', 'RF', 'CPI', 'UMP', 'EBP', 'CAPE', 'TSP', 'BAA', 'AAA', 'CPI_rolling', 'shadow_spread']]
    # .assign(BAAS = lambda x: x['BAA'] - x['TSP'] - x['RF']*12)
    .assign(BAAS = lambda x: x['BAA'] - x['AAA'])
    .drop(columns=['BAA', 'AAA'])
    .dropna(subset=['RF', 'CPI', 'UMP', 'EBP', 'CAPE', 'TSP', 'shadow_spread'])
    )


# %%
################################################################
#                          Factors
################################################################
ff3 = (
    pd.read_excel("Raw Data/Factors/ff3.xlsx")
    .set_axis(['Date', 'Mkt', 'SMB', 'HML', 'RF'], axis = 1)
    .assign(Date = lambda x: pd.to_datetime(x['Date'].astype(str) + '01') + MonthEnd(0),
            Mkt = lambda x: x['Mkt'] + x['RF'])
)

ff5 = (
    pd.read_csv("Raw Data/Factors/F-F_Research_Data_5_Factors_2x3.csv")
    .set_axis(['Date', 'Mkt', 'SMB', 'HML', 'RMW', 'CMA', 'RF'], axis = 1)
    .assign(Date = lambda x: pd.to_datetime(x['Date'].astype(str) + '01') + MonthEnd(0),
            Mkt = lambda x: x['Mkt'] + x['RF'])
)

bond_return = (
    pd.read_csv("Raw Data/Factors/CRSP bond returns.csv")
    .assign(date = lambda x: pd.to_datetime(x['MCALDT']) + MonthEnd(0)
            )
    )

bond_return_very_long = (
    bond_return.query("TTERMLBL == 'Fama Maturity Portfolios - > 120 Month'")
    .rename(columns={'TMEWRETD':'rate_very_long'})
    )

bond_return_long = (
    bond_return.query("TTERMLBL == 'Fama Maturity Portfolios - > 60 and <= 120 Month'")
    .rename(columns={'TMEWRETD':'rate_long'})
    )

bond_return_short = (
    bond_return.query("TTERMLBL == 'Fama Maturity Portfolios - <= 6 Month'")
    .rename(columns={'TMEWRETD':'rate_short'})
    )


bond_term_spread = (
    pd.merge(bond_return_long[['date', 'rate_long']], RF[['Date', 'RF']], left_on='date', right_on='Date', how='inner')
    .assign(RF_ret = lambda x: x['RF'].shift(1))
    .assign(term_spread = lambda x: (x['rate_long']*100 - x['RF_ret']) )
    .loc[:, ['date', 'term_spread']]
    )

ICE_15_plus = (
    pd.read_csv("Raw Data/Factors/BAMLCC8A015PYTRIV.csv")
    .set_axis(['Date', 'ICE_15'], axis = 1)
    .assign(ICE_15 = lambda x: pd.to_numeric(x['ICE_15'], errors='coerce'))
    .dropna()
    .assign(Date = lambda x: pd.to_datetime(x['Date']))
    .set_index('Date')
    .resample('M').last()
    .assign(ICE_15_ret = lambda x: x['ICE_15'].pct_change())
    .reset_index() 
    .assign(Date = lambda x: x['Date'] + MonthEnd(0)) 
    .dropna()
)

default_factor = (
    pd.merge(bond_return_very_long[['date', 'rate_very_long']],
             ICE_15_plus[['Date', 'ICE_15_ret']], left_on='date', right_on='Date', how='inner')
    .assign(DEF = lambda x: (x['ICE_15_ret'] - x['rate_very_long'])*100)
)

factors = (
    pd.merge(ff5, 
             bond_term_spread, 
             how='inner',
             left_on=['Date'],
             right_on=['date'])
    .merge(default_factor[['Date', 'DEF']], how='left', on='Date')
    .drop(columns = ['date',  'RF'])
    )

if convert_to_real is True:
    factors = (
        pd.merge(factors, instruments.loc[:, ['Date', 'CPI']], how = 'inner', on='Date')
        .assign(Mkt = lambda x: x['Mkt'] - x['CPI'])
        .drop(columns=['CPI'])
        .dropna()
    )


# %%
################################################################
#                          Portfolios
################################################################
# CCM
from CCM import generate_CCM
# load portfolios
from portfolio_construction import generate_portfolio

# generate CCM dataset
generate_CCM(main_path, username)

# generate portfolios
output_portfolio, output_portfolio_27 = generate_portfolio([0, 0.3, 0.7, 1], True)


output_portfolio = pd.merge(output_portfolio, 
                    instruments.loc[:, ['Date', 'CPI']], 
                    how='inner', 
                    left_on='date', 
                    right_on='Date')


output_portfolio_27 = pd.merge(output_portfolio_27, 
                    instruments.loc[:, ['Date', 'CPI']], 
                    how='inner', 
                    left_on='date', 
                    right_on='Date')

# %%
if convert_to_real is True:
    for dataset in [output_portfolio, output_portfolio_27]:
        
        for col in dataset.columns:
            if (col != 'date') and (col != 'CPI') and (col != 'Date'):
                dataset[col] = dataset[col] - dataset['CPI']
                
# %%

for dataset in [output_portfolio, output_portfolio_27]:
    dataset.drop(columns=['Date', 'CPI'], inplace = True)


# %%

if convert_to_real is True:
    name_tag = "Real"
elif convert_to_real is False:
    name_tag = "Nominal"

dates = (
    set(instruments['Date'])
    .intersection(factors['Date'])
    .intersection(output_portfolio['date'])
)

instruments\
    .query("Date in @dates")\
    .to_csv("Input/Instruments.csv", index=False)
    
factors\
    .query("Date in @dates")\
    .to_csv("Input/Factors_{0}.csv".format(name_tag), index=False)
    
    
output_portfolio\
    .query("date in @dates")\
    .dropna(axis = 1)\
    .to_csv("Input/FF5_plus_Industry_Portfolios_{0}.csv".format(name_tag), index=False)
    
    
output_portfolio_27\
    .query("date in @dates")\
    .dropna(axis = 1)\
    .to_csv("Input/27_plus_Industry_Portfolios_{0}.csv".format(name_tag), index=False)


    
# %%
################################################################
#                          Monetary Shocks
################################################################
from monetary_shocks import generate_MS
generate_MS(main_path)




################################################################
#                          Consumption
################################################################
from consumption import generate_consumption
generate_consumption(main_path)