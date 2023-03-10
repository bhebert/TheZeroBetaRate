
# %%
import pandas as pd
import numpy as np
import os



def generate_shadow_spread(main_path):
    process_path = os.path.join(main_path, "Processed Data")
    if not os.path.exists(process_path):
        os.mkdir(process_path)
    
    GBW = pd.read_csv(os.path.join(main_path, 'Raw Data/Shadow Rate/feds200628.csv'), skiprows=9)
    GBW['Date'] = pd.to_datetime(GBW['Date'])

    # %%
    # compute yield curve for different maturities
    new_var = []
    for t in range(1, 13):
        var_name = f"y{str(t)}m"
        new_var.append(var_name)
        GBW[var_name] = GBW['BETA0'] + GBW['BETA1']*np.exp(-t/12/GBW['TAU1']) + GBW['BETA2']*(t/12/GBW['TAU1'])*np.exp(-t/12/GBW['TAU1']) \
            + GBW['BETA3']*(t/12/GBW['TAU2'])*np.exp(-t/12/GBW['TAU2'])

    yc = GBW[['Date'] + new_var]

    output = yc.resample('M', on='Date')['y1m'].agg(['last']).reset_index().set_axis(['Date', 'y1m'], axis=1)
    # output = yc.resample('Q', on='Date')['y3m'].agg(['last']).reset_index().set_axis(['Date', 'y3m'], axis=1)
    output['Date'] = output['Date'] + pd.tseries.offsets.MonthEnd(0)

    # %%
    TB1M = pd.read_excel(os.path.join(main_path, "Raw Data/Factors/ff3.xlsx"))
    TB1M = TB1M.set_axis(['Date', '', '', '','RF'], axis=1)
    TB1M['Date'] = pd.to_datetime(TB1M['Date'].astype(str) + '01') - pd.tseries.offsets.MonthEnd(1)

    # %%
    output = pd.merge(output[['Date', 'y1m']],
                    TB1M[['Date', 'RF']],
                    on='Date', how='inner')

    output['shadow_spread'] = output['y1m']/12 - output['RF']
    output[['Date', 'shadow_spread']].to_csv(os.path.join(process_path, 'shadow_spread.csv'))


# output.to_csv(os.path.join(main_path, "Processed Data/yc1m.csv"))


# %%
