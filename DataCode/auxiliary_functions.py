import pandas as pd
from dateutil.relativedelta import *
from pandas.tseries.offsets import *
import numpy as np
import os
from os.path import join, exists
from os import mkdir
from datetime import timedelta
from calendar import THURSDAY





def generate_consumption(main_path):

    raw_path = join(main_path,  "Raw Data/Consumption")
    output_path = join(main_path, "Input")


    # %% reading function

    def read_data(relative_path:str, cols:str, srows:int, nrows:int, start_year:int, end_year:int, month_or_quarter:str, end_m_or_q:int, rescale:bool):
        """Perform data import and initial cleaning.

        Args:
            relative_path (str): path to the xlsx file.
            cols (str): columns to read
            srows (int): skiprows
            nrows (int): nrows
            start_year (int):  
            end_year (int): 
            month_or_quarter (str): data is in monthly or quarterly frequency
            end_m_or_q (int): last recorded month or quarter
            rescale (bool): if the variables need to be rescaled.
        """
        
        df = pd.read_excel(join(raw_path, relative_path), usecols=cols, skiprows=srows, nrows=nrows)
        col_list = ['variable']
        for year in np.arange(start_year, end_year+1):
            if month_or_quarter == "Q":
                for quarter in np.arange(1, 5):
                    if year != end_year:
                        col_list.append("{0}Q{1}".format(year, quarter))
                    else:
                        if quarter <= end_m_or_q:
                            col_list.append("{0}Q{1}".format(year, quarter))
                        else:
                            pass
            elif month_or_quarter == "M":
                for month in np.arange(1, 13):
                    if year != end_year:
                        col_list.append("{0}M{1}".format(year, month))
                    else:
                        if month <= end_m_or_q:
                            col_list.append("{0}M{1}".format(year, month))
                        else:
                            pass

        df = df.set_axis(col_list, axis = 1)
        df['variable'] = df['variable'].str.strip()
        df['variable'] = df['variable'].str.replace('\d+', '', regex=True)
        df = df.set_index('variable').T.reset_index()

        if month_or_quarter == 'Q':
            df = (
                df
                .assign(index = lambda x: pd.PeriodIndex(x['index'], freq='Q').to_timestamp() + MonthEnd(3))
                .rename(columns = {'index':'date'})
            ) 
        elif month_or_quarter == "M":
            df = (
                df
                .assign(index = lambda x: x['index'].str.replace('M', '-') + '-01')
                .assign(index = lambda x: pd.to_datetime(x['index']) + MonthEnd())
                .rename(columns = {'index':'date'})
            )
            
        for col in df.columns:
            if col != 'date':
                df[col] = pd.to_numeric(df[col], errors='coerce')
                
                if rescale is True:
                    # original data is annualized.
                    if month_or_quarter == "Q":
                        df[col] = df[col] / 4
                    elif month_or_quarter == "M":
                        df[col] = df[col] / 12
        
        return df

    # Read Data

    # Consumption Monthly
    con_m = read_data("NIPA_table_2_8_5_M.xlsx", "B:ACL", 6, 11, 1959, 2022, 'M', 8, True)  

    # Price index monthly
    pi_m = read_data("NIPA_table_2_8_4_M.xlsx", "B:ACL", 6, 11, 1959, 2022, 'M', 8, False )

    # Monthly population
    pop_m = read_data("NIPA_table_2_6_M.xlsx", "B:ACN", 6, 43, 1959, 2022, 'M', 10, False)


    ## Series
    # Hall (1988)
    hall_m = (
        pd.merge(
            con_m[['date', 'Nondurable goods']], 
            pi_m[['date', 'Nondurable goods']], 
            how = 'left', on = 'date'
            )
        .merge(pop_m[['date', 'Population (midperiod, thousands)']],
            how='left', on='date')
        .assign(non_dur_real = lambda x: x['Nondurable goods_x'] / x['Nondurable goods_y'] * 100,
                non_dur_real_pc = lambda x: x['non_dur_real'] / x['Population (midperiod, thousands)'] * 1000)
        .loc[:, ['date', 'non_dur_real', 'non_dur_real_pc']]
        .assign(cg_nd_r = lambda x: x['non_dur_real'].diff()/x['non_dur_real'].shift(1) * 100,
                cg_nd_r_pc = lambda x: x['non_dur_real_pc'].diff()/x['non_dur_real_pc'].shift(1) * 100)
        .dropna()
    )

    # Jagannathan and Wang (2007)
    jw_m = (
        pd.merge(
            con_m[['date', 'Nondurable goods', 'Services']], 
            pi_m[['date', 'Nondurable goods', 'Services']], 
            how = 'left', on = 'date'
            )
        .merge(pop_m[['date', 'Population (midperiod, thousands)']],
            how='left', on=['date'])
        .assign(non_dur_real = lambda x: x['Nondurable goods_x'] / x['Nondurable goods_y'] * 100,
                service_real = lambda x: x['Services_x'] / x['Services_y'] * 100,
                nd_sv_real = lambda x: x['non_dur_real'] + x['service_real'],
                nd_sv_real_pc = lambda x: x['nd_sv_real'] / x['Population (midperiod, thousands)'] * 1000)
        .loc[:, ['date', 'nd_sv_real', 'nd_sv_real_pc']]
        .assign(cg_nd_sv_r = lambda x: x['nd_sv_real'].diff()/x['nd_sv_real'].shift(1) * 100 ,
                cg_nd_sv_r_pc = lambda x: x['nd_sv_real_pc'].diff()/x['nd_sv_real_pc'].shift(1) * 100 )
        .dropna()
    )

    # Output Datasets
    output_list = ['hall_m', 'jw_m']
    # 
    if not exists(join(output_path, 'Consumption')):
        mkdir(join(output_path, 'Consumption'))
    
    for df in output_list:
        eval(df).to_csv(join(output_path, 'Consumption', df + '.csv'), index=False)




def generate_shadow_spread(main_path):
    process_path = os.path.join(main_path, "Processed Data")
    if not os.path.exists(process_path):
        os.mkdir(process_path)
    
    GBW = pd.read_csv(os.path.join(main_path, 'Raw Data/Shadow Rate/feds200628.csv'), skiprows=9)
    GBW['Date'] = pd.to_datetime(GBW['Date'])

    # compute yield curve for different maturities
    new_var = []
    
    def func1(GBW, t):
        new_series = GBW['BETA0'] + GBW['BETA1']*(1-np.exp(-t/12/GBW['TAU1']))/(t/12/GBW['TAU1']) + GBW['BETA2']*((1-np.exp(-t/12/GBW['TAU1']))/(t/12/GBW['TAU1']) - np.exp(-t/12/GBW['TAU1']))
        return new_series
    
    def func2(GBW, t):
        new_series = GBW['BETA0'] + GBW['BETA1']*(1-np.exp(-t/12/GBW['TAU1']))/(t/12/GBW['TAU1']) + GBW['BETA2']*((1-np.exp(-t/12/GBW['TAU1']))/(t/12/GBW['TAU1']) - np.exp(-t/12/GBW['TAU1'])) \
            + GBW['BETA3']*((1-np.exp(-t/12/GBW['TAU2']))/(t/12/GBW['TAU2']) - np.exp(-t/12/GBW['TAU2']))
        return new_series
    
    for t in range(1, 24):
        var_name = f"y{str(t)}m"
        new_var.append(var_name)
        GBW[var_name] = np.where(GBW['TAU2']<-999, func1(GBW, t), func2(GBW, t))

    yc = GBW[['Date'] + new_var]
    
    
    def last_thursday(dt):
        offset = (dt.weekday() - THURSDAY) % 7
        last_thursday = dt - timedelta(days=offset)

        return last_thursday

    yc = pd.concat(
        [yc, yc.apply(lambda row: last_thursday(row['Date']), axis=1).rename('last_thursday')],
        axis=1
    )
    yc_q = (
        yc
        .groupby('last_thursday')[new_var].mean()
        .resample('M').mean()
        .reset_index()
        .rename(columns={'last_thursday':'Date'})
        )

    TB3M = pd.read_csv(os.path.join(main_path, "Raw Data/Shadow Rate/DTB3.csv"))
    TB3M = TB3M.rename(columns={'DATE':'Date', 'DTB3':'t3m'})
    TB3M = TB3M.assign(
        Date = lambda x: pd.to_datetime(x['Date']),
        t3m = lambda x: pd.to_numeric(x['t3m'], errors='coerce')
    )
    TB3M = pd.concat(
        [TB3M, TB3M.apply(lambda row: last_thursday(row['Date']), axis=1).rename('last_thursday')],
        axis=1
    )
    
    TB3M_q = (
        TB3M
        .groupby('last_thursday')['t3m'].mean()
        .resample('M').mean()
        .reset_index()
        .rename(columns={'last_thursday':'Date'})
        )
    
    output = pd.merge(yc_q[['Date', 'y3m']],
                    TB3M_q[['Date', 't3m']],
                    on='Date', how='inner')

    output['shadow_spread'] = output['y3m']/12 - output['t3m']/12
    output[['Date', 'shadow_spread']].to_csv(os.path.join(process_path, 'shadow_spread.csv'))
    
    
    
    
def generate_MS(main_path):

    data_path = join(main_path, "Raw Data", "Monetary_shocks")
    output_path = join(main_path, "Input")

    # Nakamura and Steinsson Monetary Policy Shock
    ns = (
        pd.read_excel(join(data_path, "NS.xlsx"), sheet_name='shocks')
        .assign(Date = lambda x: pd.to_datetime(x['fomc']) + MonthEnd(0))
        .drop(columns = ['fomc'])
        .set_index(['Date'])
        .asfreq(freq='1M', fill_value=0)
        .reset_index()
    )
    
    #  Romer and Romer Shock
    rr = (
        pd.read_stata(join(data_path,"RR_monetary_shock_monthly.dta"))
        .assign(Date = lambda x: pd.to_datetime(x['date'] + MonthEnd(0))) 
        .drop(columns = ['date'])
        .set_index(['Date'])
        .asfreq(freq='1M', fill_value=0)
        .reset_index()
    )

    ns.to_csv(join(output_path, "NS_shocks.csv"), index=False)
    rr.to_csv(join(output_path, "RR_shocks.csv"), index=False)
    
    
    
    
def DP_ratio(main_path):
    """Compute dividend price ratio using crsp indices. 

    Args:
        main_path (str): 

    Returns:
        DataFrame: DataFrame of dividend price ratio
    """
    raw_path = join(main_path, "Raw Data")
    msi = (
        pd.read_csv(join(raw_path, "CCM", "msi.csv"), index_col=0)
        .assign(date = lambda x: pd.to_datetime(x['date']) + MonthEnd(0))
        .assign(A = lambda x: x['vwretd']+1,
                B = lambda x: x['vwretx']+1)
        )
    # .assign(p_p12 = lambda x: x['B'].rolling(12).apply(np.prod))

    dates = msi['date'].unique()
    outs = []
    for date in dates:
        start_dt = date - MonthEnd(11)
        dts = msi.query("date >= @start_dt and date <= @date")
        
        if len(dts) < 12:
            continue

        temp1 = dts['B'].prod()
        temp2 = 0
        for i in range(12):

            if i == 0:
                temp2 = temp2 + (dts.iloc[i]['A'] - dts.iloc[i]['B'])
            else:
                # print(len(dts.iloc[0:i]))
                temp2 = temp2 + (dts.iloc[i]['A'] - dts.iloc[i]['B']) * dts.iloc[0:i]['B'].prod()
                
        outdf = pd.DataFrame([date, temp2/temp1]).T
        outs.append(outdf)
        
    output_df = pd.concat(outs, ignore_index=True)
    output_df.columns = ['Date', 'DP_ratio']
    
    return output_df





def generate_value_spread(main_path):
    raw_path = join(main_path, "Raw Data")
    
    sbv_ret = (
        pd.read_excel(join(raw_path, "Predefined Portfolios", "6_Portfolios_2x3.xlsx"), skiprows=15, nrows=1159)
        .rename(columns = {'Unnamed: 0':'date'})
        .assign(date = lambda x: pd.to_datetime(x['date'].astype(str)+'01') + MonthEnd(0))
    )
    
    sbv_log_bm = (
        pd.read_excel(join(raw_path, "Predefined Portfolios", "6_Portfolios_2x3.xlsx"), skiprows=4871, nrows=1159)
        .rename(columns = {'Unnamed: 0':'date'})
        .assign(date = lambda x: pd.to_datetime(x['date'].astype(str)+'01') + MonthEnd(0))
        .assign(value_spread = lambda x: np.log(x['SMALL HiBM']) - np.log(x['SMALL LoBM']))
        .assign(jdate = lambda x: x['date']-MonthEnd(5),
                jyear = lambda x: x['jdate'].dt.year)
        .query("date.dt.month == 6")
    )
    
    merged = pd.merge(
        (
            sbv_ret
            .assign(jdate = lambda x: x['date']-MonthEnd(5),
                    jyear = lambda x: x['jdate'].dt.year)
            .loc[:, ['date', 'jdate', 'jyear', 'SMALL LoBM', 'SMALL HiBM']]
            ),
        sbv_log_bm[['jyear', 'value_spread']], on='jyear', how='inner')
    for col in ['SMALL LoBM', 'SMALL HiBM']:
        merged[col] = np.log(merged[col]/100 + 1)
        
    merged.loc[merged['date'].dt.month == 6, ['SMALL LoBM']]=0
    merged.loc[merged['date'].dt.month == 6, ['SMALL HiBM']]=0

    merged = (
        merged
        .assign(S_LBM = lambda x: x.groupby(['jyear'])['SMALL LoBM'].cumsum(),
                S_HBM = lambda x: x.groupby(['jyear'])['SMALL HiBM'].cumsum())
        .assign(value_spread = lambda x: x['value_spread'] + x['S_LBM'] - x['S_HBM'])
        .rename(columns={'date':'Date'})
    )
    
    return merged[['Date', 'value_spread']]
    