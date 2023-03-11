
# %%

import pandas as pd
import numpy as np
# pd.set_option('display.min_rows', 50)
from pandas.tseries.offsets import MonthEnd
from os.path import join, exists
from os import mkdir



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




    # %% Read Data

    # Consumption Monthly
    con_m = read_data("NIPA_table_2_8_5_M.xlsx", "B:ACL", 6, 11, 1959, 2022, 'M', 8, True)  

    # Price index monthly
    pi_m = read_data("NIPA_table_2_8_4_M.xlsx", "B:ACL", 6, 11, 1959, 2022, 'M', 8, False )

    # Monthly population
    pop_m = read_data("NIPA_table_2_6_M.xlsx", "B:ACN", 6, 43, 1959, 2022, 'M', 10, False)

    # %%

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






    # %%

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



    # %% Output Datasets
    output_list = ['hall_m', 'jw_m']


    
    # %%
    if not exists(join(output_path, 'Consumption')):
        mkdir(join(output_path, 'Consumption'))
    
    for df in output_list:
        eval(df).to_csv(join(output_path, 'Consumption', df + '.csv'), index=False)

# %%
