
# %%

import pandas as pd
import os

import numpy as np

# handle time
import datetime as dt
from dateutil.relativedelta import *
from pandas.tseries.offsets import *


import wrds


with open('path_variables.md', 'r') as f:
    main_path = f.readline().split('"')[1]
    username = f.readline().split('"')[1]
    
    
# %%
ccm_path = os.path.join(main_path, "Raw Data", "CCM")
if not os.path.exists(ccm_path):
    os.makedirs(ccm_path)


# %%
###################
# Connect to WRDS #
###################
conn=wrds.Connection(wrds_username=username)



# %%
###################
#    Compustat    #
###################
comp = conn.raw_sql("""
                    select gvkey, datadate, fyear, indfmt, consol, popsrc, datafmt,
                    tic, curcd, fyr, at, cogs, dvt, ib, itcb, oibdp,
                    ppegt, ppent, pstk, pstkl, pstkrv, revt, sale, seq,
                    tie, txdb, txditc, xint, xsga, costat, fic, mkvalt
                    from comp.funda
                    where indfmt='INDL' 
                    and datafmt='STD'
                    and popsrc='D'
                    and consol='C'
                    and datadate >= '01/01/1959'
                    """, date_cols=['datadate'])

comp.to_csv(os.path.join(ccm_path, "FF5_comp.csv"), index=False)


# %%
###################
#      CRSP       #
###################
crsp_m = conn.raw_sql("""
                      select permno, permco, date,
                      ret, retx, shrout, prc
                      from crsp.msf 
                      where date <= '12/31/2020'
                      """, date_cols=['date']) 

# crsp_m = crsp_m.drop(columns=['shrcd', 'exchcd'])
crsp_m.to_csv(os.path.join(ccm_path, "crsp_m.csv"), index=True)

# %%
mseall = conn.raw_sql("""
                      select date, permno, permco, divamt, shrcd, exchcd, nameendt,
                      naics
                      from crsp.mseall
                      where exchcd between 1 and 3
                      """,  date_cols=['date'])

mseall.to_csv(os.path.join(ccm_path, "mseall.csv"), index=True)

# %%
dlret = conn.raw_sql("""
                     select permno, dlret, dlstdt 
                     from crsp.msedelist
                     """, date_cols=['dlstdt'])

dlret.to_csv(os.path.join(ccm_path, "dlret.csv"), index=True)

# %%
ccm_link = conn.raw_sql("""
                   select gvkey, lpermno, linktype, linkprim, 
                   linkdt, linkenddt
                   from crsp.ccmxpf_lnkhist
                   """)
ccm_link = (
    ccm_link.query("~ lpermno.isnull()")
    .query("linktype != 'NR'")
    .query("linktype != 'NP'")
    .query("linktype != 'NU'")
    .assign(linkenddt = lambda x: pd.to_datetime(x['linkenddt']))
    .rename(columns={'lpermno':'permno'})
)

ccm_link.to_csv(os.path.join(ccm_path, "ccm_link.csv"), index=True)



# %%
msi = conn.raw_sql("""
                   select date, vwretd, vwretx
                   from crsp.msi
                   """)

msi.to_csv(os.path.join(ccm_path, "msi.csv"), index=True)
# %%
