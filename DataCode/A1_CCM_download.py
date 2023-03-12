
# %%

import pandas as pd
import os

import numpy as np

# handle time
import datetime as dt
from dateutil.relativedelta import *
from pandas.tseries.offsets import *


import wrds


main_path = "/Users/qitong/Library/CloudStorage/Dropbox/Projects/TheZeroBetaRate"


ccm_path = os.path.join(main_path, "Raw Data", "CCM")
if not os.path.exists(ccm_path):
    os.makedirs(ccm_path)


# %%
###################
# Connect to WRDS #
###################
conn=wrds.Connection(wrds_username='qitongwang')



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
                      select a.permno, a.permco, a.date,
                      a.ret, a.retx, a.shrout, a.prc, b.shrcd, b.exchcd, b.naics
                      from crsp.msf as a
                      left join crsp.msenames as b
                      on a.permno=b.permno
                      and b.namedt<=a.date
                      and a.date<=b.nameendt
                      where a.date between '01/01/1959' and '12/31/2021'
                      and b.exchcd between 1 and 3 
                      """, date_cols=['date']) 

crsp_m = crsp_m.drop(columns=['shrcd', 'exchcd'])
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
