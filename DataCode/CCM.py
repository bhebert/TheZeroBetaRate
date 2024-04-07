
# %%
import pandas as pd
from os.path import join, exists
from os import makedirs
import numpy as np

# handle time
import datetime as dt
from dateutil.relativedelta import *
from pandas.tseries.offsets import *


def generate_CCM(main_path):
    """Generate the Compustat, CRSP merged dataset that will be used in later construction.

    Args:
        main_path (string): The main directory.
    """

    raw_path = join(main_path, "Raw Data")
    process_path = join(main_path, "Processed Data")
    if not exists(process_path):
        makedirs(process_path)

    # %%
    '''
    Compustat
    '''

    comp = pd.read_csv(join(raw_path, "CCM", "FF5_comp.csv"), engine='pyarrow')

    comp['datadate'] = pd.to_datetime(comp['datadate'].astype(str))
    comp['year']=comp['datadate'].dt.year
    comp['gvkey'] = comp['gvkey'].astype(int)
    comp['fyear'] = np.where(comp['fyear'].isna(), comp['year'], comp['fyear'])
    comp['fyear'] = comp['fyear'].astype(int)
    comp['fyr'] = np.where(comp['fyr'].isna(), 12, comp['fyr'])
    comp['fyr'] = comp['fyr'].astype(int)

    # create preferrerd stock
    comp['ps']=np.where(comp['pstkrv'].isnull(), comp['pstkl'], comp['pstkrv'])
    comp['ps']=np.where(comp['ps'].isnull(),comp['pstk'], comp['ps'])
    comp['ps']=np.where(comp['ps'].isnull(),0,comp['ps'])
    comp['txditc']=np.where(comp['txditc'].isnull(), comp['txdb'] + comp['itcb'], comp['txditc'])
    comp['txditc']=comp['txditc'].fillna(0)

    # fill in nan shareholders equity
    # comp['seq'] = np.where(comp['seq'].isna(), comp['ceq'] + comp['ps'], comp['seq'])
    # comp['seq'] = np.where((comp['seq'].isna()) & (comp['ceq'].isna()), comp['at'] - comp['lt'], comp['seq'])

    # create book equity
    comp['be']=comp['seq']+comp['txditc']-comp['ps']
    comp['be']=np.where(comp['be']>0, comp['be'], np.nan)

    # Deal with duplicates by gvkey, year
    comp['month'] = comp['datadate'].dt.month
    max_month = comp.groupby(['gvkey', 'year'])['month'].max().reset_index().rename(columns={'month':'max_m'})
    comp = pd.merge(comp, max_month, how='left', on=['gvkey', 'year'])
    comp = comp[comp['month']==comp['max_m']]
    comp = comp.drop(columns=['max_m'], axis=1)

    # convert to dec timing
    comp['ddate'] = comp['datadate'] + YearEnd(0)


    # # deal with gaps in the time series
    comp = comp.set_index('ddate')
    comp = comp.groupby('gvkey').resample('YE').last()
    comp = comp.drop(['gvkey'], axis=1).reset_index()


    # profitability
    comp = comp.assign(
        op = lambda x: (x['revt'] - x['cogs'] - x['xint'] - x['xsga'])/x['be']
    )

    # investment
    comp = comp.assign(
        lag_at = lambda x: x.groupby(['gvkey'])['at'].shift(1),
        inv = lambda x: (x['at']/x['lag_at'] - 1) * 100
    )


    # %%
    '''
    CRSP
    '''

    mseall = pd.read_csv(join(raw_path, "CCM", 'mseall.csv'), index_col=0)
    mseall['date'] = pd.to_datetime(mseall['date'])
    mseall['ddate'] = mseall['date'] + MonthEnd(0)
    mseall['nameendt'] = pd.to_datetime(mseall['nameendt'])
    mseall['nameendt'] = mseall['nameendt'].fillna(pd.to_datetime('today'))


    crsp_m = pd.read_csv(join(raw_path, "CCM", "crsp_m.csv"), index_col=0, engine='pyarrow')
    # for col in crsp_m.columns:
    #     crsp_m = crsp_m.rename(columns={col:col.lower()})
    crsp_m['date'] = pd.to_datetime(crsp_m['date'])

    # change variable format to int
    crsp_m[['permco','permno']]=crsp_m[['permco','permno']].astype(int)

    # Line up date to be end of month (Dec Timing)
    crsp_m['ddate']=crsp_m['date']+MonthEnd(0) # the name refers to december date


    ###################
    # Use Delisted Information #
    ###################
    dlret = pd.read_csv(join(raw_path, "CCM", 'dlret.csv'), index_col=0)

    dlret.permno=dlret.permno.astype(int)
    #dlret['dlstdt']=pd.to_datetime(dlret['dlstdt'])
    dlret['ddate']=pd.to_datetime(dlret['dlstdt'])+MonthEnd(0)

    crsp = pd.merge(crsp_m, dlret, how='left',on=['permno','ddate'])
    crsp['dlret']=crsp['dlret'].fillna(0)
    crsp['ret']=crsp['ret'].fillna(0)

    # retadj factors in the delisting returns
    crsp['retadj']=(1+crsp['ret'])*(1+crsp['dlret'])-1

    # deal with negative prices
    crsp['prc'] = np.where(crsp['prc'] < 0, np.absolute(crsp['prc']), crsp['prc'])


    # calculate market equity
    crsp['me']=crsp['prc'].abs()*crsp['shrout']
    crsp=crsp.drop(['dlret','dlstdt','prc','shrout'], axis=1)
    crsp=crsp.sort_values(by=['ddate','permco','me'])



    ### Aggregate Market Cap ###
    # sum of me across different permno belonging to same permco a given date
    crsp_summe = crsp.groupby(['ddate','permco'])['me'].sum().reset_index()

    # largest mktcap within a permco/date
    crsp_maxme = crsp.groupby(['ddate','permco'])['me'].max().reset_index()

    # join by ddate/maxme to find the permno
    crsp1=pd.merge(crsp, crsp_maxme, how='inner', on=['ddate','permco','me'])

    # drop me column and replace with the sum me
    crsp1=crsp1.drop(['me'], axis=1)

    # join with sum of me to get the correct market cap info
    crsp2=pd.merge(crsp1, crsp_summe, how='inner', on=['ddate','permco'])

    # sort by permno and date and also drop duplicates
    crsp2=crsp2.sort_values(by=['permno','ddate']).drop_duplicates()



    # introduce June Timing
    crsp2['jdate'] = crsp2['ddate'] + MonthEnd(+6)

    # timing variable
    crsp2['dmonth'] = crsp2['ddate'].dt.month
    crsp2['jmonth'] = crsp2['jdate'].dt.month
    crsp2['dyear'] = crsp2['ddate'].dt.year
    crsp2['jyear'] = crsp2['jdate'].dt.year

    decme=crsp2[crsp2['dmonth']==12].copy() # market equity at the end of Dec
    decme=decme[['permno','ddate','me','dyear']].rename(columns={'me':'dec_me'})
    decme['dyear'] = decme['dyear'] + 1
    juneme = crsp2[crsp2['jmonth']==12].copy() # market equity at the end of June
    juneme = juneme[['permno','ddate','me','dyear']].rename(columns={'me':'june_me'})
    juneme['dyear'] = juneme['dyear'] + 1



    # compute weight used in later computation
    crsp2['1+retx']=1+crsp2['retx']
    crsp2=crsp2.sort_values(by=['permno','ddate'])

    # cumret by stock
    crsp2['cumretx']=crsp2.groupby(['permno','jyear'])['1+retx'].cumprod()

    # lag cumret
    crsp2['lcumretx']=crsp2.groupby(['permno'])['cumretx'].shift(1)

    # lag market cap
    crsp2['lme']=crsp2.groupby(['permno'])['me'].shift(1)

    # if first permno then use me/(1+retx) to replace the missing value
    crsp2['count']=crsp2.groupby(['permno']).cumcount()
    crsp2['lme']=np.where(crsp2['count']==0, crsp2['me']/crsp2['1+retx'], crsp2['lme'])


    # baseline me
    mebase=crsp2[crsp2['jmonth']==1][['permno','jyear', 'lme']].rename(columns={'lme':'mebase'})

    # merge result back together
    crsp3=pd.merge(crsp2, mebase, how='left', on=['permno','jyear'])
    # here assume the market equity grow as the return when compute the weight. 
    crsp3['wt']=np.where(crsp3['jmonth']==1, crsp3['lme'], crsp3['mebase']*crsp3['lcumretx'])

    # market equity at the end of dec
    crsp3 = pd.merge(crsp3, decme[['permno','dyear','dec_me']], how='inner', on=['permno', 'dyear'])
    crsp3 = pd.merge(crsp3, juneme[['permno','dyear','june_me']], how='inner', on=['permno', 'dyear'])


    crsp4 = pd.merge(crsp3, mseall[['permno', 'permco', 'shrcd', 'exchcd', 'ddate']], how='left', on=['permno', 'permco', 'ddate'])

    crsp4 = crsp4[['permno', 'permco', 'ddate', 'ret', 'retadj', 'me', 'dyear', 'wt', 'dec_me', 'june_me','shrcd', 'exchcd']].drop_duplicates()





    # %%
    '''
    CCM
    '''
    ccm = pd.read_csv(join(raw_path, "CCM", "ccm_link.csv"), index_col=0)

    # if linkenddt is missing then set to today date 
    ccm['linkenddt']=ccm['linkenddt'].fillna(pd.to_datetime('today'))
    ccm['linkdt'] = pd.to_datetime(ccm['linkdt'])
    ccm['linkenddt'] = pd.to_datetime(ccm['linkenddt'])

    ccm1 = pd.merge(comp, ccm, how='left', on=['gvkey'])

    ccm1['ddate'] = ccm1['datadate'] + YearEnd(0)
    ccm1['dyear'] = ccm1['ddate'].dt.year
    ccm1['jdate'] = ccm1['ddate'] + MonthEnd(6)
    ccm1['jyear'] = ccm1['ddate'].dt.year

    # set link date bounds
    ccm2=ccm1[(ccm1['jdate']>=ccm1['linkdt'])&(ccm1['jdate']<=ccm1['linkenddt'])].copy()
    ccm2 = ccm2[['gvkey', 'tic', 'permno', 'jyear', 'fyear', 'fyr', 'be', 'op', 'inv']]
    # this is key when doing June timing because we want to use last dec's financial information for this June.
    ccm2['jyear'] = ccm2['jyear'] + 1


    # merge with the stock information 
    ccm3 = pd.merge(crsp4, ccm2, how='left', left_on=['permno', 'dyear'], right_on=['permno', 'jyear'])
    ccm3 = ccm3[ccm3['dyear']>=1960]

    # deal with zero market equity
    # here I am using the BE information from Dec, and the ME information from Dec as well. 
    # This is essentially the same as the implementation by Qingyi (Freda) Song Drechsler in her example. 
    ccm3['bm'] = np.where(ccm3['dec_me']>0, ccm3['be']*1000/ccm3['dec_me'], np.nan)

    ccm3 = ccm3.drop_duplicates(subset=['permno', 'ddate'])


    ccm3.to_csv(join(process_path, 'ccm.csv'))

