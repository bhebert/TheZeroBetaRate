# %%
import pandas as pd
from os.path import join

# handle time
from dateutil.relativedelta import *
from pandas.tseries.offsets import *


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

        
# %%
