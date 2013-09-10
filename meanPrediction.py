import numpy as np
import pandas as pd
import datetime

QUARTER_HOUR = 900000000

def strToDatetime64(ix, fmt='%m/%d/%y %H:%M'):
    """Return np.array of type datetime64 for a time series index."""

    ix64 = np.array([ np.datetime64(datetime.datetime.strptime(x, fmt)) for x in ix])
    return ix64

# Each Prediction_DateTime can come from multiple predictions. We simply
# average the predictions from all runtimes.
elec = pd.read_csv('345_Pred_Sorted.csv')
elecMean = elec.groupby(['Prediction_DateTime']).mean()
elecMean.Prediction_Value

#print elecMean

pred = pd.read_csv('345_Elec_Cut.csv')
predSum = pred.groupby('ROUNDED').sum()

#print predSum

# Not all of the indicies for either will appear in the other. (More
# predictions than electric datetimes.)
# So something about registration...
# - We will systemically drop some timestamps? This is annoying.
# - Join the two tables!

predSum.index  = strToDatetime64(predSum.index).astype(np.int64)  / QUARTER_HOUR
elecMean.index = strToDatetime64(elecMean.index).astype(np.int64) / QUARTER_HOUR

outerJoin = pd.merge(elecMean, predSum, left_index=True, right_index=True, how='outer')
outerJoin.sort_index(inplace=True)

outerDiffs = np.diff(np.array(outerJoin.index))
assert np.all(outerDiffs > 0)
print outerDiffs[np.flatnonzero(outerDiffs != outerDiffs[0])]

outerJoin.to_csv('345_Outer_Join.csv')

innerJoin = pd.merge(elecMean, predSum, left_index=True, right_index=True, how='inner')
innerJoin.sort_index(inplace=True)

innerDiffs = np.diff(np.array(innerJoin.index))
assert np.all(innerDiffs > 0)
print innerDiffs[np.flatnonzero(innerDiffs != innerDiffs[0])]

innerJoin.to_csv('345_Inner_Join.csv')

# There are LONG stretches of time when one of the other series (but usually
# the prediction) is NaN. (Discrete time model is a bit problematic.)
#
# NOT FOR THE BAYESIAN MODELS!
#
# Fewer problems on recent dates.

