import pandas
c = pandas.read_csv('345_Pred.csv')
c
c['Run_DateTime']
c
c.columns
c.columns['Run_DateTime']
c
c
c.Prediction_DateTime
c.index
c
c.Prediction_DateTime
c.Prediction_DateTime[100:200]
c
cmean = c.groupby(['Prediction_DateTime'])
cmean
cmean = c.groupby(['Prediction_DateTime']).mean()
cmean
c
cmean.Prediction_Value
cmean.Prediction_Value
cmean
cmean.Lower_Bound_95
cmean
e = pandas.read_csv('345_Elec.csv')
e
esum = e.groupby(['EQUIPMENT_NO']).sum()
esum
esum = e.groupby(['EQUIPMENT_NO', 'TIMESTAMP']).sum()
esum
e
e.TIMESTAMP
e.TIMESTAMP.round()
e
e.TIMESTAMP
e.TIMESTAMP.round(0)
e.TIMESTAMP
np.datetime64(e.TIMESTAMP)
e.TIMESTAMP
e.TIMESTAMP[0]
e.TIMESTAMP[1]
dt = [ np.datetime64(t) for t in e.TIMESTAMP ]
dt
dt[0-]
dt[0]
e
dt[0].round()
%whos
e
dt = pandas.DatetimeIndex(e.TIMESTAMP)
dt
dt[0]
dt[1]
dt[2]
dt[23]
dt.round?
dt.round()
dt
dt.astype(np.int64)
dt[0]
dt[0].value
%history
dt
dt.astype(np.int64)
dt.astype(np.int64) // 10 ** 9
eTimeSteps = dt.astype(np.int64) // 10 ** 9
eTimeSteps
eTimeSteps[1] - eTimeSteps[0]
np.diff(eTimeSteps)
np.flatnonzero(np.diff(eTimeSteps) != 900)
%whos
ee = pandas.read_csv('345_Elec_Sorted.xls')
ee
ee = pandas.read_csv('345_Elec_Sorted.csv')
ee
ee.TIMESTAMP
eeSum = ee.groupby(['EQUIPMENT_NO', 'TIMESTAMP']).mean()
eeSum
eeSum[0]
eeSum[0]
%whos
eeSum
eeSum.VALUE
eeSum = ee.groupby(['EQUIPMENT_NO', 'TIMESTAMP']).sum()
ee
eeSum
ee.TIMESTAMP
ee.TIMESTAMP[0]
ee.TIMESTAMP[1]
ee.TIMESTAMP[2]
ee.TIMESTAMP[3]
ee.TIMESTAMP[4]
ee.TIMESTAMP[5]
eeSum = ee.groupby(['TIMESTAMP']).sum()
ee
eeSum
eeSum.VALUE
pp = pandas.read_csv('345_Pred_Sorted.csv')
pp
pp.Prediction_DateTime
ppMean = pp.groupby(['Prediction_DateTime'])
ppMean = pp.groupby(['Prediction_DateTime']).mean()
ppMean
eeSum
eeSum
eeSum.value
eeSum.VALUE
eeSum.select?
eeSum.columns?
eeSum
%history
eeSum = ee.groupby(['EQUIPMENT_NO']).sum()
eeSum
eeSum = ee.groupby(['EQUIPMENT_NO', 'TIMESTAMP']).sum()
eeSum
eeSum.VALUE
eeSum = ee.groupby(['TIMESTAMP', 'EQUIPMENT_NO']).sum()
eeSum
eeSum[0]
eeSum.VALUE
eeSum = ee.groupby(['TIMESTAMP']).sum()
eeSum
eeSum.EQUIPMENT_NO
eeSum = ee.groupby(['EQUIPMENT_NO']).sum()
eeSum
eeSum
eeSum = ee.groupby(['TIMESTAMP', 'EQUIPMENT_NO']).sum()
eeSum
eeSum.VALUE
eeSum.VALUE.sum()
eeSum.VALUE??
eeSum.VALUE.sum(0)
eeSum.VALUE.sum(1)
eeSum.VALUE
ee.groupby(['TIMESTAMP', 'EQUIPMENT_NO']).groups()
eeSum.VALUE.sum(level=1)
eeSum.VALUE.sum(level=0)
eeSum.VALUE
eeSum.VALUE.sum(level=1)
eeSum.VALUE.sum(level=0)
electric_measurement = eeSum.VALUE.sum(level=0)
ppMean
ppMean.Prediction_Value
ee
ee.VALUE
ee.TIMESTAMP
%whos
electric_measurement
electric_measurement.last
electric_measurement.last()
electric_measurement.last?
electric_measurement.last('1M')
electric_measurement.last(1)
electric_measurement
electric_measurement[-1]
electric_measurement[-2]
electric_measurement.index
ee
ee.index
ee.TIMESTAMP
ee.TIMESTAMP[-1]
ee.TIMESTAMP.last
ee.TIMESTAMP.last()
ee.TIMESTAMP.last(1)
ee.TIMESTAMP
%whos
electric_measurement
ee
ee.groupby(['TIMESTAMP'])
ee.groupby(['TIMESTAMP']).sum(0
)
ee.groupby(['TIMESTAMP']).sum()
ee
ee.groupby(['EQUIPMENT_NO']).sum()
ee.groupby(['timestamp']).sum()
ee.groupby(['TIMESTAMP']).sum()
ee
ee.TIMESTAMP
ts = ee.TIMESTAMP
ts
ts?
ts.in
ts.index
ts.index[5]
ts.get(100)
ts
ts.get(347451)
ee.groupby(['TIMESTAMP'])
ee.groupby(['TIMESTAMP'])?
ee.groupby(['TIMESTAMP']).groups()
grp = ee.groupby(['TIMESTAMP'])
grp
grp?
grp.ngroups
grp.last
grp.last()
ee
eeI = ee
eeI.index = eeI.TIMESTAMP
eeI
eeI.groupby('TIMESTAMP').sum()
ee
ee.TIMESTAMP
ppMean
ppMean.Prediction_Value
%save
%save?
%history
%history?
%history -f ppMean.py
