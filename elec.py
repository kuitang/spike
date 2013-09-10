import pandas
pp = pandas.read_csv('345_Elec_Sorted.csv')
pp
pp.ROUNDED == pp.TIMESTAMP
np.all(pp.ROUNDED == pp.TIMESTAMP)
np.flatnonzero(pp.ROUNDED != pp.TIMESTAMP)
pp.get(189658)
pp.get?
pp
pp.get(0)
pp.get(1)
pp.index
pp?
pp.data
pp[1]
pp.__getitem__?
pp[-1]
pp.last_valid_index
pp.last_valid_index()
pp.get(347452)
a = pp.get(347452)
a
a
%history
rounds = np.flatnonzero(pp.ROUNDED != pp.TIMESTAMP)
rounds
pp?
pp.loc
pp.ix[189658]
pp.ix[rounds]
pp.groupby(['ROUNDED']).sum(0
)
pp.groupby(['ROUNDED']).sum()
pp
ppSum0 = pp.groupby(['ROUNDED']).sum()
ppSum0.ix[0]
ppSum0.ix[1]
ppSum0.ix[3]
ppSum0.ix[4]
ppSum0.ix[-1]
ppSum1 = pp.groupby(['ROUNDED', 'EQUIPMENT_NO']).sum()
ppSum1
ppSum1.ix[0]
ppSum1.ix[1]
pp
ppSum1.ix[0:100]
ppSum1.ix[-100:-1]
ppSum1
ppSum
pp
pp.ix[-1]
pp.ix[0]
pp.ix[pp.last_valid_index]
pp.ix[pp.last_valid_index()]
pp.ix[pp.last_valid_index() - 1]
pp.ix[pp.last_valid_index() - 2]
pp.ix[pp.last_valid_index() - 3]
pp.ix[pp.last_valid_index() - 4]
ppCut = pandas.read_csv('345_Elec_Sorted.csv')
ppCut
ppCut = pandas.read_csv('345_Elec_Cut.csv')
ppCut
ppCut.groupby('ROUNDED')
ppCut.groupby('ROUNDED').sum()
ppBetter = ppCut.groupby('ROUNDED').sum()
ppBetter
ppBetter.ix[2]
ppBetter.ix[0]
ppBetter.ix[0]
%history -f elec.py
