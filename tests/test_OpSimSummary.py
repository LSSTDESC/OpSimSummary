import OpSimSummary.summarize_opsim as oss
import pandas as pd

def test_fieldIds():
    df = pd.read_hdf('/Users/rbiswas/data/LSST/OpSimData/storage.h5', key='table')
    so = oss.SummaryOpsim(df)
    assert isinstance(so.fieldIds, list)
    assert all(list(isinstance(fieldID, int) for fieldID in so.fieldIds))
