from __future__ import absolute_import
import opsimsummary as oss

def test_versionName():
    """
    The only objective of this test is to check that this is running,
    """
    version = oss.__VERSION__
    print("test_versionName")
    print("Running version ", version)
    return None

