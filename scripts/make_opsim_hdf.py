#!/usr/bin/env python
"""
Script to create an example hdf file from a OpSim database file. The example
hdf file is used for tests.
"""
import opsimsummary as oss
import os
import argparse

example_dir = os.path.join(oss.__path__[0], 'example_data')
dbName = os.path.join(example_dir, 'enigma_1189_micro.db')
default_outfile = '../opsimsummary/example_data/enigma_1189_micro.hdf'

parser = argparse.ArgumentParser(description='Write out the relevant data in'
                                 ' an OpSim sqlite file to an hdf format,'
                                 'appropriate for use with OpSimOut')
parser.add_argument('--OpSimDBPath', type=str, default=None,
                    help='absolute path to the OpSim sqlite database, defaults'
                    ' to None, which does this for `exampe_dir/enigma_1189_micro.dn')

parser.add_argument('--outFile', type=str, default=default_outfile,
                    help='absolute path to output file,'
                    'defaults to {}'.format(default_outfile))
args = parser.parse_args()

opout = oss.OpSimOutput.fromOpSimDB(dbname=args.OpSimDBPath, subset='_all')
opout.writeOpSimHDF(args.outFile)
