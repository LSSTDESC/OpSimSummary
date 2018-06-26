#!/usr/bin/env bash
gunzip -c opsimsummary/example_data/healpixels_micro.db.gz > opsimsummary/example_data/healpixels_micro.db
python setup/generate_requirements.py
python setup.py install --user
