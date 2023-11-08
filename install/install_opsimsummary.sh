#!/usr/bin/env bash
gunzip -c opsimsummary/example_data/healpixels_micro.db.gz > opsimsummary/example_data/healpixels_micro.db
python3 setup/generate_requirements.py
python3 setup.py install --user
