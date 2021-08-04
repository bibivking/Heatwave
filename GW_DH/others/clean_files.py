#!/usr/bin/env python

"""
Clean files: add pandas readable timestamp

That's all folks.
"""
__author__ = "Martin De Kauwe"
__version__ = "1.0 (24.06.2021)"
__email__ = "mdekauwe@gmail.com"


import pandas as pd
from datetime import date, timedelta
import sys
import numpy as np

def clean_files(fname, start_yr):

    df = pd.read_csv(fname, header=None)
    df.columns = ["day", "temp"]

    # Add correct timestamps
    start = date(start_yr,1,1)
    dates = []
    for i in range(len(df)):
        # need the -1 or the first day is 17 not 16
        delta = timedelta(int(df.day[i]) - 1)
        offset = start + delta
        dates.append(offset)

    df["dates"] = dates
    #df = df.set_index('dates')
    df.index = pd.to_datetime(df.dates, format = '%Y/%m/%d', utc=False)
    df["year"] = df.index.year

    ofname = fname.replace("txt", "csv")
    df.to_csv(ofname, index=True)

if __name__ == "__main__":

    start_yr = 2000
    fn = "EF_GW_rawdata_4_Python.txt"
    clean_files(fn, start_yr)

    fn = "EF_FD_rawdata_4_Python.txt"
    clean_files(fn, start_yr)
