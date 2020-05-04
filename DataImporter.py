# This is a slightly updated version of the DataAid import system
# The data comes in several formats
# getGalaxyData handles Little Things and Sparc
# getXueSofue handles the XueSofue MilkyWay
# Both return Pandas DataFrames

import pandas as pd
from numpy import sqrt
import io
from functools import reduce

# This is admittedly a little convoluted
# Using StringIO to "preprocess" what we give to Pandas
# Probably would be clearer to just process directly
def getGalaxyData(filename):
    """Reads data in either the Little Things or Sparc format
    
    Data format is:
    - Optional header line (starts with a #)
    - Column names (also with a #)
    - Units (also with a #)
    - tab separated rows of values"""
    
    def processLine(line):
        # If it starts with a #, special handling
        # Check if it contains "rad", meaning it's the column names
        if line.strip()[0] == '#':
            if "rad" in line.lower():
                # Pop off the #
                return line.strip()[1:].strip() + '\n'
            else:
                return "\n"
        else:
            return line.strip() + '\n'
        
    buffer = io.StringIO()
    
    with open(filename) as f:
        for line in f:
            buffer.write(processLine(line))
            
    buffer.seek(0)
            
    return pd.read_csv(buffer, sep='\t')

def getXueSofue(filename):
    """Reads in the XueSofue data format
    
    Data format is:
    - Several header lines starting with #
    - No column names
    - Tab separated data, two columns, radius and vLum"""
    
    return pd.read_csv(filename, sep='\t', comment='#', names=['rad', 'vLum'])

def calcVLum(galaxy, cols=['Vgas', 'Vdisk', 'Vbul'], inPlace=True):
    """This computes vLum column for the galaxy DataFrame.
    
    Computes sum of cols squared. If inPlace=True, creates a new column
    in the galaxy DataFrame. If inPlace=False, returns a DataSeries"""
    
    vLum = sqrt(reduce(lambda x,y: x+y, [galaxy[col]**2 for col in cols]))
    if inPlace:
        galaxy['vLum'] = vLum
        return None
    else:
        return vLum