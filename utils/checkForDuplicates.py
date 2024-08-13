### this file checks for duplicate files
### to run:
### python3 checkForDuplicates.py
### Note that DataAid.py does not load duplicate files, so dupes
### within a single directory are fine

from os import listdir
from os.path import join
from itertools import chain
from collections import Counter, defaultdict

# Helper functions from DataAid.py and DataImport.py
# I specifically want hashFile
import DataAid
import DataImporter

# These are the directories used in the analysis,
# Model_MarcusPaz_10_2_23.ipynb
directoryList = []
directoryList.append("data/Sparc/Rotmod_LTG/")
directoryList.append("data/Sparc/SparcSubset135/")
directoryList.append("data/Sparc/TrainingSet/")
directoryList.append("data/little-data-things/data/")
directoryList.append("data/LCMFits/data/")
directoryList.append("data/XueSofue/")
directoryList.append("data/McGaugh/")

# get a list of all the files in all the directories
fileList = []
for directory in directoryList:
    thisFileList = [join(directory,x) for x in listdir(directory) if x[0] != '.']
    fileList.append(thisFileList)

# flatten fileList
# and then get the hash of each file in the fileList
fileList = list(chain.from_iterable(fileList))
hashList = [DataAid.hashFile(x) for x in fileList]

def duplicates(lst):
    cnt= Counter(lst)
    return [key for key in cnt.keys() if cnt[key]> 1]

# returns a default dict that looks like
# {hash1: [index0, index1], hash2: [index2, index3, index4], ...}
def duplicates_indices(lst):
    dup, ind= duplicates(lst), defaultdict(list)
    for i, v in enumerate(lst):
        if v in dup: ind[v].append(i)
    return ind

duplicatesDict = duplicates_indices(hashList)

for _, duplicatesIndexList in duplicatesDict.items():
    print("The following files are duplicates:")
    for index in duplicatesIndexList:
        print(fileList[index])
    
    print("\n")