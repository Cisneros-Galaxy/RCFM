from os import listdir
from os.path import join
import hashlib

# Python program to find the SHA-1 message digest of a file
def hashFile(filename):
   """"This function returns the SHA-1 hash
   of the file passed into it"""

   # make a hash object
   h = hashlib.sha1()

   # open file for reading in binary mode
   with open(filename,'rb') as file:

       # loop till the end of the file
       chunk = 0
       while chunk != b'':
           # read only 1024 bytes at a time
           chunk = file.read(1024)
           h.update(chunk)

   # return the hex representation of digest
   return h.hexdigest()

# Goes into the specified folder for Galaxies and returns
#  all unique files in the specified relative folder
def getFiles(dirRelPath):
  fileList = [x for x in listdir(dirRelPath) if x[0] != '.']
  hashList = [hashFile(join(dirRelPath,x)) for x in fileList]

  filteredFileList = []
  filteredHashList = []
  for thisFile, thisHash in zip(fileList, hashList):
     if thisHash not in filteredHashList:
        filteredHashList.append(thisHash)
        filteredFileList.append(thisFile)

  filteredFileList = fileList

  return filteredFileList

# Given a line of data from a file, convert each string value
#  (which are sepperated by tabs) to floats
def dataLineToValues(line):
  stringValues = line.split("\t")
  floatValues = list(map(lambda x: float(x), stringValues))
  return floatValues

# Returns all the galaxy data in the given folder
def GetGalaxyData(galaxiesDirRelPath):
  files = getFiles(galaxiesDirRelPath)

  # Dictionary "galaxy name"->[list of data points]
  data = {}

  for fileName in files:
    file = open(galaxiesDirRelPath + fileName, "r")
    # Convert the file string to an array of each line
    lines = file.read().splitlines()
    # Cut off the first 3 lines which are just comments
    dataLines = lines[3:]
    # Convert each line to a list of data points
    parsedLines = list(map(dataLineToValues, dataLines))
    
    data[fileName[:-4]] = parsedLines
    file.close()

  return data
