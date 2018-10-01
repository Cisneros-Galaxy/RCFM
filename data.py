from os import listdir

class Data:
  GALAXIES_DIR_REL_PATH = "data/Rotmod_LTG/"

  def getGalaxyFiles(self, galaxiesDirRelPath=GALAXIES_DIR_REL_PATH):
    return listdir(galaxiesDirRelPath)

  def dataLineToValues(self, line):
    stringValues = line.split("\t")
    floatValues = list(map(lambda x: float(x), stringValues))
    return floatValues

  def getGalaxyData(self, galaxiesDirRelPath=GALAXIES_DIR_REL_PATH):
    files = self.getGalaxyFiles(galaxiesDirRelPath)

    # Dictionary "galaxy name"->[list of data points]
    data = {}

    for fileName in files:
      file = open(galaxiesDirRelPath + fileName, "r")
      # Convert the file string to an array of each line
      lines = file.read().splitlines()
      # Cut off the first 3 lines which are just comments
      dataLines = lines[3:]
      # Convert each line to a list of data points
      parsedLines = list(map(self.dataLineToValues, dataLines))
      
      data[fileName] = parsedLines
      file.close()

    return data