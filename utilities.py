def readTextFile(filePath):
    with open(filePath, "r") as f:
        return "".join([l.strip() for l in f.readLines()])
    
def writeTextFile(filePath, seq, mode='w'):
    with open(filePath, mode) as f:
        f.write(seq + "\n")

def readFastaFile(filePath):
    with open(filePath, "r") as f:
        FASTAFile = [l.strip() for l in f.readLines()]

    FASTADict = {}
    FASTALabel = ""

    for line in FASTAFile:
        if '>' in line:
            FASTALabel = line
            FASTADict[FASTALabel] = ""
        else:
            FASTADict[FASTALabel] += line
            
    return FASTADict
