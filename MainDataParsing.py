import sys
import glob
import numpy as np

DataFilePath = "./Data"

def GetFileList():
    pat = DataFilePath + "/*.dat"
    files = glob.glob(pat)
    return(files)

def DataReader(fileName):
    with open(fileName, 'rb') as inputFile:
        lines = inputFile.readlines()
    lines = [line.decode('windows-1252').replace('\r\n','\n') for line in lines]
    return(lines)

def DataParserAmp(rawData, n):
    parameters = {item.replace(" ","").rstrip('\n').split('=')[0]:float(item.replace(" ","").rstrip('\n').split('=')[1]) for item in rawData[1:6+n]}
    tau_data = {item.replace(" ","").rstrip('\n').split(":")[0]:float(item.replace(" ","").rstrip('\n').rstrip("ns").split(":")[1]) for item in rawData[8+n].split(";")}
    chiSquare = float(rawData[7+n].replace(": ",":").rstrip("\n").split(":")[1])
    data = rawData[10+n:len(rawData)]
    dataOut = []
    for line in data:
        line = line.rstrip('\n')
        line = line.replace(" ","")
        dataOut.append(line)

    return(parameters, chiSquare, tau_data, dataOut)

def Grouper(DataObj):
    CellTypes = list(set([myData.cellLine for myData in DataObj]))
    Output = {item:[] for item in CellTypes}
    for cell in DataObj:
        Output[cell.cellLine].append(np.average([cell.tau_data["tau_amp"], cell.tau_data["tau_int"]]))

    print("Cell\tAverage\tStd")
    for item in Output:
        print("%s\t%s\t%s"%(item,np.average(Output[item]),np.std(Output[item])))

def writePhotons(DataObj):
    with open("PhotonOutput.txt",'w') as outputFile:
        for cell in DataObj:
            for line in cell.dataOut:
                if line.startswith("Time[ns]"):
                    pass
                else:
                    outputFile.write(cell.filename.rstrip('.dat') + '\t' + line + '\n')

class Data:

    def __init__(self, file):
        self.filename = file
        self.rawData = DataReader(self.filename)
        self.cellLine = "_".join(file.split("/")[len(file.split("/"))-1].split("_")[0:4])
        if "Ampl. 2" in ','.join(self.rawData):
            self.n = 2
            self.parameters, self.chiSquare, self.tau_data, self.dataOut = DataParserAmp(self.rawData, self.n)
        else:
            self.n = 0
            self.parameters, self.chiSquare, self.tau_data, self.dataOut = DataParserAmp(self.rawData, self.n)

if __name__=="__main__":
    DataFiles = GetFileList()

    ROIs = []
    for myFile in DataFiles:
        ROIs.append(Data(myFile))

    Grouper(ROIs)

    writePhotons(ROIs)
