import sys
import glob
import numpy as np
from scipy import stats

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
    # print(CellTypes)
    Output = {item:[] for item in CellTypes}
    for cell in DataObj:
        Output[cell.cellLine].append(np.average([cell.tau_data["tau_amp"], cell.tau_data["tau_int"]]))

    print("\nCell\tAverage\tStd")
    for item in Output:
        print("%s\t%s\t%s"%(item,np.average(Output[item]),np.std(Output[item])))

    Output = {"SI_BS":[],"WT_BS":[],"WT_BS_10uM_ola_cath":[],"SI_BS_10uM_ola_cath":[],"SI_BS_H2O2_cath":[], "WT_BS_H2O2_cath":[]}
    print("\nROI\ttau")
    with open("tau_distribution.txt", "w") as outTau:
        for cell in DataObj:
            if "BS" in cell.cellLine:
                print("%s\t%s"%(cell.cellLine, np.average([cell.tau_data["tau_amp"], cell.tau_data["tau_int"]])))
                outTau.write(cell.cellLine + "\t" + repr(np.average([cell.tau_data["tau_amp"], cell.tau_data["tau_int"]])) + "\n")
                Output[cell.cellLine].append(np.average([cell.tau_data["tau_amp"], cell.tau_data["tau_int"]]))

    print("\ntau t-tests...")
    print("SI_BS vs. WT_BS")
    print(stats.ttest_ind(Output["SI_BS"],Output["WT_BS"]))
    print("SI_BS_10uM_ola_cath vs. WT_BS_10uM_ola_cath")
    print(stats.ttest_ind(Output["SI_BS_10uM_ola_cath"], Output["WT_BS_10uM_ola_cath"]))
    print("SI_BS_H2O2_cath vs. WT_BS_H2O2_cath")
    print(stats.ttest_ind(Output["SI_BS_H2O2_cath"], Output["WT_BS_H2O2_cath"]))

    Output = {"SI_BS": [], "WT_BS": [], "WT_BS_10uM_ola_cath": [], "SI_BS_10uM_ola_cath": [], "SI_BS_H2O2_cath": [], "WT_BS_H2O2_cath":[]}
    print("\nROI\tfc")
    with open("fc_distribution.txt","w") as outFc:
        for cell in DataObj:
            if "BS" in cell.cellLine:
                fc = cell.parameters["Ampl.1"] / (cell.parameters["Ampl.1"] + cell.parameters["Ampl.2"])
                print("%s\t%s"%(cell.cellLine, fc))
                outFc.write(cell.cellLine+"\t"+repr(fc)+"\n")
                Output[cell.cellLine].append(fc)

    print("\nfc t-tests...")
    print("SI_BS vs. WT_BS")
    print(stats.ttest_ind(Output["SI_BS"], Output["WT_BS"]))
    print("SI_BS_10uM_ola_cath vs. WT_BS_10uM_ola_cath")
    print(stats.ttest_ind(Output["SI_BS_10uM_ola_cath"], Output["WT_BS_10uM_ola_cath"]))
    print("SI_BS_H2O2_cath vs. WT_BS_H2O2_cath")
    print(stats.ttest_ind(Output["SI_BS_H2O2_cath"], Output["WT_BS_H2O2_cath"]))

class Data:

    def __init__(self, file):
        self.filename = file
        self.rawData = DataReader(self.filename)
        self.cellLine = file.split("/")[len(file.split("/"))-1].split("_cell")[0]
        print(self.cellLine)
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

    # writePhotons(ROIs)
