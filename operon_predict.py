#////////////////////////////////////////////////////////////////////////////////////////////////////
#/// \file operon-predict.py
#/// \brief A python program built to find the operons from genomes
#/// 
#///
#//  Author: Divya Singhal
#////////////////////////////////////////////////////////////////////////////////////////////////////

#import the required packages
import csv
import numpy as np
import pandas as pan

#////////////////////////////////////////////////////////////////////////////////////////////////////
#/// \FileToDataframes
#/// \brief function is created to load the given data in text format and reading into dataframes.
#/// \returns value i.e. dataframe with start and stop location.
#////////////////////////////////////////////////////////////////////////////////////////////////////

def FileToDataframes(filenames):
    table = list()
    with open(filenames) as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t')
        for row in spamreader:
            table.append(row)

    tableICare = list()
    for i in range(3,len(table)):
        tableICare.append(table[i])
    
    E_Coli = pan.DataFrame(tableICare,columns=table[2])
    split_E_Coli = pan.DataFrame(E_Coli['Location'].str.split('.',0).tolist(),
                                   columns = ['start','garbage','stop'])
    temp_E_Coli  =split_E_Coli.join(E_Coli)
    E_Coli = temp_E_Coli.drop('garbage',1)
    return E_Coli

#////////////////////////////////////////////////////////////////////////////////////////////////////
#/// \Task2FiletoDataframe
#/// \brief function is created to load the given data in text format and reading into dataframes.
#/// \returns value i.e. a dataframe with mutiple columns of required inofrmation.
#////////////////////////////////////////////////////////////////////////////////////////////////////

def Task2FiletoDataframes(filename):
    table = list()
    with open(filename) as csvfile:
        spamreader = csv.reader(csvfile, delimiter='\t')
        for row in spamreader:
            table.append(row)

    tableICare = list()
    for i in range(1,len(table)):
        tableICare.append(table[i])
    
    Task_2= pan.DataFrame(tableICare,columns=['Contig','img','CDS','start','stop','garbage','Strand','garbage1','ID'])
    return Task_2

#////////////////////////////////////////////////////////////////////////////////////////////////////
#/// \HandleResult
#/// \brief function is created to index the input data.
#/// \returns value i.e. list of new final output
#////////////////////////////////////////////////////////////////////////////////////////////////////

#HandleResult function created to get the desired output in a particular format.
resultList = list()
def HandleResult(index,operon,distance,HandleTable,taskNum):
    if taskNum==1:
        resultList.append([	operon, 
    					index, 
    					HandleTable['Location'][index], 
                        distance, 
                        HandleTable['Strand'][index],
                        HandleTable['Gene'][index], 
                        HandleTable['Gene'][index+1]])
    elif taskNum==2:
        resultList.append([ operon, 
                        index, 
                        HandleTable['start'][index] +".." +HandleTable['stop'][index],
                        distance, 
                        HandleTable['Strand'][index],
                        HandleTable['Contig'][index],
                        HandleTable['Contig'][index+1]])
    return

#////////////////////////////////////////////////////////////////////////////////////////////////////
#/// \findOperons
#/// \brief function is created to predict the operons from the genomes.
#/// \returns value i.e. identified operons with their locations and genes.
#////////////////////////////////////////////////////////////////////////////////////////////////////

def findOperons(HandleTable,taskNum):
    rowNum = 0
    operon = 1
    global resultList
    resultList=[]
    while(rowNum < (len(HandleTable)-1)):
        flag_50P = False
    
        while((int(HandleTable['start'][rowNum+1])-int(HandleTable['stop'][rowNum])) < 50
            and HandleTable['Strand'][rowNum] == HandleTable['Strand'][rowNum+1]):
                distance = int(HandleTable['start'][rowNum+1])-int(HandleTable['stop'][rowNum])
                if flag_50P==False:
            	    HandleResult(rowNum,operon,distance,HandleTable,taskNum)
                HandleResult(rowNum+1,operon,distance,HandleTable,taskNum)
                flag_50P = True
                rowNum += 1
            
        if flag_50P==True:
            operon += 1
        
        rowNum += 1

    intDistance1 = pan.DataFrame(resultList,columns=["Operon", "RowNum", "Location", "Distance",
                                                "Strand", "Gene1", "Gene2"])
    return intDistance1

#////////////////////////////////////////////////////////////////////////////////////////////////////
#/// \utputTableToFile
#/// \brief function is created to save the result in csv format.
#/// \returns value i.e. csv file.
#////////////////////////////////////////////////////////////////////////////////////////////////////
#outputTableToFile function created to save the output data in a csv file format.
def outputTableToFile(result,filename):
    result.to_csv(filename,sep='\t',encoding='utf-8')
    return

print ("###Be patient doing homework...... processing 5 files..###")

df_file1=FileToDataframes("E_coli_K12_MG1655.ptt")
df_file2=FileToDataframes("Synechocystis_PCC6803_uid159873.ptt")
df_file3=FileToDataframes("B_subtilis_168.ptt")
df_file4=FileToDataframes("Halobacterium_NRC1.ptt")
df_file5=Task2FiletoDataframes("2088090036.gff")

task1=1
task2=2

result_file1=findOperons(df_file1,task1)
result_file2=findOperons(df_file2,task1)
result_file3=findOperons(df_file3,task1)
result_file4=findOperons(df_file4,task1)
result_file5=findOperons(df_file5,task2)

outputTableToFile(result_file1,"result_E_coli_K12_MG1655.txt")
outputTableToFile(result_file2,"result_Synechocystis_PCC6803_uid159873.txt")
outputTableToFile(result_file3,"result_B_subtilis_168.txt")
outputTableToFile(result_file4,"result_Halobacterium_NRC1.txt")
outputTableToFile(result_file5,"result_2088090036.txt")

print ("#####output files for task1(4-files) and task2(1-file) has been created####")
print ("###Please check the folder###")
