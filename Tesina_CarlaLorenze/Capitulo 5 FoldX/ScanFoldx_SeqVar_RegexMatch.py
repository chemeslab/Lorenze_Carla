#! /usr/bin/python3
import argparse
from datetime import date
from shutil import ReadError
import pandas as pd
from Bio.SeqUtils import seq1
from Bio import SeqIO
import os
import sys
import numpy as np
import re
# import subprocess
# from ast import BinOp
# from tracemalloc import start

# ###############################################################################
# # This script runs a FoldX Matrix on the input sequences.
# # Across all over the sequence. Returns the min value for streteches equal or longer than five residues.
# # The input is a multifasta file
# # and the matrix file and its id
# # 
# # The output is written into separate files in the specified output directory.
# ###############################################################################


def loadMatrix (fileIN:str) -> dict :
    """
    Reads the Matrix and loads it to a Dictionary
    :param FileIN: Location of the matrix file to load
    :return: Dictionary in a form of {Position -> { Amino_Acid }}
    """
    print(f'Loading Matrix...')

    try:
        fh = open(fileIN,"r")
    except ReadError as e:
        print(f'Cannot read the file {fileIN}\n')
        print(e)
   
    df = pd.read_csv(fh,sep="\t")
    fh.close()
    df = df.iloc[0:,1:21]
    df.columns = [seq1(res) for res in df.columns[0:21]]
    matrix = df.T.to_dict()
    
    return(matrix)

def scanningSequences(matrix,id,fileIN,subSeqLength,regexMat):

    #print("Scanning sequences...")    
    records = list(SeqIO.parse(fileIN, "fasta"))
    nSeq = len(records)
    
    counter = 1
    winSize = len(matrix)
   # print(f'Matrix Size: {winSize}')
    results=pd.DataFrame()

    for rec in records:
       # print(f'{counter} de {nSeq} secuencias')
        counter +=1

        seqID = rec.id  # first record
        sequence = rec.seq  # last record

        sequences = []
        starts = list()
        ends = list()
        values = []

#        for i in range(0,len(sequence)+winSize-1):

        for i in range(winSize-regexMat,len(sequence)+winSize-1-(regexMat+subSeqLength-2)):

            subSeq = [pos for pos in range(i-winSize,i,1) if pos >= 0 and pos<len(sequence)]
            #print(subSeq)

            #subSeq = [pos for pos in range(i-winSize+1,i+1,1) if pos >= 0 and pos<len(sequence)]
            ss = min(subSeq)
            ee = max(subSeq)+1
            sequences.append(str(sequence[ss:ee]))
            starts.append(ss+1)
            ends.append(ee)

            if(len(subSeq)<winSize and ee<len(sequence)):
                iniMat = winSize-len(subSeq)
            else:
                iniMat = 0
            
            finMat = iniMat+len(subSeq)

            posMat = list(range(iniMat,finMat,1))

            vals = []

            for k in range(0,len(subSeq)): # Calcula el valor de FoldX para esa secuencia de péptido 
                res = sequence[subSeq[k]]
                if res in matrix[posMat[k]].keys():
                    vals.append(matrix[posMat[k]][res])
                else:
                    vals.append(matrix[posMat[k]]["A"])

            values.append(round(sum(vals),2))

        data = {'Start':starts,'End':ends,'FoldX':values,'Sub-Sequence':sequences}
        data = pd.DataFrame(data)
        data['SeqID'] = seqID
        data['Sequence'] = str(sequence)
        data['Min'] = False
        data['Matrix'] = id
        data['FoldXperRes'] = 0.0
        data['RelativeMin'] = False

        for n,sec in enumerate(data['Sub-Sequence']):
            data.loc[n,'FoldXperRes'] = round(data.loc[n,'FoldX']/len(sec),2)
        
        vipCol = f'MoreThan{subSeqLength}'

        data[vipCol] = data['Sub-Sequence'].apply(lambda x: len(x)>=subSeqLength)               

        indMin = data.loc[data[vipCol] == True,'FoldX'].nsmallest(2).index.tolist()
   
        data.at[indMin[0],'Min'] = True
        data.at[indMin[1],'Min'] = True

        minPerRes = data.loc[data[vipCol] == True,'FoldXperRes'].nsmallest(2).index.tolist()
        # if len(minPerRes) != 2:
        #     print(minPerRes)
        #     exit()
        data.at[minPerRes[0],'RelativeMin'] = True
        data.at[minPerRes[1],'RelativeMin'] = True
        
        data = data[['SeqID', 'Sequence', 'Matrix','Start', 'End', 'FoldX','Min','Sub-Sequence','FoldXperRes','RelativeMin',vipCol]]

        results = pd.concat([results,data])

    return results


def searchMotif(df:pd.DataFrame,matLen:int,regex:str, regexMat:int) -> pd.DataFrame:
    """
    Looks for the motif of interest in the data frame.
    :param df: Data frame with subsequences and values obtained from the matrix.
    :param matLen: Int size of the matrix used to scan the sequences.
    :param regex: regular expression that represents the motif
    :param regexMat: Position in the matrix where the regex match.
    :return: input DataFrame with the positions of the canonical motif and other added.
    """
    df['RegexMatch'] = False
    df['RegexPattern'] = ''

    for idx, subSeq in enumerate(df['Sub-Sequence'].values):
        #print("Entre al For")

        start = df.loc[idx,'Start']
        end = df.loc[idx,'End']

        if end>(matLen-regexMat) and start==1:
            regex = '.'+ regex

        df.loc[idx,'RegexPattern'] = regex
        pattern = re.compile(pattern=regex) #creo pattern como objeto, que es la expresion regular que debe matchear con la regex que le doy.

        
        m = pattern.match(string=subSeq) # busca el match entre la subeq y el pattern que cree mas arriba

        if m:
            df.iloc[idx,df.columns.get_loc('RegexMatch')] = True

    return(df)


def saveFiles(suffix:str,day:str,seqID:str,matID:str,matFile:str,seqFile:str,data:pd.DataFrame):

    """
    Takes the data and saves it to the output files.
    :param suffix: String defined by the user to identify the job.
    :param day: date
    :param seqID: String that states the name of the sequence
    :param matID: String to identify the matrix used to scan the sequences
    :param matFile: String matrix file used to scan the sequences
    :param seqFile: String fasta file with the sequences to scan
    :param data: Dataframe to save
    :param min: min or max. It's going to search for a min value or a max value.
    :return: None
    """
    fileOut = f'{day}_{matID}_{seqID}_{suffix}.tsv'
    print(fileOut)

    if os.path.exists(fileOut):
        print(f'The file {fileOut} already exists\n')
        #exit()

    #print('Writing file...')
    with open(fileOut, 'w') as f:
        f.write(f'# Script: {sys.argv[0]}\n')
        f.write(f'# Sequence Files: {seqFile}\n')
        f.write(f'# Matrix File: {matFile}\n')
        f.write(f'# This file contains the result of the sequence {seqID} scanned with the matrix {matID}\n')
        data.to_csv(f, index=None, sep='\t')



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This script scans sequences in a fasta file with matrices')
    
    parser.add_argument('--dirFasta',
                        metavar='.fasta',
                        type=str,
                        nargs=1,
                        help='Sequence Dir: FASTA directory',
                        dest='dirFasta',
                        action='store',
                        default='/home/juliana/Documents/2023_E1A_Linker_Evolution/05_Alignment_with_newSequences/20230426_E1A.Alignment.Full.142s.v013.fasta',
                        required=False)
    parser.add_argument('--dirMat',
                        metavar="directory",
                        type=str,
                        nargs=1,
                        help='Directory where the matrix is stored',
                        dest='matrixDir',
                        default='/home/juliana/matrices_FoldX',
                        action='store',
                        required=False)
    parser.add_argument('--matrix',
                        metavar='my_matrix',
                        type=str,
                        nargs=1,
                        help='Matrix file to use for scanning. The matrix is 21 columns where the first one is the positions and the other 20 are the amino acids. Each row is a position.',
                        dest='matrixFile',
                        default='2019_TOP12_sum_cryst_mat.txt',
                        action='store',
                        required=False)
    parser.add_argument('--matrixID',
                        metavar='my_matrix_id',
                        type=str,
                        nargs=1,
                        help='Id to identify the matrix. It is included in the name of the output files',
                        dest='idMat',
                        default='Test',
                        action='store',
                        required=False)
    parser.add_argument('--suffix',
                        metavar='descriptor',
                        type=str,
                        nargs=1,
                        help='A sufix for the job',
                        dest='suffix',
                        action='store',
                        default='FoldX',
                        required=False)
    parser.add_argument('--regularExpression',
                        metavar='Regex',
                        type=str,
                        nargs=1,
                        help='It is the regular expression to match',
                        dest='regex',
                        default='Regex',
                        action='store',
                        required=False)
    parser.add_argument('--regexInMat',
                    metavar='posInMat',
                    type=int,
                    nargs=1,
                    help='It is the position of the matrix where the regex starts',
                    dest='regexMat',
                    default='',
                    action='store',
                    required=False)
    parser.add_argument('--seqLength',
                metavar='subSeqLength',
                type=int,
                nargs=1,
                help='Threshold length of the sub-sequence that could contain the regex',
                dest='subSeqLength',
                default='',
                action='store',
                required=False)


    args = parser.parse_args()

    fastaDir = args.dirFasta[0]
    #fastaFile = args.fileFasta[0]
    matrixDir = args.matrixDir[0]
    matrixFile = args.matrixFile[0]
    slim = args.regex[0]
    regexMat = args.regexMat[0]-1
    subSeqLength = args.subSeqLength[0]


    try:
        matID = args.idMat[0]
    except TypeError:
        matID = ''

    try:
        suffix = args.suffix[0]
    except TypeError:
        suffix = ''
    
    today = date.today()
    hoy = today.strftime("%Y%m%d")

    if matrixDir.endswith('\\'):
        matrixFileFull = matrixDir + matrixFile
    else:
        matrixFileFull = matrixDir + '\\' + matrixFile

    matrix = loadMatrix(fileIN=matrixFileFull)
    winSize = len(matrix)

    for file in os.listdir(fastaDir):
        archivoFASTA = fastaDir + file

        for record in SeqIO.parse(archivoFASTA, 'fasta'):
            seqID = record.id
            sequence = str(record.seq)
            dfValues = scanningSequences(matrix=matrix,id=matID,fileIN = archivoFASTA,regexMat = regexMat,subSeqLength= subSeqLength)
            dfValues = searchMotif(df=dfValues,matLen=winSize,regex = slim, regexMat = regexMat)
            saveFiles(suffix=suffix,day=hoy,seqID=seqID,matID=matID,matFile=matrixFileFull,seqFile=archivoFASTA,data=dfValues)


