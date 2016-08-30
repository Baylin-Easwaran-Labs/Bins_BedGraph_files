import pandas as pd
import numpy as np
import re
import os
import glob
import multiprocessing
import logging
import argparse
import ntpath

######################################################################

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.FileHandler('log_bedGraph_rebin.txt')
handler.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
handler.setFormatter(formatter)
    # add the handlers to the logger
logger.addHandler(handler)
######################################################################

######################################################################
def parserMain():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c','--cpus', help = "Number of CPUs", default = 1,  type = int)
    parser.add_argument('-i','--input', help = "Input Folder Name")
    parser.add_argument('-o', '--output', help = "Output Folder Name", default= 'Output')
    parser.add_argument('-s', '--Bin_size', help = "Bin size (Intiger)", default=200, type = int)
    parser.add_argument('-t', '--Name_trimmer', help = "If this number is for example 5, it removes 5 digit fron the end of input file name",default=0, type = int)  
    args = parser.parse_args()
        
    return args.cpus, args.input, args.output, args.Bin_size, args.Name_trimmer
######################################################################
def MakeDirectory(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
######################################################################
def SaveDFtoCSV(directory,output_name, df):
    MakeDirectory(directory)
    path = os.path.join(directory, str(output_name))
    df.to_csv(path, sep=',', index = False)
######################################################################
def FindAlltheFilesInaFolder(folder_name):  
    path = os.path.join (folder_name, '*')
    dir_list = glob.glob(path)
    new_list = []
    for name in dir_list:
        if os.path.isfile(name):
            new_list.append(name)
    return new_list
######################################################################
def MultiProcessingPreparationByChromosomesOneInput(df):
    chrom_list = list(set(df['CHR']))
    array_chrom = []
    array_chrom_data = []
    for chrom in chrom_list:
        if 'LongLink' not in chrom:
            df_c = df[df['CHR'] == chrom]
            array_chrom_data.append(df_c)
            array_chrom.append(chrom)
    return array_chrom, array_chrom_data

######################################################################
def MultiprocessingModuleTwoInputsOneStable(function, df_chrom_array, df_array1, bin_size):
    results = []
    results_async = []
    
    if len(df_array1) > 0  and len(df_chrom_array) > 0:

        for chrom, df1_c in zip( df_chrom_array, df_array1):
            res = pool.apply_async(function, (chrom, df1_c, bin_size))
            results_async.append(res)

        results=[r.get() for r in results_async]

    return results
######################################################################
def FinalDataFrameReconstruction(df_array):
    df_final = pd.DataFrame()
    for df1 in df_array:
        df_final = pd.concat([df_final, df1])
    print('Final DataFrame:\n{0}\n'.format(df_final.head()))
    return df_final
######################################################################
def MakingBins(chrom, df, bin_size):
    print(chrom)
    
    df_start = list(df['Start'])
    df_end = list(df['End'])
    df_value = list(df['Value'])
    
    array_b = []
    array_v = []
    
    for start,end, value in zip(df_start, df_end, df_value):
       
        loop_end = end//bin_size
        loop_start = start//bin_size
        
        dif_pos = loop_end - loop_start
    
        if dif_pos == 0:
            value1 = (end - start)*value
          
            if len(array_b) > 0 and array_b[-1] == loop_start:
                array_v[-1] = array_v[-1] + value1
            elif len(array_b)  == 0 or array_b[-1] != loop_start:
                array_v.append(value1)
                array_b.append(loop_start)
    
        elif dif_pos == 1:
            pos = loop_end * bin_size
            n1 = pos - start
            n2 = end - pos
            v1 = n1*value
            v2 = n2*value
            
            if len(array_b) > 0 and array_b[-1] == loop_start:
                array_v[-1] = array_v[-1] + v1
                array_v.append(v2)
                array_b.append(loop_end)
            else:
                array_v.append(v1)
                array_v.append(v2)
                array_b.append(loop_start)
                array_b.append(loop_end)
                
        elif dif_pos > 1:
            pos = (loop_start + 1)*bin_size
            n1 = pos - start
            v1 = n1*value
            #print(n1, value, v1)
            if len(array_b) > 0 and array_b[-1] == loop_start:
                array_v[-1] = array_v[-1] + v1
            else:
                array_v.append(v1)
                array_b.append(loop_start)
        
            i = loop_start + 1
            while i < loop_end:
                #print(value)
                array_v.append(bin_size*(value))
                array_b.append(i)
                i = i+1
        
            pos1 = loop_end * bin_size
            n2 = end - pos1
            v2 = n2 * value
            #print(n2, value, v2)
            array_b.append(loop_end)
            array_v.append(v2)
    
    CHR_a = [chrom for x in range(len(array_v))]
    df_v = pd.DataFrame(array_v)
    df_b = pd.DataFrame(array_b)
    df_c = pd.DataFrame(CHR_a)
    df_final = pd.concat([df_c, df_b, df_v], axis = 1)
    df_final.columns = ['CHR','Bin','Value']
    
    return df_final
############################################################################
if __name__ == '__main__':

    number_of_cpus, folder_name, output, Bin_size, Name_trimmer = parserMain()
    print('\nNumber of CPUs: {0}\n'.format(number_of_cpus))
    logger.info('Number of CPUs: {0}'.format(number_of_cpus))
    pool = multiprocessing.Pool(number_of_cpus)

    ##############################################################
    file_paths = FindAlltheFilesInaFolder(folder_name)

    for tarfile in file_paths:
        print('######################################\n')
        print('\nFile: {0} \n\n'.format(tarfile))
        logger.info(tarfile)
        name = ntpath.basename(tarfile)
        name1 = name[:-Name_trimmer]
        name2 = name1+'.csv'
        #print(name2)
        
        try:
            if tarfile[:-2] == 'gz':
                df1 = pd.read_csv(tarfile, compression='gzip', sep='\t',  header=None)
            else:
                df1 = pd.read_csv(tarfile, sep='\t',  header=None)
        except:
            print('Error.. Skiping file: {0}'.format(tarfile))
            logger.info('----Error----')
            logger.info(tarfile)
            logger.info('---Error---')
            continue

        df1.columns = ['CHR','Start','End','Value']
        print(df1.head())
        print('\nStart Analysis..\n')
        print('\nMultiprocessing Preparation ..\n')
        array_chrom, array_chrom_data = MultiProcessingPreparationByChromosomesOneInput(df1)
        print('\nMain Module ..\n')
        results = MultiprocessingModuleTwoInputsOneStable(MakingBins, array_chrom, array_chrom_data, Bin_size)
        print('\nDataFrame Reconstruction .. \n')
        df_final = FinalDataFrameReconstruction(results)
        print('\nDivide by Bin Size..\n')
        df_final['test'] = df_final['Value']/Bin_size
        df_final1 = df_final.drop('Value', axis = 1)
        df_final1.columns = ['CHR', 'Bin', 'Value']

        print('\nSaving ..\n')
        SaveDFtoCSV(output,name2, df_final1)
        print('Next')
        logger.info('Done')
    print('Done')


