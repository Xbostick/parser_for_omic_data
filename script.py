import argparse
import pandas as pd
import numpy as np
import dask.dataframe as dd
import time
import dask
from dask.distributed import Client
import multiprocessing
import os

#defines
PATH = "./files"
directory = os.fsencode(PATH)
NCORES = 2
NWORKERS = 4

        #cmd line module init
cmd_line = argparse.ArgumentParser(description='Script for parsing and saving omic data')
# add fields to parser
cmd_line.add_argument(
    '--id',
    '-i',
    type=str,
    default=None,
    help='Experimental ID'
)

cmd_line.add_argument(
    '--assembly',
    '-g',
    type=str,
    default=None,
    help='Genome assembly'
)

cmd_line.add_argument(
    '--antigen_class',
    '-b',
    type=str,
    default=None,
    help='Antigen class'
)

cmd_line.add_argument(
    '--antigen',
    '-a',
    type=str,
    default=None,
    help='Antigen'
)

cmd_line.add_argument(
    '--cell_type',
    '-t',
    type=str,
    default=None,
    help='Cell type class'
)

cmd_line.add_argument(
    '--cell',
    '-c',
    type=str,
    default=None,
    help='Cell type'
)

cmd_line.add_argument(
    '--file',
    '-f',
    type=str,
    default="mm9.50.bed",
    help='File name from documentation'
)

def get_matching_experimnet_part(
            df, 
            num,
            options
        ):
    # num - chunk number
    # options - list of arguments
    part = df.partitions[num]
    
    if options[0] != '':
        tmp = options[0].split(',')
        part = part.loc[part['id'].isin(tmp)]
    if options[1] != '':
        tmp = options[1].split(',')
        part = part.loc[part['Genome assembly'].isin(tmp)]
    if options[2] != '':
        tmp = options[2].split(',')
        part = part.loc[part['Antigen class'].isin(tmp)]
    if options[3] != '':
        tmp = options[3].split(',')
        part = part.loc[part['Antigen'].isin(tmp)]
    if options[4] != '':
        tmp = options[4].split(',')
        part = part.loc[part['Cell type class'].isin(tmp)]
    if options[5] != '':
        tmp = options[5].split(',')
        part = part.loc[part['Cell type'].isin(tmp)]
    
    part = part.compute()
    return list(part.id)


def create_matching_expirement_list(
            que, 
            filename, 
            options
        ): #    Return expirement names list
    
    process_list = []
    matching_experiments = []
    
    df = dd.read_csv(
                    filename,
                    sep = ',', 
                    names = ['id', 'Genome assembly', 'Antigen class', 'Antigen', 'Cell type class', 'Cell type'],
                    blocksize = '10mb'
                )
    #TODO Можно сделать мапом. Будет красивше
    # process_list.map(             
    #         get_matching_experimnet_part, 
    #         range(df.npartitions), 
    #         [options] * df.npartitions
    #     )
    
    for part in range(df.npartitions):
       process_list.append(que.submit(
                                get_matching_experimnet_part, 
                                df, 
                                part, 
                                options
                                ))
    
    for process in process_list:
        matching_experiments.extend(process.result())

    return matching_experiments

def add_sorted_bed_2_file(
            filename,
            df,
            num,
            matching_experiments
        ):
    part = df.partitions[num]

    part = part.loc[part['id'].isin(matching_experiments)]
    part = part.compute()
    part.to_csv(filename, mode = 'a')
    return num

def create_sorted_bed_file(
        que,
        filename,
        matching_experiments
    ):

    path_2_sorted_file = PATH + "filtred_"+filename+".csv"
    process_list = []

    df = dd.read_csv(   
                PATH + filename,
                sep = "\t", 
                names = ["chr", 'begin', 'end', 'id', 'score'],
                blocksize = '100mb'
                )
    
    open(path_2_sorted_file, mode = 'w').close()  # Creating empty .csv for editing
    
    for part in range(df.npartitions):
        process_list.append(que.submit(
                                add_sorted_bed_2_file,  
                                path_2_sorted_file,
                                df,
                                part,
                                matching_experiments
                                )) 
    #TODO progress bar
    [process.result() for process in process_list]



    

if __name__ == '__main__':
    
    que = Client(n_workers=NCORES, threads_per_worker=NWORKERS)
    args = cmd_line.parse_args()

    options = {
        #Parse arguments from cmd line to special dict
        "id"            :   args.id,
        "assembly"      :   args.assembly,
        "antigen_class" :   args.antigen_class,
        "antigen"       :   args.antigen,
        "cell_type"     :   args.cell_type,
        "cell"          :   args.cell
    }

    matching_experiments = create_matching_expirement_list(que, args.file, options)
    create_sorted_bed_file(que, args.file, matching_experiments)

    que.shutdown()
    
    
    
