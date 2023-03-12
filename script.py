import argparse
import pandas as pd
import numpy as np
import dask.dataframe as dd
import time
import dask
from dask.distributed import Client
from dask.diagnostics import ProgressBar
from dask.distributed import progress

import multiprocessing
import os

#defines
FILE_PATH = "/mnt/storage/avoitetskii/chip_atlas/"
PRIVATE_PATH = "~/private_omicON.txt"
DOC_LIST = {
    "hg38"  :   "allPeaks_light.hg38.05.bed",
    "mm10"  :   "allPeaks_light.mm10.05.bed"
    }
directory = os.fsencode(FILE_PATH)
NCORES = 2
NWORKERS = 4

#cmd line module init
cmd_line = argparse.ArgumentParser(description='Script for parsing and saving omic data')
# add fields to parser
cmd_line.add_argument(
    '--id',
    '-i',
    type=str,
    default="SRX8026866",
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

cmd_line.add_argument(
    '--verbose',
    '-v',
    type=int,
    default=1,
    help='verbose operation'
)

def get_matching_experimnet_part(
            df, 
            num,
            options
        ):
    # num - chunk number
    # options - list of arguments
    part = df.partitions[num]
    
    for key in options.keys():
        if options[key]:
            tmp = options[key].split(',')
            part = part.loc[part[key].isin(tmp)]
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
                    FILE_PATH + filename,
                    sep = '\t', 
                    names = ['id', 'Genome assembly', 'Antigen class', 'Antigen', 'Cell type class', 'Cell type'],
                    usecols=range(6),
                    blocksize = '10mb'
                )
    if args.verbose:
        print("Find file " +  FILE_PATH + filename)

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

    path_2_sorted_file = FILE_PATH + "filtred_"+filename+".csv"
    process_list = []

    df = dd.read_csv(   
                FILE_PATH + filename,
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
    print(f"Your file creating. You can see progress here:\n http://{parse_private()['public']}:{parse_private()['port']}/status")
    a = [process.result() for process in process_list]
    progress(a)
    

def parse_private():
    d = {}
    with open(PRIVATE_PATH) as f:
        for line in f:
            (key, val) = line.split()
            d[str(key)] = val
    return(d)   

if __name__ == '__main__':
    
    que = Client(n_workers=NCORES, threads_per_worker=NWORKERS)
    args = cmd_line.parse_args()

    options = {
        #Parse arguments from cmd line to special dict
        "id"                :   args.id,
        "Genome assembly"   :   args.assembly,
        "Antigen class"     :   args.antigen_class,
        "Antigen"           :   args.antigen,
        "Cell type class"   :   args.cell_type,
        "Cell type'"        :   args.cell
    }

    print("Succes parse arguments!")

    matching_experiments = create_matching_expirement_list(que, "part.txt", options)
    if args.verbose:
        print(f"Was finded {len(matching_experiments)} results:\n " + str(matching_experiments))

    create_sorted_bed_file(que, args.file, matching_experiments)

    que.shutdown()
    
    
    
