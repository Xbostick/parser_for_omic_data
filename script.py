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

def compute_part(df, id, num, name):
  part = df.partitions[num]
  part = part.loc[part.id == id]
  part = part.compute()
  part.to_csv("./files/BIGKEK_"+name+".csv",mode = 'a')
  print(num)

def parse(
        #TODO args
        id,
        filename
    ):

    df = dd.read_csv(   
                "./files/"+filename,
                sep = "\t", 
                names = ["chr", 'begin', 'end', 'id', 'score'],
                blocksize = '100mb'
                )
    open("BIGKEK_"+filename+".csv", mode = 'w').close()
    comp_list = []
    for i in range(df.npartitions):
        comp_list.append(que.submit(compute_part, df, id, i, filename)) 
    #TODO progress bar
    [future.result() for future in comp_list]
    

    

if __name__ == '__main__':
    que = Client(n_workers=NCORES, threads_per_worker=NWORKERS)
    args = cmd_line.parse_args()
    print(args)
    print(args.id)
    

    if args.assembly == None:
        for file in os.listdir(directory):
            filename = os.fsdecode(file)
            parse(args.id,filename)
            que.shutdown()
    else:
        #TODO list for assembly markers
        pass
    
    
