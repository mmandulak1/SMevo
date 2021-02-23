import SMpara
import os
from mpi4py import MPI
import numpy as np

world_comm = MPI.COMM_WORLD
world_rank = world_comm.Get_rank()
world_size = world_comm.Get_size()
row_comm = MPI.Comm.Split(world_comm,world_rank/SMpara.N,world_rank)
row_rank = row_comm.Get_rank()
col_comm = MPI.Comm.Split(world_comm, world_rank%SMpara.N,world_rank)
col_rank = col_comm.Get_rank()

SMpara.initializeParaMatrix()
SMpara.setPrefsFromFilePara("preflist.txt")
if world_rank == 0:
    indiv = SMpara.Individual(True)
    indiv.print_pairs()
    SMpara.printMatrix()
else:
    indiv = None
indiv = world_comm.bcast(indiv, root=0)
nodeinfo = SMpara.ParaInfo(col_rank, row_rank)
SMpara.setMatchingPairsPara(indiv,nodeinfo, world_comm, row_comm, col_comm)
fitness = SMpara.calc_fitness_para(nodeinfo,world_comm)
indiv.fitness_val = fitness

done = False
temper = 0.01
firstSM = False
generations = 0
while not done:
    if indiv.fitness_val == 0:
        if world_rank == 0:
            print("Stable matching on generation " + str(generations))
            indiv.print_pairs()
        break
    if generations % 10 == 0 and world_rank == 0:
        print(str(generations) + "/" + str(SMpara.GENERATIONS))
    if generations == SMpara.GENERATIONS-1:
        if world_rank == 0:
            print("Max generation number reached")
            indiv.print_pairs()
        break
    indiv = SMpara.do_gen_temp_para(indiv,nodeinfo,world_comm,row_comm,col_comm, temper)
    generations += 1
