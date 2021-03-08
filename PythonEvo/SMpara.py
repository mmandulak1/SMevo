import numpy as np
import random
import copy
import sys
from mpi4py import MPI
POP_SIZE = 1
GENERATIONS = 10000

N = int(sys.argv[1])
matrix = np.empty((N,N), dtype=object)

class ParaInfo:
    def __init__(self,x,y,matching=False):
        self.pair = [x,y]
        self.matching = matching
        self.pairVal = [matrix[x][y].pair[0],matrix[x][y].pair[1]]
        self.rowMatchingVal = 0
        self.colMatchingVal = 0
        self.rowMatching = [x,0]
        self.colMatching = [0,y]

class Pairing:
    def __init__(self,x,y):
        self.pair = [x,y]

class Individual:
    def __init__(self, init=False):
        self.fitness_val = -1
        self.matchingPairs = np.empty(N, dtype=object)
        self.matchingPairValues = np.empty(N, dtype=object)
        for i in range(0,N):
            self.matchingPairs[i] = Pairing(-1,-1)
            self.matchingPairValues[i] = Pairing(-1,-1)
        if init:
            self.setRandomMatchingPairs()
            self.sort_matching_pairs()
            self.setMatchingPairValues()
            self.fitness_val = N * N
    def setRandomMatchingPairs(self):
        for i in range(0,2):
            value = 0
            while value != N:
                coord = random.randint(1,N)-1
                if self.matchingPairs[coord].pair[i] is -1:
                    self.matchingPairs[coord].pair[i] = value
                    value += 1
    def print_pairs(self):
        for i in range(0,N):
            print("(" + str(self.matchingPairs[i].pair[0]) + "," + str(self.matchingPairs[i].pair[1]))
        print(str(self.fitness_val))
    def print_pairVals(self):
        print("PairVals:")
        for i in range(0,N):
            print("(" + str(self.matchingPairValues[i].pair[0]) + "," + str(self.matchingPairValues[i].pair[1]))
        print("")
    def print_pairs_sorted(self):
        for i in range(0,N):
            for j in range(0,N):
                if self.matchingPairs[j].pair[0] == i:
                    print("(" + str(self.matchingPairs[j].pair[0]) + "," + str(self.matchingPairs[j].pair[1]))
        print(str(self.fitness_val))
    def sort_matching_pairs(self):
        for i in range(0,N):
            for j in range(0,N):
                if self.matchingPairs[j].pair[0] == i:
                    temp = self.matchingPairs[i]
                    self.matchingPairs[i] = self.matchingPairs[j]
                    self.matchingPairs[j] = temp
    def setMatchingPairValues(self):
        for i in range(0,N):
            self.matchingPairValues[i].pair[0] = matrix[self.matchingPairs[i].pair[0]][self.matchingPairs[i].pair[1]].pair[0]
            self.matchingPairValues[i].pair[1] = matrix[self.matchingPairs[i].pair[0]][self.matchingPairs[i].pair[1]].pair[1]


def initializeMatrix():
    for i in range(0,N):
        for j in range(0,N):
            matrix[i][j] = Pairing(-1,-1)

def initializeParaMatrix():
    for i in range(0,N):
        for j in range(0,N):
            matrix[i][j] = Pairing(-1,-1)
def printMatrix():
    print("")
    for i in range(0,N):
        for j in range(0,N):
            print("(" + str(matrix[i][j].pair[0]) + "," + str(matrix[i][j].pair[1]))
            #print("(" + str(matrix[i][j]),end=") ")
        print("")
    print("")

def printMatchingMatrix(Individual):
    print("")
    for i in range(0,N):
        for j in range(0,N):
            if Individual[i].pair[1] == j:
                print("1")
            else:
                print("0")
        print("")
    print("")

def findBest(population):
    lowI = population[0]
    lowFit = population[0].fitness_val
    for i in range(1, POP_SIZE):
        if population[i].fitness_val < lowFit:
            lowFit = population[i].fitness_val
            lowI = population[i]
    return lowI

def setRandomPrefs():
    for i in range(0,N):
        value = 1
        while value != N+1:
            coord = random.randint(1,N)-1
            if matrix[i][coord].pair[0] is -1:
                matrix[i][coord].pair[0] = value
                value += 1
        value = 1
        while value != N+1:
            coord = random.randint(1,N)-1
            if matrix[coord][i].pair[1] is -1:
                matrix[coord][i].pair[1] = value
                value += 1
def setPrefsFromFilePara(filename):
    inpfile = open(filename,"r")
    temp = ""
    count = 0
    row = 0
    col = 0
    while True:
        oneinp = inpfile.read(1)
        if oneinp == ',':
            matrix[row][col].pair[0] = int(temp)
            temp = ""
        elif oneinp == ' ':
            matrix[row][col].pair[1] = int(temp)
            temp = ""
            col += 1
        elif oneinp == '\n' or oneinp == '':
            matrix[row][col].pair[1] = int(temp)
            temp = ""
            row += 1
            col = 0
            count += 1
            if count == N:
                break
        else:
            temp = temp + oneinp
    inpfile.close()


def calc_fitness(Individual):
    stablePairs = 0
    checkedPairs = []
    for i in range(0,N):
        row=Individual.matchingPairs[i].pair[0]
        for j in range(0,N):
            if matrix[row][j].pair[0] >= Individual.matchingPairValues[i].pair[0]:
                stablePairs += 1
                #print(str([row,j]))
                checkedPairs.append([row,j])
    #print("done male")
    for i in range(0,N):
        col=Individual.matchingPairs[i].pair[1]
        for j in range(0,N):
            if matrix[j][col].pair[1] >= Individual.matchingPairValues[i].pair[1] and not findPairing(checkedPairs,[j,col]):
                #print(str([j,col]))
                stablePairs += 1
    return (N*N)-stablePairs

def calc_fitness_para(pairinfo, world):
    #print("ok")
    #Each matching pair broadcasts its left value to its row
    if pairinfo.rowMatchingVal > pairinfo.pairVal[0] and pairinfo.colMatchingVal > pairinfo.pairVal[1]:
        stability = 1
    else:
        stability = 0
    stability = world.reduce(stability,op=MPI.SUM,root=0)
    stability = world.bcast(stability,root=0)
    #if world.Get_rank() == 0:
        #print("hello" + str(stability))
    return stability

    #if col.Get_rank() == Individual.matchingPairs[]
    #row.Reduce([data,MPI.DOUBLE],[globaldata,MPI.DOUBLE], op=MPI.SUM,root=0)

def setMatchingPairsPara(Individual,pairinfo, world, row, col):

    if Individual.matchingPairs[col.Get_rank()].pair[1] == row.Get_rank():
        pairinfo.matching = True
        #print ("Set matching pair: " + str(col.Get_rank()) + "," + str(row.Get_rank()))
        matchingColRank = pairinfo.pair[0]
        matchingRowRank = pairinfo.pair[1]
        matchingColVal = pairinfo.pairVal[1]
        matchingRowVal = pairinfo.pairVal[0]
    else:
        matchingColRank = 0
        matchingRowRank = 0
        matchingColVal = 0
        matchingRowVal = 0
    matchingRowRank = row.allreduce(matchingRowRank, op=MPI.SUM)
    matchingColRank = col.allreduce(matchingColRank, op=MPI.SUM)
    matchingRowVal = row.allreduce(matchingRowVal, op=MPI.SUM)
    matchingColVal = col.allreduce(matchingColVal, op=MPI.SUM)
    pairinfo.rowMatching = [col.Get_rank(), matchingRowRank]
    pairinfo.colMatching = [matchingColRank, row.Get_rank()]
    pairinfo.colMatchingVal = matchingColVal
    pairinfo.rowMatchingVal = matchingRowVal

def findPairing(arr, pair):
    for i in range(0, len(arr)):
        if pair[0] == arr[i][0] and pair[1] == arr[i][1]:
            return True
    return False

def isMatchingPair(matchingPairs, pair):
    match = False
    for i in range(0,N):
        if matchingPairs[i].pair[0]==pair[0] and matchingPairs[i].pair[1]==pair[1]:
            match = True
    return match

def mutation(Individual):
    #print("Stuck on individual" + str(Individual))
    tempIndividual = copy.deepcopy(Individual)
    XY = random.randint(0,1)
    coord1 = random.randint(1,N)-1
    coord2 = random.randint(1,N)-1
    while True:
        if coord2 == coord1:
            coord2 = random.randint(1,N)-1
        else:
            break
    temp = tempIndividual.matchingPairs[coord1].pair[XY]
    tempIndividual.matchingPairs[coord1].pair[XY] = tempIndividual.matchingPairs[coord2].pair[XY]
    tempIndividual.matchingPairs[coord2].pair[XY] = temp
    Individual.setMatchingPairValues()
    #print("")
    tempIndividual.fitness_val = calc_fitness(tempIndividual)
    #tempIndividual.print_pairs()
    return tempIndividual

def mutationPara(Individual, world):
    if world.Get_rank() == 0:
        tempIndividual = copy.deepcopy(Individual)
        XY = random.randint(0,1)
        coord1 = random.randint(1,N)-1
        coord2 = random.randint(1,N)-1
        while True:
            if coord2 == coord1:
                coord2 = random.randint(1,N)-1
            else:
                break
        temp = tempIndividual.matchingPairs[coord1].pair[XY]
        tempIndividual.matchingPairs[coord1].pair[XY] = tempIndividual.matchingPairs[coord2].pair[XY]
        tempIndividual.matchingPairs[coord2].pair[XY] = temp
    else:
        tempIndividual = None
    tempIndividual = world.bcast(tempIndividual,root=0)
    return tempIndividual


def do_gen_temp_para(Individual, pairinfo, world, row, col, temper):
    while True:
        if Individual.fitness_val == 0:
            return Individual
        tempIndividual = mutationPara(Individual,world)
        tempIndividual.sort_matching_pairs()
        tempParaInfo = copy.deepcopy(pairinfo)
        setMatchingPairsPara(tempIndividual,tempParaInfo,world,row,col)
        newfit = calc_fitness_para(tempParaInfo,world)
        tempIndividual.fitness_val = newfit
        if newfit <= Individual.fitness_val:
            return tempIndividual
        else:
            probability = pow(2.7, -1*((tempIndividual.fitness_val-Individual.fitness_val)/temper))
            #print("Prob: " + str(probability) + " temper: " + str(temper))
            if probability > 0.5:
                return tempIndividual
            temper = raiseTemp(temper)



def do_gen(population):
    mutatedSM = []
    for i in range(0,POP_SIZE):
        while True:
            tempIndividual = mutation(population[i])
            newfit = tempIndividual.fitness_val
            if population[i].fitness_val == 0:
                mutatedSM.append(population[i])
                break
            elif newfit <= population[i].fitness_val:
                mutatedSM.append(tempIndividual)
                break
            else:
                mutatedSM.append(population[i])
                break
    return mutatedSM


def do_gen_temp(population, temper):
    mutatedSM = []
    for i in range(0,POP_SIZE):
        while True:
            tempIndividual = mutation(population[i])
            newfit = tempIndividual.fitness_val
            if population[i].fitness_val == 0:
                mutatedSM.append(population[i])
                break
            elif newfit <= population[i].fitness_val:
                mutatedSM.append(tempIndividual)
                break
            else:
                probability = pow(2.7, -1*((tempIndividual.fitness_val-population[i].fitness_val)/temper))
                #print("Prob: " + str(probability) + " temper: " + str(temper))
                if probability > 0.5:
                    #print("Accepted")
                    mutatedSM.append(tempIndividual)
                    break
                temper = raiseTemp(temper)
    return mutatedSM

def lowerTemp(temper):
    return (90 * temper)/100;

def raiseTemp(temper):
    return (temper*100)/90

def checkFirstSM(population, generation):
    for i in range(0,POP_SIZE):
        if population[i].fitness_val == 0:
            print("First stable matching from individual " + str(i) + " on generation " + str(generation))
            population[i].print_pairs_sorted()
            return True
    return False

def checkLastSM(population):
    last = True
    for i in range(0,POP_SIZE):
        if population[i].fitness_val != 0:
            last = False
    return last
