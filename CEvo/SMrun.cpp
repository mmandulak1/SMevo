#include <iostream>
#include <mpi.h>
#include <fstream>
#include <cmath>
#include <vector>
#include <time.h>
#include <stdlib.h>
#include <sstream>
#include <string>
using namespace std;

int GENERATION_MAX = 10000;



class Pairing{
    public:
        int left;
        int right;
        Pairing(int x, int y){left = x; right = y;}
        Pairing(){left=-1; right=-1;}
        Pairing(const Pairing &p){ left = p.left; right = p.right;}
};

class Individual{
    public:
        int fitness_val;
        vector<Pairing> matchingPairs;
        vector<Pairing> matchingPairValues;
        void setRandomMatchingPairs(int N)
        {
            int value = 0;
            int coord;
            for(int i=0;i<2;i++)
            {
                value = 0;
                while(value!=N)
                {
                    coord = rand()%N;
                    if(i==0){
                        if(this->matchingPairs[coord].left == -1)
                        {
                            this->matchingPairs[coord].left = value;
                            value++;
                        }
                    }
                    else{
                        if(this->matchingPairs[coord].right == -1)
                        {
                            this->matchingPairs[coord].right = value;
                            value++;
                        }
                    }
                }
            }
        }
        void sortMatchingPairs(int N){
            vector<Pairing> temp;
            for(int i=0;i<N;i++){
                for(int j=0;j<N;j++){
                    if(this->matchingPairs[j].left==i){
                        temp.push_back(this->matchingPairs[j]);
                    }
                }
            }
            this->matchingPairs = temp;
        }
        void setMatchingPairValues(vector<vector<Pairing>> &matrix, int N){
            for(int i=0;i<N;i++)
            {
                this->matchingPairValues[i] = matrix[this->matchingPairs[i].left][this->matchingPairs[i].right];
            }
        }
        Individual(bool init, vector< vector<Pairing> > &matrix, int N){
            fitness_val = N*N;
            Pairing temp;
            Pairing tempVal;
            for(int i=0;i<N;i++){
                temp = Pairing(-1,-1);
                tempVal = Pairing(-1,-1);
                matchingPairs.push_back(temp);
                matchingPairValues.push_back(tempVal);
            }
            if(init){
                this->setRandomMatchingPairs(N);
                this->sortMatchingPairs(N);
                this->setMatchingPairValues(matrix,N);
            }
        }
        Individual(const Individual &i, int N)
        {
            fitness_val = i.fitness_val;
            Pairing temp;
            for(int j=0;j<N;j++){
                temp = i.matchingPairs[j];
                matchingPairs.push_back(temp);
                temp = i.matchingPairValues[j];
                matchingPairValues.push_back(temp);
            }
        }
        Individual(int N){
            fitness_val = -1;
            Pairing temp;
            Pairing tempVal;
            for(int i=0;i<N;i++){
                temp = Pairing(-1,-1);
                tempVal = Pairing(-1,-1);
                matchingPairs.push_back(temp);
                matchingPairValues.push_back(tempVal);
            }
        }
        ~Individual(){
            vector<Pairing>().swap(matchingPairs);
            vector<Pairing>().swap(matchingPairValues);
        }
};
class ParaInfo{
    public:
        Pairing pair;
        Pairing pairVal;
        int rowMatchingVal;
        int colMatchingVal;
        int rowMatchingPos;
        int colMatchingPos;
    ParaInfo(int x, int y, vector<vector<Pairing>> matrix)
    {
        pair.left = x;
        pair.right = y;
        pairVal.left = matrix[x][y].left;
        pairVal.right = matrix[x][y].right;
    }
    ParaInfo(const ParaInfo &p){
        rowMatchingVal = p.rowMatchingVal;
        colMatchingVal = p.colMatchingVal;
        rowMatchingPos = p.rowMatchingPos;
        colMatchingPos = p.colMatchingPos;
        pair.left = p.pair.left;
        pair.right = p.pair.right;
        pairVal.left = p.pairVal.left;
        pairVal.right = p.pairVal.right;
    }
    
};
ostream& operator<<(ostream& os, const ParaInfo& p){
        return os << "My pair: (" << p.pair.left << "," << p.pair.right << "):" << p.pairVal.left << "," << p.pairVal.right << endl 
        << "Row matching: " << p.pair.left << "," <<p.rowMatchingPos << ":" << p.rowMatchingVal << endl 
        << "Col matching: " << p.colMatchingPos << "," << p.pair.right << ":" << p.colMatchingVal << endl;
    }
ostream& operator<<(ostream& os, const Individual& p){
        ostringstream oss;
        ostringstream oss2;
        for(int i=0;i<p.matchingPairs.size();i++)
        {
            oss << "(" << p.matchingPairs[i].left << "," << p.matchingPairs[i].right << ") ";
        }
        for(int i=0;i<p.matchingPairValues.size();i++)
        {
            oss2 << "(" << p.matchingPairValues[i].left << "," << p.matchingPairValues[i].right << ") ";
        }
        oss2 << p.fitness_val;
        string temp = oss.str();
        string temp2 = oss2.str();
        return os << temp << endl << temp2 << endl;
    }


void InitializeMatrix(vector<vector<Pairing>> &matrix, int N)
{
    vector<Pairing> v1;
    Pairing temp(-1,-1);
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            v1.push_back(temp);
        }
        matrix.push_back(v1);
    }
}
void PrintMatrix(vector<vector<Pairing>> &matrix, int N)
{
    cout << endl;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            cout << "(" << matrix[i][j].left << "," << matrix[i][j].right << ") ";
        }
        cout << endl;
    }
    cout << endl;
}
void ReadPrefsFromFile(vector<vector<Pairing>> &matrix, ifstream &inpfile, int N)
{
    string inp;
    int num;
    int rowI = 0;
    int colI = 0;
    char c;
    string value = "";
    while(getline(inpfile,inp))
    {
        for(int k=0;k<inp.length();k++)
        {
            c = inp[k];
            if(c!=' ' && c!= '\n' && c!=',')
                value += c;
            else if(c==',')
            {
                num = stoi(value);
                matrix[rowI][colI].left=num;
                value = "";
            }
            else if (c==' ' || c=='\n')
            {
                num = stoi(value);
                matrix[rowI][colI].right=num;
                value = "";
                colI++;
            }
        }
        rowI++;
        colI = 0;
    }
}
void ParseIndivBcast(Individual* indiv, int world_rank, int N)
{
    int temp, temp2, temp3, temp4;
    for(int i=0;i<N;i++)
    {
        temp = (*indiv).matchingPairs[i].left;
        //cout << "error 1" << endl;
        MPI_Bcast(&temp,1,MPI_INT,0,MPI_COMM_WORLD);
        temp2 = (*indiv).matchingPairs[i].right;
        //cout << "error 2" << endl;
        MPI_Bcast(&temp2,1,MPI_INT,0,MPI_COMM_WORLD);
        temp3 = (*indiv).matchingPairValues[i].left;
        //cout << "error 3" << endl;
        MPI_Bcast(&temp3,1,MPI_INT,0,MPI_COMM_WORLD);
        temp4 = (*indiv).matchingPairValues[i].right;
        MPI_Bcast(&temp4,1,MPI_INT,0,MPI_COMM_WORLD);
        if(world_rank!=0)
        {
            (*indiv).matchingPairs[i].left = temp;
            (*indiv).matchingPairs[i].right = temp2;
            (*indiv).matchingPairValues[i].left = temp3;
            (*indiv).matchingPairValues[i].right = temp4;

        }
    }
}

void setMatchingPairsPara(Individual *indiv, ParaInfo &nodeInfo, MPI_Comm row_comm, MPI_Comm col_comm, int row_rank, int col_rank)
{
    int matchingColRank=0, matchingRowRank=0, matchingColVal=0, matchingRowVal=0;
    if ((*indiv).matchingPairs[col_rank].right == row_rank)
    {
        matchingColRank = nodeInfo.pair.left;
        matchingRowRank = nodeInfo.pair.right;
        matchingColVal = nodeInfo.pairVal.right;
        matchingRowVal = nodeInfo.pairVal.left;
    }
    //int finalColRank, finalRowRank, finalRowVal, finalColVal = -1;
    MPI_Allreduce(&matchingColRank,&matchingColRank,1,MPI_INT,MPI_SUM,col_comm);
    MPI_Allreduce(&matchingRowRank,&matchingRowRank,1,MPI_INT,MPI_SUM,row_comm);
    MPI_Allreduce(&matchingRowVal,&matchingRowVal,1,MPI_INT,MPI_SUM,row_comm);
    MPI_Allreduce(&matchingColVal,&matchingColVal,1,MPI_INT,MPI_SUM,col_comm);
    nodeInfo.rowMatchingVal = matchingRowVal;
    nodeInfo.colMatchingVal = matchingColVal;
    nodeInfo.rowMatchingPos = matchingRowRank;
    nodeInfo.colMatchingPos = matchingColRank;
}
int calcFitnessPara(ParaInfo &nodeInfo)
{
    int stability = 0;
    if(nodeInfo.rowMatchingVal > nodeInfo.pairVal.left && nodeInfo.colMatchingVal > nodeInfo.pairVal.right)
        stability = 1;
    MPI_Allreduce(&stability,&stability,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    return stability;
}
Individual* mutationPara(Individual *indiv, int world_rank, int N){
    Individual* tempIndiv; 
    int LoR, coord1, coord2, temp;
    tempIndiv = new Individual(*indiv, N);
    if(world_rank==0){
        //cout << *indiv;
        //cout << *tempIndiv;
        LoR = rand()%2;
        coord1 = rand()%N;
        coord2 = rand()%N;
        while(true){
            if(coord2 == coord1)
                coord2 = rand()%N;
            else
                break;
        }
        if(LoR==0){
            temp = tempIndiv->matchingPairs[coord1].left;
            tempIndiv->matchingPairs[coord1].left = tempIndiv->matchingPairs[coord2].left;
            tempIndiv->matchingPairs[coord2].left = temp;
        }
        else{
            temp = tempIndiv->matchingPairs[coord1].right;
            tempIndiv->matchingPairs[coord1].right = tempIndiv->matchingPairs[coord2].right;
            tempIndiv->matchingPairs[coord2].right = temp;
        }
        //cout << *indiv;
        //cout << *tempIndiv;
    }   
    ParseIndivBcast(tempIndiv,world_rank,N);
    return tempIndiv;
}
float raiseTemp(float temper)
{
    return (temper*100)/90;
}
Individual* doGenTemperPara(Individual *indiv, ParaInfo &nodeInfo, int world_rank, int row_rank, int col_rank, MPI_Comm row_comm, MPI_Comm col_comm, float temper, int N)
{
    Individual* tempIndiv;
    ParaInfo* tempParaInfo;
    int newfit;
    float probability;
    while(true){
        if((*indiv).fitness_val==0)
            return indiv;
        tempIndiv = mutationPara(indiv, world_rank, N);
        (*tempIndiv).sortMatchingPairs(N);
        tempParaInfo = new ParaInfo(nodeInfo);
        setMatchingPairsPara(tempIndiv, *tempParaInfo, row_comm, col_comm, row_rank, col_rank);
        newfit = calcFitnessPara(*tempParaInfo);
        tempIndiv->fitness_val = newfit;
        if(newfit <= indiv->fitness_val){
            delete indiv;
            return tempIndiv;
        }
        else
        {
            probability=pow(2.7, -1*((newfit-indiv->fitness_val)/temper));
            if(probability > 0.5){
                delete indiv;
                return tempIndiv;
            }
            temper = raiseTemp(temper);
        }
    }
}


int main(int argc, char* argv[]){
    srand(time(NULL));
    
    if(argc!=3)
    {
        cout << "Invalid command line arguments (filename, N)" << endl;
        return 1;
    }
    int N = stoi(argv[2]);
    MPI_Init(NULL,NULL);
    int world_size, world_rank, row_rank, col_rank, row_color, col_color, row_size, col_size;
    MPI_Comm row_comm, col_comm;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    col_color = world_rank % N;
    row_color = world_rank / N;

    MPI_Comm_split(MPI_COMM_WORLD, col_color, world_rank, &col_comm);
    MPI_Comm_split(MPI_COMM_WORLD, row_color, world_rank, &row_comm);
    MPI_Comm_rank(col_comm, &col_rank);
    MPI_Comm_size(col_comm, &col_size);
    MPI_Comm_rank(row_comm, &row_rank);
    MPI_Comm_size(row_comm, &row_size);

    //cout << "hello from proc " << world_rank << "/" << world_size << endl;

    
    Individual* indiv;
    ifstream inpfile;
    inpfile.open(argv[1]);
    vector<vector<Pairing>> matrix;
    InitializeMatrix(matrix, N);
    ReadPrefsFromFile(matrix,inpfile,N);
    if(world_rank==0)
        PrintMatrix(matrix,N);
    if(world_rank==0){
        indiv = new Individual(true,matrix, N);
        /*for(int i=0;i<N;i++){
            cout << indiv->matchingPairs[i].left << " " << indiv->matchingPairs[i].right << ":" << indiv->matchingPairValues[i].left << " " << indiv->matchingPairValues[i].right << endl;
        }*/
    }
    else
        indiv = new Individual(false,matrix,N); 
        //for(int i=0;i<N;i++){
        //    cout << indiv.matchingPairs[i].left << " " << indiv.matchingPairs[i].right << ":" << indiv.matchingPairValues[i].left << " " << indiv.matchingPairValues[i].right << endl;
        //}
    ParseIndivBcast(indiv,world_rank,N);
    /*if(world_rank==2){
        for(int i=0;i<N;i++){
            cout << indiv->matchingPairs[i].left << " " << indiv->matchingPairs[i].right << ":" << indiv->matchingPairValues[i].left << " " << indiv->matchingPairValues[i].right << endl;
        }
    }*/
    ParaInfo nodeInfo(col_rank, row_rank, matrix);
    setMatchingPairsPara(indiv, nodeInfo, row_comm, col_comm, row_rank, col_rank);
    /*if(world_rank==2){
        cout << nodeInfo;
    }*/
    indiv->fitness_val = calcFitnessPara(nodeInfo);

    bool done = false;
    float temper = 0.01;
    bool firstSM = false;
    int generations = 0;
    while(!done){
        if(indiv->fitness_val==0)
        {
            if(world_rank==0)
            {
                cout << "Stable matching on generation " << generations << endl;
                for(int i=0;i<N;i++)
                    cout << "(" << indiv->matchingPairs[i].left << "," << indiv->matchingPairs[i].right << ") ";
                cout << endl;
            }
            MPI_Finalize();
            delete indiv;
            vector<vector<Pairing>>().swap(matrix);
            inpfile.close();
            return generations;
        }
        if(generations == GENERATION_MAX-1)
        {
            if(world_rank==0){
                cout << "Max generations reached" << endl;
                for(int i=0;i<N;i++)
                    cout << indiv->matchingPairs[i].left << " " << indiv->matchingPairs[i].right << " ";
                cout << endl;
            }
            MPI_Finalize();
            delete indiv;
            vector<vector<Pairing>>().swap(matrix);
            inpfile.close();
            return -1;
        }
        indiv = doGenTemperPara(indiv, nodeInfo,world_rank,row_rank,col_rank,row_comm,col_comm,temper, N);
        generations++;
    }
    
    //tempo = mutationPara(indiv, world_rank, N);
    //cout << "ok" << endl;
    inpfile.close();
    delete indiv;
    vector<vector<Pairing>>().swap(matrix);
    MPI_Finalize();
    return 0;
}
