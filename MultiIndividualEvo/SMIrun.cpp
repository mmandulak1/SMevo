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

int GENERATION_MAX = 100000;
bool PRINT_GENERATIONS = false;


class Pairing{
    public:
        int left;
        int right;
        bool checked = false;
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
Individual mutation(Individual indiv, int N){
    Individual tempIndiv(indiv, N);
    int LoR, coord1, coord2, temp;
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
        temp = tempIndiv.matchingPairs[coord1].left;
        tempIndiv.matchingPairs[coord1].left = tempIndiv.matchingPairs[coord2].left;
        tempIndiv.matchingPairs[coord2].left = temp;
    }
    else{
        temp = tempIndiv.matchingPairs[coord1].right;
        tempIndiv.matchingPairs[coord1].right = tempIndiv.matchingPairs[coord2].right;
        tempIndiv.matchingPairs[coord2].right = temp;
    }
    //cout << *indiv;
    //cout << *tempIndiv;
    return tempIndiv;
}
float raiseTemp(float temper)
{
    return (temper*100)/90;
}
void resetMatrixCheck(vector<vector<Pairing>> &matrix, int N){
    for(int i=0;i<N;i++)
        for(int j=0;j<N;j++)
            matrix[i][j].checked=false;
}

int calcFitness(Individual indiv,vector<vector<Pairing>> &matrix, int N){
    resetMatrixCheck(matrix,N);
    int fitness=(N*N)-N;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            if(indiv.matchingPairValues[i].left < matrix[i][j].left){
                fitness--;
                matrix[i][j].checked=true;
            }
        }
    }
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            for(int k=0;k<N;k++){
                if(indiv.matchingPairs[k].right==i)
                    if(indiv.matchingPairValues[k].right<matrix[j][i].right && !matrix[j][i].checked)
                    {
                        fitness--;
                    }
            }
        }
    }
    return fitness;
    
} 
Individual doGenTemper(Individual indiv,vector<vector<Pairing>> &matrix,float temper, int N)
{
    int newfit;
    float probability;
    while(true){
        if(indiv.fitness_val==0)
            return indiv;
        Individual tempIndiv = mutation(indiv, N);
        tempIndiv.sortMatchingPairs(N);
        tempIndiv.setMatchingPairValues(matrix, N);
        newfit = calcFitness(indiv, matrix, N);
        tempIndiv.fitness_val = newfit;
        if(newfit <= indiv.fitness_val){
            return tempIndiv;
        }
        else
        {
            probability=pow(2.7, -1*((newfit-indiv.fitness_val)/temper));
            if(probability > 0.5){
                return tempIndiv;
            }
            temper = raiseTemp(temper);
        }
    }
}

int main(int argc, char* argv[]){
    if(argc!=3)
    {
        cout << "Invalid command line arguments (filename, N)" << endl;
        return 1;
    }
    int N = stoi(argv[2]);
    MPI_Init(NULL,NULL);
    int world_size = 0; 
    int world_rank = 0;
    MPI_Comm row_comm, col_comm;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    

    //cout << "hello from proc " << world_rank << "/" << world_size << endl;
    srand(time(NULL)+world_rank);
    
    //Individual* indiv;
    ifstream inpfile;
    inpfile.open(argv[1]);
    vector<vector<Pairing>> matrix;
    InitializeMatrix(matrix, N);
    ReadPrefsFromFile(matrix,inpfile,N);
    //if(world_rank==0)
    //   PrintMatrix(matrix,N);
    Individual indiv(true,matrix, N);
    indiv.fitness_val = calcFitness(indiv, matrix, N);
    //ParaInfo nodeInfo(col_rank, row_rank, matrix);
    //setMatchingPairsPara(indiv, nodeInfo, row_comm, col_comm, row_rank, col_rank);
    
    /*if(world_rank==2){
        cout << nodeInfo;
    }*/
    //indiv->fitness_val = calcFitnessPara(nodeInfo);

    bool done = false;
    float temper = 0.01;
    bool firstSM = false;
    int generations = 0;
    int foundSM = 0;
    while(!done){
        //cout << *indiv << endl;
        if(indiv.fitness_val==0){
            foundSM += 1;
        }
        MPI_Allreduce(&foundSM, &foundSM, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if(foundSM !=0)
        {
	    if(PRINT_GENERATIONS){
	    	MPI_Allreduce(&generations, &generations, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
            	if(world_rank==0){
                	cout << generations;
                	cout << endl;
            	}
	    }
            //for(int i=0;i<N;i++)
            //    cout << "(" << indiv->matchingPairs[i].left << "," << indiv->matchingPairs[i].right << ") ";
            	MPI_Finalize();
            	vector<vector<Pairing>>().swap(matrix);
            	inpfile.close();
            	return 0;
        }
        
        if(generations == GENERATION_MAX-1)
        {
	    if(PRINT_GENERATIONS){
            	if(world_rank == 0){
                	cout << "-1";
                	cout << endl;
            	}
	    }
            //for(int i=0;i<N;i++)
            //    cout << indiv->matchingPairs[i].left << " " << indiv->matchingPairs[i].right << " ";
            MPI_Finalize();
            vector<vector<Pairing>>().swap(matrix);
            inpfile.close();
            return 0;
        }
        //if(world_rank==0)
        //    cout <<"Generation "<< generations << ":" <<  indiv;
        indiv = doGenTemper(indiv, matrix, temper, N);
	for(int i=0;i<N;i++)
		cout << indiv.matchingPairs[i].left << "," << indiv.matchingPairs[i].right << " ";
	cout << endl;
        generations++;
    }
    
    inpfile.close();
    vector<vector<Pairing>>().swap(matrix);
    MPI_Finalize();
    return 0;
}
