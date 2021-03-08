#include <iostream>
#include <cmath>
#include <vector>
#include <time.h>
#include <stdlib.h>
using namespace std;

int main(int argc, char* argv[])
{
    if(argc!=2)
    {
        cout << "Please enter a value for n" << endl;
        return 0;
    }
    int N = atoi(argv[1]);
    int leftMatrix[N][N];
    int rightMatrix[N][N];
    int numberL;
    int numberR;
    bool rFlag = false;
    bool rerollFlag = true;
    srand(time(NULL));
    numberL = rand()%N;

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            rightMatrix[i][j]=-1;
            leftMatrix[i][j]=-1;
        }
    }

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            while(rerollFlag)
            {
                numberR = rand()%N+1;
                for(int k=0;k<N;k++)
                {
                    if(numberR == rightMatrix[i][k])
                        rFlag=true;
                }
                if(!rFlag)
                {
                    rightMatrix[i][j]=numberR;
                    rerollFlag=false;
                }
                rFlag = false;
            }
            rerollFlag = true;
            while(rerollFlag)
            {
                numberL = rand()%N+1;
                for(int k=0;k<N;k++)
                {
                    if(numberL == leftMatrix[k][j])
                        rFlag=true;
                }
                if(!rFlag)
                {
                    leftMatrix[i][j]=numberL;
                    rerollFlag=false;
                }
                rFlag = false;
            }
            rerollFlag = true;
        }
    }
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            cout << rightMatrix[i][j] << "," << leftMatrix[i][j] << " ";
        }
        cout << endl;
    }

    return 0;
}
