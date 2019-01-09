#include<iostream>
#include<stdlib.h>
#include<time.h>
#include<math.h>
using namespace std; 

const int MAX = 100; 
int n = 64;

void luDecomposition(float **mat, int n) 
{ 
	int k,i,j;
	for(k=0;k<n;k++){
		//Division step
		for(i=k+1;i<n;i++){
			mat[i][k]=mat[i][k]/mat[k][k];
		}
		//Elimination step
		for(i=k+1;i<n;i++){
			for(j=k+1;j<n;j++){
				mat[i][j]=mat[i][j]-mat[i][k]*mat[k][j];
			}
		}
	}
}

void display(float **mat, int n){
	int i,j;
	for (int i = 0; i < n; i++) { 
        	for (int j = 0; j < n; j++){ 
            		cout << mat[i][j] << "\t"; 
		}  
  	cout << endl; 
    } 
}

int main() 
{ 
	clock_t time_req = clock();
    /*float mat[][MAX] = {{ 2, -1, -2,2,3 ,4}, 
                        {-4, 6, 3 ,2,3,4}, 
                       	{-4, -2, 8 ,2,3,4},
                        {2,2,2,2,3,4},
						{3,3,3,3,3,4},
			            {4,4,4,4,4,4,}}; */
	int i , j;
	srand(time(NULL));
	float **mat = (float **)malloc(sizeof(float *)*n);
	for (i=0; i<n; i++) 
        mat[i] = (float *)malloc(n * sizeof(float));
    for(i=0;i<n;i++){
  		for(j=0;j<n;j++){
		  		mat[i][j] = rand()%10;
		  	}
		}
    luDecomposition(mat, n); 
    display(mat, n);
	time_req = (clock() - time_req)*pow(10,6);
	printf( "Using Sequential LU decomposition: %f microseconds\n",(float)time_req/CLOCKS_PER_SEC);
    return 0; 
} 
