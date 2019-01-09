#include<iostream>
#include<omp.h>
#include<cmath>
#include<stdlib.h>
#include<time.h>
using namespace std; 

int producers_done[100];
const int MAX = 100;
int n = 6;
int p = 6;//Number of threads
//Defining structures 
struct record_s {
	int val;
	int prod;
	struct record_s* next_p;
};
struct buf_list {
	struct record_s* head_p;
	struct record_s* tail_p;
};

struct buf_list buff[100];

struct record_s* Create_record(int my_rank, int data){
	struct record_s* h =new record_s;
	h->val=data;
	h->prod=my_rank;
	h->next_p=NULL;
	return h;
}

void Enqueue(int my_rank, struct record_s* rec_p){
	if (buff[my_rank].tail_p == NULL) {//When zero records in queue
		buff[my_rank].head_p = rec_p;
		buff[my_rank].tail_p = rec_p;
	} 
	else {
		buff[my_rank].tail_p->next_p = rec_p;
		buff[my_rank].tail_p = rec_p;
	}
}

void Put(int my_rank, int data) {
	struct record_s* rec_p;
	rec_p = Create_record(my_rank, data);
	#pragma omp critical(queue)
	{
	Enqueue(my_rank, rec_p);
	}
	#pragma omp critical(done)
	producers_done[my_rank]++;
}


struct record_s* Dequeue(int myrank) {
	struct record_s* rec_p;
	if (buff[myrank].head_p == NULL) {//No record in queue
		return NULL;
	}
	else if (buff[myrank].head_p == buff[myrank].tail_p) {//One record in queue
		rec_p = buff[myrank].head_p;
		buff[myrank].head_p = buff[myrank].tail_p = NULL;
	} 
	else {//Multiple record in queue
		rec_p = buff[myrank].head_p;
		buff[myrank].head_p = buff[myrank].head_p->next_p;
	}
	return rec_p;
}
 

int Get(int my_rank) {
	struct record_s* rec_p;
	int data;
	while (producers_done[my_rank] < 1 || buff[my_rank].head_p != NULL) {
		#pragma omp critical (queue)
		{
		rec_p = Dequeue(my_rank);
		}
		if (rec_p != NULL) {
			data = rec_p->val;
			free(rec_p);
			return data;
		}
	}
}


  
void luDecomposition(float a[][MAX], int n) 
{ 
	int k,i,j,row,bsize;
    omp_set_num_threads(p);
	#pragma omp parallel private(k, i, j, row) shared(a) 
	{
		int my_rank = omp_get_thread_num();
	        bsize =ceil(n*1.0/p);
		//printf("yes %d\n",bsize);
		if (my_rank != 0) {
			for (k = (my_rank * bsize) ; k<(my_rank * bsize) + bsize && k<n ;k++) {
				row = Get(my_rank);
				Put(my_rank + 1, row+bsize);

				//Division step
				for (i = row+1 ;i<n ;i++) {
					a[i][row] = a[i][row]/a[row][row];
				}
				
				//Elimination step
				for (i = row+1; i<n ;i++) {
					for (j = row+1; j< n ;j++) {
						a[i][j] -= (a[i][row] * a[row][j]);
					}
				}
			}
		}

		else {
			//printf("%d\n",my_rank);
			for (k = (my_rank * bsize); k<(my_rank * bsize) + bsize;k++) {
				Put(my_rank + 1, k+bsize);

				//Division step
				for (i = k + 1;i<n;i++) {
					a[i][k] = a[i][k]/a[k][k];
				}

				//Elimination step
				for (i = k + 1 ;i<n;i++) {
					for (j = k + 1; j<n;j++) {
						a[i][j] -= (a[i][k] * a[k][j]);
					}
				}
			}
		}
	}
}

void display(float a[][MAX], int n){
	int i1,j1;
	for (i1 = 0; i1 < n; i1++) { 
        	for (j1 = 0; j1 < n; j1++){ 
            		cout << a[i1][j1] << "\t"; 
		}  
  	cout << endl; 
    } 
}

int main() 
{ 
	clock_t time_req = clock();
	float mat[][MAX] = {{ 2, -1, -2,2,3 ,4}, 
                        {-4, 6, 3 ,2,3,4}, 
                       	{-4, -2, 8 ,2,3,4},
                        {2,2,2,2,3,4},
						{3,3,3,3,3,4},
			            {4,4,4,4,4,4,}}; 
	/*int i , j;
	srand(time(NULL));
	float **mat = (float **)malloc(sizeof(float *)*n);
	for (i=0; i<n; i++) 
	         mat[i] = (float *)malloc(n * sizeof(float));	
 	for(i=0;i<n;i++){
  		for(j=0;j<n;j++){
		  	mat[i][j] = rand()%10;
		  }
	  }*/
	int x=0;
	for(x=0;x<n+1;x++){
		producers_done[x]=0;
	}
	for(x=0;x<n+1;x++){
		buff[x].tail_p=NULL;
	}
    luDecomposition(mat, n); 
    display(mat, n);
	time_req = (clock() - time_req)*pow(10,6);
	printf( "Using pipelining LU decomposition: %f microseconds\n",(float)time_req/CLOCKS_PER_SEC);
    return 0; 
} 
