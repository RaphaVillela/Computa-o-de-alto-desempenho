#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

#define ITERATIONS 100
#define MATRIX_MAX 2000

#define MATRIX_LENGHT1 100
#define MATRIX_LENGHT2 200
#define MATRIX_LENGHT3 500
#define MATRIX_LENGHT4 1000
#define MATRIX_LENGHT5 2000

#define DIVISION_PARTICLES1 1000
#define DIVISION_PARTICLES2 100
#define DIVISION_PARTICLES3 10

typedef struct
{
  int left;
  int right;
  int down;
  int up;
}Node;

void alocate_int (int *** matrix){
	
	*matrix = (int **)calloc(MATRIX_MAX, sizeof(int *));
	
	for(int i = 0; i < MATRIX_MAX; i++){
		
		(*matrix)[i] = (int *)calloc(MATRIX_MAX, sizeof(int));
	}
}

//Realease all memory space used
void realease_int (int *** matrix){	

	for(int i = 0; i < MATRIX_MAX; i++){
		
		free((*matrix)[i]);
	}
	free(*matrix);
}

void alocate_int_mpi (int *** matrix, int max){
	
	*matrix = (int **)calloc(MATRIX_MAX, sizeof(int *));
	
	for(int i = 0; i < max; i++){
		
		(*matrix)[i] = (int *)calloc(MATRIX_MAX, sizeof(int));
	}
}

//Realease all memory space used
void realease_int_mpi (int *** matrix, int max){	

	for(int i = 0; i < max; i++){
		
		free((*matrix)[i]);
	}
	free(*matrix);
}

//Alocate memory space for all analysis
void alocate (Node *** matrix){
	
	*matrix = (Node **)calloc(MATRIX_MAX, sizeof(Node *));
	
	for(int i = 0; i < MATRIX_MAX; i++){
		
		(*matrix)[i] = (Node *)calloc(MATRIX_MAX, sizeof(Node));
	}
}

//Realease all memory space used
void realease (Node ** matrix){	

	for(int i = 0; i < MATRIX_MAX; i++){
		free(matrix[i]);
	}
	free(matrix);
}

int get_chunk(int size, int lenght){
	return lenght/size;
}

void particle_move_mpi(Node** curr, Node** next, int max, int rank, int size){

	int ** moves = NULL;
	alocate_int_mpi(&moves, max);

	int extra = max % size;
	int chunk = get_chunk(size, max);

	int start_row =  rank * chunk + (rank - 1 < extra ? rank : extra);
   	int end_row = start_row + chunk + (rank < extra ? 1:0);

	for(int iterations = 0; iterations < ITERATIONS; iterations++){

		for(int i= start_row; i < end_row; i++){ //Put the moves on matrix moves 

			for(int j = 0; j < max; j = j + 4){
				//Verify the edges to change the direction
				if((j == 1) && (curr[i][j].left == 1)){
					curr[i][j].left = 0;
					curr[i][j].right = 1;
				}if((j == max-2) && (curr[i][j].right == 1)){
					curr[i][j].right = 0;
					curr[i][j].left = 1;
				}if((i == max-2) && (curr[i][j].down == 1)){
					curr[i][j].down = 0;
					curr[i][j].up = 1;
				}if((i == 1) && (curr[i][j].up == 1)){
					curr[i][j].up = 0;
					curr[i][j].down = 1;
				}

				if(curr[i][j].left == 1){//Left
					(moves[i][j-1])++;
				}if(curr[i][j].right == 1){//Right
					(moves[i][j+1])++;
				}if(curr[i][j].down == 1){//Down
					(moves[i+1][j])++;
				}if(curr[i][j].up == 1){//Up
					(moves[i-1][j])++;
				}
			}
		}

		for(int i=start_row; i < end_row; i++){//Verify the moves and modify next
			for(int j = 0; j < max; j++){

				if(moves[i][j] != 0){
					if(moves[i][j] == 1){ // 1 particule on node

						if((curr[i][j+1].left) == 1){
							next[i][j].left = 1;
						}
						else if((curr[i][j-1].right) == 1){
							next[i][j].right = 1;
						}
						else if((curr[i-1][j].down) == 1){
							next[i][j].down = 1;
						}
						else if((curr[i+1][j].up) == 1){
							next[i][j].up = 1;
						}
					}else if(moves[i][j] == 2){ //2 particules on node
						if((curr[i][j+1].left) == 1 && ((curr[i][j-1].right) == 1)){
							next[i][j].down = 1;
							next[i][j].up = 1;
						}else if((curr[i+1][j].up) == 1 && ((curr[i-1][j].down) == 1)){		
							next[i][j].left = 1;
							next[i][j].right = 1;
						}else{
							next[i][j].left = curr[i][j+1].left;
							next[i][j].right = curr[i][j-1].right;
							next[i][j].down = curr[i-1][j].down;
							next[i][j].up = curr[i+1][j].up;
						}
					}else if(moves[i][j] == 3){ // 3 particules on node
							next[i][j].left = curr[i][j+1].left;
							next[i][j].right = curr[i][j-1].right;
							next[i][j].down = curr[i-1][j].down;
							next[i][j].up = curr[i+1][j].up;
					}else if(moves[i][j] == 4){ // 4 particules on node
						next[i][j].left = 1; 
						next[i][j].right = 1;
						next[i][j].down = 1; 
						next[i][j].up = 1;
					}
				}
			}
		}

		for(int i = start_row; i < end_row; i++){ //Copy the next matrix to curr and put next to 0 again
			for(int j = 0; j < max; j++){

				 // Copy values from next to curr
				curr[i][j].left = next[i][j].left;
				curr[i][j].right = next[i][j].right;
				curr[i][j].down = next[i][j].down;
				curr[i][j].up = next[i][j].up;

				// Reset next values to 0
				next[i][j].left = 0;
				next[i][j].right = 0;
				next[i][j].down = 0;
				next[i][j].up = 0;
				moves[i][j] = 0;

				if( i == max-1){
					curr[i][j].left = 0;
					curr[i][j].right = 0;
					curr[i][j].down = 0;
					curr[i][j].up = 0;
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	realease_int_mpi(&moves, max);
}

void particle_move(Node** curr, Node** next, int max){

	int ** moves = NULL;
	alocate_int(&moves);

	for(int iterations = 0; iterations < ITERATIONS; iterations++){
		for(int i=0; i < max; i++){ //Put the moves on matrix moves 
			for(int j = 0; j < max; j++){
				//Verify the edges to change the direction
				if((j == 1) && (curr[i][j].left == 1)){
					curr[i][j].left = 0;
					curr[i][j].right = 1;
				}if((j == max-2) && (curr[i][j].right == 1)){
					curr[i][j].right = 0;
					curr[i][j].left = 1;
				}if((i == max-2) && (curr[i][j].down == 1)){
					curr[i][j].down = 0;
					curr[i][j].up = 1;
				}if((i == 1) && (curr[i][j].up == 1)){
					curr[i][j].up = 0;
					curr[i][j].down = 1;
				}

				if(curr[i][j].left == 1){//Left
					(moves[i][j-1])++;
				}if(curr[i][j].right == 1){//Right
					(moves[i][j+1])++;
				}if(curr[i][j].down == 1){//Down
					(moves[i+1][j])++;
				}if(curr[i][j].up == 1){//Up
					(moves[i-1][j])++;
				}
			}
		}

		for(int i=0; i < max; i++){//Verify the moves and modify next
			for(int j = 0; j < max; j++){
				
					if(moves[i][j] != 0){
						if(moves[i][j] == 1){ // 1 particule on node

						if((curr[i][j+1].left) == 1){
							next[i][j].left = 1;
						}
						else if((curr[i][j-1].right) == 1){
							next[i][j].right = 1;
						}
						else if((curr[i-1][j].down) == 1){
							next[i][j].down = 1;
						}
						else if((curr[i+1][j].up) == 1){
							next[i][j].up = 1;
						}
					}else if(moves[i][j] == 2){ //2 particules on node
						if((curr[i][j+1].left) == 1 && ((curr[i][j-1].right) == 1)){
							next[i][j].down = 1;
							next[i][j].up = 1;
						}else if((curr[i+1][j].up) == 1 && ((curr[i-1][j].down) == 1)){		
							next[i][j].left = 1;
							next[i][j].right = 1;
						}else{
							next[i][j].left = curr[i][j+1].left;
							next[i][j].right = curr[i][j-1].right;
							next[i][j].down = curr[i-1][j].down;
							next[i][j].up = curr[i+1][j].up;
						}
					}else if(moves[i][j] == 3){ // 3 particules on node
							next[i][j].left = curr[i][j+1].left;
							next[i][j].right = curr[i][j-1].right;
							next[i][j].down = curr[i-1][j].down;
							next[i][j].up = curr[i+1][j].up;
					}else if(moves[i][j] == 4){ // 4 particules on node
						next[i][j].left = 1; 
						next[i][j].right = 1;
						next[i][j].down = 1; 
						next[i][j].up = 1;
					}
				}
			}
		}

		for(int i=0; i < max; i++){ //Copy the next matrix to curr and put next to 0 again
			for(int j = 0; j < max; j++){

				 // Copy values from next to curr
				curr[i][j].left = next[i][j].left;
				curr[i][j].right = next[i][j].right;
				curr[i][j].down = next[i][j].down;
				curr[i][j].up = next[i][j].up;

				// Reset next values to 0
				next[i][j].left = 0;
				next[i][j].right = 0;
				next[i][j].down = 0;
				next[i][j].up = 0;
				moves[i][j] = 0;

				if( i == max-1){
					curr[i][j].left = 0;
					curr[i][j].right = 0;
					curr[i][j].down = 0;
					curr[i][j].up = 0;
				}
			}
		}
	}
	realease_int(&moves);
}

//Put 0 where will be used in matrix
void reset_matrix(Node ** matrix, int matrix_size){

	for(int i = 0; i < matrix_size; i++){
		for(int j = 0; j < matrix_size; j++){
			matrix[i][j].left = 0;
			matrix[i][j].right = 0;
			matrix[i][j].down = 0;
			matrix[i][j].up = 0;
		}
	}
}

//Add particles in the area of use
void particle_in(Node ** matrix, Node** matrix_mpi, int matrix_size, int division){

	int particles = (pow(matrix_size, 2))/division;
	int position = 0;

	int row, column;

	for(int i = 0; i < particles; i++){

		while(1){
			row = 1+(rand() % (matrix_size-2));
			column = 1+(rand() % (matrix_size-2));
			position = 1 + (rand() % 4);

			if(position == 1){
				if(matrix[row][column].left == 1){
					continue;
				}else{
					matrix[row][column].left = 1;
					matrix_mpi[row][column].left = 1;
					break;
				}
			}else if(position == 2){
				if(matrix[row][column].right == 1){
					continue;
				}else{
					matrix[row][column].right = 1;
					matrix_mpi[row][column].right = 1;
					break;
				}
			}else if(position == 3){
				if(matrix[row][column].down == 1){
					continue;
				}else{
					matrix[row][column].down = 1;
					matrix_mpi[row][column].down = 1;
					break;
				}
			}else if(position == 4){
				if(matrix[row][column].up == 1){
					continue;
				}else{
					matrix[row][column].up = 1;
					matrix_mpi[row][column].up = 1;
					break;
				}
			}
		}		
	}
}

int main(int argc, char *argv[]){

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

	//Two dimensional arrays
	Node ** curr = NULL;
	Node ** curr_mpi = NULL;
	Node ** next = NULL;
	
	alocate(&curr);
	alocate(&curr_mpi);
	alocate(&next);

	srand(time(NULL));

	//Matriz 100 Particula n^2/1000
	particle_in(curr, curr_mpi, MATRIX_LENGHT1, DIVISION_PARTICLES1);

	double T0 = MPI_Wtime();

	particle_move(curr, next, MATRIX_LENGHT1);
	double T1 = MPI_Wtime();
	double seq = T1 - T0;

	T0 = MPI_Wtime();

	particle_move_mpi(curr_mpi, next, MATRIX_LENGHT1, rank, size);

	T1 = MPI_Wtime();
	double par = T1 - T0;
	
	if(rank == 0)
		printf("Tamanho da matriz: 100  N de particulas: n^2/1000 \nSequencial: %lf   MPI:%lf SpeedUp: %lf\n\n", seq, par, seq/par);

	//Matriz 100 Particula n^2/100
	particle_in(curr, curr_mpi, MATRIX_LENGHT1, DIVISION_PARTICLES2);

	T0 = MPI_Wtime();

	particle_move(curr, next, MATRIX_LENGHT1);
	T1 = MPI_Wtime();
	seq = T1 - T0;

	T0 = MPI_Wtime();

	particle_move_mpi(curr_mpi, next, MATRIX_LENGHT1, rank, size);

	T1 = MPI_Wtime();
	par = T1 - T0;
	
	if(rank == 0)
		printf("Tamanho da matriz: 100  N de particulas: n^2/100 \nSequencial: %lf   MPI:%lf SpeedUp: %lf\n\n", seq, par, seq/par);

	//Matriz 100 Particula n^2/10
	particle_in(curr, curr_mpi, MATRIX_LENGHT1, DIVISION_PARTICLES3);

	T0 = MPI_Wtime();

	particle_move(curr, next, MATRIX_LENGHT1);
	T1 = MPI_Wtime();
	seq = T1 - T0;

	T0 = MPI_Wtime();

	particle_move_mpi(curr_mpi, next, MATRIX_LENGHT1, rank, size);

	T1 = MPI_Wtime();
	par = T1 - T0;
	
	if(rank == 0)
		printf("Tamanho da matriz: 100  N de particulas: n^2/10 \nSequencial: %lf   MPI:%lf SpeedUp: %lf\n\n", seq, par, seq/par);

	//Matriz 200 Particula n^2/1000
	particle_in(curr, curr_mpi, MATRIX_LENGHT2, DIVISION_PARTICLES1);

	T0 = MPI_Wtime();

	particle_move(curr, next, MATRIX_LENGHT2);
	T1 = MPI_Wtime();
	seq = T1 - T0;

	T0 = MPI_Wtime();

	particle_move_mpi(curr_mpi, next, MATRIX_LENGHT2, rank, size);

	T1 = MPI_Wtime();
	par = T1 - T0;
	
	if(rank == 0)
		printf("Tamanho da matriz: 200  N de particulas: n^2/1000 \nSequencial: %lf   MPI:%lf SpeedUp: %lf\n\n", seq, par, seq/par);

	//Matriz 200 Particula n^2/100
	particle_in(curr, curr_mpi, MATRIX_LENGHT2, DIVISION_PARTICLES2);

	T0 = MPI_Wtime();

	particle_move(curr, next, MATRIX_LENGHT2);
	T1 = MPI_Wtime();
	seq = T1 - T0;

	T0 = MPI_Wtime();

	particle_move_mpi(curr_mpi, next, MATRIX_LENGHT2, rank, size);

	T1 = MPI_Wtime();
	par = T1 - T0;
	
	if(rank == 0)
		printf("Tamanho da matriz: 200  N de particulas: n^2/100 \nSequencial: %lf   MPI:%lf SpeedUp: %lf\n\n", seq, par, seq/par);

	//Matriz 200 Particula n^2/10
	particle_in(curr, curr_mpi, MATRIX_LENGHT2, DIVISION_PARTICLES3);

	T0 = MPI_Wtime();

	particle_move(curr, next, MATRIX_LENGHT2);
	T1 = MPI_Wtime();
	seq = T1 - T0;

	T0 = MPI_Wtime();

	particle_move_mpi(curr_mpi, next, MATRIX_LENGHT2, rank, size);

	T1 = MPI_Wtime();
	par = T1 - T0;
	
	if(rank == 0)
		printf("Tamanho da matriz: 200  N de particulas: n^2/10 \nSequencial: %lf   MPI:%lf SpeedUp: %lf\n\n", seq, par, seq/par);

	//Matriz 500 Particula n^2/1000
	particle_in(curr, curr_mpi, MATRIX_LENGHT3, DIVISION_PARTICLES1);

	T0 = MPI_Wtime();

	particle_move(curr, next, MATRIX_LENGHT3);
	T1 = MPI_Wtime();
	seq = T1 - T0;

	T0 = MPI_Wtime();

	particle_move_mpi(curr_mpi, next, MATRIX_LENGHT3, rank, size);

	T1 = MPI_Wtime();
	par = T1 - T0;
	
	if(rank == 0)
		printf("Tamanho da matriz: 500  N de particulas: n^2/1000 \nSequencial: %lf   MPI:%lf SpeedUp: %lf\n\n", seq, par, seq/par);

	//Matriz 500 Particula n^2/100
	particle_in(curr, curr_mpi, MATRIX_LENGHT3, DIVISION_PARTICLES2);

	T0 = MPI_Wtime();

	particle_move(curr, next, MATRIX_LENGHT3);
	T1 = MPI_Wtime();
	seq = T1 - T0;

	T0 = MPI_Wtime();

	particle_move_mpi(curr_mpi, next, MATRIX_LENGHT3, rank, size);

	T1 = MPI_Wtime();
	par = T1 - T0;
	
	if(rank == 0)
		printf("Tamanho da matriz: 500  N de particulas: n^2/100 \nSequencial: %lf   MPI:%lf SpeedUp: %lf\n\n", seq, par, seq/par);

	//Matriz 500 Particula n^2/10
	particle_in(curr, curr_mpi, MATRIX_LENGHT3, DIVISION_PARTICLES3);

	T0 = MPI_Wtime();

	particle_move(curr, next, MATRIX_LENGHT3);
	T1 = MPI_Wtime();
	seq = T1 - T0;

	T0 = MPI_Wtime();

	particle_move_mpi(curr_mpi, next, MATRIX_LENGHT3, rank, size);

	T1 = MPI_Wtime();
	par = T1 - T0;
	
	if(rank == 0)
		printf("Tamanho da matriz: 500  N de particulas: n^2/10 \nSequencial: %lf   MPI:%lf SpeedUp: %lf\n\n", seq, par, seq/par);

	//Matriz 1000 Particula n^2/1000
	particle_in(curr, curr_mpi, MATRIX_LENGHT4, DIVISION_PARTICLES1);

	T0 = MPI_Wtime();

	particle_move(curr, next, MATRIX_LENGHT4);
	T1 = MPI_Wtime();
	seq = T1 - T0;

	T0 = MPI_Wtime();

	particle_move_mpi(curr_mpi, next, MATRIX_LENGHT4, rank, size);

	T1 = MPI_Wtime();
	par = T1 - T0;
	
	if(rank == 0)
		printf("Tamanho da matriz: 1000  N de particulas: n^2/1000 \nSequencial: %lf   MPI:%lf SpeedUp: %lf\n\n", seq, par, seq/par);

	//Matriz 1000 Particula n^2/100
	particle_in(curr, curr_mpi, MATRIX_LENGHT4, DIVISION_PARTICLES2);

	T0 = MPI_Wtime();

	particle_move(curr, next, MATRIX_LENGHT4);
	T1 = MPI_Wtime();
	seq = T1 - T0;

	T0 = MPI_Wtime();

	particle_move_mpi(curr_mpi, next, MATRIX_LENGHT4, rank, size);

	T1 = MPI_Wtime();
	par = T1 - T0;
	
	if(rank == 0)
		printf("Tamanho da matriz: 1000  N de particulas: n^2/100 \nSequencial: %lf   MPI:%lf SpeedUp: %lf\n\n", seq, par, seq/par);

	//Matriz 1000 Particula n^2/10
	particle_in(curr, curr_mpi, MATRIX_LENGHT4, DIVISION_PARTICLES3);

	T0 = MPI_Wtime();

	particle_move(curr, next, MATRIX_LENGHT4);
	T1 = MPI_Wtime();
	seq = T1 - T0;

	T0 = MPI_Wtime();

	particle_move_mpi(curr_mpi, next, MATRIX_LENGHT4, rank, size);

	T1 = MPI_Wtime();
	par = T1 - T0;
	
	if(rank == 0)
		printf("Tamanho da matriz: 1000  N de particulas: n^2/10 \nSequencial: %lf   MPI:%lf SpeedUp: %lf\n\n", seq, par, seq/par);

	//Matriz 2000 Particula n^2/1000
	particle_in(curr, curr_mpi, MATRIX_LENGHT5, DIVISION_PARTICLES1);

	T0 = MPI_Wtime();

	particle_move(curr, next, MATRIX_LENGHT5);
	T1 = MPI_Wtime();
	seq = T1 - T0;

	T0 = MPI_Wtime();

	particle_move_mpi(curr_mpi, next, MATRIX_LENGHT5, rank, size);

	T1 = MPI_Wtime();
	par = T1 - T0;
	
	if(rank == 0)
		printf("Tamanho da matriz: 2000  N de particulas: n^2/1000 \nSequencial: %lf   MPI:%lf SpeedUp: %lf\n\n", seq, par, seq/par);

	//Matriz 2000 Particula n^2/100
	particle_in(curr, curr_mpi, MATRIX_LENGHT5, DIVISION_PARTICLES2);

	T0 = MPI_Wtime();

	particle_move(curr, next, MATRIX_LENGHT5);
	T1 = MPI_Wtime();
	seq = T1 - T0;

	T0 = MPI_Wtime();

	particle_move_mpi(curr_mpi, next, MATRIX_LENGHT5, rank, size);

	T1 = MPI_Wtime();
	par = T1 - T0;
	
	if(rank == 0)
		printf("Tamanho da matriz: 2000  N de particulas: n^2/100 \nSequencial: %lf   MPI:%lf SpeedUp: %lf\n\n", seq, par, seq/par);

	//Matriz 2000 Particula n^2/10
	particle_in(curr, curr_mpi, MATRIX_LENGHT5, DIVISION_PARTICLES3);

	T0 = MPI_Wtime();

	particle_move(curr, next, MATRIX_LENGHT5);
	T1 = MPI_Wtime();
	seq = T1 - T0;

	T0 = MPI_Wtime();

	particle_move_mpi(curr_mpi, next, MATRIX_LENGHT5, rank, size);

	T1 = MPI_Wtime();
	par = T1 - T0;
	
	if(rank == 0)
		printf("Tamanho da matriz: 2000  N de particulas: n^2/10 \nSequencial: %lf   MPI:%lf  SpeedUp: %lf\n\n", seq, par, seq/par);

	realease(curr);
	realease(next);
	realease(curr_mpi);

    MPI_Finalize();

	return 0;
}