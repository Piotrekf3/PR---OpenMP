

// Celem tego programu jest prezentacja pomiaru i analizy 
//efektywnosci programu za pomocš  CodeAnalyst(tm).
// Implementacja mnożenia macierzy jest realizowana za pomoca typowego 
// algorytmu podręcznikowego. 
#include <stdio.h>
#include <time.h>
#include <windows.h>
#include "omp.h"
#include <iostream>
#include <iomanip>

using namespace std;

#define USE_MULTIPLE_THREADS true
#define MAXTHREADS 128
int NumThreads;
double start;

static const int ROWS = 600;     // liczba wierszy macierzy
static const int COLUMNS = 600;  // lizba kolumn macierzy

float matrix_a[ROWS][COLUMNS];    // lewy operand 
float matrix_b[ROWS][COLUMNS];    // prawy operand
float matrix_r[ROWS][COLUMNS];   // wynik
float matrix_seq[ROWS][COLUMNS];
float matrix_local[ROWS][COLUMNS][4];

FILE *result_file;

void printMatrix(const float matrix[ROWS][COLUMNS])
{
	for (int i = 0; i < ROWS; i++) {
		for (int j = 0; j < COLUMNS; j++)
			cout <<matrix[i][j]<< " ";
		cout << endl;
	}
}

void initialize_matrices()
{
	// zdefiniowanie zawarosci poczatkowej macierzy
	//#pragma omp parallel for 
	for (int i = 0; i < ROWS; i++) {
		for (int j = 0; j < COLUMNS; j++) {
			matrix_a[i][j] = (float)rand() / RAND_MAX;
			matrix_b[i][j] = (float)rand() / RAND_MAX;
			matrix_r[i][j] = 0.0;
		}
	}
}

void initialize_matricesZ()
{
	// zdefiniowanie zawarosci poczatkowej macierzy
#pragma omp parallel for 
	for (int i = 0; i < ROWS; i++) {
		for (int j = 0; j < COLUMNS; j++) {
			matrix_r[i][j] = 0.0;
		}
	}
}

void initialize_matrices_local()
{
	// zdefiniowanie zawarosci poczatkowej macierzy
#pragma omp parallel for 
	for (int k = 0; k < 4;k++)
	for (int i = 0; i < ROWS; i++) {
		for (int j = 0; j < COLUMNS; j++) {
			matrix_local[i][j][k] = 0.0;
		}
	}
}

void initialize_matricesSeq()
{
	// zdefiniowanie zawarosci poczatkowej macierzy
#pragma omp parallel for 
	for (int i = 0; i < ROWS; i++) {
		for (int j = 0; j < COLUMNS; j++) {
			matrix_seq[i][j] = 0.0;
		}
	}
}

void print_result()
{
	// wydruk wyniku
	for (int i = 0; i < ROWS; i++) {
		for (int j = 0; j < COLUMNS; j++) {
			fprintf(result_file, "%6.4f ", matrix_r[i][j]);
		}
		fprintf(result_file, "\n");
	}
}

void multiply_matrices_IJK()
{
	// mnozenie macierzy 
#pragma omp parallel for 
	for (int i = 0; i < ROWS; i++) {
		for (int j = 0; j < COLUMNS; j++) {
			float sum = 0.0;
			for (int k = 0; k < COLUMNS; k++) {
				sum = sum + matrix_a[i][k] * matrix_b[k][j];
			}
			matrix_r[i][j] = sum;
		}
	}
}

void multiply_matrices_IKJ()
{
	// mnozenie macierzy 
#pragma omp parallel for 
	for (int i = 0; i < ROWS; i++)
	for (int k = 0; k < COLUMNS; k++)
	for (int j = 0; j < COLUMNS; j++)
		matrix_r[i][j] += matrix_a[i][k] * matrix_b[k][j];

}

void multiply_matrices_JIK()
{
	// mnozenie macierzy 
#pragma omp parallel for 

	for (int j = 0; j < COLUMNS; j++) {
		for (int i = 0; i < ROWS; i++) {
			float sum = 0.0;
			for (int k = 0; k < COLUMNS; k++) {
				sum = sum + matrix_a[i][k] * matrix_b[k][j];
			}
			matrix_r[i][j] = sum;
		}
	}
}
void multiply_matrices_JKI()
{
	// mnozenie macierzy 
#pragma omp parallel for 
	for (int j = 0; j < COLUMNS; j++)
	for (int k = 0; k < COLUMNS; k++)
	for (int i = 0; i < ROWS; i++)
		matrix_r[i][j] += matrix_a[i][k] * matrix_b[k][j];

}

void multiply_matrices_KJI()
{
	// mnozenie macierzy 
	
#pragma omp parallel
	{
		int id = omp_get_thread_num();
#pragma omp for 
		for (int k = 0; k < COLUMNS; k++)
		for (int j = 0; j < COLUMNS; j++)
		for (int i = 0; i < ROWS; i++)
			matrix_local[i][j][id] += matrix_a[i][k] * matrix_b[k][j];
	}
	for (int k = 0; k < 4; k++)
	for (int j = 0; j < COLUMNS; j++)
	for (int i = 0; i < ROWS; i++)
	{
		matrix_r[i][j] += matrix_local[i][j][k];
	}
	/*
#pragma omp parallel for 
		for (int k = 0; k < COLUMNS; k++)
		{
			for (int j = 0; j < COLUMNS; j++)
			for (int i = 0; i < ROWS; i++)
#pragma omp atomic
				matrix_r[i][j] += matrix_a[i][k] * matrix_b[k][j];
		}*/
}

bool compareMatrices(float matrix_1[ROWS][COLUMNS], float matrix_2[ROWS][COLUMNS])
{
	for (int j = 0; j < COLUMNS; j++)
		for (int i = 0; i < ROWS; i++)
		{
			if (!(abs(matrix_1[i][j] - matrix_2[i][j])<0.001))
			{
				cout << setprecision(9) << matrix_1[i][j] << endl;
				cout << setprecision(9) << matrix_2[i][j] << endl;
				return false;
			}
		}
	return true;
}

void multiply_matrices_IJK_seq()
{
	// mnozenie macierzy  
	for (int k = 0; k < COLUMNS; k++)
	for (int j = 0; j < COLUMNS; j++)
	for (int i = 0; i < ROWS; i++)
		matrix_seq[i][j] += matrix_a[i][k] * matrix_b[k][j];

}

void multiply_matrices_6loops()
{
	// mnozenie macierzy 
	int r = 100;
	for (int i = 0; i < ROWS; i += r)
	for (int j = 0; j < COLUMNS; j += r)
	for (int k = 0; k < COLUMNS; k+=r)
#pragma omp parallel for 
	for (int ii = i; ii < i + r; ii++)
	for (int kk = k; kk < k + r; kk++)
	for (int jj = j; jj < j + r; jj++)
		matrix_r[ii][jj] += matrix_a[ii][kk] * matrix_b[kk][jj];

}

void print_elapsed_time()
{
	double elapsed;
	double resolution;

	// wyznaczenie i zapisanie czasu przetwarzania
	elapsed = (double)clock() / CLK_TCK;
	resolution = 1.0 / CLK_TCK;
	printf("Czas: %8.4f sec \n",
		elapsed - start);
	fprintf(result_file,
		"Czas wykonania programu: %8.4f sec (%6.4f sec rozdzielczosc pomiaru)\n",
		elapsed - start, resolution);
}

int main(int argc, char* argv[])
{
	srand(time(NULL));
	//	 start = (double) clock() / CLK_TCK ;
	if ((result_file = fopen("classic.txt", "a")) == NULL) {
		fprintf(stderr, "nie mozna otworzyc pliku wyniku \n");
		perror("classic");
		return(EXIT_FAILURE);
	}


	//Determine the number of threads to use
	if (USE_MULTIPLE_THREADS) {
		SYSTEM_INFO SysInfo;
		GetSystemInfo(&SysInfo);
		NumThreads = SysInfo.dwNumberOfProcessors;
		if (NumThreads > MAXTHREADS)
			NumThreads = MAXTHREADS;
	}
	else
		NumThreads = 1;
	fprintf(result_file, "Klasyczny algorytm mnozenia macierzy, liczba watkow %d \n", NumThreads);
	printf("liczba watkow  = %d\n\n", NumThreads);

	

	//for (int i = 0; i < 10; i++)
	//{
		initialize_matrices();
		initialize_matricesZ();
		initialize_matrices_local();
		start = (double)clock() / CLK_TCK;
		multiply_matrices_KJI();
		//printMatrix(matrix_r);
		printf("KJI ");
		print_elapsed_time();

		/*initialize_matricesSeq();
		start = (double)clock() / CLK_TCK;
		multiply_matrices_IJK_seq();
		//printMatrix(matrix_seq);
		printf("sekw ");
		print_elapsed_time();*/


		/*if (!compareMatrices(matrix_r, matrix_seq))
		{
			cout << "Nie" << endl;
			break;
		}*/




		/*initialize_matricesZ();
		start = (double)clock() / CLK_TCK;
		multiply_matrices_6loops();
		//printMatrix(matrix_r);
		printf("6 loops ");
		print_elapsed_time();*/

		/*initialize_matricesSeq();
		start = (double)clock() / CLK_TCK;
		multiply_matrices_KJI_seq();
		//printMatrix(matrix_seq);
		printf("sekw ");
		print_elapsed_time();
		if (!compareMatrices(matrix_r, matrix_seq))
		{
			cout << "Nie" << endl;
			break;
		}*/


	//}

	/*initialize_matrices();
	initialize_matricesZ();
	initialize_matrices_local();
	start = (double)clock() / CLK_TCK;
	multiply_matrices_KJI();
	//printMatrix(matrix_r);
	printf("KJI ");
	print_elapsed_time();*/

	fclose(result_file);

	return(0);
}