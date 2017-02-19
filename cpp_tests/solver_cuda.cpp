#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern"C" {
 void write_test();
 void read_matrix(int **, int **, double **, double **, int *, int *);
}

void read_values(char* filename, int** ia, int** ja, double** a, double** loads, int* neq, int* nnz){
    FILE* aux_file;
    char b[1024];

    aux_file = fopen("convert.aux", "w");
    sprintf(b, filename);
    fprintf(aux_file, "%s\n", b);
    fclose(aux_file);   

    printf("Indo para o Fortran\n");
    read_matrix(ia, ja, a, loads, neq, nnz);

    #pragma omp parallel
    {
        #pragma omp for  
        for(int i=0; i < *neq+1; i++){
            (*ia)[i] = (*ia)[i] - 1;
        }

        #pragma omp for  
        for(int i=0; i < *nnz; i++){
            (*ja)[i] = (*ja)[i] - 1;
        }
    }

}


int main(int argc, char *argv[]){
   printf("Hello World!\n");

   int *ia, *ja, neq, nnz;
   double *a, *loads;


   read_values(" /home/mateus/Solvers/Testcases/InputLS_PCG_25600els.aux", &ia, &ja, &a, &loads, &neq, &nnz);


   printf("NNZ = %d\t NEQ = %d\n", neq, nnz);
   printf(" ia = %10d%10d%10d ... %10d%10d%10d\n", ia[0], ia[1], ia[2], ia[neq-2], ia[neq-1], ia[neq]);
   printf(" ja = %10d%10d%10d ... %10d%10d%10d\n", ja[0], ja[1], ja[2], ja[nnz-3], ja[nnz-2], ja[nnz-1]);
   printf("  a = %10.2e%10.2e%10.2e ... %10.2e%10.2e%10.2e\n",  a[0],  a[1],  a[2],  a[nnz-3],  a[nnz-2],  a[nnz-1]);


   return 0;
}