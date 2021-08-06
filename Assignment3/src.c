#include <stdio.h>

#include "mpi.h"

#include <time.h>

#include <stdlib.h>

#include<string.h>

int main(int argc, char * argv[]) {

    if (argc != 2) {
        printf("to run: mpiexec -np #process ./src filename.csv\n");
        exit(2);
    }

    int world_size, world_rank;
    double stime, etime;
    double ** matrix;

    MPI_Init( & argc, & argv);
    MPI_Comm_size(MPI_COMM_WORLD, & world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, & world_rank);

    int row = 0, col = 0;

    //process 0 reading the file sequentially
    if (world_rank == 0) {

        FILE * fp = fopen(argv[1], "r");
        char * line = NULL;
        size_t len = 0;
        ssize_t read;
        //reading header line and calculating #cols
        read = getline( & line, & len, fp);
        for (int i = 0; i < read; i++) {
            if (line[i] == ',') {
                col++;
            }
        }
        col--; // #commas + 1 and ignoring 1st two cols;

        while (read = getline( & line, & len, fp) != -1) {
            row++;
        }
        //setting fp to point to begin of file
        fseek(fp, 0, SEEK_SET);
        read = getline( & line, & len, fp);

        //2d matrix to hold data from csv file
        matrix = (double ** ) malloc(row * sizeof(double * ));
        matrix[0] = (double * ) malloc(row * col * sizeof(double));
        for (int i = 1; i < row; i++) {
            matrix[i] = matrix[i - 1] + col;
        }

        long i = 0, j = 0;

        while (read = getline( & line, & len, fp) != -1) {
            //convert stirng into doubles and store in a 2d array
            char * p = line;
            double val;
            //next two line for skipping lat and long value
            strtod(p, & p);
            p++;
            strtod(p, & p);
            p++;
            while ( * p) {
                val = strtod(p, & p);
                matrix[i][j] = val;
                j++;
                p++;
            }
            j = 0;
            i++;
        }

        fclose(fp);
        if (line)
            free(line);
    }
    //timer started afer reading the whole file
    stime = MPI_Wtime();
    MPI_Status status;
    //broadcasting rows and cols to every process
    MPI_Bcast( & row, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast( & col, 1, MPI_INT, 0, MPI_COMM_WORLD);

    //making portion for each process
    int portion = row / world_size;
    //last process receive more than other in case of any remainder
    if (row % world_size && world_rank == world_size - 1) {
        portion += row % world_size;
    }

    //array to store portion of data
    double ** recvmatrix = (double ** ) malloc(portion * sizeof(double * ));
    recvmatrix[0] = (double * ) malloc(portion * col * sizeof(double));
    for (int i = 1; i < portion; i++) {
        recvmatrix[i] = recvmatrix[i - 1] + col;
    }

    //root process does a scatterv to all processes
    int send_count[world_size], send_disp[world_size];
    if (!world_rank) {
        send_count[0] = portion * col;
        send_disp[0] = 0;

        for (int i = 1; i < world_size - 1; i++) {
            send_count[i] = portion * col;
            send_disp[i] = send_count[i - 1] + send_disp[i - 1];
        }
        if (world_size - 1) {
            send_count[world_size - 1] = (portion + row % world_size) * col;
            send_disp[world_size - 1] = send_count[world_size - 2] + send_disp[world_size - 2];

        }
    }

    MPI_Scatterv( & matrix[0][0], send_count, send_disp, MPI_DOUBLE, recvmatrix[0], portion * col, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //freeing matrix
    if (world_rank == 0) {
        free(matrix[0]);
        free(matrix);
    }
    //1d array to store min temp of each portion by each process
    double yearly_min_small[col];

    for (int i = 0; i < col; i++) {
        double mini = 999;
        for (int j = 0; j < portion; j++) {
            if (recvmatrix[j][i] < mini)
                mini = recvmatrix[j][i];
        }
        yearly_min_small[i] = mini;
    }

    //freeing recvmatrix
    free(recvmatrix[0]);
    free(recvmatrix);

    //1d for year wise min temp across all location, for each year
    double yearly_min_final[col];
    MPI_Reduce(yearly_min_small, yearly_min_final, col, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    //1st line of output
    if (!world_rank) {
        for (int i = 0; i < col; i++)
            printf("%f,", yearly_min_final[i]);
        printf("\n");
    }
    //calculating global min temp
    double global_min = 999;
    if (!world_rank) {
        for (int i = 0; i < col; i++) {
            if (yearly_min_final[i] < global_min)
                global_min = yearly_min_final[i];
        }
        //2nd line of output
        printf("%f\n", global_min);
    }

    //timer ends
    etime = MPI_Wtime();
    double ttime = etime - stime;
    double maxtime;
    //max time of all processes
    MPI_Reduce( & ttime, & maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    //3rd line of output
    if (world_rank == 0)
        printf("%f\n", maxtime);

    MPI_Finalize();
    return 0;
}
