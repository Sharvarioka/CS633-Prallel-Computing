#include<stdio.h>

#include<stdio.h>

#include<string.h>

#include "mpi.h"

#include<math.h>

#include <stdlib.h>

#include <time.h>


void multiplesend(FILE ** fp, double ** matrix, double ** wrapper, int row, int iteration, int myrank, int north, int south, int east, int west, int N) {

    double stime, etime, ttime;
    int col = row;
    int wrow = row + 2;
    int wcol = col + 2;

    //initializing matrix with random values
    srand(time(NULL));
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            matrix[i][j] = rand() % 202102152200 + 1;
        }
    }

    // initializing border with -1
    for (int i = 0; i < wrow; i++) {
        for (int j = 0; j < wcol; j++) {
            if (i == 0 || j == 0 || i == row + 1 || j == col + 1)
                wrapper[i][j] = -1;
        }
    }

    stime = MPI_Wtime();

    MPI_Request reqs[8];

    while (iteration--) {

        //copy original array to wrapper array
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                wrapper[(i + 1)][(j + 1)] = matrix[i][j];
            }
        }

        //receive from north 
        for (int i = 1; i <= col; i++) {
            MPI_Irecv( & wrapper[0][i], 1, MPI_DOUBLE, north, north, MPI_COMM_WORLD, & reqs[4]);
        }
        //receive from south 
        for (int i = 1; i <= col; i++) {
            MPI_Irecv( & wrapper[(row + 1)][i], 1, MPI_DOUBLE, south, south, MPI_COMM_WORLD, & reqs[5]);
        }
        //receive from east 
        for (int i = 1; i <= row; i++) {
            MPI_Irecv( & wrapper[i][col + 1], 1, MPI_DOUBLE, east, east, MPI_COMM_WORLD, & reqs[6]);
        }
        //receive from west
        for (int i = 1; i <= row; i++) {
            MPI_Irecv( & wrapper[i][0], 1, MPI_DOUBLE, west, west, MPI_COMM_WORLD, & reqs[7]);
        }

        //send north values
        for (int i = 0; i < col; i++) {
            MPI_Isend( & matrix[0][i], 1, MPI_DOUBLE, north, myrank, MPI_COMM_WORLD, & reqs[0]);
        }

        //send south values
        for (int i = 0; i < col; i++) {
            MPI_Isend( & matrix[(col - 1)][i], 1, MPI_DOUBLE, south, myrank, MPI_COMM_WORLD, & reqs[1]);
        }

        //send east values
        for (int i = 0; i < row; i++) {
            MPI_Isend( & matrix[i][col - 1], 1, MPI_DOUBLE, east, myrank, MPI_COMM_WORLD, & reqs[2]);
        }

        //send west values

        for (int i = 0; i < row; i++) {
            MPI_Isend( & matrix[i][0], 1, MPI_DOUBLE, west, myrank, MPI_COMM_WORLD, & reqs[3]);
        }

        MPI_Waitall(8, reqs, MPI_STATUSES_IGNORE);

        //now doing stencil computation
        for (int i = 1; i <= row; i++) {
            for (int j = 1; j <= col; j++) {
                double c = 0;
                if (wrapper[(i - 1)][j] != -1) c++;
                if (wrapper[(i + 1)][j] != -1) c++;
                if (wrapper[i][j - 1] != -1) c++;
                if (wrapper[i][j + 1] != -1) c++;
                matrix[i-1][j-1] = (wrapper[(i - 1)][j] + wrapper[(i + 1)][j] + wrapper[i][j - 1] + wrapper[i][j + 1]) / c;
            }
        }

    }

    etime = MPI_Wtime();

    ttime = etime - stime;

    double maxtime;
    MPI_Reduce( & ttime, & maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (myrank == 0) {
        printf("%lf\n", maxtime);
        fprintf( * fp, "\nMultiple sends,%d,%lf", N, log2(maxtime));
    }

}

void ddtsend(FILE ** fp, double ** matrix, double ** wrapper, int row, int iteration, int myrank, int north, int south, int east, int west, int N) {

    double stime, etime, ttime;
    int col = row;
    int wrow = row + 2;
    int wcol = col + 2;

    //initializing matrix with random values
    srand(time(NULL));
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            matrix[i][j] = rand() % 202102152200 + 1;
        }
    }

    // initializing border with -1
    for (int i = 0; i < wrow; i++) {
        for (int j = 0; j < wcol; j++) {
            if (i == 0 || j == 0 || i == row + 1 || j == col + 1)
                wrapper[i][j] = -1;
        }
    }

    MPI_Request reqs[8];
    //create derived data types

    //north_south
    MPI_Datatype north_south;
    MPI_Type_contiguous(row, MPI_DOUBLE, & north_south);
    MPI_Type_commit( & north_south);

    //east_west
    MPI_Datatype east_west;
    MPI_Type_vector(row, 1, row, MPI_DOUBLE, & east_west);
    MPI_Type_commit( & east_west);

    //new datatype for wrapper matix
    MPI_Datatype weast_west;
    MPI_Type_vector(row, 1, wrow, MPI_DOUBLE, & weast_west);
    MPI_Type_commit( & weast_west);
    stime = MPI_Wtime();

    while (iteration--) {

        //copy remaining data to wrapper
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                wrapper[i + 1][j + 1] = matrix[i][j];
            }
        }

        //sending data to north
        MPI_Isend( & matrix[0][0], 1, north_south, north, myrank, MPI_COMM_WORLD, & reqs[0]);
        //send south
        MPI_Isend( & matrix[(row - 1)][0], 1, north_south, south, myrank, MPI_COMM_WORLD, & reqs[1]);
        //send east
        MPI_Isend( & matrix[0][(col - 1)], 1, east_west, east, myrank, MPI_COMM_WORLD, & reqs[2]);
        //send west
        MPI_Isend( & matrix[0][0], 1, east_west, west, myrank, MPI_COMM_WORLD, & reqs[3]);

        //recieve from north
        MPI_Irecv( & wrapper[0][1], 1, north_south, north, north, MPI_COMM_WORLD, & reqs[4]);
        //receive from south
        MPI_Irecv( & wrapper[(row + 1)][1], 1, north_south, south, south, MPI_COMM_WORLD, & reqs[5]);
        //receive from east
        MPI_Irecv( & wrapper[1][(col + 1)], 1, weast_west, east, east, MPI_COMM_WORLD, & reqs[6]);
        //receive from west
        MPI_Irecv( & wrapper[1][0], 1, weast_west, west, west, MPI_COMM_WORLD, & reqs[7]);

        MPI_Waitall(8, reqs, MPI_STATUSES_IGNORE);

        //computation
        for (int i = 1; i <= row; i++) {
            for (int j = 1; j <= col; j++) {
                double c = 0;
                if (wrapper[(i - 1)][j] != -1) c++;
                if (wrapper[(i + 1)][j] != -1) c++;
                if (wrapper[i][j - 1] != -1) c++;
                if (wrapper[i][j + 1] != -1) c++;
                matrix[i-1][j-1] = (wrapper[(i - 1)][j] + wrapper[(i + 1)][j] + wrapper[i][j - 1] + wrapper[i][j + 1]) / c;
            }
        }

    }

    etime = MPI_Wtime();

    ttime = etime - stime;
    MPI_Type_free( & north_south);
    MPI_Type_free( & east_west);
    MPI_Type_free( & weast_west);
    double maxtime;
    MPI_Reduce( & ttime, & maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (myrank == 0) {

        printf("%lf\n", maxtime);
        fprintf( * fp, "\nDerived,%d,%lf", N, log2(maxtime));
    }

}

void packed(FILE ** fp, double ** matrix, double ** wrapper, int row, int iteration, int myrank, int north, int south, int east, int west, int N) {

    MPI_Status status;
    double stime, etime, ttime;
    int col = row;
    int wrow = row + 2;
    int wcol = col + 2;
    //initializing matrix with random values
    srand(time(NULL));
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            matrix[i][j] = rand() % 202102152200 + 1;
        }
    }

    // initializing border with -1
    for (int i = 0; i < wrow; i++) {
        for (int j = 0; j < wcol; j++) {
            if (i == 0 || j == 0 || i == row + 1 || j == col + 1)
                wrapper[i][j] = -1;
        }
    }

    //buffers for packing and unpacking
    double * send_north = (double * ) malloc(row * sizeof(double));
    double * send_south = (double * ) malloc(row * sizeof(double));
    double * send_east = (double * ) malloc(row * sizeof(double));
    double * send_west = (double * ) malloc(row * sizeof(double));

    double * recv_north = (double * ) malloc(row * sizeof(double));
    double * recv_south = (double * ) malloc(row * sizeof(double));
    double * recv_east = (double * ) malloc(row * sizeof(double));
    double * recv_west = (double * ) malloc(row * sizeof(double));

    int position = 0;
    stime = MPI_Wtime();

    MPI_Request reqs[8];

    //iterating iteration times computation and communication
    while (iteration--) {

        //copy original array to wrapper array
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                wrapper[i + 1][j + 1] = matrix[i][j];
            }
        }

        //packing 1st row, sending, receiving and unpacking
        for (int i = 0; i < col; i++)
            MPI_Pack( & matrix[0][i], 1, MPI_DOUBLE, send_north, row * sizeof(double), & position, MPI_COMM_WORLD);
        MPI_Send(send_north, position, MPI_PACKED, north, myrank, MPI_COMM_WORLD);
        MPI_Recv(recv_north, row * sizeof(double), MPI_PACKED, north, north, MPI_COMM_WORLD, & status);
        position = 0;
        for (int i = 1; i <= col; i++)
            MPI_Unpack(recv_north, row * sizeof(double), & position, & wrapper[0][i], 1, MPI_DOUBLE, MPI_COMM_WORLD);

        //packing last row, sending, receiving and unpacking
        position = 0;
        for (int i = 0; i < col; i++)
            MPI_Pack( & matrix[row - 1][i], 1, MPI_DOUBLE, send_south, row * sizeof(double), & position, MPI_COMM_WORLD);
        MPI_Send(send_south, position, MPI_PACKED, south, myrank, MPI_COMM_WORLD);
        MPI_Recv(recv_south, row * sizeof(double), MPI_PACKED, south, south, MPI_COMM_WORLD, & status);
        position = 0;
        for (int i = 1; i <= col; i++)
            MPI_Unpack(recv_south, row * sizeof(double), & position, & wrapper[row + 1][i], 1, MPI_DOUBLE, MPI_COMM_WORLD);

        ////packing last column, sending, receiving and unpacking
        position = 0;
        for (int i = 0; i < col; i++)
            MPI_Pack( & matrix[i][col - 1], 1, MPI_DOUBLE, send_east, row * sizeof(double), & position, MPI_COMM_WORLD);
        MPI_Send(send_east, position, MPI_PACKED, east, myrank, MPI_COMM_WORLD);
        MPI_Recv(recv_east, row * sizeof(double), MPI_PACKED, east, east, MPI_COMM_WORLD, & status);
        position = 0;
        for (int i = 1; i <= col; i++)
            MPI_Unpack(recv_east, row * sizeof(double), & position, & wrapper[i][col + 1], 1, MPI_DOUBLE, MPI_COMM_WORLD);

        //packing 1st column, sending, receiving and unpacking
        position = 0;
        for (int i = 0; i < col; i++)
            MPI_Pack( & matrix[i][0], 1, MPI_DOUBLE, send_west, row * sizeof(double), & position, MPI_COMM_WORLD);
        MPI_Send(send_west, position, MPI_PACKED, west, myrank, MPI_COMM_WORLD);
        MPI_Recv(recv_west, row * sizeof(double), MPI_PACKED, west, west, MPI_COMM_WORLD, & status);
        position = 0;
        for (int i = 1; i <= col; i++)
            MPI_Unpack(recv_west, row * sizeof(double), & position, & wrapper[i][0], 1, MPI_DOUBLE, MPI_COMM_WORLD);
        position = 0;

        //stencil computation

        for (int i = 1; i <= row; i++) {
            for (int j = 1; j <= col; j++) {
                double c = 0;
                if (wrapper[(i - 1)][j] != -1) c++;
                if (wrapper[(i + 1)][j] != -1) c++;
                if (wrapper[i][j - 1] != -1) c++;
                if (wrapper[i][j + 1] != -1) c++;
               matrix[i-1][j-1] = (wrapper[(i - 1)][j] + wrapper[(i + 1)][j] + wrapper[i][j - 1] + wrapper[i][j + 1]) / c;
            }
        }

    }

    etime = MPI_Wtime();

    ttime = etime - stime;
    free(send_north);
    free(send_south);
    free(send_west);
    free(send_east);
    free(recv_north);
    free(recv_south);
    free(recv_east);
    free(recv_west);
    double maxtime;
    MPI_Reduce( & ttime, & maxtime, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (myrank == 0) {

        printf("%lf\n", maxtime);

        fprintf( * fp, "\nPacked,%d,%lf", N, log2(maxtime));
    }

}

int main(int argc, char * argv[]) {

    MPI_Init( & argc, & argv);

    if (argc < 3) {
        printf("3 input parameters needed\n");
        exit(2);
    }

    int cpus, myrank;
    MPI_Comm_size(MPI_COMM_WORLD, & cpus);
    MPI_Comm_rank(MPI_COMM_WORLD, & myrank);
    MPI_Status status;

    int N = atoi(argv[1]);
    int iteration = atoi(argv[2]);

    int row = sqrt(N);
    int col = row;

    int processes = cpus;
    int process_grid_row = sqrt(processes);
    int process_grid_col = process_grid_row;

    // compute all neighbours 
    int north = myrank - process_grid_row;
    if (north < 0) north = MPI_PROC_NULL;

    int south = myrank + process_grid_row;
    if (south >= processes) south = MPI_PROC_NULL;

    int east = myrank + 1;
    if ((myrank / process_grid_row) != (east / process_grid_row)) east = MPI_PROC_NULL;

    int west = myrank - 1;
    if (west < 0 || ((myrank / process_grid_row) != (west / process_grid_row))) west = MPI_PROC_NULL;

    // matrix
    double ** matrix = (double ** ) malloc(row * sizeof(double * ));
    matrix[0] = (double * ) malloc(row * col * sizeof(double));
    for (int i = 1; i < row; i++) {
        matrix[i] = matrix[i - 1] + col;
    }

    int wrow = row + 2;
    int wcol = col + 2;

    //wrapper matrix for easy computation
    double ** wrapper = (double ** ) malloc(wrow * sizeof(double * ));
    wrapper[0] = (double * ) malloc(wrow * wcol * sizeof(double));
    for (int i = 1; i < wrow; i++) {
        wrapper[i] = wrapper[i - 1] + wcol;
    }

    //datafiles
    char buf[12];
    FILE * fp;
    snprintf(buf, 12, "plot%d.csv", processes);
    fp = fopen(buf, "a");

    //calling each method
    multiplesend( & fp, matrix, wrapper, row, iteration, myrank, north, south, east, west, N);
    packed( & fp, matrix, wrapper, row, iteration, myrank, north, south, east, west, N);
    ddtsend( & fp, matrix, wrapper, row, iteration, myrank, north, south, east, west, N);

    MPI_Finalize();

    return 0;

}
