#include<stdio.h>

#include<stdio.h>

#include<string.h>

#include "mpi.h"

#include<math.h>

#include <stdlib.h>

#include <time.h>


void MPI_Bcast_default(FILE ** fp1, double * buf, int myrank, int noOfEle) {
    double ttime = 0;
    int itr = 5;
    while (itr--) {
        double sTime = MPI_Wtime();
        MPI_Bcast(buf, noOfEle, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        double eTime = MPI_Wtime();

        ttime += eTime - sTime;
    }
    if (myrank == 0)
        printf("default_Bcast: %lf\n", ttime / 5);
    fprintf( * fp1, "\ndefault_Bcast,%lf", ttime / 5);
}

void MPI_Bcast_optimized(FILE ** fp1, double * buf, int myrank, int noOfEle, int ppn, MPI_Comm first_comm, MPI_Comm second_comm, double total_comm_time) {
    double ttime = 0;
    int itr = 5;
    while (itr--) {
        double stime = MPI_Wtime();
        if (myrank % ppn == 0)
            MPI_Bcast(buf, noOfEle, MPI_DOUBLE, 0, second_comm);
        MPI_Bcast(buf, noOfEle, MPI_DOUBLE, 0, first_comm);
        double etime = MPI_Wtime();
        ttime += etime - stime;
    }
    if (myrank == 0)
        printf("optimized_Bcast: %lf\n", ttime / 5 + total_comm_time);
    fprintf( * fp1, "\noptimized_Bcast,%lf", ttime / 5 + total_comm_time);
}

void MPI_Gather_default(FILE ** fp2, double * buf, int myrank, int nodes, int noOfEle) {
    double ttime = 0;
    int itr = 5;
    while (itr--) {
        double * recvMessage = (double * ) malloc(noOfEle * nodes * sizeof(double));
        double sTime = MPI_Wtime();
        MPI_Gather(buf, noOfEle, MPI_DOUBLE, recvMessage, noOfEle, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        double eTime = MPI_Wtime();
        ttime += eTime - sTime;
    }
    if (myrank == 0) {
        printf("default_Gather: %lf\n", ttime / 5);
        fprintf( * fp2, "\ndefault_Gather,%lf", ttime / 5);
    }

}

void MPI_Gather_optimized(FILE ** fp2, double * buf, int myrank, int nodes, int noOfEle, MPI_Comm first_comm, MPI_Comm second_comm, int ppn, int P, double total_comm_time) {
    double ttime = 0;
    int itr = 5;
    while (itr--) {
        double * recvfirstcomm = (double * ) malloc(noOfEle * sizeof(double) * ppn);
        double stime = MPI_Wtime();
        MPI_Gather(buf, noOfEle, MPI_DOUBLE, recvfirstcomm, noOfEle, MPI_DOUBLE, 0, first_comm);
        double * recvsecondcomm = (double * ) malloc(noOfEle * sizeof(double) * ppn * P);
        if (myrank % ppn == 0)
            MPI_Gather(recvfirstcomm, noOfEle * ppn, MPI_DOUBLE, recvsecondcomm, noOfEle * ppn, MPI_DOUBLE, 0, second_comm);
        double etime = MPI_Wtime();
        ttime += etime - stime;
    }
    if (myrank == 0) {
        printf("optimized_Gather: %lf\n", ttime / 5 + total_comm_time);
        fprintf( * fp2, "\noptimized_Gather,%lf", ttime / 5 + total_comm_time);
    }
}

void MPI_Reduce_default(FILE ** fp3, double * buf, int myrank, int nodes, int noOfEle) {
    double ttime = 0;
    int itr = 5;
    while (itr--) {
        double maxval;
        double sTime = MPI_Wtime();
        MPI_Reduce( & buf, & maxval, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        double eTime = MPI_Wtime();

        ttime += eTime - sTime;
    }
    if (myrank == 0) {
        printf("default_Reduce: %lf\n", ttime / 5);
        fprintf( * fp3, "\ndefault_Reduce,%lf", ttime / 5);
    }
}

void MPI_Reduce_optimized(FILE ** fp3, double * buf, int myrank, int nodes, int noOfEle, MPI_Comm first_comm, MPI_Comm second_comm, int ppn, int P, double total_comm_time) {
    double ttime = 0;
    int itr = 5;
    while (itr--) {
        double firstmax, secondmax;
        double stime = MPI_Wtime();
        MPI_Reduce(buf, & firstmax, 1, MPI_DOUBLE, MPI_MAX, 0, first_comm);
        if (myrank % ppn == 0)
            MPI_Reduce( & firstmax, & secondmax, 1, MPI_DOUBLE, MPI_MAX, 0, second_comm);
        double etime = MPI_Wtime();
        ttime += etime - stime;
    }
    if (myrank == 0) {
        printf("optimized_Reduce: %lf\n", ttime / 5 + total_comm_time);
        fprintf( * fp3, "\noptimized_Reduce,%lf", ttime / 5 + total_comm_time);
    }
}

void MPI_Alltoallv_default(FILE ** f4, int myrank, int P, int ppn, int noOfEle){
  double ttime = 0;
  int itr = 5;
  while(itr--){
    int size = P*ppn;
    //send count array
    int* sendcount = (int*)malloc(size*sizeof(int));
    //displacement array 
    int* senddisp = (int*)malloc(size*sizeof(int));

    int* recvcount = (int*)malloc(size*sizeof(int));
    int* recvdisp = (int*)malloc(size*sizeof(int));
    //initialize ramdomly sendcount array
    int r = rand()%noOfEle + 1;
    for(int i = 0;i<size;i++){
      sendcount[i] = i*r;
      recvcount[i] = myrank*r;
      
    }
     senddisp[0] = 0;
    for(int i=1;i<size;i++){
      senddisp[i] = senddisp[i-1] + sendcount[i-1];
    }
    recvdisp[0] = 0;
    //initialize send displacement array
    for(int i=1;i<size;i++){
      recvdisp[i] = recvdisp[i-1] + recvcount[i-1];
    }
    
    int sendbuflen = 0;
    //calculating send buf array
    for(int i=0;i<size;i++)
      sendbuflen += sendcount[i];

    double* sendbuf = (double*)malloc(sendbuflen*sizeof(double));
    //initializing send buf array

    for(int i = 0;i<sendbuflen;i++){
      sendbuf[i] = rand() % 2021 + 1;
    }
    
    int recvbuflen = 0;
    for(int i=0;i<size;i++){
      recvbuflen += recvcount[i];
    }

        //receive buffer
    double* recvbuf = (double*)malloc(recvbuflen*sizeof(double));

    double stime = MPI_Wtime();
    MPI_Alltoallv(sendbuf, sendcount, senddisp, MPI_DOUBLE, recvbuf, recvcount, recvdisp, MPI_DOUBLE, MPI_COMM_WORLD);
    double etime = MPI_Wtime();
    ttime += etime - stime;
    free(sendcount);free(senddisp);free(recvcount);free(recvdisp);free(sendbuf);free(recvbuf);
  }
  if(myrank == 0){
    printf("default_Alltoallv: %lf\n", ttime / 5);
    //fprintf( * fp4, "\ndefault_Alltoallv,%lf", ttime / 5);
  }

}

int main(int argc, char * argv[]) {
    MPI_Init( & argc, & argv);
    int len = 7;
    char name[MPI_MAX_PROCESSOR_NAME];
    int P = atoi(argv[2]);
    int ppn = atoi(argv[3]);
    if (argc < 4) {
        printf("4 parameters are needed\n");
        exit(0);
    }

    FILE * fp1, * fp2, * fp3, * fp4;
    fp1 = fopen("plot_Bcast.csv", "a");

    fp2 = fopen("plot_Gather.csv", "a");

    fp3 = fopen("plot_Reduce.csv", "a");

    fp4 = fopen("plot_Alltoallv.csv", "a");

    int i = 0;
    MPI_Comm first_comm, second_comm;

    int cpus, myrank;
    MPI_Comm_size(MPI_COMM_WORLD, & cpus);
    MPI_Comm_rank(MPI_COMM_WORLD, & myrank);

    MPI_Status status;

    int newrank, secondnewrank;
    int colour = myrank / ppn;
    int first_comm_size, second_comm_size;

    double comm_start = MPI_Wtime(); //timing the comm split time 
    //making subcommunicator of all leader of each subcommunicator
    MPI_Comm_split(MPI_COMM_WORLD, myrank % ppn == 0, myrank, & second_comm);
    MPI_Comm_rank(second_comm, & secondnewrank);
    MPI_Comm_size(second_comm, & second_comm_size);

    //making subcommunicator of all ranks in a node
    MPI_Comm_split(MPI_COMM_WORLD, colour, myrank, & first_comm);
    MPI_Comm_rank(first_comm, & newrank);
    MPI_Comm_size(first_comm, & first_comm_size);
    double comm_end = MPI_Wtime();
    double total_comm_time = comm_end - comm_start;

    int D_inKB = atoi(argv[1]);
    int iteration = 5;
    int nodes = cpus;
    int noOfEle = D_inKB * 128;
    double buf[D_inKB * 128];
    //randomly initializing the buffer
    for (int i = 0; i < noOfEle; i++) {
        buf[i] = rand() % 2545 + 1;
    }

    MPI_Bcast_default( & fp1, buf, myrank, noOfEle);
    MPI_Bcast_optimized( & fp1, buf, myrank, noOfEle, ppn, first_comm, second_comm, total_comm_time);
    MPI_Gather_default( & fp2, buf, myrank, cpus, noOfEle);
    MPI_Gather_optimized( & fp2, buf, myrank, cpus, noOfEle, first_comm, second_comm, ppn, P, total_comm_time);
    MPI_Reduce_default( & fp3, buf, myrank, cpus, noOfEle);
    MPI_Reduce_optimized( & fp3, buf, myrank, cpus, noOfEle, first_comm, second_comm, ppn, P, total_comm_time);
    MPI_Alltoallv_default( & fp4, myrank, P, ppn, noOfEle);

    MPI_Finalize();
    return 0;
}
