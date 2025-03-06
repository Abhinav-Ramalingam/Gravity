#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define eps 1e-3

int main(int ac, char * av[]) {
    if(ac != 6){
        printf("Correct Format: %s N inputfile timesteps deltaT graphics\n",av[0]);
    }

    const int N = atoi(av[1]); 
    const double G = 100.0/N;

    const char * inputfile = av[2];
    const int timesteps = atoi(av[3]);
    const double deltaT = atof(av[4]);

    clock_t start;
    start = clock();

    //Phase 1: File to Datastucture
    int dsize = N * sizeof(double);
    int dssize = (dsize << 1) + (dsize << 2);
    double * particles = (double *) malloc(dssize);
    
    int inputfd=open(inputfile,O_RDONLY);
    if(inputfd == -1) {printf("File open error!\n"); return 1;}
    
    int inputsize = read(inputfd, particles, dssize);
    close(inputfd);

    printf("Time for File to DS(s): %lf\n", (double) (clock() - start) / CLOCKS_PER_SEC);
    start = clock();

    //Phase 2: Processing
    int instant = 0, i = 0, j = 0;
    double * accel = (double *) calloc(N << 1, sizeof(double));
    while(instant < timesteps){
        for(i = 0; i < N; i++){
            int row2 = i << 1;
            int row = row2 + (i << 2);
            double mi = particles[row + 2];
            double x = particles[row], y = particles[row + 1];
            double v_x = particles[row + 3], v_y = particles[row + 4];
            double a_x = 0.0, a_y = 0.0;

            for(j = i + 1; j < N; j++){
                int col2 = j << 1;
                int col = col2 + (j << 2);
                double dx = x - particles[col];
                double dy = y - particles[col + 1];
                double mj = particles[col + 2];
                double rij = sqrt(dx * dx + dy * dy);
                double rije = rij + eps;
                double radius = rije * rije * rije;
                double force = G / radius;
                
                //Optimization: Newton's Third Law
                double fx = force * dx;
                double fy = force * dy;
                a_x -= fx * mj;
                a_y -= fy * mj;
                accel[col2] += fx * mi;
                accel[col2 + 1] += fy * mi;
            }

            accel[row2] += a_x;
            accel[row2 + 1] += a_y;

            particles[row + 4] += deltaT * accel[row2 + 1];
            particles[row + 3] += deltaT * accel[row2];
            particles[row] += deltaT * particles[row + 3];
            particles[row + 1] += deltaT * particles[row + 4];

        }


        //reset accel
        memset(accel, 0, N * sizeof(double) << 1); 

        instant++;
    }

    printf("Time for Processing(s): %lf\n", (double) (clock() - start) / CLOCKS_PER_SEC);
    start = clock();

    //Phase 3: Datastucture to File
    int outputfd = open("results.gal",O_RDWR | O_CREAT, 0666);
    if(outputfd == -1) {printf("File create error!\n"); return 1;}

    int outputsize = write(outputfd, particles, dssize);
    close(outputfd);

    free(particles);
    free(accel);

    printf("Time for DS to File(s): %lf\n", (double) (clock() - start) / CLOCKS_PER_SEC);

    return 0;
}