#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define eps 1e-3

typedef struct {
    double *x, *y, *mass, *vx, *vy;
} ParticleData;

int main(int ac, char *av[]) {
    if(ac != 6){
        printf("Correct Format: %s N inputfile timesteps deltaT graphics\n", av[0]);
        return 1;
    }

    const int N = atoi(av[1]); 
    const double G = 100.0 / N;
    const char *inputfile = av[2];
    const int timesteps = atoi(av[3]);
    const double deltaT = atof(av[4]);

    clock_t start;
    start = clock();

    // Phase 1: File to DataStructure (SoA)
    ParticleData particles;
    particles.x = (double *) malloc(N * sizeof(double));
    particles.y = (double *) malloc(N * sizeof(double));
    particles.mass = (double *) malloc(N * sizeof(double));
    particles.vx = (double *) malloc(N * sizeof(double));
    particles.vy = (double *) malloc(N * sizeof(double));

    if (!particles.x || !particles.y || !particles.mass || !particles.vx || !particles.vy) {
        printf("Memory allocation error!\n");
        return 1;
    }

    int inputfd = open(inputfile, O_RDONLY);
    if(inputfd == -1) {
        printf("File open error!\n");
        return 1;
    }

    read(inputfd, particles.x, N * sizeof(double));
    read(inputfd, particles.y, N * sizeof(double));
    read(inputfd, particles.mass, N * sizeof(double));
    read(inputfd, particles.vx, N * sizeof(double));
    read(inputfd, particles.vy, N * sizeof(double));

    close(inputfd);

    printf("Time for File to DS(s): %lf\n", (double) (clock() - start) / CLOCKS_PER_SEC);
    start = clock();

    // Phase 2: Processing
    int instant = 0, i = 0, j = 0;
    while (instant < timesteps) {
        for (i = 0; i < N; i++) {
            double mi = particles.mass[i];
            double dmi = deltaT * mi;
            double x = particles.x[i], y = particles.y[i];
            double a_x = 0.0, a_y = 0.0;

            for (j = i + 1; j < N; j++) {
                double mj = particles.mass[j];
                double dx = particles.x[j] - x;
                double dy = particles.y[j] - y;
                double rij2 = dx * dx + dy * dy;
                double radius = (__builtin_sqrt(rij2) + eps);
                double radius_cubed = radius * radius * radius;
                double force = G / radius_cubed;
                
                double fx = force * dx;
                double fy = force * dy;

                // Calculate the acceleration directly for particle i
                a_x += fx * mj;
                a_y += fy * mj;

                // Apply the force to particle j (Newton's Third Law)
                particles.vx[j] -= dmi * fx;
                particles.vy[j] -= dmi * fy;
            }

            // Update velocities and positions directly
            particles.vx[i] += deltaT * a_x;
            particles.vy[i] += deltaT * a_y;
            particles.x[i] += deltaT * particles.vx[i];
            particles.y[i] += deltaT * particles.vy[i];
        }
        instant++;
    }

    printf("Time for Processing(s): %lf\n", (double) (clock() - start) / CLOCKS_PER_SEC);
    start = clock();

    // Phase 3: DataStructure to File
    int outputfd = open("results.gal", O_RDWR | O_CREAT, 0666);
    if(outputfd == -1) {
        printf("File create error!\n");
        return 1;
    }

    write(outputfd, particles.x, N * sizeof(double));
    write(outputfd, particles.y, N * sizeof(double));
    write(outputfd, particles.mass, N * sizeof(double));
    write(outputfd, particles.vx, N * sizeof(double));
    write(outputfd, particles.vy, N * sizeof(double));

    close(outputfd);

    free(particles.x);
    free(particles.y);
    free(particles.mass);
    free(particles.vx);
    free(particles.vy);

    printf("Time for DS to File(s): %lf\n", (double) (clock() - start) / CLOCKS_PER_SEC);

    return 0;
}
