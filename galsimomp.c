#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#define eps 1e-3

typedef struct {
    double *x, *y, *mass, *vx, *vy, *brightness;
} ParticleData;

int main(int argc, char *argv[]) {
    if (argc != 7) {
        printf("Usage: %s N inputfile timesteps deltaT graphics n_threads\n", argv[0]);
        return -1;
    }

    const int N = atoi(argv[1]);
    const char *inputfile = argv[2];
    const int timesteps = atoi(argv[3]);
    const double deltaT = atof(argv[4]);
    const int n_threads = atoi(argv[6]);

    omp_set_num_threads(n_threads);

    double start_time = omp_get_wtime();

    ParticleData particles;
    particles.x = (double *)malloc(N * sizeof(double));
    particles.y = (double *)malloc(N * sizeof(double));
    particles.mass = (double *)malloc(N * sizeof(double));
    particles.vx = (double *)malloc(N * sizeof(double));
    particles.vy = (double *)malloc(N * sizeof(double));
    particles.brightness = (double *)malloc(N * sizeof(double));

    if (!particles.x || !particles.y || !particles.mass || !particles.vx || !particles.vy || !particles.brightness) {
        printf("Memory allocation error!\n");
        return -1;
    }

    FILE *inputfd = fopen(inputfile, "rb");
    if (!inputfd) {
        printf("File open error!\n");
        return -1;
    }

    for (int i = 0; i < N; i++) {
        fread(&particles.x[i], sizeof(double), 1, inputfd);
        fread(&particles.y[i], sizeof(double), 1, inputfd);
        fread(&particles.mass[i], sizeof(double), 1, inputfd);
        fread(&particles.vx[i], sizeof(double), 1, inputfd);
        fread(&particles.vy[i], sizeof(double), 1, inputfd);
        fread(&particles.brightness[i], sizeof(double), 1, inputfd);
    }
    fclose(inputfd);

    printf("Time for File to DS(s): %lf\n", (double)(omp_get_wtime() - start_time));
    start_time = omp_get_wtime();

    double *vx_shared = (double *)calloc(N, sizeof(double));
    double *vy_shared = (double *)calloc(N, sizeof(double));

    for (int ts = 0; ts < timesteps; ts++) {
        double G = 100.0 / N;

        // Compute forces using OpenMP with reduction on shared velocity changes
        #pragma omp parallel for schedule(static, 8) reduction(+:vx_shared[:N], vy_shared[:N])
        for (int i = 0; i < N; i++) {
            double mi = particles.mass[i];
            double dmi = deltaT * mi;
            double xi = particles.x[i], yi = particles.y[i];
            double a_x = 0.0, a_y = 0.0;

            for (int j = i + 1; j < N; j++) {
                double dx = particles.x[j] - xi;
                double dy = particles.y[j] - yi;
                double rij2 = dx * dx + dy * dy;
                double radius = sqrt(rij2) + eps;
                double force = G / (radius * radius * radius);
                double fx = force * dx;
                double fy = force * dy;

                a_x += fx * particles.mass[j];
                a_y += fy * particles.mass[j];

                vx_shared[j] -= dmi * fx;
                vy_shared[j] -= dmi * fy;
            }

            vx_shared[i] += deltaT * a_x;
            vy_shared[i] += deltaT * a_y;
        }

        // Apply velocity updates in parallel
        #pragma omp parallel for
        for (int i = 0; i < N; i++) {
            particles.vx[i] += vx_shared[i];
            particles.vy[i] += vy_shared[i];
        }

        // Update positions in parallel
        #pragma omp parallel for
        for (int i = 0; i < N; i++) {
            particles.x[i] += deltaT * particles.vx[i];
            particles.y[i] += deltaT * particles.vy[i];
        }

        // Reset shared velocity arrays
        memset(vx_shared, 0, N * sizeof(double));
        memset(vy_shared, 0, N * sizeof(double));
    }

    printf("Time for Processing(s): %lf\n", (double)(omp_get_wtime() - start_time));
    start_time = omp_get_wtime();

    FILE *outputfd = fopen("resultomp.gal", "wb");
    if (!outputfd) {
        printf("File create error!\n");
        return -1;
    }

    for (int i = 0; i < N; i++) {
        fwrite(&particles.x[i], sizeof(double), 1, outputfd);
        fwrite(&particles.y[i], sizeof(double), 1, outputfd);
        fwrite(&particles.mass[i], sizeof(double), 1, outputfd);
        fwrite(&particles.vx[i], sizeof(double), 1, outputfd);
        fwrite(&particles.vy[i], sizeof(double), 1, outputfd);
        fwrite(&particles.brightness[i], sizeof(double), 1, outputfd);
    }
    fclose(outputfd);

    free(particles.x);
    free(particles.y);
    free(particles.mass);
    free(particles.vx);
    free(particles.vy);
    free(particles.brightness);
    free(vx_shared);
    free(vy_shared);

    printf("Time for DS to File(s): %lf\n", (double)(omp_get_wtime() - start_time));

    return 0;
}
