#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <pthread.h>
#include <omp.h>

#define eps 1e-3

typedef struct {
    double *x, *y, *mass, *vx, *vy, *brightness;
} ParticleData;

typedef struct {
    ParticleData *particles;
    double deltaT;
    int N, start, end;
} ThreadData;

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

void *compute_forces(void *arg) {
    ThreadData *data = (ThreadData *)arg;
    double deltaT = data->deltaT;
    int N = data->N;
    double G = 100.0 / N;
    int start = data->start, end = data->end;

    double *x = data->particles->x;
    double *y = data->particles->y;
    double *mass = data->particles->mass;
    double *vx = data->particles->vx;
    double *vy = data->particles->vy;

    double *local_vx = (double *)calloc(N, sizeof(double));
    double *local_vy = (double *)calloc(N, sizeof(double));

    if (!local_vx || !local_vy) {
        printf("Memory allocation error in thread!\n");
        return NULL;
    }

    for (int i = start; i < end; i++) {
        printf(".");
        double mi = mass[i];
        double dmi = deltaT * mi;
        double xi = x[i], yi = y[i];
        double a_x = 0.0, a_y = 0.0;

        for (int j = i + 1; j < N; j++) {
            double dx = x[j] - xi;
            double dy = y[j] - yi;
            double rij2 = dx * dx + dy * dy;
            double radius = sqrt(rij2) + eps;
            double force = G / (radius * radius * radius);
            double fx = force * dx;
            double fy = force * dy;

           
            a_x += fx * mass[j];
            a_y += fy * mass[j];


            local_vx[j] -= dmi * fx;
            local_vy[j] -= dmi * fy;
        }
        

        local_vx[i] += deltaT * a_x;
        local_vy[i] += deltaT * a_y;

        pthread_mutex_lock(&mutex);
        printf("\n");
        vx[i] += local_vx[i];
        vy[i] += local_vy[i];
        x[i] += deltaT * vx[i];
        y[i] += deltaT * vy[i];
        pthread_mutex_unlock(&mutex);
    }

    free(local_vx);
    free(local_vy);
    return NULL;
}

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

    clock_t start = omp_get_wtime();

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

    printf("Time for File to DS(s): %lf\n", (double)(omp_get_wtime() - start));
    start = omp_get_wtime();

    pthread_t threads[n_threads];
    ThreadData thread_data[n_threads];

    int rows = 2 * N / (n_threads * (n_threads + 1));
    int end = 0;

    int ts;
    for (ts = 0; ts < timesteps; ts++) {
        for (int t = 0; t < n_threads; t++) {
            int rows_t = rows * (t + 1);
            int start = end;
            end = (t == n_threads - 1) ? N : start + rows_t;

            thread_data[t].N = N;
            thread_data[t].particles = &particles;
            thread_data[t].deltaT = deltaT;
            thread_data[t].start = start;
            thread_data[t].end = end;

            pthread_create(&threads[t], NULL, compute_forces, &thread_data[t]);
        }

        for (int t = 0; t < n_threads; t++) {
            pthread_join(threads[t], NULL);
        }
    }

    printf("Current timestep is: %d\n", ts);

    printf("Time for Processing(s): %lf\n", (double)(omp_get_wtime() - start));
    start = omp_get_wtime();


    FILE *outputfd = fopen("result.gal", "wb");
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

    printf("Freeing particles.x at address: %p\n", particles.x);
    free(particles.x);
    printf("Freeing particles.y at address: %p\n", particles.y);
    free(particles.y);
    printf("Freeing particles.mass at address: %p\n", particles.mass);
    free(particles.mass);
    free(particles.vx);
    free(particles.vy);
    free(particles.brightness);


    printf("Time for DS to File(s): %lf\n", (double)(omp_get_wtime() - start));

    return 0;
}
