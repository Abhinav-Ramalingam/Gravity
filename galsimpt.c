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
    double *vx_sp;
    double *vy_sp;
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
    double *vx_sp = data->vx_sp;
    double *vy_sp = data->vy_sp;

    double *local_vx = (double *)calloc(N, sizeof(double));
    double *local_vy = (double *)calloc(N, sizeof(double));

    if (!local_vx || !local_vy) {
        printf("Memory allocation error in thread!\n");
        return NULL;
    }

    for (int i = start; i < end; i++) {
         
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

            // Update local velocities for particle j
            local_vx[j] -= dmi * fx;
            local_vy[j] -= dmi * fy;
        }

        // Update local velocities for particle i
        local_vx[i] += deltaT * a_x;
        local_vy[i] += deltaT * a_y;
    }


    // Safely add local velocities to shared velocity arrays
    pthread_mutex_lock(&mutex);
    for (int i = 0; i < N; i++) {
        vx_sp[i] += local_vx[i];
        vy_sp[i] += local_vy[i];
    }
    pthread_mutex_unlock(&mutex);

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

    double *vy_shared = (double *)calloc(N, sizeof(double));
    double *vx_shared = (double *)calloc(N, sizeof(double));

    pthread_t threads[n_threads];
    ThreadData thread_data[n_threads];

    int rows = 2 * N / (n_threads * (n_threads + 1));
    int end = 0;

    for (int ts = 0; ts < timesteps; ts++) {
        end = 0;
        for (int t = 0; t < n_threads; t++) {
            int rows_t = rows * (t + 1);
            int start = end;
            end = (t == n_threads - 1) ? N : start + rows_t;

            //printf("%d %d\n",start, end);

            thread_data[t].N = N;
            thread_data[t].particles = &particles;
            thread_data[t].deltaT = deltaT;
            thread_data[t].start = start;
            thread_data[t].end = end;
            thread_data[t].vx_sp = vx_shared;
            thread_data[t].vy_sp = vy_shared;

            pthread_create(&threads[t], NULL, compute_forces, &thread_data[t]);
        }
        //printf("-------------------\n");

        for (int t = 0; t < n_threads; t++) {
            pthread_join(threads[t], NULL);
        }

        // After all threads finish, apply the shared velocity to the global velocity
        for (int i = 0; i < N; i++) {
            particles.vx[i] += vx_shared[i];
            particles.vy[i] += vy_shared[i];
            particles.x[i] += deltaT * particles.vx[i];
            particles.y[i] += deltaT * particles.vy[i];
        }

        // Reset shared velocities for the next timestep
        memset(vx_shared, 0, N * sizeof(double));
        memset(vy_shared, 0, N * sizeof(double));
    }

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

    if (particles.x) free(particles.x);
    if (particles.y) free(particles.y);
    if (particles.mass) free(particles.mass);
    free(particles.vx);
    free(particles.vy);
    free(particles.brightness);

    printf("Time for DS to File(s): %lf\n", (double)(omp_get_wtime() - start));

    return 0;
}
