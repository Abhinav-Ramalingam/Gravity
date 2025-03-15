#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

// Constants
#define G(N) (100.0 / (N))
#define EPSILON_0 1e-3

// Particle structure
typedef struct
{
    double x, y;
    double mass;
    double vx, vy;
    double brightness;
} Particle;

// Function prototypes
void read_input(char *filename, Particle **particles, int *N);
void write_output(char *filename, Particle *particles, int N);
void gravitational_force(Particle p1, Particle p2, double G, double *fx, double *fy);
void gravitational_force_net(Particle *particles, int N, int i, double G, double *net_fx, double *net_fy);
void acceleration(double net_fx, double net_fy, double mass, double *ax, double *ay);
void velocity(double *vx, double *vy, double ax, double ay, double delta_t);
void position(double *x, double *y, double vx, double vy, double delta_t);
static double get_wall_seconds();  

// Main function
int main(int argc, char *argv[])
{
    if (argc != 6)
    {
        fprintf(stderr, "Usage: %s N filename nsteps delta_t graphics\n", argv[0]);
        return 1;
    }

    int N = atoi(argv[1]);
    char *input_file = argv[2];
    int num_steps = atoi(argv[3]);
    double delta_t = atof(argv[4]);
    int graphics_flag = atoi(argv[5]);

    Particle *particles = (Particle *)malloc(N * sizeof(Particle));
    Particle *new_particles = (Particle *)malloc(N * sizeof(Particle));
    if (particles == NULL || new_particles == NULL)
    {
        fprintf(stderr, "Memory allocation failed!\n");
        return 1;
    }

    // Read input file (particle data)
    read_input(input_file, &particles, &N);

    // Start the wall clock timer
    double start_time = get_wall_seconds();

    // Calculate net force, acceleration, velocity, and position on every particle (for each time step)
    for (int t = 0; t < num_steps; t++)
    {
        for (int i = 0; i < N; i++)
        {
            double net_fx = 0.0, net_fy = 0.0;
            double ax = 0.0, ay = 0.0;
            double vx = particles[i].vx, vy = particles[i].vy;
            double x = particles[i].x, y = particles[i].y;

            // Calculate net force on particle i
            gravitational_force_net(particles, N, i, G(N), &net_fx, &net_fy);

            // Calculate acceleration on particle i
            acceleration(net_fx, net_fy, particles[i].mass, &ax, &ay);

            // Update velocity of particle i
            velocity(&vx, &vy, ax, ay, delta_t);

            // Update position of particle i
            position(&x, &y, vx, vy, delta_t);

            // Store the updated position and velocity back to the particle
            new_particles[i].vx = vx;
            new_particles[i].vy = vy;
            new_particles[i].x = x;
            new_particles[i].y = y;
        }

        // Update particle's new position and velocities
        for (int i = 0; i < N; i++)
        {
          particles[i].vx = new_particles[i].vx;
          particles[i].vy = new_particles[i].vy;
          particles[i].x  = new_particles[i].x;
          particles[i].y  = new_particles[i].y;
        }
    }

    // Stop the wall clock timer
    double end_time = get_wall_seconds();
    printf("Elapsed time: %f seconds\n", end_time - start_time);

    // Write output to file
    write_output("resultno.gal", particles, N);

    free(particles);

    return 0;
}

// Function definitions

// FUNCTION: Read the Particle data from the input file (binary)
void read_input(char *filename, Particle **particles, int *N)
{
    FILE *file = fopen(filename, "rb");
    if (file == NULL)
    {
        fprintf(stderr, "Error opening file %s\n", filename);
        exit(1);
    }

    fread(*particles, sizeof(Particle), *N, file);
    fclose(file);
}

// FUNCTION: Write the final Particle data to the output file (binary)
void write_output(char *filename, Particle *particles, int N)
{
    FILE *file = fopen(filename, "wb");
    if (file == NULL)
    {
        fprintf(stderr, "Error opening file %s\n", filename);
        exit(1);
    }

    fwrite(particles, sizeof(Particle), N, file);
    fclose(file);
}

// FUNCTION: Compute the gravitational force between two particles
void gravitational_force(Particle p1, Particle p2, double G, double *fx, double *fy)
{
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    double r2 = pow(dx, 2) + pow(dy, 2);
    double r = sqrt(r2);

    double denominator = pow((r + EPSILON_0), 3);
    double force_magnitude = G * p1.mass * p2.mass / denominator;

    *fx = -force_magnitude * dx;
    *fy = -force_magnitude * dy;
}

// FUNCTION: Compute the net gravitational force on particle i by all other particles
void gravitational_force_net(Particle *particles, int N, int i, double G, double *net_fx, double *net_fy)
{
    *net_fx = 0.0;
    *net_fy = 0.0;

    for (int j = 0; j < N; j++)
    {
        if (j != i)
        {
            double fx, fy;
            gravitational_force(particles[i], particles[j], G, &fx, &fy);
            *net_fx += fx;
            *net_fy += fy;
        }
    }
}

// AUXILLIARY FUNCTIONS
void acceleration(double net_fx, double net_fy, double mass, double *ax, double *ay)
{
    *ax = net_fx / mass;
    *ay = net_fy / mass;
}

void velocity(double *vx, double *vy, double ax, double ay, double delta_t)
{
    *vx += ax * delta_t;
    *vy += ay * delta_t;
}

void position(double *x, double *y, double vx, double vy, double delta_t)
{
    *x += vx * delta_t;
    *y += vy * delta_t;
}

// Utility FUNCTION: Get wall clock time in seconds
static double get_wall_seconds()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
    return seconds;
}
