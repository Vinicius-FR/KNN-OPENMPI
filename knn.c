#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>

// Definir a estrutura de coordenadas da partícula
typedef struct {
    double x, y, z;
} Particle;

// Função para calcular a distância entre duas partículas
double distance(const Particle* p1, const Particle* p2) {
    double dx = p1->x - p2->x;
    double dy = p1->y - p2->y;
    double dz = p1->z - p2->z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

// Função para encontrar os k vizinhos mais próximos de uma partícula
void find_nearest_neighbors(const Particle* particles, int num_particles, int k, int particle_index, int* neighbors) {
    double* distances = (double*)malloc(num_particles * sizeof(double));
    if (distances == NULL) {
        fprintf(stderr, "Erro: Falha na alocação de memória para distâncias\n");
        exit(1);
    }

    // Calcular distâncias entre a partícula de referência e todas as outras
    for (int i = 0; i < num_particles; i++) {
        if (i != particle_index) {
            distances[i] = distance(&particles[particle_index], &particles[i]);
        } else {
            distances[i] = INFINITY; // Exclui a própria partícula
        }
    }

    // Encontrar os k menores índices das distâncias
    for (int i = 0; i < k; i++) {
        double min_distance = INFINITY;
        int min_index = -1;

        for (int j = 0; j < num_particles; j++) {
            if (distances[j] < min_distance) {
                min_distance = distances[j];
                min_index = j;
            }
        }

        neighbors[i] = min_index;
        distances[min_index] = INFINITY; // Marcar a distância como infinito para evitar recontagem
    }

    free(distances);
}

int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Uso: %s <arquivo_de_entrada> <k> <arquivo_de_saida>\n", argv[0]);
        return 1;
    }

    char *input_filename = argv[1];
    int k = atoi(argv[2]);
    char *output_filename = argv[3];

    // Abrir o arquivo de entrada
    FILE *input_file = fopen(input_filename, "r");
    if (input_file == NULL) {
        fprintf(stderr, "Erro: Não foi possível abrir o arquivo de entrada.\n");
        return 1;
    }

    int num_particles = 0; // Inicializa o número de partículas
    int current_char = EOF;
    int previous_char = EOF;

    while ((current_char = fgetc(input_file)) != EOF) {
        if (current_char == '\n') {
            num_particles++;
        }
        previous_char = current_char;
    }

    if (previous_char != '\n' && previous_char != EOF) {
        // Contar a última linha se ela não terminar com '\n'
        num_particles++;
    }

    if (ferror(input_file)) {
        perror("Erro de leitura"); // Trata erros de leitura, se houverem
    }

    // Voltar ao início do arquivo para ler as coordenadas das partículas
    rewind(input_file);

    // Alocação de memória para armazenar as partículas
    Particle *particles = (Particle*)malloc(num_particles * sizeof(Particle));
    if (particles == NULL) {
        fprintf(stderr, "Erro: Falha na alocação de memória para as partículas.\n");
        fclose(input_file);
        return 1;
    }

    // Ler as coordenadas das partículas
    for (int i = 0; i < num_particles; i++) {
        if (fscanf(input_file, "%lf %lf %lf", &particles[i].x, &particles[i].y, &particles[i].z) != 3) {
            fprintf(stderr, "Erro: Falha ao ler as coordenadas da partícula %d no arquivo de entrada.\n", i);
            fclose(input_file);
            free(particles);
            return 1;
        }

    }

    fclose(input_file);

    // Abrir o arquivo de saída para escrever os vizinhos mais próximos
    FILE *output_file = fopen(output_filename, "w");
    if (output_file == NULL) {
        fprintf(stderr, "Erro: Não foi possível criar o arquivo de saída.\n");
        free(particles);
        return 1;
    }
    
    double elapsed = 0;
    struct timeval t1, t2; // Iniciar a contagem de tempo

    // Encontrar os k vizinhos mais próximos para cada partícula
    for (int i = 0; i < num_particles; i++) {
        gettimeofday(&t1, NULL);
        int *neighbors = (int*)malloc(k * sizeof(int));
        if (neighbors == NULL) {
            fprintf(stderr, "Erro: Falha na alocação de memória para vizinhos.\n");
            fclose(output_file);
            free(particles);
            return 1;
        }

        find_nearest_neighbors(particles, num_particles, k, i, neighbors);

        gettimeofday(&t2, NULL); // Parar a contagem

        double dt = (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec) / 1e6;
        elapsed += dt;

        // Escrever os k vizinhos no arquivo de saída
        for (int j = 0; j < k; j++) {
            fprintf(output_file, "%d ", neighbors[j]);
        }
        fprintf(output_file, "\n");

        free(neighbors);
    }

    printf("%.12f\n", elapsed);

    fclose(output_file);
    free(particles);

    return 0;
}