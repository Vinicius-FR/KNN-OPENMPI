#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// OBS: Para corrigir as observações feitas e tornar o código mais otimizado,
// alterei a lógica base de resolução. Ao invés de utilizar somente listas,
// optei por organizar as partículas em uma árvore binária. Portanto, a parale-
// lização foi feita baseada nesse novo modelo de resolução.

// Definir a estrutura de coordenadas da partícula
typedef struct {
    double x, y, z;
    int index;
} Particle;

// Definir a estrutura de um nó da árvore de partículas
typedef struct Node {
    Particle particle;
    struct Node* left;
    struct Node* right;
} Node;

// Função para calcular a distância entre duas partículas
double distance(const Particle p1, const Particle p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    double dz = p1.z - p2.z;
    return sqrt(dx * dx + dy * dy + dz * dz);
}

// Função para inserir uma partícula na árvore
Node* insert(Node* root, Particle particle, int depth) {
    if (root == NULL) {
        Node* new_node = (Node*)malloc(sizeof(Node));
        new_node->particle = particle;
        new_node->left = new_node->right = NULL;
        return new_node;
    }

    int cd = depth % 3; // Alternar entre as coordenadas x e y e z

    // Ponto à esquerda
    if ((cd == 0 && particle.x < root->particle.x) || (cd == 1 && particle.y < root->particle.y) || (cd == 2 && particle.z < root->particle.z))  {
        root->left = insert(root->left, particle, depth + 1);
    // Ponto à direita
    } else {
        root->right = insert(root->right, particle, depth + 1);
    }

    return root;
}

// Função para encontrar os k vizinhos mais próximos de uma partícula
void find_nearest_neighbors(Node* root, Particle particle, int depth, int k, Node* nearest[], double distances[], int* count) {
    if (root == NULL) {
        return;
    }
	// Calcular distância
    double dist = distance(root->particle, particle);
    int cd = depth % 3;

    // Garante que a distância é calculada apenas uma vez para cada par único de partículas
    if (particle.index < root->particle.index) {
        dist = distance(root->particle, particle);
    } else if (particle.index > root->particle.index) {
        dist = distance(particle, root->particle);
    }

    // Atualizar a lista de vizinhos mais próximos se a distância for menor que a maior distância encontrada.
    if (*count < k || dist < distances[*count - 1]) { 
        if (*count < k) { // Enquanto não tiver achado todos as k partículas mais próximas
            (*count)++;
        }

        int i = *count - 1;

        while (i > 0 && distances[i - 1] > dist) {
            distances[i] = distances[i - 1];
            nearest[i] = nearest[i - 1];
            i--;
        }

        distances[i] = dist;
        nearest[i] = root;
    }

    // Identificar qual lado da árvore explorar primeiro
    if ((cd == 0 && particle.x < root->particle.x) || (cd == 1 && particle.y < root->particle.y) || (cd == 2 && particle.z < root->particle.z) ) {
        find_nearest_neighbors(root->left, particle, depth + 1, k, nearest, distances, count);

        if (*count < k || (cd == 0 && particle.x + distances[*count - 1] >= root->particle.x) ||
            (cd == 1 && particle.y + distances[*count - 1] >= root->particle.y) || (cd == 2 && particle.z + distances[*count - 1] >= root->particle.z)) {
            find_nearest_neighbors(root->right, particle, depth + 1, k, nearest, distances, count);
        }
    } else {
        find_nearest_neighbors(root->right, particle, depth + 1, k, nearest, distances, count);

        if (*count < k || (cd == 0 && particle.x - distances[*count - 1] <= root->particle.x) ||
            (cd == 1 && particle.y - distances[*count - 1] <= root->particle.y) || (cd == 2 && particle.z - distances[*count - 1] <= root->particle.z)) {
            find_nearest_neighbors(root->left, particle, depth + 1, k, nearest, distances, count);
        }
    }
}

int main(int argc, char* argv[]) {

	if (argc != 4) {
        fprintf(stderr, "Uso: %s <arquivo_de_entrada> <k> <arquivo_de_saida>\n", argv[0]);
        return 1;
    }

    char *input_filename = argv[1];
    int k = atoi(argv[2]) + 1;
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
        particles[i].index = i;
        if (fscanf(input_file, "%lf %lf %lf", &particles[i].x, &particles[i].y, &particles[i].z) != 3) {
            fprintf(stderr, "Erro: Falha ao ler as coordenadas da partícula %d no arquivo de entrada.\n", i);
            fclose(input_file);
            free(particles);
            return 1;
        }

    }

    fclose(input_file);

    MPI_Init(&argc, &argv);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int particles_per_process = num_particles / world_size;
    int remainder_particles = num_particles % world_size;

    // Distribuir as partículas entre os processos
    int start_index = world_rank * particles_per_process;
    int end_index = (world_rank + 1) * particles_per_process;

    if (world_rank == world_size - 1) {
        // Distribuir as partículas restantes para o último processo
        end_index += remainder_particles;
    }

    // Cada processo lê e armazena suas partículas
    Particle *local_particles = (Particle*)malloc(particles_per_process * sizeof(Particle));
    for (int i = start_index; i < end_index; i++) {
        local_particles[i - start_index] = particles[i];
    }

    // Sincronizar todos os processos após a leitura
    MPI_Barrier(MPI_COMM_WORLD);

    Node* root = NULL;

    // Início do timer
    clock_t start = clock();

    // Construção da árvore paralelizada
    // Loop paralelizado com MPI para inserir as partículas
    for (int i = 0; i < particles_per_process; i++) {
        root = insert(root, local_particles[i], 0);
    }

    // Sincronizar todos os processos após a construção da árvore
    MPI_Barrier(MPI_COMM_WORLD);

    int **neighbors_indexes = (int **)malloc(particles_per_process * sizeof(int*));
    for(int i = 0; i < particles_per_process; i++) neighbors_indexes[i] = (int *)malloc(k * sizeof(int));

    // Loop paralelizado com MPI para calcular os vizinhos mais próximos
    for (int i = 0; i < particles_per_process; i++) {
        Node* nearest[k];
        double distance[k];
        int count = 0;

        Particle particle = local_particles[i];
        find_nearest_neighbors(root, particle, 0, k, nearest, distance, &count);

        for (int j = 1; j < k; j++) {
            neighbors_indexes[i][j-1] = nearest[j]->particle.index;
        }
    }

    clock_t end = clock();

    double tempoExecucao = ((double) (end - start)) / CLOCKS_PER_SEC;

    // Sincronizar todos os processos após o cálculo dos vizinhos
    MPI_Barrier(MPI_COMM_WORLD);

    // Juntar os resultados de todos os processos no processo mestre (rank 0)
    if (world_rank == 0) {
        int **all_neighbors_indexes = (int **)malloc(num_particles * sizeof(int*));
        for(int i = 0; i < num_particles; i++) all_neighbors_indexes[i] = (int *)malloc(k * sizeof(int));

        // Copiar os resultados do processo mestre
        for (int i = 0; i < particles_per_process; i++) {
            for (int j = 0; j < k-1; j++) {
                all_neighbors_indexes[i][j] = neighbors_indexes[i][j];
            }
        }

        // Receber resultados dos outros processos
        for (int i = 1; i < world_size; i++) {
            int recv_particles_per_process = particles_per_process;
            if (i == world_size - 1) {
                recv_particles_per_process += remainder_particles;
            }

            for (int j = 0; j < recv_particles_per_process; j++) {
                MPI_Recv(all_neighbors_indexes[i * particles_per_process + j], k-1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        // Imprimir ou salvar os resultados
        // Abrir o arquivo de saída para escrever os vizinhos mais próximos
        FILE *output_file = fopen(output_filename, "w");
        if (output_file == NULL) {
            fprintf(stderr, "Erro: Não foi possível criar o arquivo de saída.\n");
            free(particles);
            return 1;
        }

        // Escrevendo arquivo de saída
        for (int i = 0; i < num_particles; i++) {
            for(int j = 0; j < k-1; j++){
                fprintf(output_file, "%d", neighbors_indexes[i][j]);

                // Se não está printando o último ponto, adiciona um espaço após o índice
                if(j < k-1){
                    fprintf(output_file, " ");
                }
            }
            fprintf(output_file, "\n");
        }

        fclose(output_file);

        // Imprimir tempo de execução
        printf("%lf\n", tempoExecucao);

        // Liberar memória
        for (int i = 0; i < num_particles; i++) {
            free(all_neighbors_indexes[i]);
        }
        free(all_neighbors_indexes);
    } else {
        // Enviar resultados para o processo mestre
        for (int i = 0; i < particles_per_process; i++) {
            MPI_Send(neighbors_indexes[i], k-1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
    }

    // Liberar memória e finalizar MPI
    free(local_particles);
    for (int i = 0; i < particles_per_process; i++) {
        free(neighbors_indexes[i]);
    }
    free(neighbors_indexes);

    MPI_Finalize();

    return 0;
}
