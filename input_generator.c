#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Function to generate a random bit sequence
void generate_random_sequence(char *sequence, int length) {
    for (int i = 0; i < length; i++) {
        sequence[i] = (rand() % 2) + '0';
    }
    sequence[length] = '\0';
}

// Function to flip one random bit in a sequence
void flip_random_bit(char *sequence, int length) {
    int bit_to_flip = rand() % length;
    sequence[bit_to_flip] = (sequence[bit_to_flip] == '0') ? '1' : '0';
}

// Function to generate and write sequences to a file
void generate_and_write_sequences(int num_seq, int seq_len, const char *filename, int percentage) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        perror("Failed to open file");
        exit(EXIT_FAILURE);
    }

    // Calculate the number of sequences to flip
    int num_to_flip = (percentage > 0) ? (num_seq * percentage / 100) : 0;
    int num_original = num_seq - num_to_flip;

    fprintf(file, "%d %d\n", seq_len, num_seq);

    char **sequences = (char **)malloc(num_original * sizeof(char *));
    if (sequences == NULL) {
        perror("Failed to allocate memory for sequences");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    // Generate and write the original sequences
    for (int i = 0; i < num_original; i++) {
        sequences[i] = (char *)malloc((seq_len + 1) * sizeof(char));
        if (sequences[i] == NULL) {
            perror("Failed to allocate memory for a sequence");
            for (int j = 0; j < i; j++) {
                free(sequences[j]);
            }
            free(sequences);
            fclose(file);
            exit(EXIT_FAILURE);
        }

        generate_random_sequence(sequences[i], seq_len);
        fprintf(file, "%s\n", sequences[i]);
    }

    // Generate flipped sequences from the original ones
    for (int i = 0; i < num_to_flip; i++) {
        int seq_index = rand() % num_original;
        char flipped_sequence[seq_len + 1];

        // Copy the original sequence
        for (int j = 0; j < seq_len; j++) {
            flipped_sequence[j] = sequences[seq_index][j];
        }
        flipped_sequence[seq_len] = '\0';

        // Flip one bit in the sequence
        flip_random_bit(flipped_sequence, seq_len);
        fprintf(file, "%s\n", flipped_sequence);
    }

    // Free allocated memory
    for (int i = 0; i < num_original; i++) {
        free(sequences[i]);
    }
    free(sequences);
    fclose(file);

    printf("%d sequences of length %d written to %s\n", num_seq, seq_len, filename);
}

// Generate num_seq sequences of length seq_len and write them to a file
// If a percentage is given, first generate num_seq - num_seq * percentage / 100 sequences 
// and then flip one random bit in num_seq * percentage / 100 sequences from the original ones
// this ensures that there will be at least num_seq * percentage / 100 pairs with Hamming distance 1
int main(int argc, char *argv[]) {
    if (argc < 3 || argc > 4) {
        fprintf(stderr, "Usage: %s <sequence_length> <number_of_sequences> [percentage]\n", argv[0]);
        return EXIT_FAILURE;
    }

    int seq_len = atoi(argv[1]);
    int num_seq = atoi(argv[2]);
    int percentage = 0;

    if (argc == 4) {
        percentage = atoi(argv[3]);
        if (percentage < 0 || percentage > 100) {
            fprintf(stderr, "Percentage must be between 0 and 100.\n");
            return EXIT_FAILURE;
        }
    }

    if (seq_len <= 0 || num_seq <= 0) {
        fprintf(stderr, "Both sequence length and number of sequences must be positive integers.\n");
        return EXIT_FAILURE;
    }

    srand(time(NULL));
    const char *filename = "input_sequences.txt";

    clock_t start_time = clock();
    generate_and_write_sequences(num_seq, seq_len, filename, percentage);
    clock_t end_time = clock();

    double time_taken = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    printf("Time taken to generate and write sequences: %f seconds\n", time_taken);

    return 0;
}
