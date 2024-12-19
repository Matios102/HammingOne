#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cuda_runtime.h>
#include <chrono>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#define CUDA_CHECK(call)                                                                                                  \
    {                                                                                                                     \
        cudaError_t err = call;                                                                                           \
        if (err != cudaSuccess)                                                                                           \
        {                                                                                                                 \
            std::cerr << "CUDA error in " << __FILE__ << ":" << __LINE__ << ": " << cudaGetErrorString(err) << std::endl; \
            exit(EXIT_FAILURE);                                                                                           \
        }                                                                                                                 \
    }

#define BitsInInt 32

// Srructure representing a node in the trie
struct Node {
    int *bit;  // bit value of the node
    int *bit0; // index of the child node with bit value 0
    int *bit1; // index of the child node with bit value 1
};


// Function to load the trie and the sequence vector from a file
void load_trie(const char *filename, std::vector<int> &bit, std::vector<int> &bit0, std::vector<int> &bit1, std::vector<int> &sequences, int &sequence_length, int &num_sequences) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: could not open file " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    // Add the root node which bit value does not matter
    bit.push_back(-1);
    bit0.push_back(-1);
    bit1.push_back(-1);


    // the first line represents sequence length and number of sequences
    std::string line;
    if (!std::getline(file, line)) {
        std::cerr << "Error: could not read trie" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::istringstream iss(line);
    if (!(iss >> sequence_length >> num_sequences)) {
        std::cerr << "Error: could not parse sequence length and number of sequences" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Seqences will be represented as 32 bit integers
    // num_ints is the number of integers needed to represent the sequence
    int num_ints = (sequence_length + BitsInInt - 1) / BitsInInt;

    // padding_zeros is the number of zeros needed to fill the last integer
    int padding_zeros = (num_ints * BitsInInt) - sequence_length;


    // Read the sequences
    for (int i = 0; i < num_sequences; i++) {
        std::string sequence;
        if (!(file >> sequence)) {
            std::cerr << "Error: could not read sequence " << i << std::endl;
            break;
        }


        if (sequence.length() != static_cast<std::string::size_type>(sequence_length)) {
            std::cerr << "Error: invalid sequence length (" << sequence.length()
                      << ") for sequence " << i << std::endl;
            continue;
        }

        // Convert the sequence to integers and store them in the sequences vector
        for (int j = 0; j < num_ints; j++) {
            int chunk = 0;
            int valid_bits = BitsInInt;

            if (j == num_ints - 1) {
                valid_bits = BitsInInt - padding_zeros;
            }

            for (int k = 0; k < valid_bits; k++) {
                char bit_char = sequence[j * BitsInInt + k];
                int bit = (bit_char == '1') ? 1 : 0;
                chunk |= bit << (valid_bits - 1 - k);
            }
            // Add the chunk to the sequences vector
            sequences.push_back(chunk);
        }

        int current = 0;

        // Build the trie
        for (int bit_pos = 0; bit_pos < sequence_length; ++bit_pos) {

            // Get the bit value at the current bit position
            char bit_char = sequence[bit_pos];
            int current_bit = (bit_char == '1') ? 1 : 0;

            // Get the index of the child node with the current bit
            int childIdx = (current_bit == 0) ? bit0[current] : bit1[current];

            // If the child does not exist, create it
            if (childIdx == -1) {
                childIdx = bit.size();
                bit.push_back(current_bit);
                bit0.push_back(-1);
                bit1.push_back(-1);

                if (current_bit == 0) {
                    bit0[current] = childIdx;
                } else {
                    bit1[current] = childIdx;
                }
            }

            current = childIdx;            
        }

    }
    file.close();
}

// Helper function to print the trie
void print_tire(std::vector<int> &bit, std::vector<int> &bit0, std::vector<int> &bit1) {
    for (size_t i = 0; i < bit.size(); i++) {
        std::cout << "Bit " << bit[i] << std::endl;
        std::cout << "  Bit0 " << bit0[i] << std::endl;
        std::cout << "  Bit1 " << bit1[i] << std::endl;
    }
}

// Kernel to find the pairs
__global__ void find_pairs_kernel(int *d_sequences, const int *d_bit, const int *d_bit0, const int *d_bit1, int l, int n, int *d_sequence_idx, int *d_flipped_bit_position, int *d_pairs_found)
{
    // Get the sequence index
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // If the sequence index is out of bounds, return
    if (idx >= n)
        return;

    // Number of integers needed to represent the sequence and the number of zeros needed to fill the last integer
    int num_ints = (l + BitsInInt - 1) / BitsInInt;
    int padding_zeros = (num_ints * BitsInInt) - l;
    
    // Get the sequence
    int *seq = d_sequences + idx * num_ints;

    // Iterate over the sequence chunks
    for (int chunk_idx = 0; chunk_idx < num_ints; ++chunk_idx) {
        // Number of valid bits in the chunk regarding the padding zeros
        int valid_bits = BitsInInt;

        // If it is the last chunk, calculate the number of valid bits
        if (chunk_idx == num_ints - 1) {
            valid_bits = BitsInInt - padding_zeros;
        }

        // Get the chunk
        int chunk = seq[chunk_idx];

        // Iterate over the bits in the chunk flipping them one by one and checking if the sequence exists in the trie
        for (int bit_pos = 0; bit_pos < valid_bits; ++bit_pos) {
            // Flip the bit at the current position
            unsigned int flipped_chunk = chunk ^ (1 << (valid_bits - 1 - bit_pos));

            // Flag to check if all chunks exist
            bool all_chunks_exist = true;
            
            int current = 0;

            // Iterate over the chunks and check if they exist in the trie
            for (int check_chunk_idx = 0; check_chunk_idx < num_ints; ++check_chunk_idx) {
                unsigned int check_chunk = 0;
                int check_bits = BitsInInt;

                if(check_chunk_idx == num_ints - 1) {
                    check_bits = BitsInInt - padding_zeros;
                }

                // Get the chunk
                if (check_chunk_idx == chunk_idx) {
                    check_chunk = flipped_chunk;
                } else {
                    check_chunk = seq[check_chunk_idx];
                }

                // Traverse the trie
                for (int k = 0; k < check_bits; ++k) {
                    // Get the bit at the current position and move to the child node
                    int bit = (check_chunk >> (check_bits - 1 - k)) & 1;
                    current = (bit == 0) ? d_bit0[current] : d_bit1[current];
                    
                    // If the child does not exist, break
                    if (current == -1) {
                        all_chunks_exist = false;
                        break;
                    }
                }

                if (!all_chunks_exist)
                    break;
            }

            // If the sequence exists, add it to the results
            if (all_chunks_exist) {
                int pair_idx = atomicAdd(d_pairs_found, 1);
                d_sequence_idx[pair_idx] = idx;
                d_flipped_bit_position[pair_idx] = chunk_idx * BitsInInt + bit_pos;
            }
        }
    }
}


int main() {
    // Host structure of arrays
    std::vector<int> h_bit;
    std::vector<int> h_bit0;
    std::vector<int> h_bit1;
    std::vector<int> h_sequences;
    int sequence_length;
    int num_sequences;

    auto start_time_load = std::chrono::high_resolution_clock::now();
    load_trie("input_sequences.txt", h_bit, h_bit0, h_bit1, h_sequences, sequence_length, num_sequences);
    auto end_time_load = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> load_duration = end_time_load - start_time_load;
    printf("GPU: Time to load trie and sequences: %f seconds\n", load_duration.count());
    size_t trie_size = h_bit.size();

    // Device structure of arrays
    auto start_time_alloc = std::chrono::high_resolution_clock::now();
    Node d_trie;
    CUDA_CHECK(cudaMalloc(&d_trie.bit, trie_size * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_trie.bit0, trie_size * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_trie.bit1, trie_size * sizeof(int)));

    CUDA_CHECK(cudaMemcpy(d_trie.bit, h_bit.data(), trie_size * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_trie.bit0, h_bit0.data(), trie_size * sizeof(int), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(d_trie.bit1, h_bit1.data(), trie_size * sizeof(int), cudaMemcpyHostToDevice));

    // Sequence memory
    int *d_sequences;
    CUDA_CHECK(cudaMalloc(&d_sequences, h_sequences.size() * sizeof(int)));
    CUDA_CHECK(cudaMemcpy(d_sequences, h_sequences.data(), h_sequences.size() * sizeof(int), cudaMemcpyHostToDevice));

    // Device memory for storing the results
    int *d_pairs_found;
    int *d_sequence_idx;
    int *d_flipped_bit_position;
    CUDA_CHECK(cudaMalloc(&d_pairs_found, sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_sequence_idx, sequence_length * num_sequences * sizeof(int)));
    CUDA_CHECK(cudaMalloc(&d_flipped_bit_position, sequence_length * num_sequences * sizeof(int)));
    CUDA_CHECK(cudaMemset(d_pairs_found, 0, sizeof(int)));

    auto end_time_alloc = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> alloc_duration = end_time_alloc - start_time_alloc;
    printf("GPU: Time to allocate memory and copy to device: %f seconds\n", alloc_duration.count());
    
    // Launch kernel
    int block_size = 1024;
    int grid_size = (num_sequences + block_size - 1) / block_size;

    auto start_time = std::chrono::high_resolution_clock::now();
    find_pairs_kernel<<<grid_size, block_size>>>(d_sequences, d_trie.bit, d_trie.bit0, d_trie.bit1, sequence_length, num_sequences, d_sequence_idx, d_flipped_bit_position, d_pairs_found);
    CUDA_CHECK(cudaDeviceSynchronize());
    auto end_time = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> find_duration = end_time - start_time;
    printf("GPU: Time to find pairs: %f seconds\n", find_duration.count());

    // Copy results back to host
    auto start_time_copy = std::chrono::high_resolution_clock::now();
    int pairs_found;
    CUDA_CHECK(cudaMemcpy(&pairs_found, d_pairs_found, sizeof(int), cudaMemcpyDeviceToHost));
    std::vector<int> h_sequence_idx(pairs_found);
    std::vector<int> h_flipped_bit_position(pairs_found);

    CUDA_CHECK(cudaMemcpy(h_sequence_idx.data(), d_sequence_idx, pairs_found * sizeof(int), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(h_flipped_bit_position.data(), d_flipped_bit_position, pairs_found * sizeof(int), cudaMemcpyDeviceToHost));

    auto end_time_copy = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> copy_duration = end_time_copy - start_time_copy;
    printf("GPU: Time to copy results to host: %f seconds\n", copy_duration.count());
    
    // Print results
    std::ofstream outfile("gpu_output_pairs.txt");
    for (int i = 0; i < pairs_found; i++) {
       outfile << "Sequence " << h_sequence_idx[i] << " flipped at bit " << h_flipped_bit_position[i] << std::endl;
    }
    outfile.close();



    // Cleanup
    CUDA_CHECK(cudaFree(d_trie.bit));
    CUDA_CHECK(cudaFree(d_trie.bit0));
    CUDA_CHECK(cudaFree(d_trie.bit1));
    CUDA_CHECK(cudaFree(d_sequences));
    CUDA_CHECK(cudaFree(d_pairs_found));
    CUDA_CHECK(cudaFree(d_sequence_idx));
    CUDA_CHECK(cudaFree(d_flipped_bit_position));

    return 0;
}