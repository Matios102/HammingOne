#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <chrono>

#define BitsInInt 32

using namespace std;

// Structure representing a node in the trie
struct Node {
    int bit;  // Bit value of the node
    int bit0; // Index of the child node with bit value 0
    int bit1; // Index of the child node with bit value 1
    Node(int b = -1, int b0 = -1, int b1 = -1) : bit(b), bit0(b0), bit1(b1) {}
};


// Function to load the trie and the sequence vector from a file
void load_trie(const char *filename, std::vector<Node> &trie, std::vector<int> &sequences, int &sequence_length, int &num_sequences) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    // Add the root node to the trie
    trie.emplace_back();

    // Read sequence length and number of sequences
    std::string line;
    if (!std::getline(file, line)) {
        std::cerr << "Error: Could not read trie metadata" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::istringstream iss(line);
    if (!(iss >> sequence_length >> num_sequences)) {
        std::cerr << "Error: Could not parse sequence metadata" << std::endl;
        exit(EXIT_FAILURE);
    }

    // Seqences will be represented as 32 bit integers
    // num_ints is the number of integers needed to represent the sequence
    int num_ints = (sequence_length + BitsInInt - 1) / BitsInInt;

    // padding_zeros is the number of zeros needed to fill the last integer
    int padding_zeros = (num_ints * BitsInInt) - sequence_length;

    // Read the sequences
    for (int i = 0; i < num_sequences; ++i) {
        std::string sequence;
        if (!(file >> sequence)) {
            std::cerr << "Error: Could not read sequence " << i << std::endl;
            break;
        }
        if (sequence.length() != static_cast<size_t>(sequence_length)) {
            std::cerr << "Error: Invalid sequence length for sequence " << i << std::endl;
            continue;
        }

        // Convert the sequence to integers
        for (int j = 0; j < num_ints; ++j) {
            int chunk = 0;
            int valid_bits = (j == num_ints - 1) ? BitsInInt - padding_zeros : BitsInInt;
            for (int k = 0; k < valid_bits; ++k) {
                chunk |= (sequence[j * BitsInInt + k] == '1') << (valid_bits - 1 - k);
            }
            // Add the chunk to the sequence
            sequences.push_back(chunk);
        }

        // Build the trie
        int current = 0;
        for (int j = 0; j < sequence_length; ++j) {

            // Get the bit at position j
            char bit_char = sequence[j];
            int current_bit = (bit_char == '1') ? 1 : 0;

            // Get the index of the child node with the current bit
            int child_idx = (current_bit == 0) ? trie[current].bit0 : trie[current].bit1;

            // If the child node does not exist, create it
            if (child_idx == -1) {
            
                child_idx = trie.size();

                trie.emplace_back(current_bit, -1, -1);
                
                if (current_bit == 0) {
                    trie[current].bit0 = child_idx;
                } else {
                    trie[current].bit1 = child_idx;
                }

            }
            current = child_idx;
        }
    }

    file.close();
}



// Function to find pairs by flipping bits
void find_pairs(const std::vector<int> &sequences, const std::vector<Node> &trie, int sequence_length, int num_sequences, std::vector<std::pair<int, int>> &pairs) {
    // Number of integers needed to represent the sequence and the number of zeros needed to fill the last integer
    int num_ints = (sequence_length + BitsInInt - 1) / BitsInInt;
    int padding_zeros = (num_ints * BitsInInt) - sequence_length;

    // Iterate over the sequences
    for (int idx = 0; idx < num_sequences; ++idx) {
        // Get the sequence
        const int *seq = &sequences[idx * num_ints];

        // Iterate over the chunks of the sequence
        for (int chunk_idx = 0; chunk_idx < num_ints; ++chunk_idx) {
            // Number of valid bits in the chunk regarding the padding zeros
            int valid_bits = (chunk_idx == num_ints - 1) ? BitsInInt - padding_zeros : BitsInInt;
            
            // Get the chunk
            int chunk = seq[chunk_idx];

            // Iterate over the bits of the chunk flip them and check if the new sequence exists in the trie
            for (int bit_pos = 0; bit_pos < valid_bits; ++bit_pos) {
                int flipped_chunk = chunk ^ (1 << (valid_bits - 1 - bit_pos));

                bool exists = true;

                int current = 0;

                // Iterate over the chunks and check if they exist in the trie
                for (int check_chunk_idx = 0; check_chunk_idx < num_ints; ++check_chunk_idx) {

                    // Get the chunk to check
                    int check_chunk = (check_chunk_idx == chunk_idx) ? flipped_chunk : seq[check_chunk_idx];
                    int check_valid_bits = (check_chunk_idx == num_ints - 1) ? BitsInInt - padding_zeros : BitsInInt;

                    // Traverse the trie
                    for (int k = 0; k < check_valid_bits && exists; ++k) {
                        int bit = (check_chunk >> (check_valid_bits - 1 - k)) & 1;

                        current = (bit == 0) ? trie[current].bit0 : trie[current].bit1;

                        // If the node does not exist, break
                        if (current == -1) {
                            exists = false;
                            break;
                        }
                    }
                }
                
                // If the sequence exists, add the pair to the result
                if (exists) {
                    pairs.emplace_back(idx, chunk_idx * BitsInInt + bit_pos);
                }
            }
        }
    }
}

int main() {
    std::vector<Node> trie;
    std::vector<int> sequences;
    int sequence_length;
    int num_sequences;

    auto start_time_load = std::chrono::high_resolution_clock::now();
    load_trie("input_sequences.txt", trie, sequences, sequence_length, num_sequences);
    auto end_time_load = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed_time_load = end_time_load - start_time_load;
    std::cout << "CPU: Time to load trie and sequences: " << elapsed_time_load.count() << " seconds" << std::endl;

    std::vector<std::pair<int, int>> pairs;

    auto start_time = std::chrono::high_resolution_clock::now();
    find_pairs(sequences, trie, sequence_length, num_sequences, pairs);
    auto end_time = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed_time = end_time - start_time;
    std::cout << "CPU: Time to find pairs: " << elapsed_time.count() << " seconds" << std::endl;

    // Write results to file
    std::ofstream outfile("cpu_output_pairs.txt");
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open file output_pairs.txt" << std::endl;
        exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < pairs.size(); ++i) {
        int seq_idx = pairs[i].first;
        int bit_pos = pairs[i].second;
        outfile << "Sequence " << seq_idx << " flipped at bit " << bit_pos << std::endl;
    }

    outfile.close();
    return 0;
}
