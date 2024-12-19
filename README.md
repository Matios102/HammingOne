# README

## Problem Description
Hamming distance is a metric used to measure the difference between two strings of equal length. It counts the number of positions at which the corresponding symbols (usually bits) are different. In this problem, the goal is to generate binary sequences such that some pairs of sequences have a Hamming distance of exactly one, i.e., they differ in exactly one bit position.

A sequence pair with Hamming distance one can be described as follows:

- The sequences are of equal length.
- There is exactly one bit difference between them.

## Approach
The program finds Hamming distance one pairs of binary sequences by utilizing a trie data structure and CUDA for parallel computation. It first loads a set of binary sequences and builds a trie where each node corresponds to a bit in the sequences. Then, using a CUDA kernel, it flips each bit of the sequences one by one and checks if the resulting sequence exists in the trie, indicating a Hamming distance of one. If a pair is found, it records the sequence index and the bit position that was flipped. The results are then written to a file.

## Algorithm Complexity
The time complexity of the algorithm is O(n*l^2) where n is the number of sequences and l is the length of each sequence. The algorithm first builds a trie with O(n*l) complexity and then for a given sequence, it check all possible bit flips and searches tht trie to check if such a sequence exists. The trie search has a complexity of O(l) and there are l possible bit flips, so the overall complexity is O(l^2).

## Sequence Generation Function
The input generator uses the `rand()` function to generate random binary sequences. The function generates a random number between 0 and 1 and rounds it to the nearest integer to obtain a binary digit. To ensure that there exist pairs with hamming distance one the function flips a single bit in a specified percentage of sequences.

### Command Line Execution
To run the program, use the following syntax:
```bash
./input_generator <sequence_length> <number_of_sequences> [percentage]
```
- `sequence_length`: Length of each binary sequence (positive integer).
- `number_of_sequences`: Number of sequences to generate (positive integer).
- `percentage`: (Optional) Percentage of sequences to modify by flipping one bit (0-100).

Example:
```bash
./input_generator 1000 100000 10
```
This generates 100,000 sequences, each of length 1,000. 90% of the sequences are generated randomly, and the remaining 10% are created by flipping one random bit in a random sequence from the first 90%. This ensures that at least 10% of the sequences have a Hamming distance of one.

### Makefile
The provided `Makefile` simplifies compilation and execution of the program.

#### Targets
- **`all`**: Builds and runs the input generator, CPU, and GPU implementations.
- **`input`**: Compiles the input generator.
- **`run_input`**: Executes the input generator with predefined parameters.
- **`run_cpu`**: Executes the CPU-based solution.
- **`run_gpu`**: Executes the GPU-based solution.
- **`clean`**: Removes generated files.
