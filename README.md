# Sequence Alignment Program

This Python script performs sequence alignment using Needleman-Wunsch algorithm. It reads sequences from a FASTA file, computes the alignment score matrix, determines the optimal alignment path, and prints both the alignment score matrix and the final aligned sequences.

## Prerequisites

Make sure you have Python installed, along with the Biopython library for handling FASTA files. You can install Biopython with:

```bash
pip install biopython
```
## Usage

The script takes a path to a FASTA file containing two sequences. To run the script, use the following command:

python script_name.py <fasta_path>
```
<fasta_path>: Path to the FASTA file with two sequences.
```
## Example

Input FASTA File

Create a FASTA file, for example, example.fasta, with the following content:

```fasta
>Sequence1
GTCGACGCA
>Sequence2
GATTACA
```

#### Running the Program

Execute the script with the path to the example.fasta file:

python script_name.py example.fasta

#### Expected Output

The script will generate a `output.txt` file alongside the script with the score matrix, alignment path, and aligned sequences. Here is the expected output:
```
Score Matrix:
0   | -2  | -4 | -6 | -8 | -10 | -12 | -14 | -16 | -18
-2  | 1   | -1 | -3 | -5 | -7  | -9  | -11 | -13 | -15
-4  | -1  | 0  | -2 | -4 | -4  | -6  | -8  | -10 | -12
-6  | -3  | 0  | -1 | -3 | -5  | -5  | -7  | -9  | -11
-8  | -5  | -2 | -1 | -2 | -4  | -6  | -6  | -8  | -10
-10 | -7  | -4 | -3 | -2 | -1  | -3  | -5  | -7  | -7 
-12 | -9  | -6 | -3 | -4 | -3  | 0   | -2  | -4  | -6 
-14 | -11 | -8 | -5 | -4 | -3  | -2  | -1  | -3  | -3 
Alignment:
G T C G A C G C A
| . . . |     | |
G A T T A _ _ C A

```


Explanation of Output

1.	Score Matrix: The alignment score matrix shows the cumulative scores for aligning subsequences of the two sequences with gap penalties, match scores, and mismatch scores.
2.	Alignment Path: Displays the optimal alignment between the two sequences using:
    - _ to indicate gaps.
    - Matching nucleotides are shown with |
    - Mismatched nucleotides are shown with .
