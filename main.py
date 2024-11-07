import argparse
from Bio import SeqIO

# Set up argument parser
parser = argparse.ArgumentParser(description="Parse a FASTA file with two sequences.")
parser.add_argument("fasta_path", type=str, help="Path to the FASTA file")
args = parser.parse_args()

sequences = list(SeqIO.parse(args.fasta_path, "fasta"))

match_score=1
mismatch_score=-1
gap_penalty_score=-2
score_matrix = []
direction_matrix = []
path = []
populated_matrices = False

seq2 = sequences[0];
seq1 = sequences[1];


def get_score_on_pos(pos1, pos2):
    if seq1[pos1] == seq2[pos2]:
        return match_score
    else:
        return mismatch_score

def populate_matrices(seq1, seq2):

    global score_matrix, direction_matrix, populated_matrices
    matrix_width = len(seq1) + 1
    matrix_height = len(seq2) + 1

    score_matrix = [[0 for i in range(matrix_width)] for j in range(matrix_height)]
    direction_matrix = [[[] for i in range(matrix_width)] for j in range(matrix_height)]
    for i in range(matrix_width):
        score_matrix[0][i] = i * gap_penalty_score
        direction_matrix[0][i] = [0,0]

    for i in range(matrix_height):
        score_matrix[i][0] = i * gap_penalty_score
        direction_matrix[i][0] = [0,0]

    choice_translation_dict = {0 : [-1,-1], 1 : [-1,0], 2 : [0,-1]}
    for i in range(1,matrix_height):
        for j in range(1,matrix_width):
            choiceList = [score_matrix[i-1][j-1]+get_score_on_pos(j-1,i-1),
                            score_matrix[i-1][j]+gap_penalty_score,
                            score_matrix[i][j-1]+gap_penalty_score
                            ]
            score_matrix[i][j] = max(choiceList)
            direction_matrix[i][j] = choice_translation_dict[choiceList.index(max(choiceList))]
    populated_matrices = True

def find_path():
    global path, direction_matrix, populated_matrices

    if not populated_matrices :
        raise ValueError("Matrices have not been populated")
    y = len(direction_matrix) - 1
    x = len(direction_matrix[0]) - 1
    current_directions = direction_matrix[y][x]
    tmp_path = []
    while current_directions != [0,0] :
        tmp_path.insert(0,current_directions)
        y += current_directions[0]
        x += current_directions[1]
        current_directions = direction_matrix[y][x]
    path = tmp_path

def generate_alignment_output():
    global seq1, seq2, path
    pos = [0, 0]
    seq1_str = ""
    seq2_str = ""

    # Traverse through the path and build aligned sequences
    for step in path:
        pos = [p - s for p, s in zip(pos, step)]
        match step:
            case [-1, -1]:
                seq1_str += seq1[pos[1] - 1]
                seq2_str += seq2[pos[0] - 1]
            case [-1, 0]:
                seq1_str += "_"
                seq2_str += seq2[pos[0] - 1]
            case [0, -1]:
                seq1_str += seq1[pos[1] - 1]
                seq2_str += "_"

    # Build the alignment output
    alignment_output = (
        " ".join(list(seq1_str)) + "\n" +
        " ".join([ " " if i == "_" or j == "_" else ("|" if i == j else ".") for i, j in zip(seq1_str, seq2_str)]) + "\n" +
        " ".join(list(seq2_str)) + "\n"
    )
    return alignment_output


def generate_score_matrix_output(data):
    col_widths = [max(len(str(item)) for item in column) for column in zip(*data)]

    output_lines = []
    for row in data:
        formatted_row = " | ".join(f"{str(item):<{width}}" for item, width in zip(row, col_widths))
        output_lines.append(formatted_row)

    return "\n".join(output_lines)

def save_output_to_file():
    populate_matrices(seq1, seq2)
    find_path()
    score_matrix_output = generate_score_matrix_output(score_matrix)
    alignment_output = generate_alignment_output()
    with open("output.txt", "w") as f:
        f.write("Score Matrix:\n")
        f.write(score_matrix_output)
        f.write("\nAlignment:\n")
        f.write(alignment_output)

save_output_to_file()
