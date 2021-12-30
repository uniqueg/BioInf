from pathlib import Path

from src import BioInfExercises


if __name__ == '__main__':
    file_one = str(Path(__file__).parents[1] / "tests/test_input/sequences.fasta")
    file_two = str(Path(__file__).parents[1] / "tests/test_input/genome.fasta")

    ids, seqs = BioInfExercises.parse_fasta(file_one)
    print(ids)
    print(seqs)
    valid_seqs = BioInfExercises.discard_ambiguous_seqs(seqs)
    print(valid_seqs)
    BioInfExercises.nucleotide_frequencies(valid_seqs)
    dictionary = BioInfExercises.map_reads(file_one, file_two)
    print(dictionary)

    # add commands to create STAR index, run STAR and analyze the alignments (2.1.1, 2.1.2, 2.2.1)

    path = str(Path(__file__).parents[1] / "tests/test_input/Aligned.out.sam")
    BioInfExercises.convert_sam_fasta(path)

