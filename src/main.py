import BioInfExercises


if __name__ == '__main__':
    file_one = r"C:\Users\fabri\PycharmProjects\BioInf\tests\test_input\sequences.fasta"
    file_two = r"C:\Users\fabri\PycharmProjects\BioInf\tests\test_input\genome.fasta"

    tup = BioInfExercises.parse_fasta(file_two)
    BioInfExercises.discard_ambiguous_seqs(tup[1])
    BioInfExercises.nucleotide_frequencies(tup[1])
    dictionary = BioInfExercises.map_reads(file_one, file_two)

    path = r"C:\Users\fabri\PycharmProjects\BioInf\tests\test_input\Aligned.out.sam"
    BioInfExercises.convert_sam_fasta(path)
