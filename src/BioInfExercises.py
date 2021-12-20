
import re


def parse_fasta(
        path: str
) -> tuple:
    """ reads a fasta file and turns it into a list.

    Args:
        path():   the path to the fasta file.

    Returns:
        tuple:  a tuple of lists one containing the header of the fasta file,
                the other containing the sequence
    """

    header_list: list[str] = []
    seq_list: list[str] = []

    with open(path, 'r') as file:
        for line in file:
            line = line.rstrip()
            if line.startswith('>'):
                header_list.append(line)
                seq_list.append('')
            else:
                seq_list[-1] = seq_list[-1] + line

    tup = (header_list, seq_list)
    return tup


def discard_ambiguous_seqs(
        seq_list: list
) -> list:
    """ reads a string of nucleotide sequence and
        removes everything that isn't ACTG.

    Args:
        seq_list: containing nucloetide sequences.

    Returns:
        list: edited input list
    """

    for i in range(len(seq_list)):
        seq_list[i] = re.sub("[^ACGTacgt]+", "", seq_list[i])
    return seq_list


def nucleotide_frequencies(
        seq_list: list[str]
) -> None:
    """ caluculates total frequency of each nucleotide
        across all input sequences.

    Args:
        seq_list: containing nucloetide sequences.

    Returns:
        None
    """
    for i in range(len(seq_list)):
        count_a = (seq_list[i].count('a') + seq_list[i].count('A')) / len(seq_list[i])
        count_t = (seq_list[i].count('t') + seq_list[i].count('T')) / len(seq_list[i])
        count_g = (seq_list[i].count('g') + seq_list[i].count('G')) / len(seq_list[i])
        count_c = (seq_list[i].count('c') + seq_list[i].count('C')) / len(seq_list[i])

        print(f"sequence number: {i + 1}")
        print(f"A: {round(count_a, 2)}")
        print(f"T: {round(count_t, 2)}")
        print(f"G: {round(count_g, 2)}")
        print(f"C: {round(count_c, 2)}")


def map_reads(
        file_one: str,
        file_two: str
) -> dict:
    """ reads two fasta files and turns them in
        to a dictionary of dictionarys.

    Takes two FASTA files as input, the first containing short
    read sequences, and the second containing reference sequences.
    Print the nucleotide fractions for both files to the console
    and returns a dictionary of dictionaries, where the outer
    dictionary uses the names of query sequences as its keys,
    and the inner dictionary uses reference sequence names as
    keys and a list of 1-based indices indicating at which position
    in the reference sequence the sequence occurs as an exact
    substring.


    Args:
        file_one:   FASTA file containing short read sequences
        file_two:   FASTA file containing reference sequences

    Returns:
        dict:       A dictionary containing a dictionaries
                    containing seq_list
    """

    tup1 = parse_fasta(file_one)
    key1 = tup1[0]

    tup2 = parse_fasta(file_two)
    key2 = tup2[0]

    list1 = discard_ambiguous_seqs(tup1[1])
    nucleotide_frequencies(list1)

    list2 = discard_ambiguous_seqs(tup2[1])
    nucleotide_frequencies(list2)

    dicti = {}

    seq_list: list[int] = {}
    for i, word in enumerate(tup2[1]):
        seq_list[i] = word.upper().find(list1[i].upper())
        dicti[key1[i]] = {key2[i]: seq_list[i]}

    return dicti


def convert_sam_fasta(
        sam_file: str
) -> None:
    """ converts a SAM file into a fasta file
        with the same name as the SAM file.

    Args:
        sam_file: path to SAM file

    Returns:
        None
    """
    file_name = sam_file.replace('.sam', '')
    open(file_name + '.fasta', 'w').close()
    wc: int = 0
    with open(sam_file, 'r') as file:
        for i, line in enumerate(file):
            if i <= 3:
                continue
            for j, word in enumerate(line.split()):
                if j == 9:
                    wc += 1
                    with open(file_name + '.fasta', 'a') as new_file:

                        if i == 4:
                            new_file.write('>alignment ' + str(wc) + '\n')
                        else:
                            new_file.write('\n>alignment ' + str(wc) + '\n')
                        k = 0
                        for letter in word:
                            if k < 70:
                                new_file.write(letter)
                                k += 1
                            else:
                                k = 1
                                new_file.write('\n' + letter)
