from _pytest.capture import capsys
from src.BioInfExercises import parse_fasta, discard_ambiguous_seqs, nucleotide_frequencies, map_reads, \
    convert_sam_fasta


def test_parse_fasta():
    """ test for funtion parse_fasta()
    tests if it returns right headers
    """
    arr = parse_fasta(r"C:\Users\fabri\PycharmProjects\BioInf\tests\test_input\genome.fasta")
    assert arr[0][0] == '>chr1'
    assert arr[0][1] == '>chr2'
    assert arr[0][2] == '>chr3'
    assert arr[0][3] == '>chr4'


def test_discard_ambiguous_seqs():
    """ test for funtion discard_ambiguous_seqs()
    tests if function recognizes wrong letters and is able to check capitalized and lowercase letters
    """
    wrong_letters = discard_ambiguous_seqs(['TAT', 'SSSSSSSNAKE', 'GACCTATGTGTCCGGTAAC', 'GAGA'])
    assert wrong_letters == ['TAT', 'A', 'GACCTATGTGTCCGGTAAC', 'GAGA']

    capitalized_lowercase = discard_ambiguous_seqs(['ccCcT', 'GCtAtTgTTAG', 'ATAgccca', 'GooooOoAt'])
    assert capitalized_lowercase == ['ccCcT', 'GCtAtTgTTAG', 'ATAgccca', 'GAt']


def test_nucleotide_frequencies(capsys):
    """ test for funtion nucleotide_frequencies()
    tests if function prints corr calculated frequencies out
    """
    nucleotide_frequencies(['TATA'])
    captured = capsys.readouterr()
    assert captured.out == 'sequence number: 1\nA: 0.5\nT: 0.5\nG: 0.0\nC: 0.0\n'


def test_map_reads():
    """ test for funtion map_reads()
    tests if function returns correct index of the position of the substring
    """

    file_one = r"C:\Users\fabri\PycharmProjects\BioInf\tests\test_input\sequences.fasta"
    file_two = r"C:\Users\fabri\PycharmProjects\BioInf\tests\test_input\genome.fasta"
    dic = map_reads(file_one,file_two)
    assert dic['>sequence1']['>chr1'] == 758
    assert dic['>sequence2']['>chr2'] == 1421
    assert dic['>sequence4']['>chr4'] == 1454


def test_convert_sam_fasta():
    """ test for funtion sam_to_fasta_file()
    tests samples if the new created fasta file has correct format
    """
    convert_sam_fasta(r"C:\Users\fabri\PycharmProjects\BioInf\tests\test_input\Aligned.out.sam")
    with open(r"C:\Users\fabri\PycharmProjects\BioInf\tests\test_input\Aligned.out.fasta") as fh:
        data = fh.readlines()
    assert data[0][0] == '>'
    assert data[12][1:10] == 'alignment'
    assert data[13][0:4] == 'TGTG'


if __name__ == '__main__':
    test_parse_fasta()
    test_discard_ambiguous_seqs()
    test_nucleotide_frequencies(capsys)
    test_map_reads()
    test_convert_sam_fasta()
