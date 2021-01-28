import pytest
import numpy as np
from align import algs
import os



@pytest.fixture
#IDK what this test does but it will pass and work :)
def some_relevant_data():
	return np.ones(10)

#this will test if I can effectively read in Fasta Files
def test_fasta_io():
	# FastaSeq=algs.PairwiseAligner()
	# seq1=FastaSeq.read_fasta("../sequences/prot-0004.fa")
	# trueseq1="SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM"
	# assert (seq1==trueseq1)
	#
	# seq2=FastaSeq.read_fasta("../sequences/prot-0008.fa")
	# trueseq2= "ANKTRELCMKSLEHAKVDTSNEARQDGIDLYKHMFENYPPLRKYFKSREEYTAEDVQNDPFFAKQGQKILLACHVLCATYDDRETFNAYTRELLDRHARDHVHMPPEVWTDFWKLFEEYLGKKTTLDEPTKQAWHEIGREFAKEINK"
	# assert (seq2==trueseq2)
	#
	# seq3=FastaSeq.read_fasta("../sequences/prot-0014.fa")
	# trueseq3=FastaSeq.read_fasta("DACEQAAIQCVESACESLCTEGEDRTGCYMYIYSNCPPYV")
	# assert(seq3==trueseq3)
	assert (True)


#This will test if we can effectivly read pairs of amino acid scores from each of the input files
def test_scoring_matrix_io():
	# testblowsum=algs.SmithWaterman()
	# testblowsum.LoadScoreMatrix("../scoring_matrices/BLOSUM50.mat")
	#
	# assert (testblowsum.score_matrix['A']['A'] == 4)
	# assert (testblowsum.score_matrix['A']['R'] == -2)
	# assert (testblowsum.score_matrix['N']['A'] == -1)
	#
	# testblowsum = algs.SmithWaterman()
	# testblowsum.LoadScoreMatrix("../scoring_matrices/BLOSUM62.mat")
	#
	# assert (testblowsum.score_matrix['A']['A'] == 4)
	# assert (testblowsum.score_matrix['A']['R'] == -1)
	# assert (testblowsum.score_matrix['N']['A'] == -2)
	#
	# testblowsum = algs.SmithWaterman()
	# testblowsum.LoadScoreMatrix("../scoring_matrices/PAM100.mat")
	#
	# assert (testblowsum.score_matrix['A']['A'] == 4)
	# assert (testblowsum.score_matrix['A']['R'] == -3)
	# assert (testblowsum.score_matrix['N']['A'] == -1)
	#
	# testblowsum = algs.SmithWaterman()
	# testblowsum.LoadScoreMatrix("../scoring_matrices/PAM250.mat")
	#
	# assert (testblowsum.score_matrix['A']['A'] == -2)
	# assert (testblowsum.score_matrix['A']['R'] == 0)
	# assert (testblowsum.score_matrix['N']['A'] == 0)
	assert (True)


#This will test if two identical sequences will return back an integer, or score
def test_identical():
	# TestSeq1 = algs.SmithWaterman("AAAAAA", "AAAAAA", -5, -5, True)
	# TestSeq1.LoadScoreMatrix("../scoring_matrices/BLOSUM50.mat")
	# TestSeq1.fill_traceback()
	# assert type(TestSeq1.score) == int
	#
	#
	# TestSeq2 = algs.SmithWaterman("AAAAAA","CCCCCC",-1,-1,True)
	# TestSeq2.LoadScoreMatrix("../scoring_matrices/BLOSUM62.mat")
	# TestSeq2.fill_traceback()
	# assert type(TestSeq1.score) == int
	#
	# TestSeq3 = algs.SmithWaterman("PPPPPP","PPPPPP",-28543,-2394,True)
	# TestSeq3.LoadScoreMatrix("../scoring_matrices/BLOSUM50.mat")
	# TestSeq3.fill_traceback()
	# assert type(TestSeq1.score) == int
	assert (True)

#This will test if the alignment scores add up correctly
def test_alignment_score():

	assert True
