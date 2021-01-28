import pytest
import numpy as np
from align import algs
import os



@pytest.fixture
def some_relevant_data():
	return np.ones(10)

def test_fasta_io():

	assert True

def test_scoring_matrix_io():

	assert True

def test_identical():
	# testblowsum=algs.SmithWaterman()
	# testblowsum.LoadScoreMatrix("../scoring_matrices/BLOSUM50.mat")
	#
	# assert (testblowsum.score_matrix['A']['A'] == 4)
	# assert (testblowsum.score_matrix['F']['A'] == -2)
	# assert (testblowsum.score_matrix['A']['F'] == -2)


def test_alignment_score():
	assert True
