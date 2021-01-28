import numpy as np


'''
Class PairwiseAligner
This is a parent class of Smith Waterman and Needleman Wunsh, this class is used to read and clean both alignments 
sequences. You are able to call this class also but it will have no functionality if you want to run or fill a traceback 
for best alignment. You can use this class to obtain a score matrix or to read in a fastafile or clean a sequence
Parameters:
	seq1= first raw sequence that will be aligned
	seq2 = second raw sequence that will be aligned
	score_matrix = Will contain score matrix that will return amino acid score by key ['A']['A']
	pointer_mat = contains the pointer matrix for Mp,Xp,Yp. will need to unpack if desired to use
	Ix = contains scores for best alignement for seq1
	Iy = contains scores for best alignment for seq2
	M = contains best scores for overall alignments
	gap_start = gap open penatly for sequence you are trying to align
	gap_end = gap extention penalty for sequence you are trying to align
	
	
Functions:

LoadScoreMatrix(self,filename):
filename= path to score matrix file

this function can be used to read in a score matrix, this will initalize the score matrix depending on the matrix you 
input. MUST INCLUDE FILEPATH TO MATRIX AND NOT JUST THE FILENAME

def clean_seq(self,isprotein):
this will clean the sequence of the matrix, only including values that are amino acid sequences/ nucleotide sequences. 


def read_fasta(filename):
this will read in the fasta file for any function, this is not implemented in my final code but still can be used if 
you include the filepath to read it, will also return the fasta and not assign it to any function. 


'''
class PairwiseAligner:
	def __init__(self, seq1, seq2, gap_start, gap_ext,isprotein):
		#Assign fasta file path to the
		self.seq1 = seq1
		self.seq2 = seq2
		#self.read_fasta()

		self.score_matrix = 0
		self.pointer_mat = np.zeros((3,len(self.seq1) + 1, len(self.seq2) + 1), dtype="<U10")

		self.Ix = np.zeros([len(self.seq1) + 1, len(self.seq2) + 1])
		self.Iy = np.zeros([len(self.seq1) + 1, len(self.seq2) + 1])
		self.M = np.zeros([len(self.seq1) + 1, len(self.seq2) + 1])
		self.gap_start = gap_start
		self.gap_ext = gap_ext
		#clean sequence once everything is loaded.
		self.clean_seq(isprotein)


	def LoadScoreMatrix(self, filename):
		AminoAcidScore = {}

		with open(filename, 'r') as l:
			lines = [line for line in l.read().splitlines() if line[0] != "#"]

		AminoAcid = lines[0].split()
		Scores = lines[1:]

		for i, a1 in enumerate(AminoAcid):
			for j, a2 in enumerate(AminoAcid):
				if a1 not in AminoAcidScore:
					AminoAcidScore[a1] = {}
				AminoAcidScore[a1][a2] = float(Scores[i].split()[j])

		self.score_matrix = AminoAcidScore

	#Will clean function for any unwanted sequence characters, must input True/False for if the sequence passed for a
	# protein or a nucleotide sequence. Will assign the clean sequences back to self.seq1/self.seq2
	def clean_seq(self,isprotein):
		protein=['A', 'R', 'N', 'D', 'C',  'Q',  'E',  'G',  'H',  'I',  'L',  'K',
				 'M',  'F',  'P',  'S',  'T',  'W',  'Y',  'V',  'B',  'Z',  'X']
		nucleotide=['A','C','G','T']
		seq1=""
		seq2=""
		# Assign * for any character that is not allowed
		if(isprotein==True):
			sequence_parameters=protein
		else:
			sequence_parameters = nucleotide
		for char in self.seq1.upper():
			if char not in sequence_parameters:
				seq1 += "*"
			else:
				seq1 += char
		#Assign * for any character that is not allowed
		for char in self.seq2.upper():
			if char not in sequence_parameters:
				seq2 += "*"
			else:
				seq2 += char
		self.seq1=seq1
		self.seq2=seq2

	#will return a fasta sequence, MUST INPUT FILEPATH
	def read_fasta(filename):
		contents = []
		with open(filename, "r") as f:
			for line in f:
				line = line.strip('\n')
				contents.append(line)
		seq = ''.join(contents[1:])
		return (seq)

	#Parent class for this function, pass and ignore function
	def run_traceback(self):
		pass
	#Parent class for this functin, pass and ignore function
	def fill_traceback(self):
		pass
'''
class SmithWaterman:
This is a child class of PairwiseAligner that will preform local sequence alignment of two sequences. all input requirements
are the same with the only difference in this class if the two functions inputted


fill_traceback():
	This function will fill in the traceback matrix for alignment based on a local smith-waterman approach. This will fill in 
	values for M,Ix,Iy and Pointer_matrix. This will input scores based off the max of previous steps, or 0. This alg
	does not allow any negative scores for best alignment.
	
	
	returns:
	self.pointer_mat = 3 pointer matrix for each of the corresponding scores
	self.score , score of global alignment which is found in bottom left corner
	self.M,Ix,Iy - score matrix of each of the three dynamic programming iterations. 

run_traceback():
	this function will run the traceback sequence for the sequence alignment test.
	
	returns:
	self.seq1_aligned = aligned seq1
	self.seq2_aligned = align_seq2
	self.align_pattern = pattern of the best fit 
	recommended to print these three lines all at once!!


'''

class SmithWaterman(PairwiseAligner):
	def fill_traceback(self):
		#Input all pointer matrix to be zero, starting point for a local alignment
		Mp=np.zeros([len(self.seq1) + 1,len(self.seq2) + 1],dtype="<U10")
		Xp=np.zeros([len(self.seq1) + 1,len(self.seq2) + 1],dtype="<U10")
		Yp=np.zeros([len(self.seq1) + 1,len(self.seq2) + 1],dtype="<U10")

		#Traverse though each of the two sequences to input recursive scores.
		for i, a1 in enumerate(self.seq1, start=1):
			for j, a2 in enumerate(self.seq2, start=1):
				# Obtain score for curent two seq based of scoring matrix
				Score = self.score_matrix[a1][a2]
				# Create index relationship for direction matrix
				DirectionLegend = {0: "M", 1: "Ix", 2: "Iy"}
				# Find direction of current M step, if values are equal then opt to point in direction M
				if (self.M[i - 1, j - 1] == self.Ix[i - 1, j - 1] or self.M[i - 1, j - 1] == self.Iy[i - 1, j - 1]):
					idx = 0
				else:
					idx = np.argmax([self.M[i - 1, j - 1], self.Ix[i - 1, j - 1], self.Iy[i - 1, j - 1]])
				# Append direction to Mp
				Mp[i, j] = DirectionLegend[idx]
				# Fill score from direction Mp, will correspond to the direction the score is coming from, or zero if no
				# positive values
				self.M[i, j] = max(
					0,
					self.M[i - 1, j - 1]+Score ,
					self.Ix[i - 1, j - 1]+Score,
					self.Iy[i - 1, j - 1]+Score
				)
				# Check direction of seq1 best alignments, if equal then opt in direction M
				if (self.M[i, j - 1] + self.gap_start + self.gap_ext == self.Ix[i, j - 1] + self.gap_ext):
					idx = 0
				else:
					idx = np.argmax(
						[self.M[i, j - 1] + self.gap_start + self.gap_ext, self.Ix[i, j - 1] + self.gap_ext])
				# Append Direction for Xp
				Xp[i, j] = DirectionLegend[idx]
				# Fill score of seq1 will correspond to max score direction found
				self.Ix[i, j] = max(
					self.M[i, j - 1] + self.gap_start + self.gap_ext,
					self.Ix[i, j - 1] + self.gap_ext,

				)
				# Check direction of seq2 best alignments, if equal then opt in direction M
				if (self.M[i - 1, j] + self.gap_start + self.gap_ext == self.Iy[i - 1, j] + self.gap_ext):
					idx = 0
				else:
					idx = np.argmax(
						[self.M[i - 1, j] + self.gap_start + self.gap_ext,-1, self.Iy[i - 1, j] + self.gap_ext])
				# Append Direction for Yp
				Yp[i, j] = DirectionLegend[idx]
				# Fill score of seq2 will correspond to max score direction found
				self.Iy[i, j] = max(
					self.M[i - 1, j] + self.gap_start + self.gap_ext,
					self.Iy[i - 1, j] + self.gap_ext,
				)

		# will combine pointer matrix based on 0,1,2 direction matrix order and assign to self.pointer_mat
		pointer_matrix = [Mp, Xp, Yp]
		self.pointer_mat = pointer_matrix
		# Score is value found on bottom right corner of M score matrix, can now acsses score freely without run_traceback
		self.score=np.max([np.max(self.M),np.max(self.Ix),np.max(self.Iy)])


	def run_traceback(self):
		#Variables that will hold alignment traceback alignment
		seq1_align = ""
		seq2_align = ""
		align = ""

		# Mp, Xp, Yp = self.pointer_mat
		#This will unpack all of the pointer matrix into a useable interface to traverse back based on index variable
		mat_total = (self.M, self.pointer_mat[0]), (self.Ix, self.pointer_mat[1]), (self.Iy, self.pointer_mat[2])

		#Find matrix that hols the best alignment score
		max_idx=np.argmax([np.max(self.M),np.max(self.Ix),np.max(self.Iy)])

		# Assign the current matrix/ current direction to be where we found the max value to lie
		cur_mat,cur_dir=mat_total[max_idx]
		#Find index of max value that is currently inside curr_mat
		i,j = np.argwhere(cur_mat == np.max(cur_mat))[0]


		idx=0
		while cur_mat[i,j] > 0 and (i>0 and j>0):
			# If direction is based off M then we will go back one step diagonally, or one step back in both the x,y direction
			if cur_dir[i,j]=='M':
				seq1_align = self.seq1[i - 1] + seq1_align
				seq2_align = self.seq2[j - 1] + seq2_align

				if (seq1_align[idx] == seq2_align[idx]):
					align = "|" + align
				else:
					align = "*" + align

				i = i - 1
				j = j - 1
				#If we are going in direction of M, then shift back to M score/pointer matrix, index 0
				cur_mat, cur_dir = mat_total[0]
			# Check if direction is based of seq1 being the best alignment, if so then we will only go back in one
			# step in the Y direction and shift to X score/pointer matrix
			elif cur_dir[i,j]=='Ix':
				seq1_align = '-' + seq1_align
				seq2_align = self.seq2[j] + seq2_align
				align = " " + align
				j = j - 1
				# If we are going in direction of M, then shift back to M score/pointer matrix, index 1
				cur_mat, cur_dir = mat_total[1]
			# Check if direction is based of seq2 being the best alignment, if so then we will only go back in one
			# step in the X direction and shift to the Y score/pointer matrix
			elif cur_dir[i,j]=='Iy':
				seq1_align = self.seq1[i] + seq1_align
				seq2_align = '-' + seq2_align
				align = " " + align
				i = i - 1
				# If we are going in direction of Iy, then shift back to M score/pointer matrix, index 2
				cur_mat,cur_dir=mat_total[2]
			idx=idx+1

		# Reutrn back best alignment score for each of the of the two sequneces and attach 3 new variables to the object.
		self.seq1_aligned = seq1_align
		self.seq2_aligned = seq2_align
		self.align_pattern = align


'''
class NeedlemanWunsch:
This is a child class of PairwiseAligner that will preform global sequence alignment of two sequences. All input requirements
are the same with the only difference in this class if the two functions inputted


fill_traceback():
	This function will fill in the traceback matrix for alignment based on a global Needleman-Wunsch approach. This will fill in 
	values for M,Ix,Iy and Pointer_matrix. This will input scores based off the max of previous steps. This alg
	allows negative scores with -inf being the floor case for any of the edges.


	returns:
	self.pointer_mat = 3 pointer matrix for each of the corresponding scores
	self.score , score of global alignment which is found in bottom left corner
	self.M,Ix,Iy - score matrix of each of the three dynamic programming iterations. 

run_traceback():
	This function will run the traceback sequence for the sequence alignment test, based on a global alignment. This 
	will start at the bottom most right case and traverse back until we reach an edge or corner

	returns:
	self.seq1_aligned = aligned seq1
	self.seq2_aligned = align_seq2
	self.align_pattern = pattern of the best fit 
	recommended to print these three lines all at once!!


'''


class NeedlemanWunsch(PairwiseAligner):
	def fill_traceback(self):
		#Create variable that will keep trace of negative infinity instances
		Inf = (float('-inf'))
		# I will first Initialize the three matrix M,Iy,Ix and input the initial values
		Mp=np.zeros([len(self.seq1) + 1,len(self.seq2) + 1],dtype="<U10")
		Xp=np.zeros([len(self.seq1) + 1,len(self.seq2) + 1],dtype="<U10")
		Yp=np.zeros([len(self.seq1) + 1,len(self.seq2) + 1],dtype="<U10")
		Mp[0, 0] = 0
		#input all initial values found on the x axis
		for i in range(1, len(self.seq1) + 1):
			self.M[i][0] = Inf
			self.Ix[i][0] = Inf
			self.Iy[i][0] = self.gap_start + i * self.gap_ext

			Mp[i][0] = 'Iy'
			Yp[i][0] = 'Iy'
		#input all initial values for the y axis
		for i in range(1, len(self.seq2) + 1):
			self.M[0][i] = Inf
			self.Ix[0][i] = self.gap_start + i * self.gap_ext
			self.Iy[0][i] = Inf

			Mp[0][i] = 'Ix'
			Xp[0][i]= 'Ix'

		#iterate though each score
		for i, a1 in enumerate(self.seq1, start=1):
			for j, a2 in enumerate(self.seq2, start=1):
				#Obtain score for curent two seq based of scoring matrix
				Score=self.score_matrix[a1][a2]
				#Create index relationship for direction matrix
				DirectionLegend = {0: "M", 1: "Ix", 2: "Iy"}

				#Find direction of current M step, if values are equal then opt to point in direction M
				if (self.M[i - 1, j - 1]==self.Ix[i - 1, j - 1] or self.M[i - 1, j - 1]==self.Iy[i - 1, j - 1] ):
					idx=0
				else:
					idx = np.argmax([self.M[i - 1, j - 1], self.Ix[i - 1, j - 1], self.Iy[i - 1, j - 1]])
				#Append direction to Mp
				Mp[i, j] = DirectionLegend[idx]
				#Fill score from direction Mp, will correspond to the direction the score is coming from
				self.M[i, j] = Score + max(
					self.M[i - 1, j - 1],
					self.Ix[i - 1, j - 1],
					self.Iy[i - 1, j - 1]
				)
				#Check direction of seq1 best alignments, if equal then opt in direction M
				if (self.M[i, j - 1] + self.gap_start + self.gap_ext == self.Ix[i, j - 1] + self.gap_ext):
					idx=0
				else:
					idx = np.argmax([self.M[i, j - 1] + self.gap_start + self.gap_ext,self.Ix[i, j - 1] + self.gap_ext])
				# Append direction to Xp
				Xp[i, j] = DirectionLegend[idx]
				# Fill score of seq1 will correspond to max score direction found
				self.Ix[i, j] = max(
					self.M[i, j - 1] + self.gap_start + self.gap_ext,
					self.Ix[i, j - 1] + self.gap_ext,

				)
				# Check direction of seq2 best alignments, if equal then opt in direction M
				if (self.M[i - 1, j] + self.gap_start + self.gap_ext==self.Iy[i - 1, j] + self.gap_ext):
					idx=0
				else:
					idx = np.argmax([self.M[i - 1, j] + self.gap_start + self.gap_ext,Inf,self.Iy[i-1, j] + self.gap_ext])
				# Append direction to Yp
				Yp[i,j] = DirectionLegend[idx]
				# Fill score of seq2 will correspond to max score direction found
				self.Iy[i, j] = max(
					self.M[i - 1, j] + self.gap_start + self.gap_ext,
					self.Iy[i - 1, j] + self.gap_ext,
				)
		# will combine pointer matrix based on 0,1,2 direction matrix order and assign to self.pointer_mat
		pointer_matrix=[Mp,Xp,Yp]
		self.pointer_mat=pointer_matrix
		#Score is value found on bottom right corner of M score matrix, can now acsses score freely without run_traceback
		self.score = self.M[i, j]



	def run_traceback(self):

		seq1_align = ""
		seq2_align = ""
		align = ""
		#Initialize start point of the trace back, which is bottom right corner of M, also same value as self.score
		i = len(self.seq1)
		j = len(self.seq2)
		idx=0
		Mp,Xp,Yp=self.pointer_mat
		# will traverse back in the matrix until we reach zero for the index.
		while (i != 0 or j != 0):
				#If direction is based off M then we will go back one step diagonally, or one step back in both the x,y direction
				if Mp[i,j] == 'M':
					seq1_align = self.seq1[i-1] + seq1_align
					seq2_align = self.seq2[j-1] + seq2_align

					if (seq1_align[idx] == seq2_align[idx]):
						align = "|" + align
					else:
						align = "*" + align
					i = i - 1
					j = j - 1
				# Check if direction is based of seq1 being the best alignment, if so then we will only go back in one
				# step in the y direction
				elif Mp[i,j] == 'Ix':

					seq1_align = '-' + seq1_align
					seq2_align = self.seq2[j-1] + seq2_align
					align = " " + align
					j = j - 1
				# Check if direction is based of seq2 being the best alignment, if so then we will only go back in one
				# step in the x direction
				elif Mp[i,j] == 'Iy':

					seq1_align = self.seq1[i-1] + seq1_align
					seq2_align = '-' + seq2_align
					align = " " + align
					i = i - 1
				idx = idx + 1

		#Reutrn back best alignment score for each of the of the two sequneces and attach 3 new variables to the object.
		self.seq1_aligned=seq1_align
		self.seq2_aligned=seq2_align
		self.align_pattern=align

