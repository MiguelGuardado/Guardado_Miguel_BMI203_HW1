import numpy as np



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


	def clean_seq(self,isprotein):
		protein=['A', 'R', 'N', 'D', 'C',  'Q',  'E',  'G',  'H',  'I',  'L',  'K',
				 'M',  'F',  'P',  'S',  'T',  'W',  'Y',  'V',  'B',  'Z',  'X']
		nucleotide=['A','C','G','T']
		seq1=""
		seq2=""

		if(isprotein==True):
			sequence_parameters=protein
		else:
			sequence_parameters = nucleotide
		for char in self.seq1.upper():
			if char not in sequence_parameters:
				seq1 += "*"
			else:
				seq1 += char
		for char in self.seq2.upper():
			if char not in sequence_parameters:
				seq2 += "*"
			else:
				seq2 += char
		self.seq1=seq1
		self.seq2=seq2


	def run_traceback(self):
		pass

	def fill_traceback(self):
		pass


class SmithWaterman(PairwiseAligner):
	def fill_traceback(self):
		Mp=np.zeros([len(self.seq1) + 1,len(self.seq2) + 1],dtype="<U10")
		Xp=np.zeros([len(self.seq1) + 1,len(self.seq2) + 1],dtype="<U10")
		Yp=np.zeros([len(self.seq1) + 1,len(self.seq2) + 1],dtype="<U10")


		for i, a1 in enumerate(self.seq1, start=1):
			for j, a2 in enumerate(self.seq2, start=1):
				Score = self.score_matrix[a1][a2]
				DirectionLegend = {0: "M", 1: "Ix", 2: "Iy"}
				if (self.M[i - 1, j - 1] == self.Ix[i - 1, j - 1] or self.M[i - 1, j - 1] == self.Iy[i - 1, j - 1]):
					idx = 0
				else:
					idx = np.argmax([self.M[i - 1, j - 1], self.Ix[i - 1, j - 1], self.Iy[i - 1, j - 1]])

				Mp[i, j] = DirectionLegend[idx]

				self.M[i, j] = max(
					0,
					self.M[i - 1, j - 1]+Score ,
					self.Ix[i - 1, j - 1]+Score,
					self.Iy[i - 1, j - 1]+Score
				)
				if (self.M[i, j - 1] + self.gap_start + self.gap_ext == self.Ix[i, j - 1] + self.gap_ext):
					idx = 0
				else:
					idx = np.argmax(
						[self.M[i, j - 1] + self.gap_start + self.gap_ext, self.Ix[i, j - 1] + self.gap_ext])

				Xp[i, j] = DirectionLegend[idx]

				self.Ix[i, j] = max(
					self.M[i, j - 1] + self.gap_start + self.gap_ext,
					self.Ix[i, j - 1] + self.gap_ext,

				)

				if (self.M[i - 1, j] + self.gap_start + self.gap_ext == self.Iy[i - 1, j] + self.gap_ext):
					idx = 0
				else:
					idx = np.argmax(
						[self.M[i - 1, j] + self.gap_start + self.gap_ext,-1, self.Iy[i - 1, j] + self.gap_ext])

				Yp[i, j] = DirectionLegend[idx]
				self.Iy[i, j] = max(
					self.M[i - 1, j] + self.gap_start + self.gap_ext,
					self.Iy[i - 1, j] + self.gap_ext,
				)

		pointer_matrix = [Mp, Xp, Yp]
		self.pointer_mat = pointer_matrix
		self.score=np.max([np.max(self.M),np.max(self.Ix),np.max(self.Iy)])


	def run_traceback(self):

		seq1_align = ""
		seq2_align = ""
		align = ""

		# Mp, Xp, Yp = self.pointer_mat
		mat_total = (self.M, self.pointer_mat[0]), (self.Ix, self.pointer_mat[1]), (self.Iy, self.pointer_mat[2])
		max_idx=np.argmax([np.max(self.M),np.max(self.Ix),np.max(self.Iy)])

		mat_val = [self.M, self.Ix, self.Iy]
		cur_mat,cur_dir=mat_total[max_idx]
		i,j = np.argwhere(cur_mat == np.max(cur_mat))[0]


		idx=0
		while cur_mat[i,j] > 0 and (i>0 and j>0):

			if cur_dir[i,j]=='M':
				seq1_align = self.seq1[i - 1] + seq1_align
				seq2_align = self.seq2[j - 1] + seq2_align

				if (seq1_align[idx] == seq2_align[idx]):
					align = "|" + align
				else:
					align = "*" + align
				i = i - 1
				j = j - 1
				cur_mat, cur_dir = mat_total[0]

			elif cur_dir[i,j]=='Ix':
				seq1_align = '-' + seq1_align
				seq2_align = self.seq2[j] + seq2_align
				align = " " + align
				j = j - 1
				cur_mat, cur_dir = mat_total[1]

			elif cur_dir[i,j]=='Iy':
				seq1_align = self.seq1[i] + seq1_align
				seq2_align = '-' + seq2_align
				align = " " + align
				i = i - 1
				cur_mat,cur_dir=mat_total[2]
			idx=idx+1


		self.seq1_aligned = seq1_align
		self.seq2_aligned = seq2_align
		self.align_pattern = align




class NeedlemanWunsch(PairwiseAligner):
	def fill_traceback(self):
		Inf = (float('-inf'))
		# I will first Initialize the three matrix M,Iy,Ix and input the initial values
		Mp=np.zeros([len(self.seq1) + 1,len(self.seq2) + 1],dtype="<U10")
		Xp=np.zeros([len(self.seq1) + 1,len(self.seq2) + 1],dtype="<U10")
		Yp=np.zeros([len(self.seq1) + 1,len(self.seq2) + 1],dtype="<U10")
		Mp[0, 0] = 0
		for i in range(1, len(self.seq1) + 1):
			self.M[i][0] = Inf
			self.Ix[i][0] = Inf
			self.Iy[i][0] = self.gap_start + i * self.gap_ext

			Mp[i][0] = 'Iy'
			Yp[i][0] = 'Iy'

		for i in range(1, len(self.seq2) + 1):
			self.M[0][i] = Inf
			self.Ix[0][i] = self.gap_start + i * self.gap_ext
			self.Iy[0][i] = Inf

			Mp[0][i] = 'Ix'
			Xp[0][i]= 'Ix'


		for i, a1 in enumerate(self.seq1, start=1):
			for j, a2 in enumerate(self.seq2, start=1):
				Score=self.score_matrix[a1][a2]
				DirectionLegend = {0: "M", 1: "Ix", 2: "Iy"}

				if (self.M[i - 1, j - 1]==self.Ix[i - 1, j - 1] or self.M[i - 1, j - 1]==self.Iy[i - 1, j - 1] ):
					idx=0
				else:
					idx = np.argmax([self.M[i - 1, j - 1], self.Ix[i - 1, j - 1], self.Iy[i - 1, j - 1]])

				Mp[i, j] = DirectionLegend[idx]

				self.M[i, j] = Score + max(
					self.M[i - 1, j - 1],
					self.Ix[i - 1, j - 1],
					self.Iy[i - 1, j - 1]
				)
				if (self.M[i, j - 1] + self.gap_start + self.gap_ext == self.Ix[i, j - 1] + self.gap_ext):
					idx=0
				else:
					idx = np.argmax([self.M[i, j - 1] + self.gap_start + self.gap_ext,self.Ix[i, j - 1] + self.gap_ext])

				Xp[i, j] = DirectionLegend[idx]

				self.Ix[i, j] = max(
					self.M[i, j - 1] + self.gap_start + self.gap_ext,
					self.Ix[i, j - 1] + self.gap_ext,

				)

				if (self.M[i - 1, j] + self.gap_start + self.gap_ext==self.Iy[i - 1, j] + self.gap_ext):
					idx=0
				else:
					idx = np.argmax([self.M[i - 1, j] + self.gap_start + self.gap_ext,Inf,self.Iy[i-1, j] + self.gap_ext])

				Yp[i,j] = DirectionLegend[idx]
				self.Iy[i, j] = max(
					self.M[i - 1, j] + self.gap_start + self.gap_ext,
					self.Iy[i - 1, j] + self.gap_ext,
				)

		pointer_matrix=[Mp,Xp,Yp]
		self.pointer_mat=pointer_matrix
		self.score = self.M[i, j]



	def run_traceback(self):
		seq1_align = ""
		seq2_align = ""
		align = ""
		i = len(self.seq1)
		j = len(self.seq2)
		idx=0
		Mp,Xp,Yp=self.pointer_mat

		while (i != 0 or j != 0):
				if Mp[i,j] == 'M':
					seq1_align = self.seq1[i-1] + seq1_align
					seq2_align = self.seq2[j-1] + seq2_align

					if (seq1_align[idx] == seq2_align[idx]):
						align = "|" + align
					else:
						align = "*" + align
					i = i - 1
					j = j - 1

				elif Mp[i,j] == 'Ix':

					seq1_align = '-' + seq1_align
					seq2_align = self.seq2[j-1] + seq2_align
					align = " " + align
					j = j - 1

				elif Mp[i,j] == 'Iy':

					seq1_align = self.seq1[i-1] + seq1_align
					seq2_align = '-' + seq2_align
					align = " " + align
					i = i - 1
				idx = idx + 1


		self.seq1_aligned=seq1_align
		self.seq2_aligned=seq2_align
		self.align_pattern=align



def main():
	Seq1 = NeedlemanWunsch("HEAGAWGHEE", "PAWHEAE", -5, -2, True)
	Seq1.LoadScoreMatrix('../scoring_matrices/BLOSUM50.mat')
	Seq1.fill_traceback()
	Seq1.run_traceback()
	print(Seq1.seq1_aligned)
	print(Seq1.align_pattern)
	print(Seq1.seq2_aligned)



	# Seq1 = SmithWaterman("HEAGAWGHEE", "PAWHEAE", -8, -8,True)
	#
	# Seq1.LoadScoreMatrix('../scoring_matrices/BLOSUM50.mat')
	# Seq1.fill_traceback()
	# Seq1.run_traceback()




if __name__ == "__main__":
	main()
