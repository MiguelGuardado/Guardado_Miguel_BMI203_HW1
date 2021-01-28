from align import algs
import numpy as np
from sklearn import metrics
import matplotlib.pyplot as plt

#This is the function used to read in fasta, same as the one used inside algs but does not require creating an
#instant of an object.
def read_fasta(filename):
    contents = []
    with open(filename, "r") as f:
        for line in f:
            line = line.strip('\n')
            contents.append(line)
    seq = ''.join(contents[1:])
    return( seq )


#This will be used to load in PosPairs.txt/NegPairs.txt
def read_sequence_tests(filepath):
    seq_test=[]
    with open(filepath,"r") as f:
        for line in f:
            line=line.split()
            line= ["../" + l for l in line]
            seq_test.append(line)
    return (seq_test)

#This will be used to create roc plots for each of the test preformed, fpr/tpr for making the roc plot and the
# other input parameters are to create a proper file name for all the tests.
def plot_roc_curve(fpr, tpr,gap_open_cost,gap_ext_cost,Score_Mat):
    plt.plot(fpr, tpr, color='orange', label='ROC')
    plt.plot([0, 1], [0, 1], color='darkblue', linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC) Curve')
    plt.legend()
    FileName='../Images/PosNegRoc_GapOpen_'+str(gap_open_cost)+'_GapExt_'+str(gap_ext_cost)+'ScoreMat_'+str(Score_Mat)+'.png'
    plt.savefig(FileName)

    plt.close()



#This will input a pair of positive/negative sequence files to test, and require the gap open/ext cost for
# each of the tests preformed, alignment_test should only get inputted 0/1 depending if you want a
#SW test (0), or NW(1). Finally you also input the file path to the score matrix to input and read for the test.


'''
We will get outputted 4 different parameters

TotalScore= np.array(len(Pos)+len(Neg))array of the total score outputted by each sequence alignment test
threshold= integer of the threshold of each value
confusion_matrix = np.array(2,2) confusion matrix of the alignment test based on the threshold defined earlier 
auc= area under the roc curve for the current test ran


'''
def return_score(Pos, Neg, gap_open_cost, gap_ext_cost, alignment_test,score_file):
    PosScores=[]
    NegScores=[]
    score_name=score_file

    score_file="../scoring_matrices/"+score_file
    #Iterate through the entire length of pos, which should be the same length of index
    for i in range(0,len(Pos)):
        #Read in iteration of test for pos/neg at the same index i
        PosSeq1=read_fasta(Pos[i][0])
        PosSeq2=read_fasta(Pos[i][1])
        NegSeq1=read_fasta(Neg[i][0])
        NegSeq2=read_fasta(Neg[i][1])
        if(alignment_test==0):
            CurrSeqTestPos=algs.SmithWaterman(PosSeq1, PosSeq2, gap_open_cost, gap_ext_cost, True)
            CurrSeqTestNeg = algs.SmithWaterman(NegSeq1, NegSeq2, gap_open_cost, gap_ext_cost, True)
        else:
            CurrSeqTestPos = algs.NeedlemanWunsch(PosSeq1, PosSeq2, gap_open_cost, gap_ext_cost, True)
            CurrSeqTestNeg = algs.NeedlemanWunsch(NegSeq1, NegSeq2, gap_open_cost, gap_ext_cost, True)

        CurrSeqTestPos.LoadScoreMatrix(score_file)
        CurrSeqTestPos.fill_traceback()
        PosScores.append(CurrSeqTestPos.score)

        CurrSeqTestNeg.LoadScoreMatrix(score_file)
        CurrSeqTestNeg.fill_traceback()
        NegScores.append(CurrSeqTestNeg.score)

    PosScores=np.array(PosScores)
    NegScores=np.array(NegScores)
    TotalScores=np.concat=np.concatenate((PosScores,NegScores),axis=None)

    #np.savetxt("TempArray.txt", TotalScores)
    #TotalScores=np.loadtxt("TempArray.txt")
    #Get threshold based on the mean of total scores
    threshold=np.mean(TotalScores)
    confusion_matrix=np.zeros((2,2))

    y_true = np.zeros(100)
    y_true[0:50] = 1

    #create confusion matrix based on the threshold
    for score in PosScores:
        if (score > threshold):
            confusion_matrix[0,0]=confusion_matrix[0,0]+1
        else:
            confusion_matrix[0, 1] = confusion_matrix[0, 1] + 1
    for score in NegScores:
        if (score < threshold):
            confusion_matrix[1, 1] = confusion_matrix[1, 1] + 1
        else:
            confusion_matrix[1, 0] = confusion_matrix[1, 0] + 1


    fpr,tpr,threshold = metrics.roc_curve(y_true,TotalScores)
    plot_roc_curve(fpr, tpr, gap_open_cost, gap_ext_cost,score_name)

    auc=metrics.roc_auc_score(y_true,TotalScores)

    return(TotalScores,threshold,confusion_matrix,auc)











def main():
    #
    # #This code will be used for questions 1-4.

    #load in Pospair/Negpair
    Pospair=read_sequence_tests('../scoring_matrices/Pospairs.txt')
    Negpair = read_sequence_tests('../scoring_matrices/Negpairs.txt')

    #Run test for pos/neg pair with open_gap=-11 and gap_ext=-3 on a smith-waterman(0) test, using a BLOSUM50.mat
    TotalScores,threshold,confusion_matrix,auc=return_score(Pospair,Negpair,-11,-3,0,"BLOSUM50.mat")

    #Create histogram for the Totalscore
    plt.hist(TotalScores)
    plt.xlabel('Scores')
    plt.ylabel('Count')
    plt.title('Histogram of Pos/Neg scores with gap_open=-11 and gap_ended=-3')
    filepath='../Images/sw_histogram_gapopen_-11_gapext_-3.png'
    plt.savefig(filepath)
    plt.close()

    # This code will be used for questions 5-8

    #Create matrix for varying rates of gap_open/gap_close
    gap_open=np.arange(1,21)
    gap_ext=np.arange(1,6)
    auc_mat=np.zeros((20,5))

    #Iterate through rates of gap_open/gap_test to obtain auc score for each test, and append to the matrix that is indexd
    #by the test run. row=gap_open, col=gap_ext.
    for i in range(0,20):
        for j in range(0,5):
            print("Running Test on:",i,j)
            Gap_o=-gap_open[i]
            Gap_e=-gap_ext[j]
            TotalScores, threshold, confusion_matrix, auc = return_score(Pospair, Negpair, Gap_o, Gap_e, 0,"BLOSUM62.mat")

            auc_mat[i,j]=auc
            print(auc)
            print(auc_mat)

    #np.save("../scoring_matrices/auc_mat.txt",auc_mat)

    #Load in auc_mat so you dont have to rerun the test
    auc_mat=np.load("../scoring_matrices/auc_mat.txt.npy")

    i,j = np.argwhere(auc_mat == np.max(auc_mat))[0]
    #Index where the highest auroc value is found
    print(i,j)

    #Run test on the optimal gap open/ext scores for each of the 4 different scoring matricies
    ScoreMatricies=["BLOSUM50.mat","BLOSUM62.mat","PAM100.mat","PAM250.mat"]

    GlobalAuc=[]
    for score in ScoreMatricies:
        print(score)
        TotalScores, threshold, confusion_matrix, auc = return_score(Pospair, Negpair, -i, -j, 1, score)
        GlobalAuc.append(auc)

    print(GlobalAuc)






if __name__ == "__main__":
    main()