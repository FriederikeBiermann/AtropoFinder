
import sklearn
import pandas as pd
import pickle
import numpy as np
from Bio import SeqIO, pairwise2, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt

FASTAS_ALIGNED_BEFORE = True

INCLUDE_CHARGE_FEATURES = False



# enter names of your files here
filenames_alignment = ["Fasta/final_run_tryptorubins/identical_protein_groups_p450_300_450_alignment_1.fasta",
                       "Fasta/final_run_tryptorubins/identical_protein_groups_p450_300_450_alignment_2.fasta",
                       "Fasta/final_run_tryptorubins/identical_protein_groups_p450_300_450_alignment_3.fasta",
                       "Fasta/final_run_tryptorubins/identical_protein_groups_p450_300_450_alignment_4.fasta"
                       ]
FOLDERNAME_OUTPUT = "Output/Tryptorubinlike_peptides_25_01_2023/Classifiers_bigger_trainingsset/"
filename_permutations = FOLDERNAME_OUTPUT+"/permutations.txt"
filename_classifier = FOLDERNAME_OUTPUT + \
    '/RandomForest_optimized_classifier.sav'
filename_index = FOLDERNAME_OUTPUT+"/newindex.txt"
splitting_list = [["begin", 0, 92], ["sbr1", 93, 192], [
    "sbr2", 193, 275], ["core", 276, 395], ["end", 396, 400]]
fragments = ["begin", "sbr1", "sbr2", "core", "end"]

with open(filename_index, 'r') as name_fragment:
    index = [line.rstrip('\n') for line in name_fragment]

with open(filename_permutations, 'r') as name_fragment:
    permutations = [line.rstrip('\n') for line in name_fragment]


def merge_two_dicts(x, y):
    #input-> 2 dictionaries output->  merged dictionary
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z


def calculate_charge(sequence):
    # uses aa sequence as input and calculates the approximate charge of it
    AACharge = {"C": -.045, "D": -.999,  "E": -.998,
                "H": .091, "K": 1, "R": 1, "Y": -.001}
    charge = -0.002
    seqstr = str(sequence)
    seqlist = list(seqstr)
    for aa in seqlist:
        if aa in AACharge:
            charge += AACharge[aa]
    return charge


def easysequence(sequence):
    #creates a string out of the sequence file, that only states if AA is acidic (a), basic (b), polar (p), neutral/unpolar (n),aromatic (r),Cystein (s) or a Prolin (t)
    seqstr = str(sequence)
    seqlist = list(seqstr)
    easylist = []
    for i in seqlist:
        if i == 'E' or i == 'D':
            easylist = easylist+['a']
        if i == 'K' or i == 'R' or i == 'H':
            easylist = easylist+['b']
        if i == 'S' or i == 'T' or i == 'N' or i == 'Q':
            easylist = easylist+['p']
        if i == 'F' or i == 'Y' or i == 'W':
            easylist = easylist+['r']
        if i == 'C':
            easylist = easylist+['s']
        if i == 'P':
            easylist = easylist+['t']
        if i == 'G' or i == 'A' or i == 'V' or i == 'L' or i == 'I' or i == 'M':
            easylist = easylist+['n']

    seperator = ''
    easysequence = seperator.join(easylist)
    return easysequence


def indexing_reference(record):
    # index the reference sequence without ignoring gaps
    list_reference = list(str(record.seq))

    index_aa = 0
    index_mapping = []
    for index, AA in enumerate(list_reference):
        if AA != "-":
            index_aa += 1
            index_mapping.append([index_aa, index])

    return (index_mapping)


def convert_splitting_list(splitting_list, index_reference):
    #-> convert the canonic splitting list to also reflect eventual gaps in the reference sequence
    converted_splitting_list = []
    for fragment in splitting_list:
        converted_splitting_list.append(
            [fragment[0], index_reference[fragment[1]][1], index_reference[fragment[2]-1][1]])
    return converted_splitting_list


def split_alignment(alignment, fragment, fastas_aligned_before):
    # split the aligned sequences at the positions determined by the splitting list
    start = fragment[1]
    end = fragment[2]
    if fastas_aligned_before == False:
        alignment = [alignment]
    seqRecord_list_per_fragment = []
    if fragment[0] == "begin":
        start = 1
    if fragment[0] != "end":
        for record in alignment:
            subsequence = str(record.seq)[start-1:end].replace('-', '')

            seqRecord_list_per_fragment.append(
                [record.id, subsequence])
    else:
        for record in alignment:
            subsequence = str(record.seq)[start-1:].replace('-', '')
            seqRecord_list_per_fragment.append(
                                                    [record.id, subsequence])
    seqRecord_array_per_fragment = np.array(seqRecord_list_per_fragment)

    return seqRecord_array_per_fragment


def fragment_alignment(alignment, splitting_list, fastas_aligned_before):
    # create a matrix from the splitted alignment
    fragment_matrix = pd.DataFrame()
    if fastas_aligned_before == False:

        seqa = alignment[0]
        seqb = alignment[1]
        index_reference = indexing_reference(SeqRecord(Seq(seqa), id=seqa))

        converted_splitting_list = convert_splitting_list(
            splitting_list, index_reference)
        for fragment in converted_splitting_list:
            name_fragment = fragment[0]
            seqRecord_list_per_fragment = split_alignment(
                SeqRecord(Seq(seqb), id=seqb), fragment, fastas_aligned_before)

            fragment_matrix[name_fragment] = seqRecord_list_per_fragment[:, 1]
            fragment_matrix.set_index(
                pd.Index(seqRecord_list_per_fragment[:, 0]))
    else:
        for record in alignment:
            if record.id == "Reference":
                print("Found Reference")
                index_reference = indexing_reference(record)
                converted_splitting_list = convert_splitting_list(
                    splitting_list, index_reference)
                print (converted_splitting_list)
                for fragment in converted_splitting_list:
                    name_fragment = fragment[0]
                    seqRecord_list_per_fragment = split_alignment(
                        alignment, fragment, fastas_aligned_before)
                    fragment_matrix[name_fragment] = seqRecord_list_per_fragment[:, 1]
                    print(fragment_matrix)
                break
    fragment_matrix = fragment_matrix.set_index(pd.Index(seqRecord_list_per_fragment[:, 0]))
    print(fragment_matrix)
    return fragment_matrix


def featurize(fragment_matrix, permutations, fragments, include_charge_features):
    #create feature_matrix from fragment_matrix, count motifs in each fragemnt
    feature_matrix = pd.DataFrame()
    new_rows = []
    for index, row in fragment_matrix.iterrows():
        new_row = {}
        for fragment in fragments:
            sequence_fragment = row[fragment]

            easysequence_fragment = easysequence(sequence_fragment)
            for motif in permutations:
                name_column = motif+fragment
                new_row[name_column] = easysequence_fragment.count(motif)

            if include_charge_features == True:
                new_row = append_charge_features(new_row, fragment, easysequence_fragment, sequence_fragment)
        new_rows += [new_row]
    feature_matrix=feature_matrix.append(new_rows, ignore_index=True)
    if include_charge_features==True:
        feature_matrix=sum_charge_features(feature_matrix,fragments) 

    return feature_matrix


def append_charge_features(new_row, fragment, easysequence_fragment, sequence_fragment):
    #append features indicating the charge to the feature matrix
    acidic = fragment+"acidic"
    new_row = merge_two_dicts(new_row, {acidic: (
        easysequence_fragment.count("a")/(len(easysequence_fragment)+1))})
    acidic_absolute = fragment+"acidic absolute"
    new_row = merge_two_dicts(
        new_row, {acidic_absolute: (easysequence_fragment.count("a"))})
    charge_name = fragment+"charge"
    new_row = merge_two_dicts(
        new_row, {charge_name: (calculate_charge(sequence_fragment))})
    basic = fragment+"basic"
    basic_absolute = fragment+"basic absolute"
    new_row = merge_two_dicts(new_row, {basic: (
        easysequence_fragment.count("b")/(len(easysequence_fragment)+1))})
    new_row = merge_two_dicts(
        new_row, {basic_absolute: (easysequence_fragment.count("b"))})
    return new_row


def sum_charge_features(feature_matrix, fragments):
    #sum up charge features to obtain the charge of the whole protein
    chargerows = []
    acidicrows = []
    basicrows = []
    absacidicrows = []
    absbasicrows = []
    for fragment in fragments:
        chargerows.append(str(fragment)+"charge")
        acidicrows.append(str(fragment)+"acidic")
        basicrows.append(str(fragment)+"basic")
        absacidicrows.append(str(fragment)+"acidic absolute")
        absbasicrows.append(str(fragment)+"basic absolute")
    feature_matrix['complete charge'] = feature_matrix[chargerows].sum(axis=1)
    feature_matrix['mean acidic'] = feature_matrix[acidicrows].mean(axis=1)
    feature_matrix['mean basic'] = feature_matrix[basicrows].mean(axis=1)
    feature_matrix['absolute acidic'] = feature_matrix[absacidicrows].sum(
        axis=1)
    feature_matrix['absolute basic'] = feature_matrix[absbasicrows].sum(axis=1)
    return feature_matrix


# complete_feature_matrix = pd.DataFrame()
# for dataset in filenames_alignment:
#     if FASTAS_ALIGNED_BEFORE == True:
#         alignment = AlignIO.read(open(dataset), "fasta")
#         fragment_matrix = fragment_alignment(
#             alignment, splitting_list, FASTAS_ALIGNED_BEFORE)
#     if FASTAS_ALIGNED_BEFORE == False:
#         fragment_matrix = pd.DataFrame()
#         seq_record_ids = []
#         for seq_record in SeqIO.parse(dataset, "fasta"):
#             def fewgaps(x, y): return -20 - y
#             def specificgaps(x, y): return (-2 - y)
#             alignment = pairwise2.align.globalmc(
#                 alignmentfa, seq_record.seq, 1, -1, fewgaps, specificgaps)
#             fragment_matrix_for_record = fragment_alignment(
#                 alignment[0], splitting_list, FASTAS_ALIGNED_BEFORE)

#             fragment_matrix = fragment_matrix.append(
#                 fragment_matrix_for_record, ignore_index=True)
#             seq_record_ids = seq_record_ids+[seq_record.id]

#     feature_matrix = featurize(
#         fragment_matrix, permutations, fragments, INCLUDE_CHARGE_FEATURES)
#     print("Featurization completed")
#     complete_feature_matrix = complete_feature_matrix.append(
#         feature_matrix, ignore_index=True)
# feature_matrix_reindexed = complete_feature_matrix.reindex(columns=index)
# feature_matrix_reindexed.to_csv(
#     FOLDERNAME_OUTPUT + "/feature_matrix_ncbi_300_450_1.csv", index=False)
feature_matrix_reindexed = pd.read_csv(FOLDERNAME_OUTPUT + "/feature_matrix_ncbi_300_450_1.csv")
classifier = pickle.load(open(filename_classifier, 'rb'))
predictions = classifier.predict_proba(feature_matrix_reindexed)
# predicts all non Tryptorubinlike p450, if algorithm says "Tryptorubinlike"-> writes to fasta
print ("Prediction completed!")
counter = -1
predicitonsfile = pd.DataFrame(columns=['title', 'probability'])
listpositive = []
rows = []

for dataset in filenames_alignment:
    for seq_record in SeqIO.parse(dataset, "fasta"):
        counter += 1
        prediction = predictions[counter][1]
        if prediction > 0.1:
            row = {'title': seq_record.description, 'probability': prediction}
            rows += [row]
            listpositive.append(seq_record)
predicitonsfile = predicitonsfile.append(rows, ignore_index=True)
filenamelistpositive = FOLDERNAME_OUTPUT+"/listpositive_threshold_1_01.fasta"
print( len(listpositive))
SeqIO.write(listpositive, filenamelistpositive, "fasta")
pathcompletetable = FOLDERNAME_OUTPUT+"/predictiontablee_threshold_1_01.csv"
predicitonsfile.to_csv(pathcompletetable, index=False)
# plot histogramm of probabilities
probabilities = [probability[1] for probability in predictions]
plt.hist(probabilities, bins=40, edgecolor='blue')
plt.xlabel('Probability')
plt.ylabel('Sequence count')
plt.savefig(FOLDERNAME_OUTPUT + "/histogram_of_probabilities_1.png",
            format="png", dpi=1000)
plt.show()
