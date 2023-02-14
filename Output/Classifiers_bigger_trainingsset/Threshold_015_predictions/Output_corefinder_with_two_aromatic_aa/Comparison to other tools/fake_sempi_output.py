from Bio.Seq import Seq
from Bio import Entrez
from os import scandir
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import copy
import pandas as pd
INPUT_FOLDER = "Output/Tryptorubinlike_peptides_25_01_2023/Classifiers_bigger_trainingsset/Thrshold_015_predictions/genbank_files/"
whole_length_before = 0
cluster_starts = []
cluster_ends = []
dataframe = pd.DataFrame()
start = None
with scandir(INPUT_FOLDER) as folder:
    for filename in folder:
        for record in SeqIO.parse(filename.path, "gb"):
            positions = []
            for feature in record.features:
                if ("p450" in str(feature).lower() or "tryptorubin" in str(feature).lower()):
                    positions += [feature.location.start + 1,
                                   int(feature.location.end)]
            start = min (positions) + whole_length_before
            end = max(positions) + whole_length_before
            if start and end:
                cluster_starts.append(start)
                cluster_ends.append(end)
            start = None
            end = None
            whole_length_before += len(str(record.seq))
print(cluster_starts)
dataframe["Start"] = cluster_starts
dataframe["Stop"] = cluster_ends
dataframe["Type"] = "ripp"
dataframe.index.name = "Cluster"
dataframe.index +=1
dataframe.to_csv("AtropoFinderOutput.csv")