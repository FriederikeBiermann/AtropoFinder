#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  2 10:11:37 2022

@author: Friederike Biermann
"""
import os
import re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Entrez
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
import argparse
from xml.etree import ElementTree

# create a parser object
parser = argparse.ArgumentParser(
    description="Corefiner makes identifying RiPP cores easy.")

# add arguments
parser.add_argument("-i", "--input", type=str, nargs=1,
                    metavar="file_name", default=None,
                    help="Opens and reads the specified fasta file.", required=True)

parser.add_argument("-b", "--boundary", type=int, nargs=1,
                    metavar="boundary", default=3000,
                    help="The window around the bait enzyme in bp in which to look for precursor genes.")

parser.add_argument("-d", "--dynamic_core_detection", type=bool, nargs=1,
                    metavar="parameter", default=True,
                    help="If dynamic core detection is enabled, the core will annotated starting from the amino acid before tryptophan, leaving the core at least 5 amino acids long. The maximum length of the core peptide will be 12 aa.")
# parser.add_argument("-m", "--min", type=int, nargs=1,
#                     metavar="min", default=10,
#                     help="Minimum legth of precursor gene in bp.")

# parser.add_argument("-M", "--max", type=int, nargs=1,
#                     metavar="max", default=60,
#                     help="Maximum legth of precursor gene in bp.")



parser.add_argument("-e", "--email", type=str, nargs=1,
                    metavar="min", default=None,
                    help="Email account to use for NCBI entrez.", required=True)

parser.add_argument("-o", "--output", type=str, nargs=1,
                    metavar="directory_name", default="Output/",
                    help="Output directory")

# parse the arguments from standard input
args = parser.parse_args()

boundary = args.boundary
dynamic_core_prediciton = args.dynamic_core_detection
gb_output_directory = args.output[0] + "Genbank_files/"
try:
    os.mkdir(args.output[0])
    os.mkdir(gb_output_directory)
except:
    print("WARNING: output directory already existing and not empty.")
stopcodon1 = ('TGA', 'TAA', 'TAG')
startcodon1 = ('ATG', "GTG")
stopcodon2 = ("TCA", "TTA", "CTA")
startcodon2 = ("CAT", "CAC")

filename_input = args.input[0]
filename_output = args.output[0] + "info_precursor_peptides.csv"
tableofproteins = pd.DataFrame(columns=['Protein Reference', 'NCBI Reference', 'Species', "Definition", "gb_whole_genome_accession",
                               "gb_whole_genome_position", "p450 count", "p450 position", "seqencednafragment", "openreadingframes", "coreproteins"])
Entrez.email = args.email[0]
filenamereadingframes = "readingframes_only_common_binding_sites.txt"
fastalist_openframe = []
fastalist_core = []


def find_between(s, first, last):
    try:
        start = s.index(first) + len(first)
        end = s.index(last, start)
        return s[start:end]
    except ValueError:
        return ""


def find_infrontof(s, last):
    try:
        start = 0
        end = s.index(last, start)
        return s[start:end]
    except ValueError:
        return ""


def find_after(s, first):
    try:
        start = s.index(first) + len(first)
        end = len(s)
        return s[start:end]
    except ValueError:
        return ""


def find_gi_infasta_header(record_id):
    if "ref" in record_id:
        print("ref")

        gi = find_between(record_id, "gi|", "|")
        ref = find_between(record_id, "ref|", "|")
    elif "gb" in record_id:
        print("gb")

        gi = find_between(record_id, "gi|", "|")
        ref = find_between(record_id, "gb|", "|")
    else:
        gi = record_id
        ref = record_id
    return gi, ref


def format_feature(feature):
    final_feature = feature["GBFeature_intervals"][0]
    for feature_qualifier in feature['GBFeature_quals']:
        try:
            final_feature[feature_qualifier['GBQualifier_name']
                          ] = feature_qualifier['GBQualifier_value']
        except:
            final_feature[feature_qualifier['GBQualifier_name']
                          ] = feature_qualifier
    return final_feature


listreadingframe = []
with open(filenamereadingframes, 'r') as f:
    listreadingframe = [line.rstrip('\n') for line in f]
rows =[]

# actually finding in fastas
for record in SeqIO.parse(filename_input, "fasta"):
    try:
        if "Mycobacteroides abscessus" in str(record):
            continue
        gi, ref = find_gi_infasta_header(record.id)

        print(gi, ref)
        # fetch from NCBI all genomes linked to it
        handle = Entrez.elink(dbFrom="protein", db="nucleotide", id=gi)
        xml = handle.read()
        records = ElementTree.fromstring(xml)
        listgenomes = []
        try:
            for genome_record in records:
                for index in range(0, len(genome_record[2])):
                    if genome_record[2][index].tag == "Link":
                      listgenomes.append(genome_record[2][index][0].text)
        except:
            for genome_record in records:
                for index in range(0, len(genome_record[1])):
                    if genome_record[1][index].tag == "Link":
                      listgenomes.append(genome_record[1][index][0].text)
        handle.close()
        print(listgenomes)
        for genome_id in listgenomes:
            print(genome_id)
            dnahandle = Entrez.efetch(
                db="nucleotide", id=genome_id, rettype='gbwithparts', retmode="xml")
            genomes = Entrez.parse(dnahandle)

            for genome in genomes:
                definition = genome["GBSeq_definition"]
                species = genome["GBSeq_organism"]
                print(species)
                CDS_index = 0
                complete_sequence = genome['GBSeq_sequence']
                GBInterval_from = None
                # find sequence in genomes
                for index_feature, feature in enumerate(genome['GBSeq_feature-table']):
                    if feature['GBFeature_key'] == "CDS":
                        if feature['GBFeature_quals'][-1]['GBQualifier_value'] == str(record.seq).replace("-", ""):
                            GBInterval_from = int(
                                feature['GBFeature_intervals'][0]['GBInterval_from'])
                            GBInterval_to = int(
                                feature['GBFeature_intervals'][0]['GBInterval_to'])
                            index_p450 = index_feature
                            xml_features = [
                                genome['GBSeq_feature-table'][index_feature-1], feature]
                            break
                if not GBInterval_from:
                    print("protein not in", definition)
                    continue
                framebegin = GBInterval_from-boundary
                if framebegin <= int(genome['GBSeq_feature-table'][0]['GBFeature_intervals'][0]['GBInterval_from']):
                    framebegin = int(
                        genome['GBSeq_feature-table'][0]['GBFeature_intervals'][0]['GBInterval_from'])
                frameend = GBInterval_to+boundary
                if frameend >= int(genome['GBSeq_feature-table'][len(genome['GBSeq_feature-table'])-1]['GBFeature_intervals'][0]['GBInterval_to']):
                    frameend = int(
                        genome['GBSeq_feature-table'][len(genome['GBSeq_feature-table'])-1]['GBFeature_intervals'][0]['GBInterval_to'])

                # obtain range of "region" for xml output
                index_protein_begin = index_p450-1
                begin_xml_interval = GBInterval_from
                p450count = 0

                while begin_xml_interval > framebegin:
                    index_protein_begin -= 2
                    xml_features = [genome['GBSeq_feature-table'][index_protein_begin],
                                    genome['GBSeq_feature-table'][index_protein_begin+1]]+xml_features
                    if "p450" in str(genome['GBSeq_feature-table'][index_protein_begin]) or "P450" in str(genome['GBSeq_feature-table'][index_protein_begin]):
                        p450count += 1
                    if int(genome['GBSeq_feature-table'][index_protein_begin]['GBFeature_intervals'][0]['GBInterval_from']) < int(genome['GBSeq_feature-table'][index_protein_begin]['GBFeature_intervals'][0]['GBInterval_to']):
                        begin_xml_interval = int(
                            genome['GBSeq_feature-table'][index_protein_begin]['GBFeature_intervals'][0]['GBInterval_from'])
                    else:
                        begin_xml_interval = int(
                            genome['GBSeq_feature-table'][index_protein_begin]['GBFeature_intervals'][0]['GBInterval_to'])
                index_protein_end = index_p450
                end_xml_interval = GBInterval_to
                while end_xml_interval < frameend:
                    if "p450" in str(genome['GBSeq_feature-table'][index_protein_begin]) or "P450" in str(genome['GBSeq_feature-table'][index_protein_begin]):
                        p450count += 1
                    index_protein_end += 2
                    xml_features += [genome['GBSeq_feature-table'][index_protein_end-1],
                                     genome['GBSeq_feature-table'][index_protein_end]]
                    if int(genome['GBSeq_feature-table'][index_protein_end]['GBFeature_intervals'][0]['GBInterval_from']) < int(genome['GBSeq_feature-table'][index_protein_end]['GBFeature_intervals'][0]['GBInterval_to']):
                        end_xml_interval = int(
                            genome['GBSeq_feature-table'][index_protein_end]['GBFeature_intervals'][0]['GBInterval_to'])
                    else:
                        end_xml_interval = int(
                            genome['GBSeq_feature-table'][index_protein_end]['GBFeature_intervals'][0]['GBInterval_from'])

                list_feature_locations = []
                begin_xml_interval -=1
                chosensequence = complete_sequence[begin_xml_interval:end_xml_interval].upper(
                )
                filename_genbank_file = f'{gb_output_directory}{genome_id}.gb'
                genbank_seq_record = SeqRecord(Seq(chosensequence), id=genome_id, annotations={
                                               "molecule_type": "dna"}, name=genome_id+"-"+str(begin_xml_interval)+"-"+str(end_xml_interval))
                seq_feature = SeqFeature(FeatureLocation(
                    0, end_xml_interval-begin_xml_interval), type="source", qualifiers=genome['GBSeq_feature-table'][0])
                genbank_seq_record.features += [seq_feature]

                for feature in xml_features:
                    formatted_feature = format_feature(feature)
                    if int(feature['GBFeature_intervals'][0]['GBInterval_from']) < int(feature['GBFeature_intervals'][0]['GBInterval_to']):
                        if feature['GBFeature_key'] == "CDS" and "tryptorubin family RiPP precursor CDS" not in feature:
                            list_feature_locations += [[int(feature['GBFeature_intervals'][0]['GBInterval_from']), int(
                                feature['GBFeature_intervals'][0]['GBInterval_to'])]]
                        feature['GBFeature_location'] = "complement("+str(int(feature['GBFeature_intervals'][0]['GBInterval_from'])-begin_xml_interval-1)+".."+str(
                            int(feature['GBFeature_intervals'][0]['GBInterval_to'])-begin_xml_interval)+")"
                        seq_feature = SeqFeature(FeatureLocation(int(feature['GBFeature_intervals'][0]['GBInterval_from'])-begin_xml_interval -1, int(
                            feature['GBFeature_intervals'][0]['GBInterval_to'])-begin_xml_interval, strand=1), type=feature['GBFeature_key'], qualifiers=formatted_feature)
                        genbank_seq_record.features += [seq_feature]
                    else:
                        if feature['GBFeature_key'] == "CDS" and "tryptorubin family RiPP precursor CDS" not in feature:
                            list_feature_locations += [[int(feature['GBFeature_intervals'][0]['GBInterval_to']), int(
                                feature['GBFeature_intervals'][0]['GBInterval_from'])]]
                        feature['GBFeature_location'] = "complement("+str(int(feature['GBFeature_intervals'][0]['GBInterval_from'])-begin_xml_interval-1)+".."+str(
                            int(feature['GBFeature_intervals'][0]['GBInterval_to'])-begin_xml_interval)+")"
                        seq_feature = SeqFeature(FeatureLocation(int(feature['GBFeature_intervals'][0]['GBInterval_to'])-begin_xml_interval, int(
                            feature['GBFeature_intervals'][0]['GBInterval_from'])-begin_xml_interval, strand=-1), type=feature['GBFeature_key'], qualifiers=formatted_feature)
                        genbank_seq_record.features += [seq_feature]

                dnahandle.close()
                listopenframe = ""
                listcore = ""
                for frame in listreadingframe:
                    for match in re.finditer(frame, chosensequence):
                        begin = match.span()[0]+begin_xml_interval
                        end = match.span()[1]+begin_xml_interval
                        overlapping = False
                        for CDS_location in list_feature_locations:

                            if begin > CDS_location[0] and begin < CDS_location[1] or end > CDS_location[0] and end < CDS_location[1]:
                                overlapping = True
                        if overlapping == False or overlapping == True:
                            if match.group()[:3] in startcodon1:
                                dna_openframe = Seq(match.group())
                                dna_openframe_record = SeqRecord(
                                    dna_openframe, id=genome_id)
                                aa = Seq(match.group()).translate()
                                if dynamic_core_prediciton == False:
                                    coredna = Seq(match.group()[-21:])
                                    coreaa = coredna.translate()
                                    
                                    coreaa_record = SeqRecord(
                                        coreaa, id=genome_id, name=ref)
                                else:
                                    preliminary_core = aa[-13:]
                                    coreaa = re.findall(
                                        r'.W.{3,}', str(preliminary_core))
                                    if len(coreaa)>0:
                                        coreaa = Seq(coreaa[0])
                                        print(type(coreaa))
                                    else:
                                        coreaa = preliminary_core
                                    coreaa_record = SeqRecord(
                                        coreaa, id=genome_id, name=ref)

                                if "*" in aa[:-1]:
                                    pass
                                if "*" not in aa[:-1] and "*" not in coreaa[:-1]:

                                    if str(begin) not in str(listcore) and str(end) not in str(listcore):
                                        print("+", coreaa)
                                        listopenframe = listopenframe+", "+match.group()+", "+str(begin) + \
                                            ", "+str(end)+"(+)"
                                        listcore = listcore+", " + \
                                            str(coreaa)+", "+str(begin) + \
                                            ", "+str(end)+"(+)"
                                        feature_gene_openframe = SeqFeature(FeatureLocation(
                                            begin-begin_xml_interval, end-begin_xml_interval, strand=1), type="gene", id="putative_precursor")
                                        feature_CDS_openframe = SeqFeature(FeatureLocation(begin-begin_xml_interval, end-begin_xml_interval, strand=1),
                                                                           type="CDS", id="putative_precursor", qualifiers={'Product': "Putative atropopeptide precursor", 'Translation': aa})

                                        genbank_seq_record.features += [
                                            feature_CDS_openframe]
                                        genbank_seq_record.features += [
                                            feature_gene_openframe]
                                        fastalist_openframe.append(
                                            dna_openframe_record)
                                        fastalist_core.append(coreaa_record)
                            if match.group()[:3] in stopcodon2:
                                reversedna = Seq(
                                    match.group()).back_transcribe().reverse_complement()
                                reversedna = str(reversedna)
                                aa = Seq(reversedna).translate()
                                if dynamic_core_prediciton == False:
                                    coredna = reversedna[-21:]
                                    coreaa = Seq(coredna).translate()
                                    coreaa_record = SeqRecord(
                                        coreaa, id=genome_id, name=ref)
                                else:
                                    preliminary_core = aa[-13:]
                                    coreaa = re.findall(
                                        r'.W.{3,}', str(preliminary_core))
                                    if len(coreaa) > 0:
                                        coreaa = Seq(coreaa[0])
                                        coreaa_record = SeqRecord(
                                            coreaa, id=genome_id, name=ref)
                                    else:
                                        coreaa = preliminary_core
                                        coreaa_record = SeqRecord(
                                            preliminary_core, id=genome_id, name=ref)                           
                                dna_openframe_record = SeqRecord(
                                    Seq(reversedna), id=genome_id, name=ref)
                                if "*" in aa[:-1]:
                                    pass
                                if "*" not in aa[:-1] and "*" not in coreaa[:-1]:
                                    if str(begin) not in str(listcore) and str(end) not in str(listcore):
                                        print("-", coreaa)
                                        listopenframe = listopenframe+", " + \
                                            str(reversedna)+", " + \
                                            str(begin)+", "+str(end)+"(-)"

                                        feature_gene_openframe = SeqFeature(FeatureLocation(
                                            begin-begin_xml_interval, end-begin_xml_interval, strand=-1), type="gene", id="putative_precursor")
                                        feature_CDS_openframe = SeqFeature(FeatureLocation(begin-begin_xml_interval, end-begin_xml_interval, strand=-1),
                                                                           type="CDS", id="putative_precursor", qualifiers={'Product': "Putative tryptorubin precursor", 'Translation': aa})
                                        genbank_seq_record.features += [
                                            feature_CDS_openframe]
                                        genbank_seq_record.features += [
                                            feature_gene_openframe]
                                        listcore = listcore+", " + \
                                            str(coreaa)+", "+str(begin) + \
                                            ", "+str(end)+"(-)"
                                        fastalist_openframe.append(
                                            dna_openframe_record)
                                        fastalist_core.append(coreaa_record)
                position4 = str(begin_xml_interval)+".."+str(end_xml_interval)

                new_row = {'Protein Reference': ref, 'NCBI Reference': genome_id, 'Species': species, "Definition": definition, "gb_whole_genome_accession": genome_id,
                           "gb_whole_genome_position": position4, "p450 count": p450count, "seqencednafragment": chosensequence, "openreadingframes": str(listopenframe), "coreproteins": str(listcore)}
                if listcore != "":
                    SeqIO.write(genbank_seq_record,
                                filename_genbank_file, 'gb')
                rows.append(new_row)
                
                dnahandle.close()

    except:
        print("Error")
tableofproteins = tableofproteins.append(
    rows, ignore_index=True)
tableofproteins.to_csv(filename_output, index=False)
SeqIO.write(fastalist_openframe,
            args.output[0] + "open_reading_frames.fasta", 'fasta')
SeqIO.write(fastalist_core, args.output[0] + "core_peptides.fasta", 'fasta')

