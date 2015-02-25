################################################################################################################
# Alignment editing program. It removes indels and stop codons and provides an option to adjust the alignment  #
# according to user defined protein or species sequence                                                        #
#                                                                                                              #
# Copyright (C) {2014}  {Ambuj Kumar, Kimball-Brain lab group, Biology Department, University of Florida}      #
#                                                                                                              #
# This program is free software: you can redistribute it and/or modify                                         #
# it under the terms of the GNU General Public License as published by                                         #
# the Free Software Foundation, either version 3 of the License, or                                            #
# (at your option) any later version.                                                                          #
#                                                                                                              #
# This program is distributed in the hope that it will be useful,                                              #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                               #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                #
# GNU General Public License for more details.                                                                 #
#                                                                                                              #
# This program comes with ABSOLUTELY NO WARRANTY;                                                              #
# This is free software, and you are welcome to redistribute it                                                #
# under certain conditions;                                                                                    #
#                                                                                                              #
################################################################################################################


import sys
import os
import platform
import subprocess


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
from Bio.Alphabet import IUPAC
from Bio.Seq import translate
from Bio.Alphabet import SingleLetterAlphabet



def _spliter(str, num):
    '''Splits the string object'''
    return [ str[start:start+num] for start in range(0, len(str), num) ]



def _groupy(L):
    first = last = L[0]
    for n in L[1:]:
        if n - 1 == last:
            last = n
        else:
            yield first, last
            first = last = n
    yield first, last


def _translator(recordData, ign, omit, table):
    proteinSeqList = list()
    recordsFunc = recordData
    
    for i, rec in enumerate(recordsFunc):
        counter = dict()
        seqT = rec.seq.translate()
        
        for j, obj in enumerate(seqT):
            if '*' in obj:
                seqT = seqT[:j] + 'Z' + seqT[j+1:]
        
        proteinSeqList.append(SeqRecord(Seq(str(seqT), IUPAC.protein), id=rec.id, name=rec.name, description=rec.description))


    with open('translated.fas', 'w') as fp:
        SeqIO.write(proteinSeqList, fp, 'fasta')
    
    
    return recordsFunc




def _alignP(pkg, arguments=None):
    if pkg == 'muscle':
        if 'Darwin' in platform.system():
            subprocess.call("./src/muscle/muscle -in translated.fas -out tAligned.fas", shell=True)
        else:
            subprocess.call("./src/muscle/muscleLinux -in translated.fas -out tAligned.fas", shell=True)
    
    else:
        arguments = arguments.replace('[', '').replace(']', '')
        subprocess.call("./src/mafft/mafft.bat %s translated.fas > tAligned.fas" %arguments, shell=True)


def _cleanAli(recordNuc, omit, fileName):
    handleP = open('tAligned.fas', 'rU')
    records = list(SeqIO.parse(handleP, 'fasta'))
    
    store = list()
    for i, rec in enumerate(records):
        nucData = [x.seq for x in recordNuc if x.id in rec.id]
        nucSeqData = _spliter(nucData[0], 3)
        sequence = Seq("", SingleLetterAlphabet()); pos = 0
        
        #print len([x for x in rec.seq if x!="-"]), len(nucSeqData)
        for j, amino in enumerate(rec.seq):
            if amino == '-':
                sequence = sequence + Seq("---", SingleLetterAlphabet())
            elif amino == "Z":
                sequence = sequence + Seq("NNN", SingleLetterAlphabet())
                pos = pos + 1
            else:
                if pos == 0 or pos == len(nucSeqData) - 1:
                    sequence = sequence + nucSeqData[pos].strip("N")
                else:
                    sequence = sequence + nucSeqData[pos]
                
                pos = pos + 1
        
        records[i].seq = Seq(str(sequence).strip("N"), SingleLetterAlphabet())


    #for rec in records:
    #print rec.seq


    with open(fileName, 'w') as fp:
        SeqIO.write(records, fp, "fasta")


    os.remove('translated.fas')
    os.remove('tAligned.fas')


def cdsAlign(inputFile, outFile, pkg='muscle', ign=True, CT=None):
    
    
    
    codonTables = ['Ascidian Mitochondrial', 'SGC9', 'Coelenterate Mitochondrial', 'Protozoan Mitochondrial', 'Vertebrate Mitochondrial', 'Plant Plastid', 'Thraustochytrium Mitochondrial', 'Blepharisma Macronuclear', 'Mold Mitochondrial', 'Invertebrate Mitochondrial', 'Standard', 'Trematode Mitochondrial', 'Scenedesmus obliquus Mitochondrial', 'Euplotid Nuclear', 'Yeast Mitochondrial', 'Spiroplasma', 'Alternative Flatworm Mitochondrial', 'Ciliate Nuclear', 'SGC8', 'Alternative Yeast Nuclear', 'Hexamita Nuclear', 'SGC5', 'SGC4', 'SGC3', 'SGC2', 'SGC1', 'SGC0', 'Flatworm Mitochondrial', 'Dasycladacean Nuclear', 'Chlorophycean Mitochondrial', 'Mycoplasma', 'Bacterial', 'Echinoderm Mitochondrial']
    
    omit = False
    
    if CT == None:
        table = CodonTable.ambiguous_dna_by_id[1]
    elif CT != None and CT in codonTables:
        table = CodonTable.ambiguous_generic_by_name[CT]
    else:
        table = CodonTable.ambiguous_generic_by_name['Standard']

    
    handle = open(inputFile, 'rU')
    records = list(SeqIO.parse(handle, "fasta"))
    print "\nPahle\n"
    for j, rec in enumerate(records):
        #print rec.seq
        if 'TAA' in rec.seq[-3:] or 'TGA' in rec.seq[-3:] or 'TAG' in rec.seq[-3:]:
            records[j].seq = rec.seq[0:-3]


    records = _translator(records, ign, omit, table)
    _alignP(pkg)
    _cleanAli(records, omit, outFile)











