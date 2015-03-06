################################################################################################################
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


import os
import sys

fullpath = os.getcwd() + "/src/Utils"
sys.path.append(fullpath)

import time
import glob
import shutil
import argparse
import textwrap
import warnings

from src.mappy import fetcher
from src.fetchSeq import cdsImport
from src.aligner import cdsAlign
from Bio import AlignIO, SeqIO
from StringIO import StringIO
from Bio.Alphabet import IUPAC, Gapped


parser = argparse.ArgumentParser(prog='Consensus',
                                 version= 'Consensus-1.0',
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=textwrap.dedent('''\
    ----------------------------------------------------------------------------------------------------------
    \t\t\t\t Welcome to ConCat-1.0 alignment handling program
    \t\t\t ConCat-Align conducts sequence alignment using Muscle/Mafft 
    \t\t\t sequence alignment programs.
    \t\t\tWritten by Ambuj Kumar, University of Florida
    \t      ConCat sequence alignment concatenation program was designed at Kimball-Braun
    \t\tlab group to handle large amount of sequence dataset for phylogeny study.
    
    ----------------------------------------------------------------------------------------------------------
    
    '''))


parser.add_argument('-o', type=str, default=None,
                    help='Raw sequence output folder name')

parser.add_argument('-log', type=str, default=None,
                    help='Enter log file name')

parser.add_argument('-cds', type=str, default=None,
                    help='Takes gene name via file for CDS import')

parser.add_argument('-orgn', type=str, default=None,
                    help='Takes group name to extract sequence data')

parser.add_argument('-ortho', type=str, default=None,
                    help='Takes species name for orthologue search')

parser.add_argument('-pkg', type=str, default='muscle',
                    choices=['muscle', 'mafft'],
                    help='User defined program selection')

parser.add_argument('-args', type=str,
                    help='Arguments to run MAFFT. EXAMPLE: "--retree 2 --maxiterate 10"')

parser.add_argument('-argf', type=str, default=None,
                    help='Takes argument file as input')


argmnts = parser.parse_args()


if argmnts.pkg == 'mafft' and not argmnts.args:
    parser.error('-args argument is required in "mafft" mode.')


if argmnts.cds == True and argmnts.orgn == None and argmnts.ortho == None:
    parser.error('-orgn or -ortho argument is required in "-cds" mode.')



try:
    os.mkdir("Align")
except OSError:
    shutil.rmtree("Align")
    os.mkdir("Align")




def _warnfxn():
    warnings.warn("deprecated", DeprecationWarning)


with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    _warnfxn()

def _is_empty(any_structure):
    if any_structure:
        return False
    else:
        return True


def _remDuplicate(filename):
    files = glob.glob('Align/' + filename)
    for file in files:
        handle = open(file, 'rU')
        records = list(SeqIO.parse(handle, 'fasta'))
        store = list()
        print("Removing duplicate sequences from file %s" %file)
        for x in records:
            for y in [z for z in records if z.id != x.id]:
                if x.id.split('|')[0] == y.id.split('|')[0]:
                    if x.seq.count('-') > y.seq.count('-'):
                        store.append(x.id)
                    elif x.seq.count('-') <= y.seq.count('-'):
                        store.append(y.id)
    
        newRec = [x for x in records if x.id not in store]
        with open(file, 'w') as fp:
            SeqIO.write(newRec, fp, 'fasta')


def _remGeneDuplicate(filename):
    handle = open(filename, 'rU')
    records = list(SeqIO.parse(handle, "fasta"))
    newRec = list()
    store=list()
    for x in records:
        for y in records:
            if x.id == y.id:
                flag = True
                if len(x.seq) > len(y.seq) and x.id not in [z.id for z in newRec]:
                    newRec.append(x)
                elif len(x.seq) < len(y.seq) and y.id not in [z.id for z in newRec]:
                    newRec.append(y)
            elif x.id != y.id and x.id not in [z.id for z in newRec]:
                newRec.append(x)

    with open("Align/" + argmnts.orgn + ".fas", 'w') as fp:
        SeqIO.write(newRec, fp, 'fasta')


def main():
    if argmnts.cds != None:
        store = list()
        
        
        try:
            os.mkdir("Sequences")
        except OSError:
            shutil.rmtree("Sequences")
            os.mkdir("Sequences")
        
        try:
            genes = [x for x in open(argmnts.cds, 'r').readlines() if x != '' and x != '\n']
        except IOError:
            print("\n%s file not found. Using %s as gene name\n" %(argmnts.cds, argmnts.cds))
            genes = [argmnts.cds]
        
        for geneName in genes:

            warnings.filterwarnings("ignore")
            try:
                cdsImport(geneName.rstrip('\n'), argmnts.orgn, argmnts.ortho)
            except:
                print("\nConnection Failure. Rebooting program in 10 seconds\n")
                time.sleep(10)
                try:
                    cdsImport(geneName.rstrip('\n'), argmnts.orgn, argmnts.ortho)
                except:
                    print("\nConnection Failure. Rebooting program in 10 seconds\n")
                    time.sleep(10)
                    try:
                        cdsImport(geneName.rstrip('\n'), argmnts.orgn, argmnts.ortho)
                    except:
                        print("Failed to import sequences for %s" %geneName.strip("\n"))
    
    listObject_human = [x.strip('\n') for x in genes]

    fetcher(listObject_human, output=argmnts.o, logObj = "exonName.log")



if __name__ == "__main__":
    main()





