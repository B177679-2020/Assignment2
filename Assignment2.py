#!/usr/bin/python3

##FIRST STEP: GET ALL THE SEQUENCES

import os
import subprocess

#Example query to get the birds G6P protein sequences in fasta format
#ESearchBirds = 'esearch -db protein -query "txid8782[Organism] AND glucose-6-phosphatase[Protein]" | efetch -db protein -format fasta > birds_G6P.fasta'

#ESearchBText = subprocess.check_output(ESearchBirds, shell=True)

TaxonID = input("Please type in your Taxon ID (4-number digit):")
Protein = input("Please type in the protein of interest (ex.: phosphatase, kinase...)")
GeneralES = 'esearch -db protein -query "{}[Organism] AND {}[Protein]" | efetch -db protein -format fasta > {}.fasta'.format(TaxonID, Protein, TaxonID)

GeneralESout = subprocess.check_output(GeneralES, shell=True)

 
