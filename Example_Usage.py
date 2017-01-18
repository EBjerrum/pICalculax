#!/usr/bin/python
""" Example usage of the peppKalib for pI prediction of unmodified and modified peptides """

from __future__ import print_function
from peppKalib import find_pKas, pI
from rdkit import Chem
from rdkit.Chem import Draw


#Load a protein from SD file in condensed format
sdsup = Chem.SDMolSupplier('Datasets/example_mols.sdf')

def predict_show(mol):
	#Get list of identified pKa values and charge
	pkalist, charge = find_pKas(mol)
	#Predict the pI from the identified pKa values
	piPred = pI(pkalist, charge)
	#Report and Visualize
	print("Predicted pI:%0.2F"%piPred)
	Draw.ShowMol(mol, legend = "Predicted pI:%0.2F"%piPred)
	Draw.tkRoot.update()
	txt = raw_input('Press <ENTER> to continue')

# An unmodified peptide
mol = sdsup[0]
predict_show(mol)

# A peptide with a modification
mol = sdsup[1]
predict_show(mol)


#Proteax can be used for manipulation of sequence <=> condensed <=> full formula
from proteax_desktop import *
prtx = ProteaxDesktop()

mol = Chem.MolFromMolBlock(prtx.as_molfile('H-GHANYEA-OH','expansion-mode=condensed'))
predict_show(mol)

mol = Chem.MolFromMolBlock(prtx.as_molfile('H-GHANY[Gla]A-OH','expansion-mode=condensed'))
predict_show(mol)

