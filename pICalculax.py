#PepPkaLib
"""Functions for working with pKa values"""

from rdkit import Chem
from rdkit.Chem import AllChem
from numpy import unique, argmin
from numpy import abs as np_abs
from Rules import pkatable, condensedtable

def find_pKas(mol, pkatable=pkatable, condensedtable=condensedtable, debug=False, returnindices=False):
	"""This function finds pKa values in a supplied molecule by matching it up to SMARTS rules in a table
	Returns the list of pKa values and a list with the charge of the group at pH levels blow the pKa.
	Example: So histidine has a +1 charge at pH levels below the pKa of ~6"""
	chargelist = []
	pkalist = []
	dummy = Chem.MolFromSmarts('[#178]') #Used Gly instead of #0 or *, for compatibility with building blocks

	if returnindices:
		atomindices = []
		for atom in mol.GetAtoms():
			atom.SetProp('O_i', str(atom.GetIdx()))  #Save the Original Idx as a property

	#Look for condensed AA atoms with potential pKa values. Substitute with Glycine
	for s in condensedtable:
		if mol.HasSubstructMatch(s[1]):
			if debug: print "Match: %s"%(s[0]),
			#How many??
			substructmatches = mol.GetSubstructMatches(s[1])
			n = len(substructmatches)
			if debug: print "times %s"%(n)
			pkalist += [s[2]] * n
			chargelist += [s[3]] * n
			if returnindices:
				for match in substructmatches:
					atomindices += [[int(mol.GetAtomWithIdx(idx).GetProp('O_i')) for idx in match]] #Retrieve the original index of matches
			mol = AllChem.ReplaceSubstructs(mol, s[1], dummy)[0] # Delete acidic/basic groups to prevent multiple matching

	#Look for molecular fragments with known pKa's
	for s in pkatable:
		#print s[0]
		if mol.HasSubstructMatch(s[1]):
			if debug: print "Match: %s"%(s[0]),
			#How many??
			substructmatches = mol.GetSubstructMatches(s[1])
			n = len(substructmatches)
			if debug: print "times %s"%(n)
			if type(s[2])==type(float()):
				pkalist += [s[2]] * n
				chargelist += [s[3]] * n
			else: #If there's multiple pKa values found for the molecular substructure
				for i in range(len(s[2])):
					pkalist += [s[2][i]] * n
					chargelist += [s[3][i]] * n
			if returnindices:
				for match in substructmatches:
					atomindices += [[int(mol.GetAtomWithIdx(idx).GetProp('O_i')) for idx in match]] #Retrieve the original i
			mol = AllChem.ReplaceSubstructs(mol, s[1], dummy)[0]			
	if returnindices:
		return pkalist, chargelist, atomindices
	return pkalist, chargelist

def charge(ph,pkalist,chargelist):
	"""Jacob Tolborgs charge model where the charge is assigned from partial charges from all pKa values at the pH point"""
	chargesum = []
	for charge,pka in zip(chargelist, pkalist):
		#print charge, pka
		if charge == 1:
			charge = 1/(1+10**(ph-pka))
			chargesum.append(charge)
		else:
			charge = -1/(1+10**-(ph-pka))
			chargesum.append(charge)
	return sum(chargesum)

def pI(pkalist,chargelist):
	"""Uses Jacob Tolborgs charge function and calculates pI.
	If only only acidic or basic groups are found, the pI is ill defined and returned as 42 or -42"""
	#Check if Only acidic or basic groups are present.
	if len(unique(chargelist))== 0: #No pKa groups present.
		pI = 21
	elif len(unique(chargelist))== 1: # Only one type is present
		if chargelist[0] == 0: #Only acidic groups are present
			pI = -42
		elif chargelist[0] == 1: #Only basic groups are present
			pI = 42 #42 The answer to everything ;-)
	else:
		#Find pI by simulation in the pH range 0 to 0
		chargecol = []
		for i in range(0,1400):
			ph= i/100.
			chargecol.append(charge(ph,pkalist,chargelist)) #Calculate charge
		pI = argmin(np_abs(chargecol))/100. # Simply taking the smallest absolute value, and dividing the index with 100
	#print "pI %.1f"%pI
	return pI


def pred_mol(mol, sequence):
	#find pKa values and charge class
	pkalist, charge = find_pKas(mol)
	#Calculate pI
	pIpred = pI(pkalist, charge)
	print("%0.2F \t%s"%(pIpred, sequence))

def pred_fasta(fasta):
	for sequence in fasta:
		mol = Chem.MolFromFASTA(sequence)
		if mol != None:
			pred_mol(mol, sequence)
		else: print("Conversion error with %s"%sequence)

def pred_pln(pln):
	from proteax_desktop import ProteaxDesktop
	prtx = ProteaxDesktop()
	for sequence in pln:
		mol = Chem.MolFromMolBlock(prtx.as_molfile(sequence,'expansion-mode=condensed'))
		if mol != None:
			pred_mol(mol, sequence)
		else: print("Conversion error with %s"%sequence)


if __name__ == '__main__':
	#parse command line
	import argparse
	parser = argparse.ArgumentParser(description='Predict isoeletric point pI of peptides and modified peptides')
	parser.add_argument('--fasta', nargs='+', type=str, help='Predict fasta sequence')
	parser.add_argument('--pln', nargs='+', type=str, help='Predict PLN sequence (Requires Proteax Desktop)')
	args = parser.parse_args()
	# Predict lists.
	if args.fasta != None: pred_fasta(args.fasta)
	if args.pln != None: pred_pln(args.pln)
