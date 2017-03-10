#pICalculax
Isoelectric point (pI) predictor for chemically modified peptides and proteins.

For handling conversion of PLN to condensed molformat, Proteax desktop is needed: http://www.biochemfusion.com/products/proteax_desktop/
A modification database for import is found in the mods_db directory

For handling condensed molfile formats, RDKit needs to be patched. Patch can be found in the rdkit_patch directory.

Example usage can be found in the file Example_usage.py


##Example Usage interactive session
```
fasta = 'ICECREAM'

from peppKalib import find_pKas, pI
from rdkit import Chem

mol = Chem.MolFromFASTA(fasta)

#find pKa values and charge class
pkalist, charge = find_pKas(mol)
#Calculate pI
pIpred = pI(pkalist, charge)

print pIpred
```

The peptides can be loaded from a SDfile

```
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
```

With Proteax Desktop protein line notation of modified peptides can be converted to a RDKit mol object and the pI predicted
```
from proteax_desktop import *
prtx = ProteaxDesktop()

pln = 'H-GHANY[Gla]A-OH'

mol = Chem.MolFromMolBlock(prtx.as_molfile(pln,'expansion-mode=condensed'))

#find pKa values and charge class
pkalist, charge = find_pKas(mol)
#Calculate pI
pIpred = pI(pkalist, charge)

print pIpred
```

#Command line usage



