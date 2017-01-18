#RuleTables
from rdkit import Chem

#(Grimsley 2009) with extensions
#Grimsley, G. R.; Scholtz, J. M. & Pace, C. N. (2009), 'A summary of the measured pK values of the ionizable groups in folded proteins.', Protein science : a publication of the Protein Society 18, 247--251.

#Format of ruletable tuple (Name, SMARTS, pKa (evt. as list), type (evt.as list) )
#Type is 0 = acidic (charge below pKa)
#Type 1 is basic (charge below pKa)

pkatable_grimsley = [
		('Tetramethylrhodamine #3',Chem.MolFromSmiles('CN(C)C1=C-C=C2-C(C4=C(C(=O)O)C=CC=C4)=C3-C=C-C(=[N+](C)C)C=C3OC2=C1'), [-10, 1.19,3.59, 10],[1, 1,0,1]), #Guess DOI: 10.1039/C0OB01045F 
		('6-hydroxy-3H-xanthen-3-one',Chem.MolFromSmarts('c2c3ccc(O)cc3oc3cc(=O)ccc2-3'),[2.94,6.56],[1,0]),
		('Anilin', Chem.MolFromSmarts('[$([NH2]c1ccccc1)]'),4.87,1),
		('N-acetyl Hydrazine', Chem.MolFromSmarts('[$([NH2]NC(=O))]'),2.81,1),#Reaxys
		('Hydrazine', Chem.MolFromSmarts('[$([NH2]N)]'),8.1,1),#Wikiped
		('N-acetyl Piperazine',	Chem.MolFromSmarts('[N+&!H0,NX3&H0;$(N1(C)CCN(C(=O)C)CC1)]'),7.1,1), #Reaxys
		('Piperazine', Chem.MolFromSmarts('C1CNCCN1'),9.73,1),#Wikiped
		('Pyridine',Chem.MolFromSmarts('n1ccccc1'),5.25,1),#Wikiped
		('Sulfonate',Chem.MolFromSmarts('[$([OH]S(=O)=O)]'), -1.5,0),#estimated
		('Phosphonate',Chem.MolFromSmarts('[OH]P([OH])(=O)[#6]'), [1.5, 7.0],[0,-1]),#estimated
		('7H-Pyrimido[4,5-b][1,4]oxazine-2,4-diamine',Chem.MolFromSmarts('C1C=NC2=CN=C(N=C2O1)N'),3.6,1), #PKA estimated
		('Tyrosine',Chem.MolFromSmarts('[$([O-,OH]c1ccccc1)]'),10.3,0),
		('Histidine', Chem.MolFromSmarts('[CH2X4][#6X3]1:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]:[$([#7X3H+,#7X2H0+0]:[#6X3H]:[#7X3H]),$([#7X3H])]:[#6X3H]1'),6.6,1),
		('Histidine1',Chem.MolFromSmarts('[$(cC)]1cnc[nH]1'),6.6,1),
		('Histidine2',Chem.MolFromSmarts('[$(cC)]1c[nH]cn1'),6.6,1),
		('Glutamic Acid', Chem.MolFromSmarts('[$([CX3][CD2][CD2]C(N)C(=O))](=O)[OX1H0-,OX2H1]'),4.2,0),
		('Cysteine',Chem.MolFromSmarts('[$([S-,SH][CD2]C(N)C(=O))]'),6.8,0),
		('Aspartic Acid v2',Chem.MolFromSmarts('[$([CX3][CD2]C(N)C(=O))](=O)[OX1H0-,OX2H1]') ,3.5,0),
		('Arginine',Chem.MolFromSmarts('[$(N[CD2][CD2])]C(=N)N'),12.0,1),
		('Lysine',Chem.MolFromSmarts('[$([ND1][CD2][CD2][CD2][CD2]C(N)C(=O))]'),10.5,1),		
		('End Carboxylate (avg)',Chem.MolFromSmarts('[$([CX3]CN)](=O)[OX1H0-,OX2H1]'),3.3,0),
		('End Amine (avg) PseudoAtom safe', Chem.MolFromSmarts('[$([NX3;H2,H1;!$(NC=O);!$(N[#171,#172,#173,#174,#175,#176,#177,#178,#179,#180,#181,#182,#183,#184,#185,#186,#187,#188,#189,#190,#191,#192,#193])]CC=O)]'),7.7,1),
		('Acrydine',Chem.MolFromSmarts('c1c2ccccc2nc2ccccc21'),6.15,1),
		('Generic Carboxylic Acid', Chem.MolFromSmarts('[CX3](=O)[OX1H0-,OX2H1]'),4.75,0),
		('Generic tert. Amine, not amide', Chem.MolFromSmarts('[N+&!H0&NX3,NX3;$(N(C)(C)C);!$(NC=[!#6]);!$(NC#[!#6])]'),9.8,1),
		('Generic Amine, not amide', Chem.MolFromSmarts('[N+,NX3;!H0;!$(NC=[!#6]);!$(NC#[!#6])!$(N[#171,#172,#173,#174,#175,#176,#177,#178,#179,#180,#181,#182,#183,#184,#185,#186,#187,#188,#189,#190,#191,#192,#193])]'),10.,1),
	]


condensed_grimsley = [
		('Condensed Tyrosine',Chem.MolFromSmarts('[#189]'),10.3,0),
		('Condensed Histidine', Chem.MolFromSmarts('[#179]'),6.6,1),
		('Condensed Glutamic Acid', Chem.MolFromSmarts('[#177]'),4.2,0),
		('Condensed Cysteine',Chem.MolFromSmarts('[#175][SH]'),6.8,0),
		('Condensed Aspartic Acid', Chem.MolFromSmarts('[#174]'),3.5,0),
		('Condensed Arginine',Chem.MolFromSmarts('[#172]'),12.0,1),
		('Condensed Lysine',Chem.MolFromSmarts('[#182]'),10.5,1),
		('Condensed AA end amino', Chem.MolFromSmarts('[#171,#172,#173,#174,#176,#177,#178,#179,#180,#181,#182,#183,#184,#185,#186,#187,#188,#189,#190,#191,#192,#193;D1]'),7.7,1),		
		('Condensed cys end amino', Chem.MolFromSmarts('[#175;D2]'),7.7,1),		
		('Condensed AA end Carboxylic', Chem.MolFromSmarts('[OH][#171,#172,#173,#174,#175,#176,#177,#178,#179,#180,#181,#182,#183,#184,#185,#186,#187,#188,#189,#190,#191,#192,#193]'),3.3,0)
		]


pkatable = pkatable_grimsley
condensedtable = condensed_grimsley
