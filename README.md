# reassigner
This software reads pdb files (or downloads them from the web) and assigns secondary structure. These assignments are based on the synthesis of two different models: a virtual geometry model and a phi/psi2 motif model. This software can identify many aspects of secondary structure, including beta-bulges. Human readable or machine readable output is available.

At its core, there are two methods for structure identification, a virtual geometry-based method first described by King and Johnson (1) with some modifications made, and an application of (φ/ψ) 2 motifs (a metric defined by Hollingsworth et al. in 2012) (2) to determine secondary structure. The software is able to intelligently combine these two interpretations of structure in order to identify a very wide range of secondary structures, including β-bulges and non-hydrogen bonded conformations.

![Example of Reassigner vs DSSP]
(1osh320-340example.png)

	DSSP:		HHHHHHHHHH	HHHHHHHHHI	IIIIS HHHH
	Reassigner:	HHHHHHtTHH	HHHHHHHHHH	Httt-THHNH	
	Residue Type	QIALLKGSAV	EAMFLRSAEI	FNKKLPSGHS	
	Residue Number	320		330		340 


Shown here is an example of the Reassigner vs the more commonly used program DSSP. The protein in the figure is PDB code 1OSH, and the residues are 320-340. It's clear from the example that the DSSP assignments (red text) misidentify the mid-helix turn as more helix residues, whereas the Reassigner assignments (blue text) correctly identify it. Furthermore, residues 339 and 340 are labeled as "I" (a pi-helix) by DSSP, while it is clear from visual inspection that this is still an ordinary alpha-helix residue, as identified by the Reassigner.



# References:
	
1: S. M. King, W. C. Johnson, Assigning secondary structure from protein coordinate data. Proteins Struct. Funct. Bioinforma. 35, 313–320 (1999).
	
2: S. A. Hollingsworth, M. C. Lewis, D. S. Berkholz, W.-K. Wong, P. A. Karplus, (φ,ψ)2 Motifs: A Purely Conformation-Based Fine-Grained Enumeration of Protein Parts at the Two-Residue Level. J. Mol. Biol. 416, 78–93 (2012).

# Help:

	usage: reassigner.v14.py [-h] [-i FILE] [-c CODE] [-o FILE] [--legend] [--verbose]
	
	Assigns secondary structure without using hydrogen bonds. One method uses
	virtual dihedrals and bond angles, the other using phi/psi2 motifs to identify
	secondary structure.
	
	optional arguments:
		-h, --help	show this help message and exit
		-i FILE	This should be a pdb file. Do not use in combination with the "-c" option.
		-c CODE	This should be a four letter pdb code. Only use this option if you want to download 
				directly from the PDB.
		-o FILE	Output file, a table using tabs as seperators. If not specified, default is to output
				to STDOUT as human readable output.
		--legend	Option to print a legend for the Secondary Structure codes.
		--verbose	Option to print a all the output, including the behind the scenes methods 
				for structure assignment.
	
	IMPORTANT NOTES: Expect there to be strange ? residues at the end of the
	output if there are any ligands or other non-water HETATOMS present. This can
	be safely ignored.
 
# Legend:
	
	_:	Break in the chain (as detected by checking bond distance)
	-:	Unassigned trans-residue
	=:	Unassigned cis-residue
	P:	PII-Helix
	t:	Turn defined by CA to CA Distance
	N:	Non-Hydrogen Bonded Turn (Similar to 's')
	T:	Hydrogen Bonded Turn T
	E:	Extended or Beta-strand conformation
	H:	Alpha Helical Conformation
	G:	3,10 Helix
	Bb:	Beta-bulge
	U:	Pi-helical bulge
	X:	The Stig
