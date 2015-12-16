# reassigner
This software reads pdb files (or downloads them from the web) and assigns secondary structure. These assignments are based on the synthesis of two different models: a virtual geometry model and a phi/psi2 motif model. This software can identify many aspects of secondary structure, including beta-bulges. Human readable or machine readable output is available.



	usage: reassigner.v14.py [-h] [-i FILE] [-c CODE] [-o FILE] [--legend] [--verbose]
	
	Assigns secondary structure without using hydrogen bonds. One method uses
	virtual dihedrals and bond angles, the other using phi/psi2 motifs to identify
	secondary structure.
	
	optional arguments:
		-h, --help  show this help message and exit
		-i FILE     This should be a pdb file. Do not use in combination with the "-c" option.
		-c CODE     This should be a four letter pdb code. Only use this option if you want to download directly from the PDB.
		-o FILE     Output file, a table using tabs as seperators. If not specified, default is to output to STDOUT as human readable output.
		--legend    Option to print a legend for the Secondary Structure codes.
		--verbose   Option to print a all the output, including the behind the scenes methods for structure assignment.
	
	IMPORTANT NOTES: Expect there to be strange ? residues at the end of the
	output if there are any ligands or other non-water HETATOMS present. This can
	be safely ignored.
 
 
	Legend:
	
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
