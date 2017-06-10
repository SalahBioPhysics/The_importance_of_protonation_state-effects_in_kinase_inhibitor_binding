"""
Description
-----------
make_tpl_1.py is created by Salah Salah with the help of the tpl_maker.py created by Denise M Kilburg for the Gunner Lab at CUNY-City College for the
generation of ligand topology files for use in MCCE from epik output file from John Chodera
lab.  
 

Input Requirements
------------------
An epik file output *-merged.mol2 format and *-state-penalties.out. For example, Bosutinib-merged.mol2 and Bosutinib-state-penalties.out 

Example
-------
python epik2tpl.py Bosutinib

"""
#!/usr/bin/python
import sys
import re
import math
import datetime
import commands
import os
import re
import tempfile


## vdw_dict contains vdw radius of atoms as written in
## Bondi,A. van Der Waals Volumes and Radii,
## J.Phys.Chem, 63, 3, 1964 and http://periodic.lanl.gov/index.shtml
## (Los Alamos National Lab Periodic Table)

vdw_dict = {'H':1.20,
			'HN':1.20,
			'HA':1.20,
			'HB':1.20,
            'He':1.40,
            'C':1.70,
            'N':1.55,
			'B':1.92,  # Added by Salah 
            'O':1.52,
            'F':1.47,
            'NE':1.54,
            'SI':2.10,
            'P':1.80,
            'S':1.80,
            'CL':1.75,
			'Cl':1.75,
            'AR':1.88,
            'AS':1.85,
            'SE':1.90,
            'BR':1.85,
            'KR':2.02,
            'TE':2.06,
            'I':1.98,
            'XE':2.16,
            'ZN':1.39,
            'CU':1.40,
            'HG':1.55,
            'CD':1.58,
            'NI':1.63,
            'PD':1.63,
            'AU':1.66,
            'AG':1.72,
            'MG':1.73,
            'PT':1.75,
            'LI':1.82,
            'U':1.86,
            'GA':1.87,
            'PB':2.02,
            'SN':2.17,
            'NA':2.27,
            'FE':0.00,
            'K':2.75}

bond_table = {'H':'s',
			  'H.1':'s',
			  'C.3':'sp3',
			  'C.2':'sp2',
			  'C.1':'sp',
			  'C.ar':'sp2',
			  'C.cat':'unknown',
			  'N.3':'sp3',
			  'N.2':'sp2',
			  'N.1':'sp',
			  'N.ar':'sp2',
			  'N.am':'sp2',
			  'N.4':'sp3',
			  'N.pl3':'sp3',
			  'O.3':'sp3',
			  'O.2':'sp2',
			  'O.co2':'unknown',
			  'O.spc':'unknown',
			  'O.t3p':'unknown',
			  'S.3':'sp3',
			  'S.2':'sp2',
			  'S.O':'sp2',
			  'S.O2':'unknown',
			  'S.o2':'unknown',
			  'P.3':'sp3',	
			  'Cl' : 'sp',
			  'F'  : 'sp',
			  'HA' : 's',
			  'CA' : 'unknown',
			  'I'  : 'sp',
			  'Br' : 'sp',
			  }
			
bond_table_2 = {'1':'sp3',
				'2':'sp2',
				'3':'sp',
				'ar':'sp2',
				}
			
			
class Atom():
	serial = ''
	name = ''
	xyz = (0.0, 0.0, 0.0)
	SYBYL = ''
	chainID = ''
	resName = ''
	charge = 0.0
	connect = []
	atomDic = {}
	
	def __init__(self,line):
		list_of_atoms_ids = []
		if len(line) == 9:
			self.serial = int(line[0])    # the ID number of the atom at the time the file was created
			self.name =  line[1]           # the name of the atom
			self.xyz = (float(line[2]), float(line[3]), float(line[4]))      #  the x y z coordinate of the atom
			self.SYBYL = line[5]          # the SYBYL atom type for the atom (hybrid)
			self.chainID =  line[6]        # the ID number of the substructure containing the atom
			self.resName = line[7]        # the name of the substructure containing the atom
			self.charge  =  float(line[8]) # the charge associated with the atom
			self.atomDic[self.serial] = self.name
			
		if len(line) == 4 or len(line) == 5:
			list_of_atoms_ids.append(line[1])
			list_of_atoms_ids.append(line[2])
			
		self.connect = list_of_atoms_ids
		


class Tautomer():
	
	
	def __init__(self,mylist,nameofLig):
		self.mylist = mylist
		self.nameofLig = nameofLig.strip('\n')
		
		charge = 0.0
		numberofatoms = 0
		atoms_list = []
		connect_list = []
		for line in self.mylist:
			#print line
			if len(line) == 9:
				charge += float(line[8])
				numberofatoms += 1
				atoms_list.append(Atom(line))
		self.thecharge = charge
		self.numberofatoms = numberofatoms
		self.theatoms_list = atoms_list
		
		
		for line in self.mylist:
			if len(line) == 4 or len(line) == 5:
				#connect_list.append(Connect(line))
				connect_list.append(Atom(line))
		self.theconnect_list = connect_list
	
		
	
def write_comment_header(tpl,inhibitor_name):
	tpl.write('####################################\n')
	tpl.write('# Topology File for:\n')
	tpl.write('# {}\n'.format(tpl.name))
	tpl.write('# {}\n'.format(inhibitor_name))
	tpl.write('#\n')
	tpl.write('# Created on: {}\n'.format(datetime.date.today()))
	tpl.write('#\n')
	tpl.write('# Created with: make_tpl_1.py by Salah Salah\n')
	tpl.write('####################################\n')
	tpl.write('\n')
	return
	
def write_conformers(tpl,myTautomer_list_obj):
	conf_list = []
	crg_list = []
	for everyObj in myTautomer_list_obj:
		crg_list.append(int(round(everyObj.thecharge)))
	
	index_two = 0
	index_zero = 0
	index_one = 0
	index_three = 0
	for everyObj in myTautomer_list_obj:
		
		charge_ = int(round(everyObj.thecharge))
		if charge_ > 0:
			sign_ = '+'
		elif charge_ < 0:
			sign_ = '-'
		
		if charge_ == 0:
			index_zero += 1
			conf_list.append(everyObj.nameofLig + str(charge_)+ str(index_zero))
		elif abs(charge_) == 1:
			index_one += 1
			conf_list.append(everyObj.nameofLig + sign_+ str(index_one))
		elif abs(charge_) == 2:
			conf_list.append(everyObj.nameofLig + sign_+ chr(ord('a')+index_two))
			index_two += 1
		elif abs(charge_) == 3:
			conf_list.append(everyObj.nameofLig + sign_+ chr(ord('A')+index_three))
			index_three += 1
	
	tpl.write('# neural always starts with 0\n')
	tpl.write('# numberical value means the charge is 1\n')
	tpl.write('# alphabet lower case means the charge is 2\n')
	tpl.write('# alphabet upper case means the charge is 3\n')
	tpl.write('# nothing for charge of 4, this code will not work\n\n')
	tpl.write('CONFLIST {}        {}BK '.format(conf_list[0][:3],conf_list[0][:3]))
	for conf in conf_list:
		tpl.write('{} '.format(conf))
		
	tpl.write('{} '.format(conf_list[0][:3] + "DM"))
	tpl.write('\n')
	tpl.write('\n')
	return conf_list
	
def write_natom(tpl,myTautomer_list_obj,conf_list):
	#conf_list = write_conformers(tpl,myTautomer_list_obj)
	#numberofatoms = [1]
	i = 0
	template = '{0:9}{1:11}{2:1}\n'
	tpl.write(template.format('NATOM',conf_list[0][:3]+"BK", '0'))
	for conf in conf_list:
		tpl.write(template.format('NATOM',
										conf, myTautomer_list_obj[i].numberofatoms))
		i += 1
	tpl.write(template.format('NATOM',conf_list[0][:3]+"DM", '0'))
	tpl.write('\n')
	return

def write_iatom(tpl,myTautomer_list_obj,conf_list):
	
	template = '{0:9s}{1:6s}{2:>5s}{3:5d}\n'
	i = 0
	for conf in myTautomer_list_obj:
		i += 1
		count = 0
		for atom in conf.theatoms_list:
			tpl.write(template.format('IATOM', conf_list[i-1],
									atom.name, count))
			count += 1
		tpl.write('\n')
	return
	
def write_ATOMNAME(tpl,myTautomer_list_obj,conf_list):
	template = '{0:9s}{1:6s}{2:5d}{3:>4s}\n'
	i = 0
	for conf in myTautomer_list_obj:
		i += 1
		count = 0
		for atom in conf.theatoms_list:
			tpl.write(template.format('ATOMNAME', conf_list[i-1], \
                                  count,atom.name))
			count += 1
		tpl.write('\n')

	return

def write_basic_info(tpl):
	tpl.write('#1.Basic Conformer Information: name, pka, em, rxn.\n')
	tpl.write('#23456789A123456789B123456789C\n')
	return
	
def write_proton(tpl,myTautomer_list_obj,conf_list):

	
	tpl.write('\n# PROTON SECTION: PROTON means charge\n\n')
	template = '{0:9}{1:11}{2:5}\n'
	i = 0
	for conf in myTautomer_list_obj:
		tpl.write(template.format("PROTON",conf_list[i], \
                                      '0'))
    
		i += 1
	tpl.write(template.format("PROTON",conf_list[0][:3]+"DM", '0'))
	tpl.write('\n')
	
	return
	
def write_pka(tpl,myTautomer_list_obj,conf_list):
	
	tpl.write('# Solution pKa Section: pKa data from CRC Handbook of Chemistry and Physics\n')
	tpl.write('# pka is set to zero\n')
	template = '{0:9s}{1:6s}     {2:8.3f}\n'
	
	i = 0
	for conf in myTautomer_list_obj:		
		tpl.write(template.format("PKA", conf_list[i], 0.0))
		i += 1
	tpl.write(template.format("PKA", conf_list[0][:3]+"DM", 0.0))
	tpl.write('\n')
	
	return

def write_electron(tpl,myTautomer_list_obj,conf_list):
	
	tpl.write("#ELECTRON SECTION:\n")
	template = '{0:9}{1:11}{2:5}\n'
 
	for conf in conf_list:
		tpl.write(template.format("ELECTRON",conf,"0.0"))
	tpl.write(template.format("ELECTRON",conf_list[0][:3]+"DM","0.0"))
	tpl.write('\n')
	return
	
def write_EM(tpl,myTautomer_list_obj,conf_list): 
	
	template = '{0:9}{1:11}{2:5}\n'
	tpl.write("# EM SECTION:\n")
	
	for conf in conf_list:
		tpl.write(template.format("EM",conf,"0.0"))
	tpl.write(template.format("EM",conf_list[0][:3]+"DM","0.0"))
	tpl.write('\n')
	return
	
def write_RXN(tpl,myTautomer_list_obj,conf_list): 
	template = '{0:9}{1:11}{2:5}\n'
	tpl.write("# REACTION FIELD ENERGY SECTION:\n")
	
	for conf in conf_list:
		tpl.write(template.format("RXN",conf,"0.0"))
	tpl.write('\n')
	return
	
def write_con_header(tpl): 
    string = "ires" + " " + "conn "
    tpl.write("#ONNECT" + "   " + "conf" + " " + "atom" + "  " + "orbital" \
                  + "  " +string*4 + "\n")
    tpl.write("#ONNECT" +" " + "|-----|----|---------|" + "----|"*10 + "\n")
    return	

	
def write_con_section(tpl,myTautomer_list_obj,conf_list):
	i = 0

	for conf in conf_list:
		tpl.write("#  " + conf + "\n")
		write_con_header(tpl)
		template1 = '{0:9}{1:5} {2:^4} {3:^10}'
		template2 = '{0:^4} {1:^4} '
		for atom in myTautomer_list_obj[i].theatoms_list:
			tpl.write(template1.format("CONNECT",conf,atom.name,bond_table[atom.SYBYL]))
			for blahhh in myTautomer_list_obj[i].theconnect_list:
				if atom.serial == int(blahhh.connect[0]):
					#print 'atom ' + str(atom.serial) + ' is connected to atom ' + 
					tpl.write(template2.format("0", atom.atomDic[int(blahhh.connect[1])]))
				elif atom.serial == int(blahhh.connect[1]):
					tpl.write(template2.format("0", atom.atomDic[int(blahhh.connect[0])]))
			tpl.write('\n')
		i += 1
		tpl.write('\n')
	return

def write_atom_param_section(tpl,myTautomer_list_obj,conf_list,vdw_dict):
	
	dummyList = []
	for conf in myTautomer_list_obj:
		dummyList.append(conf.numberofatoms)
	
	conf_with_highest_numberofatoms = dummyList.index(max(dummyList))
	
	#a = len(atoms_list)/len(protonation_states.ligand_list)
	tpl.write("# Atom Parameters:\n")
	tpl.write("# Van Der Waals Radii. See source for reference\n")
	template = '{0:9}{1:7}{2:5} {3:7}\n'
	listof_radii = []
	for atom in myTautomer_list_obj[conf_with_highest_numberofatoms].theatoms_list:
		element = ''.join([i for i in atom.name if not i.isdigit()])[:1]
		tpl.write(template.format("RADIUS",conf_list[0][:3], \
								  atom.name,vdw_dict[element.upper()]))
			
	
	tpl.write('\n')
	return
	
def write_charge(tpl,myTautomer_list_obj,conf_list):
	
	template = '{0:9}{1:7}{2:3} {3:7}\n'
	i = 0
	for conf in conf_list:
		for atom in myTautomer_list_obj[i].theatoms_list:
			tpl.write(template.format("CHARGE",conf,atom.name,atom.charge))
		tpl.write('\n')
		i += 1
	return
	
def write_extra(tpl,conf_list,epik_State_Penalty):
	
	template = '{0:9s}{1:6s}{2:5s}{3:8.3f}\n'
	tpl.write("# EXTRA energy for tautomers:\n")
	i = 0
	for state_penalty in epik_State_Penalty:
		tpl.write(template.format("EXTRA", conf_list[i], "", float(state_penalty)))
		i += 1
	tpl.write('\n')
	return
	
def write_tpl(tpl,fname,epik_State_Penalty_file,inhibitor_name):
	#Save epik_State_Penalty into a list
	i = 1
	mylist = []
	myTautomer_list_obj = []
	linesthefile = open(fname).readlines()
	nameofLig = linesthefile[1]
	for line_ in linesthefile:
		line_list = line_.split()
		if len(line_list) != 0:	
			if len(line_list) > 3: #line_list[0] != '@<TRIPOS>ATOM':
				mylist.append(line_list)
			if line_list[0] == '@<TRIPOS>SUBSTRUCTURE':
				myTautomer_list_obj.append(Tautomer(mylist,nameofLig))
				#Tautomer(mylist).getTautomer()
				del mylist[:]

						
	# Write header into the tpl file
	write_comment_header(tpl,inhibitor_name)
	
	# Write CONFLIST section
	conf_list = write_conformers(tpl,myTautomer_list_obj)
				
	# Write block describing number of atoms present in each conformer
	write_natom(tpl,myTautomer_list_obj,conf_list)
	
	# Write block describing name of atoms present in each conformer.
	write_iatom(tpl,myTautomer_list_obj,conf_list)
	
	write_ATOMNAME(tpl,myTautomer_list_obj,conf_list)
	
	write_basic_info(tpl)
	
	write_proton(tpl,myTautomer_list_obj,conf_list)
	
	write_pka(tpl,myTautomer_list_obj,conf_list)
	
	write_electron(tpl,myTautomer_list_obj,conf_list)
	
	write_EM(tpl,myTautomer_list_obj,conf_list)
	
	write_RXN(tpl,myTautomer_list_obj,conf_list)
	
	write_con_section(tpl,myTautomer_list_obj,conf_list)
	
	write_atom_param_section(tpl,myTautomer_list_obj,conf_list,vdw_dict)
	
	write_charge(tpl,myTautomer_list_obj,conf_list)
	
	
	#Save epik_State_Penalty into a list
	with open(epik_State_Penalty_file, 'r') as f:
		epik_State_Penalty = [line.strip() for line in f]
	
	
	write_extra(tpl,conf_list,epik_State_Penalty)
	
	
	'''for conf in conf_list:
		print str(conf)
	print 'are created'
	'''
	
		
	

def main(fname_base,fname, epik_State_Penalty_file): 
	
	lines = [line .strip() for line in open(fname).readlines()]
	removeExt = fname_base.replace(".mol2", "")
	three_characters = removeExt[:3]
	paramDir = "param"
	if not os.path.exists(paramDir):
		os.makedirs(paramDir)
	#print 
	
	if len(lines[1]) > 3:
		lines[1] = three_characters.upper()
		
	inhibitor_name = removeExt
	print 'conf ' + inhibitor_name
	with open(os.path.join(paramDir, lines[1]+'.tpl'),'w') as tpl:
		write_tpl(tpl,fname,epik_State_Penalty_file,inhibitor_name)


	
			
if __name__ == "__main__":
	if len(sys.argv) < 2:
		print "Specify an input\nFor example python epik2tpl.py Bosutinib"
		sys.exit()
	fname_base = sys.argv[1]
	
	
	fname_merged = fname_base+'-merged.mol2'
	fname_state_penalties =  fname_base+'-state-penalties.out'
	
	main(fname_base,fname_merged,fname_state_penalties)
	
	
	
	
	
	
