'''
This is a code to send data to github from the mcce runs
Protonation state effects in selective kinase inhibitors

'''

#!/usr/bin/python
import sys
import os
import time
import os.path
from tempfile import mkstemp
from shutil import move
from os import remove, close
from subprocess import call # This is needed to submit jobs
from argparse import ArgumentParser
from argparse import RawTextHelpFormatter

from_ = '/home/salah/johnChodera_project_on_github/mcce-charges/mcce_runs/default/'
to_ = '/home/salah/johnChodera_project_on_github/mcce-charges/mcce_runs/data/default/'

dirs = [d for d in os.listdir(from_) if os.path.isdir(os.path.join(from_, d))]



def send_what(str_1):
	for topDir in dirs:
		dirs_2 = [d for d in os.listdir(from_+'/'+str(topDir)) if os.path.isdir(os.path.join(from_+'/'+str(topDir), d))]
		makeThisDir = to_ + topDir
		if not os.path.exists(makeThisDir):
			sys_call_2 = 'mkdir ' + makeThisDir
			os.system(sys_call_2)
		for secondDir in dirs_2:
			if str_1 == 'send_together':
				if '_together' in secondDir:
					print secondDir 
					makeThisDir = to_ + topDir+'/'+secondDir
					if not os.path.exists(makeThisDir):
						sys_call_2 = 'mkdir ' + makeThisDir
						os.system(sys_call_2)
				
					getFrom = from_ +'/'+topDir+'/'+secondDir
					sendTo = to_ + topDir+'/'+secondDir
					
					pdb_file = 'cp ' +getFrom+'/good* '+ sendTo
					head3 = 'cp ' +getFrom+'/head3.lst '+ sendTo
					sum_charge = 'cp ' +getFrom+'/sum_crg.out '+sendTo
					run_prm = 'cp ' +getFrom+'/run.prm '+sendTo
					fort38 = 'cp ' +getFrom+'/fort.38 '+sendTo
							
					os.system(pdb_file)
					os.system(head3)
					os.system(sum_charge)
					os.system(run_prm)
					os.system(fort38)
					
			elif str_1 == 'send_proteinalone':
				if '_proteinalone' in secondDir:
					print secondDir 
					makeThisDir = to_ + topDir+'/'+secondDir
					if not os.path.exists(makeThisDir):
						sys_call_2 = 'mkdir ' + makeThisDir
						os.system(sys_call_2)
				
					getFrom = from_ +'/'+topDir+'/'+secondDir
					sendTo = to_ + topDir+'/'+secondDir
					
					head3 = 'cp ' +getFrom+'/head3.lst '+ sendTo
					sum_charge = 'cp ' +getFrom+'/sum_crg.out '+sendTo
					run_prm = 'cp ' +getFrom+'/run.prm '+sendTo
					fort38 = 'cp ' +getFrom+'/fort.38 '+sendTo
							
					os.system(head3)
					os.system(sum_charge)
					os.system(run_prm)
					os.system(fort38)

			elif str_1 == 'send_inhibitoralone':
				if '_inhibitoralone' in secondDir:
					print secondDir 
					makeThisDir = to_ + topDir+'/'+secondDir
					if not os.path.exists(makeThisDir):
						sys_call_2 = 'mkdir ' + makeThisDir
						os.system(sys_call_2)
				
					getFrom = from_ +'/'+topDir+'/'+secondDir
					sendTo = to_ + topDir+'/'+secondDir
					
					head3 = 'cp ' +getFrom+'/head3.lst '+ sendTo
					sum_charge = 'cp ' +getFrom+'/sum_crg.out '+sendTo
					run_prm = 'cp ' +getFrom+'/run.prm '+sendTo
					fort38 = 'cp ' +getFrom+'/fort.38 '+sendTo
							
					os.system(head3)
					os.system(sum_charge)
					os.system(run_prm)
					os.system(fort38)
			else:
				print 'Nothing happened'
			
			
			


	

def parse_args():
	"""Parse the command line arguments and perform some validation on the
	arguments
	Returns
	-------
	args : argparse.Namespace
		The namespace containing the arguments
	"""
	parser = ArgumentParser(description='''Send data to github''',formatter_class=RawTextHelpFormatter)
	required = parser.add_argument_group('required argument')
	required.add_argument('-send', '--str_1', required=True, type=str,
                          help='''Send to github''')
	args = parser.parse_args()
	return args
	
def main():
	args = parse_args()
	send_what(args.str_1)
	
def entry_point():
    main()
			
if __name__ == '__main__':
    entry_point()			
	
