'''
MIT License

Copyright (c) 2018 Aprahamian, Kim, Lindert

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Usage: predict the protein target/ligand pair from docking scores using a combined z-score

python TargetID.py -f /path/to/docking_input_file.txt (-verbose)

Input: white space seperated file with 3 columns: target, ligand, docking score

Output: list of predicted target/ligand pairs
	Optional: avg. scores, std. deviation of scores, calculated zscores

Authors: Stephanie Kim and Melanie Aprahamian
PI: Steffen Lindert
Date: August 22, 2018
Contact: lindert.1@osu.edu
'''

import numpy as np
from operator import itemgetter
import warnings
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f',type=str,help='input file containing 3 columns, white space seperated, comprised of: target, ligand, docking score')
parser.add_argument('-verbose',help='print(more in depth details',action='store_true')
args = parser.parse_args()

#Z-Score Calculator
def zscore( type, docking_results_list, unique_targets_list, unique_ligands_list):
	with warnings.catch_warnings():
		warnings.simplefilter("ignore", category=RuntimeWarning)
		avg_score_list = []
		stddev_score_list = []
		temp_list = []
		sorted_results = sorted(docking_results_list,key=itemgetter(type))
		if type == 'target':
			unique_list = set(unique_targets_list)
		if type == 'ligand':
			unique_list = set(unique_ligands_list)
		for unique in set(unique_list):
			temp_score_list = []
			for i in range(0,len(sorted_results)):
				if sorted_results[i][type] in unique:
					temp_score_list.append(sorted_results[i]['score'])
			temp_score_list = np.array(temp_score_list)
			avg_score_list.append(np.mean(temp_score_list))
			stddev_score_list.append(np.std(temp_score_list))
			temp_list.append(unique)
		return avg_score_list, stddev_score_list, temp_list

#Read in file
def file_reader( input_file ):
	docking_results_list = []
	targets = []
	ligands = []
	with open(input_file) as input:
		for line in input:
			targets.append(line.split()[0])
			ligands.append(line.split()[1])
			docking_results_dict = {}
			docking_results_dict = {
				'target': line.split()[0],
				'ligand': line.split()[1],
				'score': float(line.split()[2]),
			}
			docking_results_list.append(docking_results_dict)
	return docking_results_list, set(targets), set(ligands)

#Calculate the combined zscores
def combined_zscore( input_list, targets, ligands ):
	avg_score_target_list, stddev_score_target_list, temp_target_list = zscore('target', input_list, targets, ligands)
	avg_score_ligand_list, stddev_score_ligand_list, temp_ligand_list = zscore('ligand', input_list, targets, ligands)
	if args.verbose:
		print('Target, Average Score, Standard Deviation:')
		for i in range(len(temp_target_list)):
			print(temp_target_list[i], '%.2f'%avg_score_target_list[i], '%.2f'%stddev_score_target_list[i])
		print(' ')
		print('Ligand, Average Score, Standard Deviation:')
		for i in range(len(temp_ligand_list)):
			print(temp_ligand_list[i], '%.2f'%avg_score_ligand_list[i], '%.2f'%stddev_score_ligand_list[i])
		print(' ')
	sorted_results_by_ligand = sorted(input_list,key=itemgetter('ligand'))
	all_ligand_list = []
	for tp in ligands:
		count = 0
		per_ligand_list = []
		for j in range(0,len(sorted_results_by_ligand)):
			per_ligand_list.append(sorted_results_by_ligand[j])
			target_index = temp_target_list.index(sorted_results_by_ligand[j]['target'])
			ligand_index = temp_ligand_list.index(sorted_results_by_ligand[j]['ligand'])
			per_ligand_list[count]['zscore'] = 0.7*((per_ligand_list[count]['score'] - avg_score_target_list[target_index])/stddev_score_target_list[target_index])+0.3*((per_ligand_list[count]['score'] - avg_score_ligand_list[ligand_index])/stddev_score_ligand_list[ligand_index])
			count += 1
		all_ligand_list = all_ligand_list + per_ligand_list
	sorted_ligand_results_by_zscore = sorted(per_ligand_list,key=itemgetter('zscore'))
	return sorted_ligand_results_by_zscore


#MAIN
if args.verbose:
	print('Calculated averages and standard deviations across both the various targets and ligands:')
docked_list, targets, ligands = file_reader(args.f)
zscore_list = combined_zscore(docked_list, targets, ligands)
if args.verbose:
	print('Predicted target/ligand pairs and combined zscore:')
else:
	print('Predicted target/ligand pairs:')
for ligand in ligands:
	for i in range(0,len(zscore_list)):
		if ligand == zscore_list[i]['ligand']:
			if args.verbose:
				print(zscore_list[i]['target'], zscore_list[i]['ligand'], '%.2f'%zscore_list[i]['zscore'])
			else:
				print(zscore_list[i]['target'], zscore_list[i]['ligand'])
			break

