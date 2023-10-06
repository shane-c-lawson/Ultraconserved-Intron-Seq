	
import os
import collections
import statistics
import subprocess

best_scores = {}
c=0
with open('500_upstream.fa', 'r') as infile1 :
	for line1 in infile1 :
		if '>' in line1 :
			gene = line1[1:].split(':')[0]
			if '_' in gene :
				gene = gene.split('_')[0]
			outtemp = open('temp.fa', 'a+')
			outtemp.write(line1)
			outtemp.close()
			continue
		outtemp = open('temp.fa', 'a+')
		outtemp.write(line1)
		outtemp.close()
		seq = line1.strip()
		cur_best_scores = {}
		self_found = False
		os.system("python3 /rbpamp/bin/RBPamp_eval.py --realm=hg38 --fasta temp.fa /rbpamp/bin/motifs_no_Ns/*.txt > out.txt")
		rbp_list = []
		err = False
		if c > 5 :
			break
		with open('out.txt', 'r') as o :
			for line2 in o :
				d = 0
				line2 = line2.strip()
				lines2 = line2.split('\t')
				if err == True :
					break
				if 'name' in line2 :
					for RBP in lines2 :
						rbp_list.append(RBP)
				else :
					for scores in lines2 :
						d += 1
						if '.' in scores and d > 2 :
							avg = float(scores)
							check = d - 1
							rbp1 = rbp_list[check]
							if '-7' not in rbp1 and 'eclip' not in rbp1 and 'rbpmap' not in rbp1 :
								continue
							if '-7' in rbp1 :
								rbp = rbp1.split('-7')[0]
							else :
								rbp = rbp1.split('_')[0]
							if "SRSF" in rbp or 'TIA' in rbp :
								if "SR" not in cur_best_scores :
									cur_best_scores["SR"] = avg
								else :
									cur_avg = cur_best_scores["SR"]
									if avg > cur_avg :
										cur_best_scores["SR"] = avg
							if "HNRNP" in rbp or 'PTBP' in rbp :
								if "HNRNP" not in cur_best_scores :
									cur_best_scores["HNRNP"] = avg
								else :
									cur_avg = cur_best_scores["HNRNP"]
									if avg > cur_avg :
										cur_best_scores["HNRNP"] = avg
							if gene == rbp :
								self_found = True
								if "Self" not in cur_best_scores :
									cur_best_scores["Self"] = avg
								else :
									cur_avg = cur_best_scores["Self"]
									if avg > cur_avg :
										cur_best_scores["Self"] = avg
							if rbp not in cur_best_scores :
								cur_best_scores[rbp] = avg
							else :
								cur_avg = cur_best_scores[rbp]
								if avg > cur_avg :
									cur_best_scores[rbp] = avg
						elif d > 2 :
							print(scores)
							err = True
							break
				anytot=0
				for rbpz in cur_best_scores :
					anytot += cur_best_scores[rbpz]
				cur_best_scores['Total Motifs'] = round(float(anytot), 5)
				for rbps in cur_best_scores :
					if rbps not in best_scores :
						best_scores[rbps] = [cur_best_scores[rbps]]
					else :
						cur_bests = best_scores[rbps]
						cur_bests.append(cur_best_scores[rbps])
						best_scores[rbps] = cur_bests
		os.remove('out.txt')
		os.remove('temp.fa')

rbp_avgsi3 = {}
sort_contsany = {}
sort_totsany = {}
sort_resany = {}
for rbpi in best_scores :
	all_avg = []
	tots = best_scores[rbpi]
	if len(tots) == 0 :
		continue
	count = 0
	for tot in tots :
		count += tot
		#if tot < 1 :
		#	continue
		if rbpi != 'SR' and rbpi != 'Self' and rbpi != 'HNRNP' and rbpi != 'Total Motifs' :
			all_avg.append(tot)
		else :
			Selfs = open('grouped_results.txt', 'a+')
			Selfs.write('Conserved NMD\t' + rbpi + '\t' + str(tot) + '\n')
			Selfs.close()
		if rbpi == 'UPF1' :
			Selfs = open('grouped_results.txt', 'a+')
			Selfs.write('Conserved NMD\t' + rbpi + '\t' + str(tot) + '\n')
			Selfs.close()
	avg = float(count / float(len(tots)))
	rbp_avgsi3[rbpi] = round(avg, 5)
	res = round(statistics.pstdev(tots), 5)
	sort_resany[rbpi] = res
	sort_contsany[rbpi] = count
	sort_totsany[rbpi] = len(tots)
sort_avgsi31 = sorted(rbp_avgsi3.items(), key=lambda x: x[1], reverse=True)
sort_avgsi3 = collections.OrderedDict(sort_avgsi31)
Selfs = open('ordered_results.txt', 'w')
for rbps in sort_avgsi3 :
	Selfs.write(rbps + '\t' + str(sort_avgsi3[rbps]) + '\t' + str(sort_contsany[rbps]) + '\t' + str(sort_totsany[rbps]) + '\t' + str(sort_resany[rbps]) + '\n')
Selfs.close()
