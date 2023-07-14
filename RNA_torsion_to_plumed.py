print("-> Remove header 2 lines and 1 footer line from .gro")
print("-> First residue should start from 1 and assumed to be continuous until the end")
print("-> Assuming just 1 LIG residue present with heavy atoms N, O, C, P")
print("-> Nucleic Acid Sugar atom number array -> listed as C1', C2', C3', C4', O4'")

import numpy as np
import pandas as pd


########### USER CHOICES#######
#-----------------------------#
save_dir = 'D:/disha/'
a = pd.read_csv(save_dir + 'test.pdb', header=None) ## processing file name
num_RNA_residues = 76


#-----------------------------#

plumed_master = []
print_master = []


def torsion():
    O3_p = []
    P = []
    O5_p = []
    C5_p = []
    C4_p = []
    C3_p = []
    
    O4_p = []
    C1_p = []
    N1 = []
    C2 = []
    N9 = []
    C4 = []

    residue_name_list = []
    previous_residue = 0
    for i in range(len(a.index)):                
        entry = a.iloc[i,0]
        current_residue = int(entry[24:27])
        if current_residue == previous_residue+1:
            residue_name_list.append(entry[19:20])
            previous_residue = current_residue
    #print(residue_name_list)
    print("Finding torsion angles ...")
    for i in range(len(a.index)):
        entry = a.iloc[i,0]
        current_residue = int(entry[24:27])
        if current_residue <= num_RNA_residues:
            if " O3' " in entry:
                O3_p.append(int(entry[7:11]))
            elif " P " in entry:
                P.append(int(entry[7:11]))
            elif " O5' " in entry:
                O5_p.append(int(entry[7:11]))
            elif " C5' " in entry:
                C5_p.append(int(entry[7:11]))
            elif " C4' " in entry:
                C4_p.append(int(entry[7:11]))
            elif " C3' " in entry:
                C3_p.append(int(entry[7:11]))
            elif " O4' " in entry:
                O4_p.append(int(entry[7:11]))
            elif " C1' " in entry:
                C1_p.append(int(entry[7:11]))
            
            if 'RC' in entry or 'RU' in entry: ## Pyrimidine
                if " N1 " in entry:
                    N1.append(int(entry[7:11]))
                elif " C2 " in entry:
                    C2.append(int(entry[7:11]))
            if 'RA' in entry or 'RG' in entry: ## Purine
                if " N9 " in entry:
                    N9.append(int(entry[7:11]))
                elif " C4 " in entry:
                    C4.append(int(entry[7:11]))
    
    chi_pyr = 0
    chi_pur = 0
    if 'RC' in residue_name_list[0] or 'RU' in residue_name_list[0]:
        chi_pyr = 1
    elif 'RA' in residue_name_list[0] or 'RG' in residue_name_list[0]:
        chi_pur = 1
        
    purine_counter = 0
    pyrimidine_counter = 0
    for i in range(num_RNA_residues-2):
        plumed_master.append('alpha_' + str(i+2) + ':  TORSION ATOMS=' + str(O3_p[i]) + ',' + str(P[i]) + ',' + str(O5_p[i+1]) + ',' + str(C5_p[i+1]) + ' NOPBC') ##alpha
        plumed_master.append('beta_' + str(i+2) + ':  TORSION ATOMS=' + str(P[i]) + ',' + str(O5_p[i+1]) + ',' + str(C5_p[i+1]) + ',' + str(C4_p[i+1]) + ' NOPBC') ##beta
        plumed_master.append('gamma_' + str(i+2) + ':  TORSION ATOMS=' + str(O5_p[i+1]) + ',' + str(C5_p[i+1]) + ',' + str(C4_p[i+1]) + ',' + str(C3_p[i+1]) + ' NOPBC') ##gamma
        plumed_master.append('delta_' + str(i+2) + ':  TORSION ATOMS=' + str(C5_p[i+1]) + ',' + str(C4_p[i+1]) + ',' + str(C3_p[i+1]) + ',' + str(O3_p[i+1]) + ' NOPBC') ##delta
        plumed_master.append('epsilon_' + str(i+2) + ':  TORSION ATOMS=' + str(C4_p[i+1]) + ',' + str(C3_p[i+1]) + ',' + str(O3_p[i+1]) + ',' + str(P[i+1]) + ' NOPBC') ##epsilon
        plumed_master.append('zeta_' + str(i+2) + ':  TORSION ATOMS=' + str(C3_p[i+1]) + ',' + str(O3_p[i+1]) + ',' + str(P[i+1]) + ',' + str(O5_p[i+2]) + ' NOPBC') ##zeta
        plumed_master.append('eta_' + str(i+2) + ':  TORSION ATOMS=' + str(C4_p[i]) + ',' + str(P[i]) + ',' + str(C4_p[i+1]) + ',' + str(P[i+1]) + ' NOPBC') ##eta
        plumed_master.append('theta_' + str(i+2) + ':  TORSION ATOMS=' + str(P[i]) + ',' + str(C4_p[i+1]) + ',' + str(P[i+1]) + ',' + str(C4_p[i+2]) + ' NOPBC') ##theta
        if residue_name_list[i+1] == 'RC' or residue_name_list[i+1] == 'RU': ## Pyrimidine
            plumed_master.append('chi_' + str(i+2) + ':  TORSION ATOMS=' + str(O4_p[i+1]) + ',' + str(C1_p[i+1]) + ',' + str(N1[pyrimidine_counter+chi_pyr]) + ',' + str(C2[pyrimidine_counter+chi_pyr]) + ' NOPBC') ##chi
            pyrimidine_counter += 1
        elif residue_name_list[i+1] == 'RA' or residue_name_list[i+1] == 'RG': ## Purine
            plumed_master.append('chi_' + str(i+2) + ':  TORSION ATOMS=' + str(O4_p[i+1]) + ',' + str(C1_p[i+1]) + ',' + str(N9[purine_counter+chi_pur]) + ',' + str(C4[purine_counter+chi_pur]) + ' NOPBC') ##chi
            purine_counter += 1
    
    plumed_master.append('\n')
    torsion = ['alpha_', 'beta_', 'gamma_', 'delta_', 'epsilon_', 'zeta_', 'eta_', 'theta_', 'chi_']
    for i in range(len(torsion)):
        temp=''
        for j in range(1,num_RNA_residues-1):
            temp = temp + torsion[i] +  str(j+1) + ','
        print_master.append("PRINT ARG=" + temp[:-1] + " STRIDE=100 FILE=plum/" + torsion[i][:-1])
    print("Torsion search complete!")
    print("-"*50)
    
torsion()

with open(save_dir + 'plumed.dat', 'w') as f:
    for line in plumed_master:
        f.write(line)
        f.write('\n')
   # f.write('c63_spec: COM ATOMS=2008,2010,2012,2016,2017,2007 NOPBC\nd_63_spec: DISTANCE ATOMS=c63_spec,clig_3\n')
   # f.write('PRINT ARG=d_63_spec   STRIDE=100 FILE=plum/d_63_spec\n')    
    for line in print_master:
        f.write(line)
        f.write('\n')