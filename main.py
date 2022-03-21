#!/usr/bin/env python
# -*- coding: utf-8 -*-

## python modules ##
####################
import numpy as np
import os
######################################################################
## Title: LibrePhsysics VASP Freature extraction (LPVFE)            ##
## Description: Program to colect some stufs on VASP Program        ##
## Version: 0.8                                                     ##
## author: Natan Moreira Regis                                      ##
## Project Name: LibrePhysics (librephysics.xyz)                    ##
## license: GPLv3                                                   ##
## maintainer: "Natan Moreira Regis"                                ##
## email: "n.m.regis@df.ufscar.br"                                  ##
## status: "Production"                                             ##
######################################################################


###########################
## Using for latex table ##
###########################
def table_header_props():
    k=0
    columns=[]
    for element in props[0]:
        if element in prop:
            columns.append(props[1][k]) 
        k+=1
    return columns

################################
## open File and Return lines ##
################################
def openfile(filename):
    file = open(filename, 'r')
    lines = file.readlines()
    return lines
    file.close()

################################
## Check if string is in file ##
################################
def CheckStringInFile(File,String):
    datafile=openfile("OUTCAR")
    for line in datafile:
        if String in line:
            return True
    return False

###########################
## Search string in list ##
###########################
def search_in_list(lista, item):
    try:
        return lista.index(item)
    except ValueError:
        return -1

####################################
## Energies of molecule and layer ##
####################################
def EnergyTotMoleLayer():
    Molecules = [['CH4',-24.03376510], ['CO',-14.90940997], ['CO2',-23.17457153],["H2",-6.75664087],['H2O',-14.26411655],['N2',-16.77991136],['N2O',-21.57599354],['NH3',-19.54729190],['NO',-12.39971547],['NO2',-18.58681627],['O2',-9.99926374],['SO2',-17.01174353]]
    Monolayer = [['Pristine 5x5',-577.66298749], ['Defective 5x5',-570.60150797]]# Old: -577.66300013 e  -570.60142894  
    return Molecules, Monolayer
  
##################
## Total Energy ##
##################
def E_tot():
    lines=openfile("OUTCAR")
    LinesEtot=[]
    for line in lines:
        if "energy  w" in line:
            LinesEtot.append(line)
    FinalEtot=LinesEtot[-1]    
    EEotEWntropy = float(FinalEtot.split()[-1])
    return EEotEWntropy
    
###################################
## Adsorption energy calculation ##
###################################
def E_ads(Mol):
    EnergyAdsorbed = E_tot()
    Molecules, Monolayer = EnergyTotMoleLayer()
    for Molecule in Molecules:
        if Mol==Molecule[0]:
            MolEnergy=Molecule[1]
    if S_defects==0:
        LayerEnergy = Monolayer[0][1]
    else:
        LayerEnergy = Monolayer[1][1]
    AdsorptionEnergy = float(EnergyAdsorbed) - float(MolEnergy) - float(LayerEnergy)
    return '{:.4f}'.format(round(AdsorptionEnergy,4))    

#####################
## gap gamma point ##
#####################
def gap_gamma():
    return 0

#####################    
## relative energy ##  
#####################
def E_rel(*args):
    for el in args:
        if type(el)!= list:
            pass
            return None
        else:
            mvalue=min(el) 
            for i in range(len(el)):
                el[i]='{:.4f}'.format(round(float(el[i])-float(mvalue),4))
            return el 

#####################################
## Geometric parameters extraction ##
#####################################
def Dist_two_vectors(v1,v2):
        return np.linalg.norm(v2-v1)

def lowest_dist_ignoring_0(Dist):
    Dist = sorted(Dist)
    Dist = min(numero for numero in Dist if numero != 0) 
    return Dist

def dist_mol_to_atom(Atom):# Using mainly for dist Mol to Mo and S
    positions = atoms_info()# Positions of all atoms
    Atoms, TotAtoms = species_info()# Atomic atoms and total number atoms 
    vMed=center_mol()# Vector of the center of the molecule
    IndexAtom=int(search_in_list(Atoms[0], Atom))
    LineInitiAtom=0# Line in the positions of begining of Atom passing in the argument 
    i=0
    while i<IndexAtom:
        LineInitiAtom+=int(Atoms[1][i])# Sum the precedent atoms
        i+=1# Loop in the list of atoms
    Dist = []
    for i in range(LineInitiAtom,int(Atoms[1][IndexAtom])):
        x1 = float(positions[i,0])
        y1 = float(positions[i,1])
        z1 = float(positions[i,2])
        if z1>10.0:# Correct the standart positions using by vasp, ie, atoms in the layers images
            z1-= float(basis_info()[2,2])
        v1 = np.matrix([x1, y1, z1])
        Dist.append(Dist_two_vectors(v1,vMed))
    return '{.:2f}'.format(round(lowest_dist_ignoring_0(Dist),2))
    
def center_mol_dir_z():
    positions = atoms_info()
    Atoms, TotAtoms = species_info()
    IndexMo=int(search_in_list(Atoms[0], "Mo"))
    IndexS=int(search_in_list(Atoms[0], "S"))
    IndexInit=int(Atoms[1][IndexMo]) + int(Atoms[1][IndexS])
    p = 0
    u = 0
    for i in range(IndexInit,len(positions)):
       u += positions[i,2]
       i +=1
       p +=1
       uMed = u/float(p)
       vMed=[0,0,uMed]
       return vMed

def dist_mol_to_atom_dir_z(Atom):
    positions = atoms_info()
    Atoms, TotAtoms = species_info()
    vMed=center_mol_dir_z()
    IndexAtom=int(search_in_list(Atoms[0], Atom))
    LineInitiAtom=0
    i=0
    while i<IndexAtom:
        LineInitiAtom+=int(Atoms[1][i])
        i+=1
    Dist = []
    for i in range(LineInitiAtom,int(Atoms[1][IndexAtom])):
        z1 = float(positions[i,2])
        if z1>6.0:
            z1-= float(basis_info()[2,2])
        v1 = np.matrix([0, 0, z1])
        Dist.append(Dist_two_vectors(v1,vMed))
    return '{:.2f}'.format(round(min(numero for numero in Dist if numero != 0),2))
    
def distance_atoms_mol():
    CONTCAR = openfile('CONTCAR')
    positions = atoms_info()
    atoms, tot_atoms = species_info()
    list_atoms=[]
    count=0
    dist_atoms=[]
    for a in range(len(atoms[0])):
        molecule=atoms[0][a]
        for b in range(int(atoms[1][a])):
            list_atoms.append(molecule)
            count+=1
    for i in range(2,len(atoms[1])):
        if len(atoms[0])==3:
            suplimit=int(int(atoms[1][i]))
        else:
            suplimit=int(int(atoms[1][i])+1)

        for k in range(0,suplimit):
            k=75+k-S_defects
            x1 = float(positions[k,0])
            y1 = float(positions[k,1])
            z1 = float(positions[k,2])
            v1 = np.matrix([x1,y1,z1])
            atom1=list_atoms[k]
            for l in range(0,suplimit):
                l=75+l-S_defects
                x2 = float(positions[l,0])
                y2 = float(positions[l,1])
                z2 = float(positions[l,2])
                v2 = np.matrix([x2,y2,z2])
                atom2=list_atoms[l]
                Dist = Dist_two_vectors(v1,v2)
                if k==l:
                    pass
                else:
                    mol_dist=atom1+'-'+atom2
                    dist_atoms.append([mol_dist,Dist])
    uniques=[]
    final_dist_atoms=[]
    for item in dist_atoms:
        if round(item[1],2) not in uniques:
            uniques.append(round(item[1],2))
            final_dist_atoms.append(item)
    tmp=[]
    for item in final_dist_atoms:
        item[1] = round(item[1],2)
        tmp.append([item[0],'{:.2f}'.format(item[1])])
    final_dist_atoms=tmp
    for i in range(len(final_dist_atoms)):
        final_dist_atoms[i]=str(": ".join(final_dist_atoms[i]))
    final_dist_atoms=' \\\ '.join(final_dist_atoms)
    final_dist_atoms='\\begin{{tabular}}{{c}} {} \\end{{tabular}}'.format(final_dist_atoms)
    return final_dist_atoms

def angle(Mol):
    Positions = atoms_info()
    a=0
    b=0
    c=0
    if len(Positions) == 77-S_defects:
        return '{:.1f}'.format(round(0.0,1))
    else:
        if Mol=='N2O':
            i = 75-S_defects
            x1 = float(Positions[i,0])
            y1 = float(Positions[i,1])
            z1 = float(Positions[i,2])

            x3 = float(Positions[i+2,0])
            y3 = float(Positions[i+2,1])
            z3 = float(Positions[i+2,2])

            x2 = float(Positions[i+1,0])
            y2 = float(Positions[i+1,1])
            z2 = float(Positions[i+1,2])
            b = np.array([x1,y1,z1])
            a = np.array([x2,y2,z2])
            c = np.array([x3,y3,z3])

        elif Mol=='CO2' or Mol=='H2O' or Mol=='NO2':
            i=75-S_defects
            x1 = float(Positions[i,0])
            y1 = float(Positions[i,1])
            z1 = float(Positions[i,2])

            x3 = float(Positions[i+2,0])
            y3 = float(Positions[i+2,1])
            z3 = float(Positions[i+2,2])

            x2 = float(Positions[i+1,0])
            y2 = float(Positions[i+1,1])
            z2 = float(Positions[i+1,2])
            b = np.array([x3,y3,z3])
            a = np.array([x2,y2,z2])
            c = np.array([x1,y1,z1])

        elif Mol=='SO2':
            i=75-S_defects
            x1 = float(Positions[i,0])
            y1 = float(Positions[i,1])
            z1 = float(Positions[i,2])

            x3 = float(Positions[i+2,0])
            y3 = float(Positions[i+2,1])
            z3 = float(Positions[i+2,2])

            x2 = float(Positions[i+1,0])
            y2 = float(Positions[i+1,1])
            z2 = float(Positions[i+1,2])
            b = np.array([x1,y1,z1])
            a = np.array([x2,y2,z2])
            c = np.array([x3,y3,z3])
 
        elif Mol=='CH4':
            i = 75-S_defects
            x1 = float(Positions[i,0])
            y1 = float(Positions[i,1])
            z1 = float(Positions[i,2])
            x3 = float(Positions[i+1,0])
            y3 = float(Positions[i+1,1])
            z3 = float(Positions[i+1,2])
            x2 = float(Positions[-1,0])
            y2 = float(Positions[-1,1])
            z2 = float(Positions[-1,2])
            b = np.array([x2,y2,z2])
            a = np.array([x1,y1,z1])
            c = np.array([x3,y3,z3])

        elif Mol=='NH3':
            i = 75-S_defects
            x1 = float(Positions[i,0])
            y1 = float(Positions[i,1])
            z1 = float(Positions[i,2])
            x3 = float(Positions[i+1,0])
            y3 = float(Positions[i+1,1])
            z3 = float(Positions[i+1,2])
            x2 = float(Positions[i+2,0])
            y2 = float(Positions[i+2,1])
            z2 = float(Positions[i+2,2])
            b = np.array([x1,y1,z1])
            a = np.array([x2,y2,z2])
            c = np.array([x3,y3,z3])
    ba = a - b
    bc = c - b
    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)
    angle = np.degrees(angle)
    return '{:.1f}'.format(round(angle,1))

def CONTCAR_figs():
    lines=openfile("CONTCAR")
    Species, TotAtoms = AtomsInfo.species_info()
    Positions=[]    
    if lines[7].split()[0]=="Direct" or lines[7].split()[0]=="Cartesian":   
        cabecalho=lines[:8]
        for i in range(8,TotAtoms+8):
            v = lines[i].split()
            x = float(v[0])
            y = float(v[1])
            z = float(v[2])
            if z>0.7:
                z-=1. #float(AtomsInfo.basis_info()[2,2])
                Positions.append('{}  {}  {} \n'.format(x,y,z))
            else:
                cabecalho=lines[:9]
            for j in range(9,TotAtoms+9):
                v = lines[j].split()
                x = float(v[0])
                y = float(v[1])
                z = float(v[2])
                if z>0.7:
                    z-= 1.#float(AtomsInfo.basis_info()[2,2])
                Positions.append('{}  {}  {} \n'.format(x,y,z))
       
    file_out=open('CONTCAR_fig','w')
    file_out.writelines(cabecalho)
    file_out.writelines(Positions)
    file_out.close()

def run_fig(Mol,Path, Root):
        figssh = str(Root)+'/'+'figs.sh'
        os.system('cp {} ./figs.sh'.format(figssh)) 
        os.system('chmod +x figs.sh')
        os.system('sh figs.sh')

        os.system('jmol -x top.spt && mv top.png /run/media/natan/home-natan/estudos/faculdade/IC/MoS2/results/allfigures/{} '.format(Mol+'_'+Path+'_top.png'))
        os.system('jmol -x side.spt && mv side.png /run/media/natan/home-natan/estudos/faculdade/IC/MoS2/results/allfigures/{} '.format(Mol+'_'+Path+'_side.png'))

################
## basis info ##
################
def basis_info():
    basis=[]
    lines = openfile("CONTCAR")
    for i in range(2,5,1):
        line = lines[i]
        line=line.split()
        xb = float(line[0])
        yb = float(line[1])
        zb = float(line[2])
        basis.append([xb,yb,zb])
    basis=np.matrix(basis)
    return basis
 
#########################   
## Atomic species info ##
#########################
def species_info():
    lines = openfile("CONTCAR")
    LineAtoms=lines[5].split()
    LineNumAtoms=lines[6].split()
    Atoms=[LineAtoms,LineNumAtoms]
    TotAtoms=0
    for j in range(len(LineNumAtoms)):
        TotAtoms+=int(LineNumAtoms[j])
    return Atoms, TotAtoms
    
######################
## Atomic positions ##
######################
def atoms_info():
    basis = basis_info()
    lines = openfile("CONTCAR")
    Species, TotAtoms = species_info()
    positions=[]   
    if lines[7].split()[0]=="Direct" or lines[7].split()[0]=="Cartesian":
        for i in range(8,TotAtoms+8):
            v = lines[i].split()
            x = float(v[0])
            y = float(v[1])
            z = float(v[2])
            positions.append([x,y,z])
    else:
        for j in range(9,TotAtoms+9):
            v = lines[j].split()
            x = float(v[0])
            y = float(v[1])
            z = float(v[2])
        positions.append([x,y,z])
    positions = np.array(positions)
    positions= np.matrix(positions)*np.matrix(basis)
    return positions

##############################
## Center adsorbed molecule ##
##############################
def center_mol():
    positions = atoms_info()
    Atoms, TotAtoms = species_info()
    IndexMo=int(search_in_list(Atoms[0], "Mo"))
    IndexS=int(search_in_list(Atoms[0], "S"))
    IndexInit=int(Atoms[1][IndexMo]) + int(Atoms[1][IndexS])
    p = 0; x = 0; y = 0; z = 0
    for i in range(IndexInit,len(positions)):
        x += positions[i,0]
        y += positions[i,1]
        z += positions[i,2]
        i +=1
        p +=1
    xMed = x/float(p)
    yMed = y/float(p)
    zMed = z/float(p)
    vMed=[xMed,yMed,zMed]
    return vMed

######################
## Latice Parameter ##
######################
def lattice_vectors():
    basis=basis_info()
    na, nb, nc = type_supercel 
    a = np.sqrt(float(basis[0,0])**2 + float(basis[0,1])**2 + float(basis[0,2])**2)
    b = np.sqrt(float(basis[1,0])**2 + float(basis[1,1])**2 + float(basis[1,2])**2)
    c = np.sqrt(float(basis[2,0])**2 + float(basis[2,1])**2 + float(basis[2,2])**2)
    
    if a/na == b/nb: 
        return a/na
    else: 
        return "erro"

def configuration(folder):
    Path = folder.replace('./','')
    Path = Path.replace('/','-')
    return Path

def colect_results(Mol, folder):
    files = [x[2] for x in os.walk('.')][0]
    if 'OUTCAR' in files:
        Check = CheckStringInFile('OUTCAR','reached required accuracy - stopping structural energy minimisation')
        if Check ==True:
            print(Mol,folder)
            aux=[]
            for i in prop:
                aux.append(eval(i))
            return aux
        else:
            return 0

def saida_tex(lista_tot):
    columns = table_header_props()
    file = open('/run/media/natan/home-natan/estudos/faculdade/IC/MoS2/results/saida.tex','w')
    file.write("\\setlength{{\\tabcolsep}}{{3pt}} \n \\footnotesize \n \\begin{{longtable}}{{|{}|}} \n \\hline \n" .format("c"*len(columns)))
    file.write('\\endfirsthead')
    file.write('\\multicolumn{{{}}}{{c}}'.format(len(columns)))
    file.write('{{\\bfseries \\tablename\ \\thetable{} -- continued from previous page}} \\\ \\hline \n \\endhead')
    file.write('\\hline \n \\multicolumn{{{}}}{{|r|}}{{Continued on next page}} \\\ \\hline \\endfoot'.format(len(columns)))
    file.write('\\hline  \\hline  \\endlastfoot \n') 
    file.write('\hline \n'+'\n'+' & '.join(columns)+' \\\ \\hline \n')
    saida=[]
    for lista in lista_tot: 
        saida.append('\\hline \\hline \n \multicolumn{{{}}}{{|c|}}{{{}}} \\\ \\hline \\hline \n'.format(len(prop),'pristine'))
        for i in range(len(lista)):
            for j in range(len(lista[i])):
                if len(lista[i][j])>3:
                    a = ' & '.join(map(str,lista[i][j]))
                    a = a +' \\\ \n'
                else:
                    a='\\hline \\hline \n \multicolumn{{{}}}{{|c|}}{{{}}} \\\ \\hline \\hline \n'.format(len(prop),lista[i][j])
                saida.append(a)
    file.writelines(saida)
    file.write('\\hline  \n \end{longtable}')
    file.close()

def navigate_dir(Molecules,root):
    print(root)
    out=[]
    os.chdir(root)
    if (recursive=='yes'):
        for Mol in Molecules:
            print(Mol)
            etot=[]
            saida=[]
            saida.append(Mol)
            os.chdir(Mol)
            Root = os.getcwd()
            folders = [x[0] for x in os.walk('.')]
            for folder in folders:
                os.chdir(folder)
                aux=colect_results(Mol, folder)
                if aux==None:
                    pass
                else:
                    if 'E_rel()' in prop:
                        index1=prop.index('E_rel()') 
                        index2=props[0].index('E_tot()')
                    saida.append(aux)  
                    etot.append(aux[index2])
                os.chdir(Root)
            os.chdir(root)
            erel=[]
            erel = E_rel(etot)
            k=0
            for i in range(len(saida)):
                if type(saida[i])==str:
                    pass
                else:
                    saida[i][index1] = erel[k]
                    k=k+1
            out.append(saida)
    return out


###################
## Main function ##
###################
def main():
    global Molecules, recursive, S_defects, type_supercel, prop,  props, figures

    ######################
    ## general Options  ##
    ######################
    
    type_supercel=[5,5,1]# Typer of supercell. Will be uset to compute the total number of atoms 
    prop=['configuration(folder)','E_tot()','E_ads(Mol)','E_rel()','dist_mol_to_atom_dir_z("S")','distance_atoms_mol()', 'angle(Mol)']
    props=[['configuration(folder)','E_tot()','E_ads(Mol)', 'E_rel()', 'gap_gamma()' ,'dist_mol_to_atom_dir_z("S")', 'dist_mol_to_atom("Mo")', 'dist_mol_to_atom("S")', 'distance_atoms_mol()', 'angle(Mol)'],['Configurations','$E_{Tot}$ (\\si{\\electronvolt})', '$E_{Ads}$ (\\si{\\electronvolt})', '$E_{Rel}$ (\\si{\\electronvolt})', '$E_{gap-\\Gamma}$ (\\si{\\electronvolt}))', 'H (\\si{\\angstrom})', 'D-\\ce{Mo}(\\si{\\angstrom})','D-\\ce{S}(\\si{\\angstrom})', 'd (\\si{\\angstrom})', '$\\gamma$ (deg)' ]]# List all avaliable proprieties and the formated name (with units) for the Tex tables
    figures='yes'# Print the figures of the positions and the table for the figures


    print('-----------------------------------------------')
    print('Programa para extração de dados do VASP')
    print('Autor: Natan M. Regis - n.m.regis@df.ufscar.br')
    print('-----------------------------------------------')
    
    out=[]

    root1='/run/media/natan/home-natan/estudos/faculdade/IC/MoS2/results/ads_pristine'
    layer='pristine'
    os.chdir(root1)
    Molecules=['H2', 'H2O', 'NH3','CH4','N2','CO','O2','N2O','CO2','NO2','SO2']# Molecules that will be calculated the proprieties, e.g. the roots directories
    recursive='yes'
    S_defects=0# Number of sulphur defects in the system. Set 0 to neither defects
    pris = navigate_dir(Molecules,root1)
   
    root2='/run/media/natan/home-natan/estudos/faculdade/IC/MoS2/results/ads_defective_far'
    os.chdir(root2)
    layer='far'
    Molecules=['H2', 'H2O', 'NH3','CH4','N2','CO','O2','N2O','CO2','NO2','SO2']# Molecules that will be calculated the proprieties, e.g. the roots directories
    recursive='yes'
    S_defects=1# Number of sulphur defects in the system. Set 0 to neither defects
    far = navigate_dir(Molecules,root2)

    root3='/run/media/natan/home-natan/estudos/faculdade/IC/MoS2/results/ads_defective_close'
    os.chdir(root3)
    layer='close'
    Molecules=['H2', 'H2O', 'NH3','CH4','N2','CO','O2','N2O','CO2','NO2','SO2']# Molecules that will be calculated the proprieties, e.g. the roots directories
    recursive='yes'
    S_defects=1# Number of sulphur defects in the system. Set 0 to neither defects
    close = navigate_dir(Molecules,root3)

    out=[pris,far,close]

    saida_tex(out)

if __name__ == "__main__":
    main();
