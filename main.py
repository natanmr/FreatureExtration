#!/usr/bin/env python
# -*- coding: utf-8 -*-

####################
## python modules ##
####################
import numpy as np
import os

######################################################################
## Title: LibrePhsysics VASP Freature extraction (LPVFE)            ##
## Description: Program to colect some stufs on VASP Program        ##
## Version: 0.3                                                     ##
## author: Natan Moreira Regis                                      ##
## Project Name: LibrePhysics (librephysics.xyz)                    ##
## license: GPLv3                                                   ##
## maintainer: "Natan Moreira Regis"                                ##
## email: "n.m.regis@df.ufscar.br"                                  ##
## status: "Production"                                             ##
######################################################################

#####################
## Control Options ##
#####################


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

##########################
## Atoms and basis info ##
##########################
class AtomsInfo:

    ## Basis information
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
        basis=np.array(basis)
        return basis
    
    ## Atomic species Info
    def species_info():
        lines = openfile("CONTCAR")
        LineAtoms=lines[5].split()
        LineNumAtoms=lines[6].split()
        Atoms=[LineAtoms,LineNumAtoms]
        TotAtoms=0
        for j in range(len(LineNumAtoms)):
            TotAtoms+=int(LineNumAtoms[j])
        return Atoms, TotAtoms
    
    ## Final atomic positions
    def atoms_info():
        basis = AtomsInfo.basis_info()
        lines = openfile("CONTCAR")
        Species, TotAtoms = AtomsInfo.species_info()
        Positions=[]   
        if lines[7].split()[0]=="Direct" or lines[7].split()[0]=="Cartesian":
            for i in range(8,TotAtoms+8):
                v = lines[i].split()
                x = float(v[0])
                y = float(v[1])
                z = float(v[2])
                Positions.append([x,y,z])
        else:
            for j in range(9,TotAtoms+9):
                v = lines[j].split()
                x = float(v[0])
                y = float(v[1])
                z = float(v[2])
                Positions.append([x,y,z])
        Positions = np.array(Positions)
        Positions= np.matrix(Positions)*np.matrix(basis)
        return Positions

##################
## Total Energy ##
##################
def enmax():
    lines=openfile("OUTCAR")
    LinesEtot=[]
    for line in lines:
        if "energy  w" in line:
            LinesEtot.append(line)
    FinalEtot=LinesEtot[-1]    
    EEotEWntropy = float(FinalEtot.split()[-1])
    return EEotEWntropy

###############
## Distances ##
###############
class Distance:
    
    def search_in_list(lista, item):
        try:
            return lista.index(item)
        except ValueError:
            return -1    
    
    def Dist_two_vectors(v1,v2):
        return np.linalg.norm(v2-v1)
    
    def lowest_dist_ignoring_0(Dist):
        Dist = sorted(Dist)
        #print(Dist)
        Dist = min(numero for numero in Dist if numero != 0) 
        return Dist
    
    def center_mol():
        positions = AtomsInfo.atoms_info()
        Atoms, TotAtoms = AtomsInfo.species_info()
        IndexMo=int(Distance.search_in_list(Atoms[0], "Mo"))
        IndexS=int(Distance.search_in_list(Atoms[0], "S"))
        IndexInit=int(Atoms[1][IndexMo]) + int(Atoms[1][IndexS])
        p = 0
        x = 0
        y = 0
        z = 0
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
        #print(vMed)
        return vMed
        
    def dist_mol_to_atom(Atom):
        positions = AtomsInfo.atoms_info()
        Atoms, TotAtoms = AtomsInfo.species_info()
        vMed=Distance.center_mol()
        IndexAtom=int(Distance.search_in_list(Atoms[0], Atom))
        LineInitiAtom=0
        i=0
        while i<IndexAtom:
            LineInitiAtom+=int(Atoms[1][i])
            i+=1
        Dist = []
        for i in range(LineInitiAtom,int(Atoms[1][IndexAtom])):
            x1 = float(positions[i,0])
            y1 = float(positions[i,1])
            z1 = float(positions[i,2])
            if z1>5.0:
                z1-= float(AtomsInfo.basis_info()[2,2])
            v1 = np.matrix([x1, y1, z1])
            Dist.append(Distance.Dist_two_vectors(v1,vMed))
        return Distance.lowest_dist_ignoring_0(Dist)
    
    def distance(direction,Atom1,Atom2):
        positions = AtomsInfo.atoms_info()
        Atoms, TotAtoms = AtomsInfo.species_info()
        IndexAtom1=int(Distance.search_in_list(Atoms[0], Atom1))
        IndexAtom2=int(Distance.search_in_list(Atoms[0], Atom2))
        i=0
        LineInitiAtom1=0
        LineInitiAtom2=0
        while i<IndexAtom1:
            LineInitiAtom1+=int(Atoms[1][i])
            i+=1
        j=0
        while j<IndexAtom2:
            LineInitiAtom2+=int(Atoms[1][j])
            j+=1        
        Dist = []
        for i in range(LineInitiAtom1,int(Atoms[1][IndexAtom1])):
            x1 = float(positions[i,0])
            y1 = float(positions[i,1])
            z1 = float(positions[i,2])
            v1 = np.matrix([x1, y1, z1])
            for j in range (LineInitiAtom2,int(Atoms[1][IndexAtom2])):
                x2 = float(positions[j,0])
                y2 = float(positions[j,1])
                z2 = float(positions[j,2])
                v2 = np.matrix([x2, y2, z2])
                
                Dist.append(Distance.Dist_two_vectors(v1,v2))
        return Distance.lowest_dist_ignoring_0(Dist)

##########################################
## Distância interatômica das moléculas ##
##########################################
def distance(Mol):
    positions = AtomsInfo.atoms_info()
    if Mol =='H2' or Mol=='N2' or Mol =='O2' or Mol=='NO' or Mol=='CO'or Mol=='SO2':
        i = 75
        xAtom1 = positions[i,0]
        yAtom1 = positions[i,1]
        zAtom1 = positions[i,2]
        
        xAtom2 = positions[i+1,0]
        yAtom2 = positions[i+1,1]
        zAtom2 = positions[i+1,2]

    elif Mol=='H2O' or Mol =='CO2':
        i = 75
        xAtom1 = positions[i,0]
        yAtom1 = positions[i,1]
        zAtom1 = positions[i,2]

        xAtom2 = positions[i+2,0]
        yAtom2 = positions[i+2,1]
        zAtom2 = positions[i+2,2]

    elif Mol=='NH3':
        i = 75
        xAtom1 = positions[i,0]
        yAtom1 = positions[i,1]
        zAtom1 = positions[i,2]

        xAtom2 = positions[i+2,0]
        yAtom2 = positions[i+2,1]
        zAtom2 = positions[i+2,2]

    elif Mol=='CH4':
        i =75
        xAtom1 = positions[i+1,0]
        yAtom1 = positions[i+1,1]
        zAtom1 = positions[i+1,2]

        xAtom2 = positions[i+4,0]
        yAtom2 = positions[i+4,1]
        zAtom2 = positions[i+4,2]

    v1 = np.matrix([xAtom1, yAtom1, zAtom1])
    v2 = np.matrix([xAtom2, yAtom2, zAtom2])

    DistanceMol = np.linalg.norm(v1-v2)
    return DistanceMol

######################
## Angulo moleculas ##
######################
def angle(Mol):
    Positions = AtomsInfo.atoms_info()
    if len(Positions) == 77:
        return 0
    else:
        if Mol =='H2O' or Mol =='CO2' or Mol=='N2O' or Mol=='SO2':
            i = 75
            x1 = float(Positions[i,0])
            y1 = float(Positions[i,1])
            z1 = float(Positions[i,2])

            x3 = float(Positions[i+1,0])
            y3 = float(Positions[i+1,1])
            z3 = float(Positions[i+1,2])

            x2 = float(Positions[i+2,0])
            y2 = float(Positions[i+2,1])
            z2 = float(Positions[i+2,2])
        elif Mol=='CH4':
            i = 75
            x1 = float(Positions[i,0])
            y1 = float(Positions[i,1])
            z1 = float(Positions[i,2])
            x3 = float(Positions[i+1,0])
            y3 = float(Positions[i+1,1])
            z3 = float(Positions[i+1,2])
            x2 = float(Positions[-1,0])
            y2 = float(Positions[-1,1])
            z2 = float(Positions[-1,2])
        elif Mol=='NH3':
            i = 75
            x1 = float(Positions[i,0])
            y1 = float(Positions[i,1])
            z1 = float(Positions[i,2])
            x3 = float(Positions[i+1,0])
            y3 = float(Positions[i+1,1])
            z3 = float(Positions[i+1,2])
            x2 = float(Positions[i+2,0])
            y2 = float(Positions[i+2,1])
            z2 = float(Positions[i+2,2])

    a = np.array([x1,y1,z1])
    b = np.array([x2,y2,z2])
    c = np.array([x3,y3,z3])

    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)
    angle = np.degrees(angle)
    return angle

############################
## Lowest value of a list ##
############################
def LineLowestValue(List,col):
    menor = List[0][col]
    LinhaMenor=0
    for i in range(len(List)):
        if List[i][col] < menor:
            menor = List[i][col]
            LinhaMenor=i
    return LinhaMenor

####################################
## Energies of molecule and layer ##
####################################
def EnergyTotMoleLayer():
    Molecules = [['CH4',-24.03376510], ['CO',-14.90940997], ['CO2',-23.17457153],["H2",-6.75664087],['H2O',-14.26411655],['N2',-16.77991136],['N2O',-21.57599354],['NH3',-19.54729190],['NO',-12.39971547],['NO2',-18.58681627],['O2',-9.99926374],['SO2',-17.01174353]]
    Monolayer = [['Pristine 5x5',-577.66298683], ['Defective 5x5', -570.79075056]] 
    return Molecules, Monolayer

###################################
## Adsorption energy calculation ##
###################################
def adsorption_energy(Mol):
    EnergyAdsorbed = enmax()
    Molecules, Monolayer = EnergyTotMoleLayer()
    for Molecule in Molecules:
        if Mol==Molecule[0]:
            MolEnergy=Molecule[1]
    LayerEnergy = Monolayer[0][1]
    AdsorptionEnergy = float(EnergyAdsorbed) - float(MolEnergy) - float(LayerEnergy)
    return AdsorptionEnergy

#########################
## Print to LaTex file ##
#########################
class PrintToTex:

    def open_file_tex():
        FileTEX = open('table.tex','w')# Exit file
        return FileTEX
    
    def table_header(FileTEX):
        # Table Header:
        FileTEX.write('\\setlength{\\tabcolsep}{3pt} \n')
        FileTEX.write('\\begin{longtable}{ccccccc} \n')
        FileTEX.write('\\hline \n')
        FileTEX.write('Configuration & $E_{Tot}$ (eV) & $E_{Ads}$ (eV) & D-\\ce{Mo} (\AA) & D-\\ce{S} (\AA) & d(\AA) & Angle ($^{\circ}$) \\\ \n')
    
    def table_line_mol(FileTEX,Mol):
        FileTEX.write('\\hline \n \\multicolumn{{7}}{{c}}{{{}}} \\\ \n \\hline \n'.format(Mol))
    
    def table_line(FileTEX,ValuestoPrint,i):
        FileTEX.write('{} & \\num{{{}}} & \\num{{{}}} &\\tablenum{{{}}} & \\tablenum{{{}}} & \\tablenum{{{}}} & \\tablenum{{{}}} \\\ \n'.format(ValuestoPrint[i][0], ValuestoPrint[i][1], ValuestoPrint[i][2], ValuestoPrint[i][3], ValuestoPrint[i][4], ValuestoPrint[i][5], ValuestoPrint[i][6]))
    
    def table_line_lowest_eads(FileTEX,ValuestoPrint,i):
        FileTEX.write('\\textbf{{ {} }} & \\num{{{}}} &\\textbf{{ \\num{{{}}} }} & \\textbf{{ \\tablenum{{{}}} }} & \\textbf{{ \\tablenum{{{}}} }} & \\textbf{{ \\tablenum{{{}}} }} & \\textbf{{ \\tablenum{{{}}} }}  \\\ \n'.format(ValuestoPrint[i][0], ValuestoPrint[i][1],ValuestoPrint[i][2],ValuestoPrint[i][3],ValuestoPrint[i][4], ValuestoPrint[i][5], ValuestoPrint[i][6] ))
            
    def table_footer(FileTEX):
        # Table Footer:
        FileTEX.write('\\hline \n')
        FileTEX.write('\\end{longtable} \n')
        
    def close_file_tex(FileTEX):
        FileTEX.close()
    
#########################################################
## Navega entre os diretórios para extrair informações ##
#########################################################
def navigate_dir():
    
    # Open file and print the header
    FileTEX = PrintToTex.open_file_tex()    
    PrintToTex.table_header(FileTEX)

    Molecules=['CH4','CO','CO2','H2','H2O','N2','NH3','NO','O2','SO2']
    #Molecules=['O2']
    RootMol = os.getcwd()
    MinusValueEads=[]
    for Mol in Molecules:
        ValuestoPrint=[]
        os.chdir(Mol)
        
        ## Print line with molecule
        PrintToTex.table_line_mol(FileTEX,Mol)
        
        Root = os.getcwd()
        Folders = [x[0] for x in os.walk('.')]

        for Folder in Folders:
            os.chdir(Folder)
            print(Folder)
            Files = [x[2] for x in os.walk('.')][0]
            if 'vasprun.xml' in Files:
                Check = CheckStringInFile('OUTCAR','reached required accuracy - stopping structural energy minimisation')
                if Check ==True:
                    
                    ## system with adsorbed molecule
                    EnergyAdsorbed = round(enmax(),8)
                    
                    ## Adsorption energy calculation 
                    AdsorptionEnergy = round(adsorption_energy(Mol),4)
                    
                    ## Interatomic distances of molecule 
                    Dist = round(distance(Mol),2)
                    
                    ## Geometric center mol to Mo and S atoms
                    DisttoMo=round(Distance.dist_mol_to_atom("Mo"),2)
                    DistoS=round(Distance.dist_mol_to_atom("S"),2)
                    
                    ## Molecular angle 
                    Angle = round(angle(Mol),2)
                    
                    ####Verificar angulo do SO2
                    #### Suspeito que tem alguns erros nos indices das distancias
                    
                    ## Testes:
                    Distance.distance("XYZ","Mo","Mo")                        
                    #####    
                    
                    ## Style stufs for the name of dir 
                    Path = Folder.replace('./','')
                    Path = Path.replace('/','-')
                    
                    ## Exit values in list to ordering 
                    ValuestoPrint.append([str(Path), float(EnergyAdsorbed), float(AdsorptionEnergy), float(DisttoMo), float(DistoS) , float(Dist), float(Angle)])
                    
                    ## List of all adsorption energies 
                    MinusValueEads.append([float(AdsorptionEnergy),str(Path),str(Mol)])                
            os.chdir(Root)
             
        ## Stuf to highlight do lowest Adsorption energy 
        LineMenorEads = LineLowestValue(ValuestoPrint,1)
        for i in range(len(ValuestoPrint)):
            if LineMenorEads == i:
                PrintToTex.table_line_lowest_eads(FileTEX,ValuestoPrint,i)
            else:
                PrintToTex.table_line(FileTEX,ValuestoPrint,i)
        os.chdir(RootMol)

    ## 2 Lasts lines of table in file and close file 
    PrintToTex.table_footer(FileTEX)  
    PrintToTex.close_file_tex(FileTEX)
    
###################
## Main function ##
###################
def main():
    navigate_dir()

if __name__ == "__main__":
    main();
