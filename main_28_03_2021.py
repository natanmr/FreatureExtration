#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 13:36:16 2021
Programa para extrair informações do vasp
@author: natan
"""

#---Modulos---#

# Definição de modulos para a execussão do programa

import xml.etree.ElementTree as ET # ler o xml
import pandas as pd # Prático para criar uns dataframe
from math import sqrt
import numpy as np
import os  # Usado para navegar entre diretorios
import sys


#---Opções Gobais---#
pd.set_option("display.max_rows", None, "display.max_columns", None)# Imprime todo dataframe


#---Pharser do xml---#
def xmlpharser():
    tree = ET.parse('vasprun.xml')
    root = tree.getroot()
    return root #Retorna a arvore do xml
    
def CheckStringInFile(File,String):
    with open(File) as f:
        datafile = f.readlines()
    found = False  
    for line in datafile:
        if String in line:
            return True
    return False  


#---Informações da base de átomos---#
def BasisInfo():
    root = xmlpharser() # Recebe a arvore do xml
    basis = root.findall("structure/crystal//varray[@name='basis']")[0]
    elements = basis.findall("v")
    basis = [ ]
    for element in elements:
        element = element.text.split()
        x = float(element[0])
        y = float(element[1])
        z = float(element[2])
        basis.append([x, y, z])
    return basis    


#---Informações dos átomos final---#
def AtomsInfo():
    basis = BasisInfo()
    root = xmlpharser() # Recebe a arvore do xml
    atoms = root.findall("structure[@name='finalpos']/varray[@name='positions']/")
    species = root.findall("atominfo/array[@name='atoms']/set/rc/c[1]")
    elements=[]
    for element in species:
        elements.append(element.text)    
    basis = np.matrix(basis)    
    position=[]
    i = 0
    for v in atoms:
        x = float(v.text.split()[0])
        y = float(v.text.split()[1])
        z = float(v.text.split()[2])
        position.append([x, y, z])
        i = i+1
    position = np.matrix(position) 
    position = position*basis
    positions = pd.DataFrame(position, columns=(['x','y','z']))
    positions['Atom'] = elements
    #print(positions)
    return positions

#---Verifica se tudo certo com ecutoff        
def Ecutoff():  
    root = xmlpharser()
    ENCUT_default = float(488.73400) 
    ENCUT = root.findall("incar/i[@name='ENCUT']")[0]
    ENCUT = float(ENCUT.text)
    if ENCUT == ENCUT_default:
        print('ENCUT ok!')
    else:
        print('Incorrect ENCUT!')
        print('erro na conta....')
        sys.exit()
    return ENCUT    

#---Enegia total---#
def Enmax():
    root = xmlpharser()
    EEotEWntropy = 0
    CalculationStep = root.findall("calculation")[-1]
    EEotEWntropy = CalculationStep.findall("energy/i[@name='e_wo_entrp']")[0]
    EEotEWntropy = float(EEotEWntropy.text)
    return EEotEWntropy

#---Parâmetro de rede---#
def a0():    
    positions = AtomsInfo()      
    dist = []     
    for i in range(0,25-1,1):
        x1 = float(positions.loc[i][0]) 
        y1 = float(positions.loc[i][1])
        z1 = float(positions.loc[i][2])
        v1 = np.matrix([x1, y1, z1])
        for i in range (0,25-1,1):
            x2 = float(positions.loc[i][0])
            y2 = float(positions.loc[i][1])
            z2 = float(positions.loc[i][2])    
            v2 = np.matrix([x2, y2, z2])
            dist.append(np.linalg.norm(v1-v2))
    dist = sorted(dist)
    dist = min(numero for numero in dist if numero != 0) # Menor distância excluindo zero
    return dist 

#---Altura das moleculas para layer---#
def Hight():
    positions = AtomsInfo()
    LenPositions = int(len(positions))
    p = 0
    x = 0
    y = 0
    z = 0
    for i in range(75,LenPositions):
        x = x + positions['x'][i]
        y = y + positions['y'][i]
        z = z + positions['z'][i]
        i = i+1
        p = p+1
    xMed = x/float(p)
    yMed = y/float(p)
    zMed = z/float(p)
    zS = positions['z'][43]
    H = np.linalg.norm(zS-zMed)
    return H

def Distance():
    positions = AtomsInfo()
    LenPositions = int(len(positions))
    
    dist = ['Dist','Value'] 
    for i in range(76,LenPositions):
        xAtom1 = positions['x'][i-1]
        yAtom1 = positions['y'][i-1]
        zAtom1 = positions['z'][i-1]
        Specie1 = positions['Atom'][i-1]
        
        xAtom2 = positions['x'][i]
        yAtom2 = positions['y'][i]
        zAtom2 = positions['z'][i]        
        Specie2 = positions['Atom'][i]
        
        i = i+1
        
        v1 = np.matrix([xAtom1, yAtom1, zAtom1])
        v2 = np.matrix([xAtom2, yAtom2, zAtom2])
        
        DistanceMol = np.linalg.norm(v1-v2)
        
        dist.append(['{}--{}'.format(Specie1,Specie2), DistanceMol])
    print(dist)
    return DistanceMol


#---Informações das moléculas que podem ser usadas para alguns cálculos---#
def InfoMol():
    Molecules = ({
    'Molecule': ['H2', 'O2', 'N2', 'CO', 'NO', 'CO2', 'NO2', 'N2O', 'H2O', 'NH3', 'SO2', 'CH4'],    
    'Energy__w': [-6.75643835, -9.99914716, -16.77964988, -14.90957902, -12.39910767, -23.17391481, -18.58684634, -21.57542293, -14.26654709, -19.53889317, -17.01184620, -24.03235987], 
    'Distance': [0.75243, 1.22470, 1.10606, 1.13830, 1.16232, 1.17272, 1.20745, [1.19360, 1.13997], 0.97084, 1.02255, 1.45386, 1.09666],
     'Angle': [None, None, None, None, None, 180., 134.0, 180, 104.0, 106.1, 119.4, 109.5]       
        } ) 
    Molecules = pd.DataFrame(Molecules)
    MoS25x5 = ({
        'System': ['Pristine 5x5', 'Defective 5x5'],
        'Energy_w':[-577.73633023, -570.79075056],
        'a0': [3.16073, 3.04771],       
        })
    MoS25x5 = pd.DataFrame(MoS25x5)
    return Molecules, MoS25x5

#---Extrai a energia de adsorção---#
def Adsorption_energy(IndexMol):
    EnergyAdsorbed = Enmax()
    Molecules, MoS25x5 = InfoMol()
    MolEnergy= Molecules.loc[0,'Energy__w']
    LayerEnergy = MoS25x5.loc[0,'Energy_w']
    AdsorptionEnergy = float(EnergyAdsorbed) - float(MolEnergy) - float(LayerEnergy)
    #print(EnergyAdsorbed, MolEnergy, LayerEnergy)
    return AdsorptionEnergy

#---Navega entre os diretórios para extrair informações---#
def NavigateDir():
    FileTEX = open('table.tex','w')
    FileTEX.write('Configuration & E_{tot} (eV) & E_{Ads} & H (\AA) & Dist & Angle  \\\ \n')
    Molecules=['CH4','CO2','H2','H2O','N2','O2']
    RootMol = os.getcwd()
    for Mol in Molecules:
        os.chdir(Mol)
        FileTEX.write('{} \\\ \n'.format(Mol))
        Root = os.getcwd()
        Folders = [x[0] for x in os.walk('.')]    
        for Folder in Folders:
            os.chdir(Folder)
            print(Folder)
            Files = [x[2] for x in os.walk('.')][0]
            if 'vasprun.xml' in Files:
                EnergyAdsorbed = Enmax() 
                Altura = Hight()
                AdsorptionEnergy=0
                Dist = 0
                Angle = 0
                FileTEX.write('{} & {} & {} & {} & {} & {} \\\ \n'.format(Folder,EnergyAdsorbed, Adsorption_energy, Altura, Dist, Angle) )
            else:
                FileTEX.write('{} & Não Convergiu & Não Convergiu & Não Convergiu & Não Convergiu  \\\ \n'.format(Folder))

            os.chdir(Root)
        os.chdir(RootMol)    
def main():    
    NavigateDir()
    
if __name__ == "__main__":
    main();  
    
    
