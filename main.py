#!/usr/bin/env python3
# -*- coding: utf-8 -*-

    #############################################
    #  Programa para extração de Informações    #
    #  do VASP                                  #
    #  Autor: Natan Moreira Regis               #
    #  Ultima Atualização: 22 jun 2021          #
    #############################################

#---Modulos---#
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

#---Pharser de qualquer arquivos---#
def openfile(filename):
    file = open(filename, 'r')
    lines = file.readlines()
    return lines

#---Checar se uma string está no arquivo---#
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

#---Verifica se tudo certo com ecutoff---#
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

#---Energia de Fermi---#
def FermiEnergy():
    lines = openfile('OUTCAR')
    Linesefermi = []
    # exibe a energia de fermi final
    for line in lines:
        if 'E-fermi' in line:
            Linesefermi.append(line)
    efermi = Linesefermi[-1].split()[2]
    #print('Fermi energy',efermi)
    return float(efermi)

#---Parâmetro de rede (a0)---#
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

#---Distância Centro Molecula atomos de Mo--->
def DistCenterMolMo():
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

    dist=[]
    for i in range(1,26):
        xMo = positions['x'][i]
        yMo = positions['y'][i]
        zMo = positions['z'][i]

        xMoltoMo = np.linalg.norm(xMo-xMed)
        yMoltoMo = np.linalg.norm(yMo-yMed)
        zMoltoMo = np.linalg.norm(zMo-zMed)

        Dist= np.sqrt(float(xMoltoMo)**2+float(yMoltoMo)**2 + float(zMoltoMo)**2)
        dist.append(Dist)

    Dist = sorted(dist)
    Dist = min(numero for numero in dist if numero != 0)

    return Dist

#---Distância Centro Molecula atomos de S--->
def DistCenterMolS():
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
    dist=[]
    for i in range(27,77):
        xS = positions['x'][i]
        yS = positions['y'][i]
        zS = positions['z'][i]
        xMoltoS = np.linalg.norm(xS-xMed)
        yMoltoS = np.linalg.norm(yS-yMed)
        zMoltoS = np.linalg.norm(zS-zMed)
        Dist = np.sqrt(float(xMoltoS)**2+float(yMoltoS)**2 + float(zMoltoS)**2)
        dist.append(Dist)
    Dist = sorted(dist)
    Dist = min(numero for numero in dist if numero != 0)
    return Dist

#---Distância interatômica das moléculas--->
def Distance(Mol):
    positions = AtomsInfo()
    LenPositions = int(len(positions))
    if Mol =='H2' or Mol=='N2' or Mol =='O2' or Mol=='NO' or Mol=='CO':
        i = 75
        xAtom1 = positions['x'][i]
        yAtom1 = positions['y'][i]
        zAtom1 = positions['z'][i]

        xAtom2 = positions['x'][i+1]
        yAtom2 = positions['y'][i+1]
        zAtom2 = positions['z'][i+1]

    elif Mol=='H2O' or Mol =='CO2':
        i = 75
        xAtom1 = positions['x'][i]
        yAtom1 = positions['y'][i]
        zAtom1 = positions['z'][i]

        xAtom2 = positions['x'][i+2]
        yAtom2 = positions['y'][i+2]
        zAtom2 = positions['z'][i+2]

    elif Mol=='NH3':
        i = 75
        xAtom1 = positions['x'][i]
        yAtom1 = positions['y'][i]
        zAtom1 = positions['z'][i]

        xAtom2 = positions['x'][i+2]
        yAtom2 = positions['y'][i+2]
        zAtom2 = positions['z'][i+2]

    elif Mol=='CH4':
        i =75
        xAtom1 = positions['x'][i+1]
        yAtom1 = positions['y'][i+1]
        zAtom1 = positions['z'][i+1]

        xAtom2 = positions['x'][i+4]
        yAtom2 = positions['y'][i+4]
        zAtom2 = positions['z'][i+4]

    v1 = np.matrix([xAtom1, yAtom1, zAtom1])
    v2 = np.matrix([xAtom2, yAtom2, zAtom2])

    DistanceMol = np.linalg.norm(v1-v2)
    return DistanceMol

#---Angulo moleculas---#
def Angulo(Mol):
    Positions = AtomsInfo()
    if len(Positions) == 77:
        return 0
    else:
        if Mol =='H2O' or Mol =='CO2' or Mol=='N2O':
            i = 75
            x1 = float(Positions.loc[i][0])
            y1 = float(Positions.loc[i][1])
            z1 = float(Positions.loc[i][2])

            x3 = float(Positions.loc[i+1][0])
            y3 = float(Positions.loc[i+1][1])
            z3 = float(Positions.loc[i+1][2])

            x2 = float(Positions.loc[i+2][0])
            y2 = float(Positions.loc[i+2][1])
            z2 = float(Positions.loc[i+2][2])
        elif Mol=='CH4':
            i = 75
            x1 = float(Positions.loc[i][0])
            y1 = float(Positions.loc[i][1])
            z1 = float(Positions.loc[i][2])
            x3 = float(Positions.loc[i+1][0])
            y3 = float(Positions.loc[i+1][1])
            z3 = float(Positions.loc[i+1][2])
            x2 = float(Positions.loc[79][0])
            y2 = float(Positions.loc[79][1])
            z2 = float(Positions.loc[79][2])
        elif Mol=='NH3':
            i = 75
            x1 = float(Positions.loc[i][0])
            y1 = float(Positions.loc[i][1])
            z1 = float(Positions.loc[i][2])
            x3 = float(Positions.loc[i+1][0])
            y3 = float(Positions.loc[i+1][1])
            z3 = float(Positions.loc[i+1][2])
            x2 = float(Positions.loc[i+2][0])
            y2 = float(Positions.loc[i+2][1])
            z2 = float(Positions.loc[i+2][2])

    a = np.array([x1,y1,z1])
    b = np.array([x2,y2,z2])
    c = np.array([x3,y3,z3])

    ba = a - b
    bc = c - b

    cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
    angle = np.arccos(cosine_angle)
    angle = np.degrees(angle)
    return angle

#---Menor valor de uma lista---#
def LineLowestValue(List,col):
    menor = List[0][col]
    LinhaMenor=0
    for i in range(len(List)):
        if List[i][col] < menor:
            menor = List[i][col]
            LinhaMenor=i
    return LinhaMenor

#---Informações das moléculas---#
def InfoMol():
    Molecules = ({
    'Molecule': ['H2', 'O2', 'N2', 'CO', 'NO', 'CO2', 'NO2', 'N2O', 'H2O', 'NH3', 'SO2', 'CH4'],
    'Energy__w': [-6.75643835,-9.99914716, -16.77964988, -14.90957902, -12.39910767, -23.17391481, -18.58684634, -21.57542293, -14.26654709, -19.53889317, -17.01184620, -24.03235987],
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
def Adsorption_energy(Mol):
    EnergyAdsorbed = Enmax()
    Molecules, MoS25x5 = InfoMol()
    if Mol=='H2':
        IndexMol = 0
    elif Mol=='O2':
        IndexMol = 1
    elif Mol=='N2':
        IndexMol = 2
    elif Mol=='CO':
        IndexMol= 3
    elif Mol=='NO':
        IndexMol= 4
    elif Mol=='CO2':
        IndexMol = 5
    elif Mol=='N2O':
        IndexMol= 7
    elif Mol=='H2O':
        IndexMol = 8
    elif Mol=='NH3':
        IndexMol=9
    elif Mol=='CH4':
        IndexMol = 11

    MolEnergy= Molecules.loc[IndexMol,'Energy__w']
    LayerEnergy = MoS25x5.loc[0,'Energy_w']
    AdsorptionEnergy = float(EnergyAdsorbed) - float(MolEnergy) - float(LayerEnergy)
    return AdsorptionEnergy

#---Navega entre os diretórios para extrair informações---#
def NavigateDir():
    FileTEX = open('table.tex','w')
    FileTEX.write('\\begin{longtable}{ccccccc} \n')
    FileTEX.write('\hline Configuration & $E_{Tot}$ (eV) & $E_{Ads}$ (eV) & D-\\ce{Mo} (\AA) & D-\\ce{S} (\AA) & d(\AA) & Angle (degrees) \\\ \n')
    Molecules=['H2','CH4','CO','CO2','H2O','N2','NH3','NO','O2']
    #Molecules=['CO']
    RootMol = os.getcwd()
    MinusValueEads=[]
    for Mol in Molecules:
        ValuestoPrint=[]
        os.chdir(Mol)
        FileTEX.write('\\hline \n \\multicolumn{{7}}{{c}}{{{}}} \\\ \n \\hline \n'.format(Mol))
        Root = os.getcwd()
        Folders = [x[0] for x in os.walk('.')]

        for Folder in Folders:
            os.chdir(Folder)
            print(Folder)
            Files = [x[2] for x in os.walk('.')][0]
            if 'vasprun.xml' in Files:
                Check = CheckStringInFile('OUTCAR','reached required accuracy - stopping structural energy minimisation')
                if Check ==True:
                    EnergyAdsorbed = round(Enmax(),8)
                    #Altura = round(Hight(),2)
                    AdsorptionEnergy = round(Adsorption_energy(Mol),4)
                    Dist = round(Distance(Mol),2)
                    DistoMo = round(DistCenterMolMo(),2)
                    DistoS = round(DistCenterMolS(),2)
                    Angle = round(Angulo(Mol),2)
                    Path = Folder.replace('./','')
                    Path = Path.replace('/','-')
                    ValuestoPrint.append([str(Path), float(EnergyAdsorbed), float(AdsorptionEnergy), float(DistoMo), float(DistoS) , float(Dist), float(Angle)])
                    MinusValueEads.append([float(AdsorptionEnergy),str(Path)])
            os.chdir(Root)
        LineMenorEads = LineLowestValue(ValuestoPrint,1)
        for i in range(len(ValuestoPrint)):
            if LineMenorEads == i:
                FileTEX.write('\\textbf{{ {} }} & \\num{{{}}} &\\textbf{{ \\num{{{}}} }} & \\textbf{{ \\tablenum{{{}}} }} & \\textbf{{ \\tablenum{{{}}} }} & \\textbf{{ \\tablenum{{{}}} }} & \\textbf{{ \\tablenum{{{}}} }}  \\\ \n'.format(ValuestoPrint[i][0], ValuestoPrint[i][1],ValuestoPrint[i][2],ValuestoPrint[i][3],ValuestoPrint[i][4], ValuestoPrint[i][5], ValuestoPrint[i][6] ))
            else:
                FileTEX.write('{} & \\num{{{}}} & \\num{{{}}} &\\tablenum{{{}}} & \\tablenum{{{}}} & \\tablenum{{{}}} & \\tablenum{{{}}} \\\ \n'.format(ValuestoPrint[i][0], ValuestoPrint[i][1], ValuestoPrint[i][2], ValuestoPrint[i][3], ValuestoPrint[i][4], ValuestoPrint[i][5], ValuestoPrint[i][6]))

        os.chdir(RootMol)

    for i in range(len(MinusValueEads)):
        MinusValue = 0
        if MinusValueEads[i][0]<0:
            if abs(MinusValueEads[i][0]) > MinusValue:
                MinusValue = abs(MinusValueEads[i][0])
                res = i
    print("menorEads:",MinusValueEads[res][1],':',MinusValueEads[res][0])

    FileTEX.write('\\hline')
    FileTEX.write('\\end{longtable} \n')

#---Função main do programa---#
def main():
    NavigateDir()

if __name__ == "__main__":
    main();
