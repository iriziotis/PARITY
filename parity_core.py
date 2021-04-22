#!/usr/bin/env python3
import sys
from subprocess import call
import os
import signal
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdchem import EditableMol
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem.rdmolfiles import MolFromMolFile
from rdkit.Chem import MACCSkeys
from rdkit.Chem.AtomPairs import Pairs
from rdkit.Chem.AtomPairs import Torsions
from rdkit.Chem import rdFMCS as MCS
from rdkit.Chem.rdFMCS import BondCompare, AtomCompare


def generate_parity(invariable, variable, faster=1):
    """ Calculates PARITY score between two molecules (SMILES or .mol).
    If one of the molecules has an R group, it has to be put in the 
    second argument (sm_2)."""

    # Check if input is a mol file or SMILES string
    try:
        mol_1 = MolFromMolFile(invariable)
    except Exception as e:
        mol_1 = Chem.MolFromSmiles(invariable)
    try:
        mol_2 = MolFromMolFile(variable)
    except Exception:
        mol_2 = Chem.MolFromSmiles(variable)

    ms = [mol_1,mol_2]

    smiles_1=Chem.MolToSmiles(mol_1)
    smiles_2=Chem.MolToSmiles(mol_2)
    
    matches_best=None
    sim_score_best=None
    markush_best=None

    #First match bond types and elements, then match any bond type and elements, 
    #then any bond type and any element to find the best score for a low similar fragment, 
    #for more similar fragments then any any will return the best match

    mcs_graph_bondtypes_elements=MCS.FindMCS(ms,bondCompare=BondCompare.CompareOrder,atomCompare=AtomCompare.CompareElements,timeout=int(10/faster),completeRingsOnly=True)
    matches_best,sim_score_best,markush_best=generate_sim_score(mol_1,mol_2,mcs_graph_bondtypes_elements.smartsString)

    mcs_graph_any_elements=MCS.FindMCS(ms,bondCompare=BondCompare.CompareAny,atomCompare=AtomCompare.CompareElements,timeout=int(10/faster),completeRingsOnly=True)
    matches,sim_score,markush=generate_sim_score(mol_1,mol_2,mcs_graph_any_elements.smartsString)
    if sim_score>sim_score_best:
        matches_best, sim_score_best, markush_best = matches, sim_score, markush

    mcs_graph_any_any=MCS.FindMCS(ms,bondCompare=BondCompare.CompareAny,atomCompare=AtomCompare.CompareAny,timeout=int(30/faster),completeRingsOnly=True)
    matches,sim_score,markush=generate_sim_score(mol_1,mol_2,mcs_graph_any_any.smartsString)
    if sim_score>sim_score_best:
        matches_best, sim_score_best, markush_best = matches, sim_score, markush

    return sim_score_best

def generate_sim_score(mol_1,mol_2,smarts):
    smiles_1=Chem.MolToSmiles(mol_1)
    smiles_2=Chem.MolToSmiles(mol_2)
    #print smiles_1,smiles_2
    best_markush=0
    if smiles_1==smiles_2:
        best_matches=mol_1.GetNumAtoms()
        best_sim_score=1.0
    elif mol_1.GetNumAtoms()==1 and smiles_1=="[H+]":
        best_matches=0
        best_sim_score=0.0
    elif mol_2.GetNumAtoms()==1 and smiles_2=="[H+]":
        best_matches=0
        best_sim_score=0.0
    elif mol_1.GetNumAtoms()==2 and smiles_1=="[HH]":
        best_matches=0
        best_sim_score=0.0
    elif mol_2.GetNumAtoms()==2 and smiles_2=="[HH]":
        best_matches=0
        best_sim_score=0.0
    else:
        if mol_1.GetNumAtoms()==1:
            smarts=Chem.MolToSmarts(mol_1)
        elif mol_2.GetNumAtoms()==1:
            smarts=Chem.MolToSmarts(mol_2)
        if smarts==None:
            best_matches=0
            best_markush=0
            best_sim_score=0.0
        else:
            matches_1=get_matches(mol_1,smarts)
            matches_2=get_matches(mol_2,smarts)
            best_sim_score=0.0
            best_matches=0
            best_markush=0
            size_1=mol_1.GetNumAtoms()
            size_2=mol_2.GetNumAtoms()

            for match_1 in matches_1:
                for match_2 in matches_2:
                    sim_score=0.0
                    matches=0
                    #need to count markush to gross up size_2
                    markush=0
                    for i in range(0,len(match_1)):
                        atom_1=mol_1.GetAtomWithIdx(match_1[i])
                        atom_2=mol_2.GetAtomWithIdx(match_2[i])
                        symbol_1=get_symbol(atom_1)
                        symbol_2=get_symbol(atom_2)
                        if symbol_2=="*":
                            markush_size=get_atoms_in_markush(mol_1,match_1,atom_1)
                            matches+=1
                            markush+=(markush_size-1)
                        elif symbol_1==symbol_2:
                            matches+=1
                    sim_score=float(matches+markush)/float(size_1 + (size_2+markush)  - (matches+markush))
                    if matches>best_matches:
                        best_sim_score=sim_score
                        best_matches=matches
                        best_markush=markush
                    elif matches==best_matches and sim_score>best_sim_score:
                        best_sim_score=sim_score
                        best_matches=matches
                        best_markush=markush
    return (best_matches+best_markush),best_sim_score,best_markush

def generate_reaction_score(l):
    matches=0
    atoms=0
    for i in range(0,len(l)):
        counts=l[i][1]
        if counts[1]!='NA':
            atoms+=counts[0]
            atoms+=counts[1]
            matches+=counts[2]
    rxn_score=0
    if (atoms-matches)==0:
        rxn_score=0
    else:
        rxn_score=float(matches)/float(atoms-matches)
    return rxn_score

def get_atoms_in_markush(mol,match,atom):
    #print "markush",match
    match_list=list(match)
    match=0
    to_visit=set([atom.GetIdx()])
    visited=set()

    while len(to_visit)>0:
        next_id=list(to_visit)[0]
        match+=1
        visited.add(next_id)
        next_atom=mol.GetAtomWithIdx(next_id)
        neighbors=next_atom.GetNeighbors()
        for neigh in neighbors:
            if neigh.GetIdx() not in match_list:
                to_visit.add(neigh.GetIdx())
                #print "to visit",to_visit      
        to_visit=to_visit.difference(visited)
        #print "to visit",to_visit
        #print match
    return match

def get_matches(mol,smarts):
    patt=Chem.MolFromSmarts(smarts)
    matches=mol.GetSubstructMatches(patt,uniquify=False)
    return matches

def get_symbol(atom):
    symbol=atom.GetSymbol()
    return symbol

if __name__ == "__main__":
    # FAD
    sm_1 = "Cc1cc2c(cc1C)N(C3=NC(=O)NC(=O)C3=N2)C[C@@H]([C@@H]([C@@H](CO[P@@](=O)(O)O[P@](=O)(O)OC[C@@H]4[C@H]([C@H]([C@@H](O4)n5cnc6c5ncnc6N)O)O)O)O)O"
    # NAD
    sm_2 = "c1cc(c[n+](c1)[C@H]2[C@@H]([C@@H]([C@H](O2)CO[P@@](=O)([O-])O[P@@](=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)n4cnc5c4ncnc5N)O)O)O)O)C(=O)N"
    print(generate_parity(sm_1, sm_2))
