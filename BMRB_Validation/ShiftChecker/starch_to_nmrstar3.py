import pynmrstar
import sys
ATOM_DICT = {
    'A': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H3'", "H4'", "H5'", "H5''"],
        'base': ["C2", "H2", "C8", "H8", "N6", "H61", "H62", "N1", "H1"],
    },

    'DA': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H2''", "H3'", "H4'", "H5'", "H5''"],
        'base': ["C2", "H2", "C8", "H8", "N6", "H61", "H62", "N1", "H1"],
    },

    'G': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H3'", "H4'", "H5'", "H5''"],
        'base': ["N1", "H1", "N2", "H21", "H22", "N3", "H3", "C8", "H8", "N7", "H7"],
    },

    'DG': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H2''", "H3'", "H4'", "H5'", "H5''"],
        'base': ["N1", "H1", "N2", "H21", "H22", "N3", "H3", "C8", "H8", "N7", "H7"],
    },

    'C': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H3'", "H4'", "H5'", "H5''"],
        'base': ["C6", "H6", "C5", "H5", "N4", "H41", "H42", "N3", "H3"],
    },

    'DC': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H2''", "H3'", "H4'", "H5'", "H5''"],
        'base': ["C6", "H6", "C5", "H5", "N4", "H41", "H42", "N3", "H3"],
    },

    'T': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H3'", "H4'", "H5'", "H5''"],
        'base': ["C6", "H6", "C5", "H5", "N3", "H3"],
    },

    'DT': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H2''", "H3'", "H4'", "H5'", "H5''"],
        'base': ["C6", "H6", "C5", "H5", "N3", "H3"],
    },

    'U': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H3'", "H4'", "H5'", "H5''"],
        'base': ["C6", "H6", "C7", "H71", "N3", "H3"],
    },

    'DU': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H2''", "H3'", "H4'", "H5'", "H5''"],
        'base': ["C6", "H6", "C7", "H71", "N3", "H3"],
    },

    'GLY': {
        'backbone': ['H', 'CA', 'HA2', 'HA3', 'N', 'C'],
    },

    'ALA': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB1', 'HB2', 'HB3'],
    },

    'VAL': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB', 'CG1', 'CG2', 'HG11', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23'],
    },

    'LEU': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG', 'HG', 'CD1', 'CD2', 'HD11', 'HD12', 'HD13', 'HD21', 'HD22', 'HD23'],
    },

    'ILE': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB', 'CG1', 'HG21', 'HG22', 'HG23', 'CD1', 'HD11', 'HD12', 'HD13', 'CG2', 'HG12',
                      'HG13'],
    },

    'MET': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CE', 'HE1', 'HE2', 'HE3'],
    },

    'PHE': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'aromatic': ['CD1', 'CD2', 'CE1', 'CE2', 'CZ','HD1', 'HD2', 'HE1', 'HE2', 'HZ'],
        'sidechain': ['CB', 'HB2', 'HB3'],
    },

    'TYR': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'aromatic': ['CD1', 'CD2', 'CE1', 'CE2', 'CZ','HD1', 'HD2', 'HE1', 'HE2'],
        'sidechain': ['CB', 'HB2', 'HB3'],
    },

    'TRP': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'aromatic': ['CD1', 'CE3', 'CZ2', 'CZ3', 'CH2', 'HD1', 'HE3', 'HZ2', 'HZ3', 'HH2', 'HE1', 'NE1'],
        'sidechain': ['CB', 'HB2', 'HB3'],
    },

    'CSE': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3'],
    },

    'CYS': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3'],
    },

    'PRO': {
        'backbone': ['HA', 'CA', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3'],
    },

    'SER': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3'],
    },

    'THR': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB', 'CG2', 'HG21', 'HG22', 'HG23'],
    },

    'ASN': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG', 'ND2', 'HD21', 'HD22'],
    },

    'GLN': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'NE2', 'HE21', 'HE22'],
    },

    'ASP': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG'],
    },

    'GLU': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD'],
    },

    'LYS': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'CE', 'HE2', 'HE3', 'NZ'],
    },

    'HIS': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'aromatic': ['CD2', 'HD2', 'ND1', 'HD1', 'NE2', 'HE2', 'CE1', 'HE1'],
        'sidechain': ['CB', 'HB2', 'HB3'],
    },

    'ARG': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        #			 'aromatic' : ['CD2','HD2','ND1','HD1','NE2','HE2','CE1','HE1'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'NE', 'HE', 'CZ', 'NH1', 'NH2',
                      'HH11', 'HH12', 'HH21', 'HH22'],
    },

}

ATOM_LIST = {}
for k in ATOM_DICT:
    ATOM_LIST[k] =[]
    for k1 in ATOM_DICT[k]:
        ATOM_LIST[k]+=ATOM_DICT[k][k1]

def read_starch_loop(fname,entry_id,list_id):
    cs_loop = pynmrstar.Loop.from_file(fname)
    cs_loop.add_missing_tags()
    cs_loop.delete_tag('_Atom_chem_shift.Entity_assembly_asym_ID')
    cs_loop.delete_tag('_Atom_chem_shift.Original_PDB_strand_ID')
    cs_loop.delete_tag('_Atom_chem_shift.Original_PDB_residue_no')
    cs_loop.delete_tag('_Atom_chem_shift.Original_PDB_residue_name')
    cs_loop.delete_tag('_Atom_chem_shift.Original_PDB_atom_name')
    columns = cs_loop.get_tag_names()
    comp_index_idx = columns.index('_Atom_chem_shift.Comp_index_ID')
    comp_idx = columns.index('_Atom_chem_shift.Comp_ID')
    atom_idx = columns.index('_Atom_chem_shift.Atom_ID')
    entry_idx = columns.index('_Atom_chem_shift.Entry_ID')
    list_idx = columns.index('_Atom_chem_shift.Assigned_chem_shift_list_ID')
    auth_seq_idx=columns.index('_Atom_chem_shift.Auth_seq_ID')
    auth_comp_idx=columns.index('_Atom_chem_shift.Auth_comp_ID')
    auth_atom_idx=columns.index('_Atom_chem_shift.Auth_atom_ID')
    j=0
    for row in cs_loop.data:
        cs_loop.data[j][comp_idx] = row[comp_idx].upper()
        if row[atom_idx] in ATOM_LIST[row[comp_idx]]:
            cs_loop.data[j][entry_idx]=entry_id
            cs_loop.data[j][list_idx]=list_id
            cs_loop.data[j][auth_seq_idx] = row[comp_index_idx]
            cs_loop.data[j][auth_comp_idx] = row[comp_idx]
            cs_loop.data[j][auth_atom_idx] = row[atom_idx]
        else:
            print (row)
        j+=1
    f=open(f'{fname.split(".str")[0]}_withalltags.str','w')
    f.write(str(cs_loop))
    f.close()
if __name__ == "__main__":
    fname = sys.argv[1]
    eid = sys.argv[2]
    lid = sys.argv[3]
    #read_starch_loop('/Users/kumaranbaskaran/bmr52506/work/data/bmr52487/work/data/test.str','52487','1')
    read_starch_loop(fname,eid,lid)