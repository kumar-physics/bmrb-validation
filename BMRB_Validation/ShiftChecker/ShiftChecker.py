import pynmrstar
import os
import csv

__version__ = '1.2'


def check_chemical_shifts(str_file=None, bmrb_id=None):
    print(bmrb_id, str_file)
    meta_data = {}
    starData = {}
    bmrb_stat = {}
    os.path.dirname(__file__) + '/lib/'
    bmrb_prot_stat = read_statistics_csv(os.path.dirname(__file__) + '/lib/bmrb_protein.csv')
    bmrb_dna_stat = read_statistics_csv(os.path.dirname(__file__) + '/lib/bmrb_dna.csv')
    bmrb_rna_stat = read_statistics_csv(os.path.dirname(__file__) + '/lib/bmrb_rna.csv')
    bmrb_stat.update(bmrb_prot_stat)
    bmrb_stat.update(bmrb_dna_stat)
    bmrb_stat.update(bmrb_rna_stat)
    if str_file is None and bmrb_id is None:
        raise ValueError('File name with full path or BMRB ID should be given')
    elif str_file is None:
        ent = pynmrstar.Entry.from_database(bmrb_id)
    else:
        ent = pynmrstar.Entry.from_file(str_file)
    entity_info = {}
    completeness = {}
    assmebly_info = {}
    assembly_size = ent.get_tag('_Assembly.Number_of_components')[0]
    assmebly_info['size'] = assembly_size
    assembly = list(zip(ent.get_tag('_Entity_assembly.ID'), ent.get_tag('_Entity_assembly.Entity_ID'),
                        ent.get_tag('_Entity_assembly.Asym_ID')))
    assmebly_info['assembly'] = assembly
    for entity_id, entity_seq, poly_type, entity_name, entity_type in zip(ent.get_tag('_Entity.ID'),
                                                                          ent.get_tag(
                                                                              '_Entity.Polymer_seq_one_letter_code'),
                                                                          ent.get_tag('_Entity.Polymer_type'),
                                                                          ent.get_tag('_Entity.Name'),
                                                                          ent.get_tag('_Entity.Type')):
        entity_info[entity_id] = (entity_type, poly_type, entity_seq.replace('\n', ''), entity_name)
        completeness[entity_id] = {
            'stereomethyl': [0, 0],
            'Total': {
                'overall': [0, 0],
                'backbone': [0, 0],
                'sidechain': [0, 0],
                'aromatic': [0, 0],
                'sugar': [0, 0],
                'base': [0, 0],
                'phosphate': [0, 0],
                'other': [0,0],
            },
            'H': {
                'overall': [0, 0],
                'backbone': [0, 0],
                'sidechain': [0, 0],
                'aromatic': [0, 0],
                'sugar': [0, 0],
                'base': [0, 0],
                'phosphate': [0, 0],
                'other': [0, 0],
            },
            'C': {
                'overall': [0, 0],
                'backbone': [0, 0],
                'sidechain': [0, 0],
                'aromatic': [0, 0],
                'sugar': [0, 0],
                'base': [0, 0],
                'phosphate': [0, 0],
                'other': [0, 0],
            },
            'N': {
                'overall': [0, 0],
                'backbone': [0, 0],
                'sidechain': [0, 0],
                'aromatic': [0, 0],
                'sugar': [0, 0],
                'base': [0, 0],
                'phosphate': [0, 0],
                'other': [0, 0],
            },
            'P': {
                'overall': [0, 0],
                'backbone': [0, 0],
                'sidechain': [0, 0],
                'aromatic': [0, 0],
                'sugar': [0, 0],
                'base': [0, 0],
                'phosphate': [0, 0],
                'other': [0, 0],
            },
        }
    streo_methyl = {}
    completeness_all = {}
    ent_list = ent.get_tag('_Entity.ID')
    con_seq_list = ent.get_loops_by_category('_Entity_comp_index')
    for entity in entity_info:
        if entity_info[entity][0] == 'polymer':
            streo_methyl[entity] = {}
            # for seq in entity_info[entity][2]:
            try:
                for seq in con_seq_list[ent_list.index(entity)].get_tag('Comp_ID'):
                    if seq in ['LEU', 'VAL']:
                        completeness[entity]['stereomethyl'][0] += 1
                    try:
                        for k in ATOM_DICT[seq]:
                            for atm in ATOM_DICT[seq][k]:
                                if 'H' in atm:
                                    completeness[entity]['H'][k][0] += 1
                                    completeness[entity]['H']['overall'][0] += 1
                                if 'C' in atm:
                                    completeness[entity]['C'][k][0] += 1
                                    completeness[entity]['C']['overall'][0] += 1
                                if 'N' in atm:
                                    completeness[entity]['N'][k][0] += 1
                                    completeness[entity]['N']['overall'][0] += 1
                                if 'P' in atm:
                                    completeness[entity]['P'][k][0] += 1
                                    completeness[entity]['P']['overall'][0] += 1
                                completeness[entity]['Total'][k][0] += 1
                                completeness[entity]['Total']['overall'][0] += 1
                    except KeyError as e:
                        print (e)
                        pass
            except IndexError:
                pass
    sample_cond_loops = ent.get_loops_by_category('Sample_condition_variable')
    sample_conditions = {}
    for sc in sample_cond_loops:
        sample_conditions[sc.get_tag('Sample_condition_list_ID')[0]] = {}
        for var, val, unt in zip(sc.get_tag('Type'), sc.get_tag('Val'), sc.get_tag('Val_units')):
            sample_conditions[sc.get_tag('Sample_condition_list_ID')[0]][var] = (val, unt)
    cs_saveframes = ent.get_saveframes_by_category('assigned_chemical_shifts')
    ent_saveframes = ent.get_saveframes_by_category('entry_information')
    meta_data['sample_conditions'] = sample_conditions
    meta_data['title'] = ent_saveframes[0].get_tag('title')[0]
    meta_data['ID'] = ent_saveframes[0].get_tag('ID')[0]
    meta_data['sub_date'] = ent_saveframes[0].get_tag('_Entry.Submission_date')[0]
    authors_loop = ent_saveframes[0].get_loop('Entry_author')
    authors_loop_tags = authors_loop.get_tag_names()
    first_name = authors_loop_tags.index('_Entry_author.Given_name')
    last_name = authors_loop_tags.index('_Entry_author.Family_name')
    middle_initital = authors_loop_tags.index('_Entry_author.Middle_initials')
    author_list = []
    sample_condition_map = {}
    for author in authors_loop.data:
        if author[middle_initital] != '.':
            author_list.append(f'{author[first_name]} {author[middle_initital]} {author[last_name]}')
        else:
            author_list.append(f'{author[first_name]} {author[last_name]}')
    meta_data['authors'] = author_list
    cs_dict = {}
    cs_error = {}
    cs_error['CS_VALUE'] = {}
    cs_error['CS_DUPLICATE'] = {}
    cs_error['CS_OUTLIER'] = {}
    row_count = {}

    for cs_saveframe in cs_saveframes:
        sf_name = cs_saveframe.get_tag('Sf_framecode')[0]
        sf_condition_id = cs_saveframe.get_tag('Sample_condition_list_ID')[0]
        try:
            cs_loop = cs_saveframe.get_loop('Atom_chem_shift')
            sf_list_id = cs_saveframe.get_tag('ID')[0]
            if sf_list_id not in completeness_all:
                completeness_all[sf_list_id] = {}
            sample_condition_map[sf_list_id] = sf_condition_id
            col_names = cs_loop.get_tag_names()
            title = (sf_name, sf_list_id)
            starData[title] = {}
            if sf_list_id not in row_count:
                row_count[sf_list_id] = 0
            if sf_list_id not in cs_error['CS_VALUE']:
                cs_error['CS_VALUE'][sf_list_id] = []
            if sf_list_id not in cs_error['CS_DUPLICATE']:
                cs_error['CS_DUPLICATE'][sf_list_id] = []
            if sf_list_id not in cs_error['CS_OUTLIER']:
                cs_error['CS_OUTLIER'][sf_list_id] = []
            cs_id = col_names.index('_Atom_chem_shift.ID')
            seq_id = col_names.index('_Atom_chem_shift.Comp_index_ID')
            asym_id = col_names.index('_Atom_chem_shift.Entity_assembly_ID')
            entity_id = col_names.index('_Atom_chem_shift.Entity_ID')
            comp_id = col_names.index('_Atom_chem_shift.Comp_ID')
            atom_id = col_names.index('_Atom_chem_shift.Atom_ID')
            cs_val_id = col_names.index('_Atom_chem_shift.Val')
            cs_err_id = col_names.index('_Atom_chem_shift.Val_err')
            amb_id = col_names.index('_Atom_chem_shift.Ambiguity_code')
            cs_list_id = col_names.index('_Atom_chem_shift.Assigned_chem_shift_list_ID')
            for cs in cs_loop.data:
                if cs[cs_list_id] not in completeness_all:
                    completeness_all[cs[cs_list_id]] = {}
                if cs[entity_id] not in completeness_all[cs[cs_list_id]]:
                    completeness_all[cs[cs_list_id]][cs[entity_id]] = {
                        'stereomethyl': [0, 0],
                        'Total': {
                            'overall': [0, 0],
                            'backbone': [0, 0],
                            'sidechain': [0, 0],
                            'aromatic': [0, 0],
                            'sugar': [0, 0],
                            'base': [0, 0],
                            'phosphate': [0, 0],
                        },
                        'H': {
                            'overall': [0, 0],
                            'backbone': [0, 0],
                            'sidechain': [0, 0],
                            'aromatic': [0, 0],
                            'sugar': [0, 0],
                            'base': [0, 0],
                            'phosphate': [0, 0],
                        },
                        'C': {
                            'overall': [0, 0],
                            'backbone': [0, 0],
                            'sidechain': [0, 0],
                            'aromatic': [0, 0],
                            'sugar': [0, 0],
                            'base': [0, 0],
                            'phosphate': [0, 0],
                        },
                        'N': {
                            'overall': [0, 0],
                            'backbone': [0, 0],
                            'sidechain': [0, 0],
                            'aromatic': [0, 0],
                            'sugar': [0, 0],
                            'base': [0, 0],
                            'phosphate': [0, 0],
                        },
                        'P': {
                            'overall': [0, 0],
                            'backbone': [0, 0],
                            'sidechain': [0, 0],
                            'aromatic': [0, 0],
                            'sugar': [0, 0],
                            'base': [0, 0],
                            'phosphate': [0, 0],
                        },
                    }
                row_count[sf_list_id] += 1
                if cs[cs_list_id] not in cs_dict:
                    cs_dict[cs[cs_list_id]] = {}
                if cs[asym_id] == '.':
                    eid = '1'  # Temp fix
                else:
                    eid = cs[entity_id]
                atom = (cs[seq_id], eid, cs[comp_id], cs[atom_id])
                if atom not in cs_dict[cs[cs_list_id]]:
                    try:
                        cs_dict[cs[cs_list_id]][atom] = (float(cs[cs_val_id]), cs[cs_err_id], cs[amb_id])
                        if cs[comp_id] in ['LEU', 'VAL']:
                            if cs[atom_id] in ['CG1', 'CG2', 'CD1', 'CD2']:
                                if (cs[seq_id], eid, cs[comp_id]) not in streo_methyl[eid]:
                                    streo_methyl[eid][(cs[seq_id], eid, cs[comp_id])] = []
                                streo_methyl[eid][(cs[seq_id], eid, cs[comp_id])].append(
                                    (float(cs[cs_val_id]), cs[cs_err_id], cs[amb_id]))
                        atm_type = get_atom_type(cs[comp_id], cs[atom_id])
                        try:
                            if cs[comp_id] in ATOM_DICT:
                                completeness_all[cs[cs_list_id]][eid]['Total']['overall'][1] += 1
                                completeness_all[cs[cs_list_id]][eid]['Total'][atm_type][1] += 1
                                if 'H' in cs[atom_id]:
                                    completeness_all[cs[cs_list_id]][eid]['H']['overall'][1] += 1
                                    completeness_all[cs[cs_list_id]][eid]['H'][atm_type][1] += 1
                                if 'C' in cs[atom_id]:
                                    completeness_all[cs[cs_list_id]][eid]['C']['overall'][1] += 1
                                    completeness_all[cs[cs_list_id]][eid]['C'][atm_type][1] += 1
                                if 'N' in cs[atom_id]:
                                    completeness_all[cs[cs_list_id]][eid]['N']['overall'][1] += 1
                                    completeness_all[cs[cs_list_id]][eid]['N'][atm_type][1] += 1
                                if 'P' in cs[atom_id]:
                                    completeness_all[cs[cs_list_id]][eid]['P']['overall'][1] += 1
                                    completeness_all[cs[cs_list_id]][eid]['P'][atm_type][1] += 1
                        except KeyError:
                            pass
                        try:
                            cs_mean = bmrb_stat[atom[2]][atom[3]][0]
                            cs_std = bmrb_stat[atom[2]][atom[3]][1]
                            if cs_mean is not None and cs_std is not None:
                                z = zscore(float(cs[cs_val_id]), cs_mean, cs_std)
                            else:
                                z = 0.0
                            starData[title][(eid, cs[seq_id], cs[comp_id], cs[atom_id])] = {
                                'value': float(cs[cs_val_id]),
                                'ambig': cs[amb_id],
                                'error': cs[cs_err_id]}
                            if abs(z) > 5.0:
                                cs_error['CS_OUTLIER'][sf_list_id].append(
                                    (atom[0], atom[1], atom[2], atom[3], cs[cs_list_id], cs[amb_id],
                                     cs[cs_val_id], z, round(cs_mean - 5.0 * cs_std, 2),
                                     round(cs_mean + 5.0 * cs_std, 2)))
                        except KeyError:
                            print("CS_TESTING-1", atom)
                    except ValueError:
                        cs_error['CS_VALUE'][sf_list_id].append(
                            (atom[0], atom[1], atom[2], atom[3], cs[cs_list_id], cs[cs_val_id], cs[cs_err_id],
                             cs[amb_id]))
                        # logging.WARNING(
                        #     '{}-{}-{}-{}-{} Can not read Chemical shift value'.format(atom[0], atom[1], atom[2], atom[3],
                        #                                                          cs[cs_list_id]))

                else:
                    cs_error['CS_DUPLICATE'][sf_list_id].append(
                        (atom[0], atom[1], atom[2], atom[3], cs[cs_list_id], cs[cs_val_id], cs[cs_err_id], cs[amb_id]))
                    # logging.WARNING('{}-{}-{}-{}-{} Duplicate chemical shifts found'.format(atom[0], atom[1], atom[2], atom[3], cs[cs_list_id]))
            meta_data['sample_conditions_map'] = sample_condition_map
            for ent1 in streo_methyl:
                for res1 in streo_methyl[ent1]:
                    if len(streo_methyl[ent1][res1]) == 2 and streo_methyl[ent1][res1][0][2] == '1' and \
                            streo_methyl[ent1][res1][1][2]:
                        try:
                            completeness_all[sf_list_id][ent1]['stereomethyl'][1] += 1
                        except KeyError:
                            pass
            for lid in completeness_all:
                for eid in completeness_all[lid]:
                    for t in completeness_all[lid][eid]:
                        if t == 'stereomethyl':
                            completeness_all[lid][eid][t][0] = completeness[eid][t][0]
                        else:
                            for st in completeness_all[lid][eid][t]:
                                completeness_all[lid][eid][t][st][0] = completeness[eid][t][st][0]
        except KeyError:
            pass
    result = {}
    cc = [0, 0]
    for lid in completeness_all:
        for aid in completeness_all[lid]:
            cc[0] += completeness_all[lid][aid]['Total']['overall'][0]
            cc[1] += completeness_all[lid][aid]['Total']['overall'][1]
    if cc[0] > 0:
        cc.append(round(float(cc[1]) / float(cc[0]) * 100.0, 1))
    else:
        cc.append(0.0)

    result['cs_dict'] = cs_dict
    result['cs_error'] = cs_error
    result['starData'] = starData
    result['meta_data'] = meta_data
    result['row_count'] = row_count
    result['completeness'] = completeness_all
    result['entity_info'] = entity_info
    result['assembly'] = assmebly_info
    result['overall_completeness'] = cc
    return result


def read_statistics_csv(fname):
    d = {}
    with open(fname) as csvfile:
        sp = csv.reader(csvfile, delimiter=',')
        header = next(sp)
        for row in sp:
            try:
                if row[header.index('comp_id')] not in d:
                    d[row[header.index('comp_id')]] = {}
                if row[header.index('atom_id')] == 'MG' and row[header.index('comp_id')] in ['ILE', 'THR']:
                    d[row[header.index('comp_id')]]['HG21'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                    d[row[header.index('comp_id')]]['HG22'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                    d[row[header.index('comp_id')]]['HG23'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                if row[header.index('atom_id')] == 'MD' and row[header.index('comp_id')] == 'ILE':
                    d[row[header.index('comp_id')]]['HD11'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                    d[row[header.index('comp_id')]]['HD12'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                    d[row[header.index('comp_id')]]['HD13'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                elif row[header.index('atom_id')] == 'MB':
                    d[row[header.index('comp_id')]]['HB1'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                    d[row[header.index('comp_id')]]['HB2'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                    d[row[header.index('comp_id')]]['HB3'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                elif row[header.index('atom_id')] == 'MD':
                    d[row[header.index('comp_id')]]['HD1'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                    d[row[header.index('comp_id')]]['HD2'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                    d[row[header.index('comp_id')]]['HD3'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                # elif row[header.index('atom_id')] == 'MG' and row[header.index('comp_id')] == 'THR':
                #     d[row[header.index('comp_id')]]['HG21'] = (
                #     float(row[header.index('avg')]), float(row[header.index('std')]))
                #     d[row[header.index('comp_id')]]['HG22'] = (
                #         float(row[header.index('avg')]), float(row[header.index('std')]))
                #     d[row[header.index('comp_id')]]['HG23'] = (
                #         float(row[header.index('avg')]), float(row[header.index('std')]))
                elif row[header.index('atom_id')] == 'MG':
                    d[row[header.index('comp_id')]]['HG1'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                    d[row[header.index('comp_id')]]['HG2'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                    d[row[header.index('comp_id')]]['HG3'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                elif row[header.index('atom_id')] == 'MD1':
                    d[row[header.index('comp_id')]]['HD11'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                    d[row[header.index('comp_id')]]['HD12'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                    d[row[header.index('comp_id')]]['HD13'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                elif row[header.index('atom_id')] == 'MD2':
                    d[row[header.index('comp_id')]]['HD21'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                    d[row[header.index('comp_id')]]['HD22'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                    d[row[header.index('comp_id')]]['HD23'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                elif row[header.index('atom_id')] == 'ME':
                    d[row[header.index('comp_id')]]['HE1'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                    d[row[header.index('comp_id')]]['HE2'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                    d[row[header.index('comp_id')]]['HE3'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                elif row[header.index('atom_id')] == 'MG1':
                    d[row[header.index('comp_id')]]['HG11'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                    d[row[header.index('comp_id')]]['HG12'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                    d[row[header.index('comp_id')]]['HG13'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                elif row[header.index('atom_id')] == 'MG2':
                    d[row[header.index('comp_id')]]['HG21'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                    d[row[header.index('comp_id')]]['HG22'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                    d[row[header.index('comp_id')]]['HG23'] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
                else:
                    d[row[header.index('comp_id')]][row[header.index('atom_id')]] = (
                        float(row[header.index('avg')]), float(row[header.index('std')]))
            except ValueError:
                d[row[header.index('comp_id')]][row[header.index('atom_id')]] = (None, None)
    return d


def zscore(x, m, s):
    return round(((x - m) / s), 2)


ATOM_DICT = {
    'A': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H3'", "H4'", "H5'", "H5''"],
        'base': ["C2", "H2", "C8", "H8", "N6", "H61", "H62", "N1", "H1"],
        'phosphate': ["P"],
        'other' : ["HO2'","C4","C6","N3","N7","N9"]
    },

    'DA': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H2''", "H3'", "H4'", "H5'", "H5''"],
        'base': ["C2", "H2", "C8", "H8", "N6", "H61", "H62", "N1", "H1"],
        'phosphate': ["P"],
        'other': ["C4","C5","C6","N3","N7","N9"]
    },

    'G': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H3'", "H4'", "H5'", "H5''"],
        'base': ["N1", "H1", "N2", "H21", "H22", "N3", "H3", "C8", "H8", "N7", "H7"],
        'phosphate': ["P"],
        'other': ["H2","H4","H5","H6","HO2","C2"]
    },

    'DG': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H2''", "H3'", "H4'", "H5'", "H5''"],
        'base': ["N1", "H1", "N2", "H21", "H22", "N3", "H3", "C8", "H8", "N7", "H7"],
        'phosphate': ["P"],
        'other': ["C2","C4","C5","C6","N9"]
    },

    'C': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H3'", "H4'", "H5'", "H5''"],
        'base': ["C6", "H6", "C5", "H5", "N4", "H41", "H42", "N3", "H3"],
        'phosphate': ['P'],
    },

    'DC': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H2''", "H3'", "H4'", "H5'", "H5''"],
        'base': ["C6", "H6", "C5", "H5", "N4", "H41", "H42", "N3", "H3"],
        'phosphate': ["P"],
        'other': ["C2","C4","N1"],
    },

    'T': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H3'", "H4'", "H5'", "H5''"],
        'base': ["C6", "H6", "C5", "H5", "N3", "H3"],
        'phosphate': ['P'],
    },

    'DT': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H2''", "H3'", "H4'", "H5'", "H5''"],
        'base': ["C6", "H6", "C5", "H5", "N3", "H3", "C7", "H71", "H72", "H73"],
        'phosphate': ["P"],
        'other': ["C2","C4","N1"],
    },

    'U': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H3'", "H4'", "H5'", "H5''"],
        'base': ["C6", "H6", "C7", "H71", "N3", "H3"],
        'phosphate': ["P"],
        'other': ["H5","HO2'","C2","C4","C5","N1"],
    },

    'DU': {
        'sugar': ["C1'", "C2'", "C3'", "C4'", "C5'", "H1'", "H2'", "H2''", "H3'", "H4'", "H5'", "H5''"],
        'base': ["C6", "H6", "C7", "H71", "N3", "H3"],
        'phosphate': ["P"],
        'other': ["H5","HO2'","C2","C4","C5","N1"],
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
        'sidechain': ['CB', 'HB', 'CG1', 'HG12', 'HG13', 'CD1', 'HD11', 'HD12', 'HD13', 'CG2', 'HG21', 'HG22',
                      'HG23', ],
    },

    'MET': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CE', 'HE1', 'HE2', 'HE3'],
    },

    'PHE': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'aromatic': ['CD1', 'CD2', 'CE1', 'CE2', 'HD1', 'HD2', 'HE1', 'HE2', 'HZ', 'CZ'],
        'sidechain': ['CB', 'HB2', 'HB3'],
    },

    'TYR': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'aromatic': ['CD1', 'CD2', 'CE1', 'CE2', 'HD1', 'HD2', 'HE1', 'HE2', 'CG', 'HH'],
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
        'sidechain': ['CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'NE', 'HE', 'CZ', 'NH1', 'NH2',
                      'HH11', 'HH12', 'HH21', 'HH22'],
    },

    'DGLY': {
        'backbone': ['H', 'CA', 'HA2', 'HA3', 'N', 'C'],
    },

    'DAL': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB1', 'HB2', 'HB3'],
    },

    'DVA': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB', 'CG1', 'CG2', 'HG11', 'HG12', 'HG13', 'HG21', 'HG22', 'HG23'],
    },

    'DLE': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG', 'HG', 'CD1', 'CD2', 'HD11', 'HD12', 'HD13', 'HD21', 'HD22', 'HD23'],
    },

    'DIL': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB', 'CG1', 'HG12', 'HG13', 'CD1', 'HD11', 'HD12', 'HD13', 'CG2', 'HG21', 'HG22',
                      'HG23', ],
    },

    'DME': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CE', 'HE1', 'HE2', 'HE3'],
    },

    'DPN': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'aromatic': ['CD1', 'CD2', 'CE1', 'CE2', 'HD1', 'HD2', 'HE1', 'HE2', 'HZ', 'CZ'],
        'sidechain': ['CB', 'HB2', 'HB3'],
    },

    'DTY': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'aromatic': ['CD1', 'CD2', 'CE1', 'CE2', 'HD1', 'HD2', 'HE1', 'HE2', 'CG', 'HH'],
        'sidechain': ['CB', 'HB2', 'HB3'],
    },

    'DTR': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'aromatic': ['CD1', 'CE3', 'CZ2', 'CZ3', 'CH2', 'HD1', 'HE3', 'HZ2', 'HZ3', 'HH2', 'HE1', 'NE1'],
        'sidechain': ['CB', 'HB2', 'HB3'],
    },

    'DCS': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3'],
    },

    'DCY': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3'],
    },

    'DPR': {
        'backbone': ['HA', 'CA', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3'],
    },

    'DSN': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3'],
    },

    'DTH': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB', 'CG2', 'HG21', 'HG22', 'HG23'],
    },

    'DSN': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG', 'ND2', 'HD21', 'HD22'],
    },

    'DGN': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'NE2', 'HE21', 'HE22'],
    },

    'DAS': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG'],
    },

    'DGL': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD'],
    },

    'DLY': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'CE', 'HE2', 'HE3', 'NZ'],
    },

    'DHI': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'aromatic': ['CD2', 'HD2', 'ND1', 'HD1', 'NE2', 'HE2', 'CE1', 'HE1'],
        'sidechain': ['CB', 'HB2', 'HB3'],
    },

    'DAR': {
        'backbone': ['H', 'HA', 'CA', 'N', 'C'],
        'sidechain': ['CB', 'HB2', 'HB3', 'CG', 'HG2', 'HG3', 'CD', 'HD2', 'HD3', 'NE', 'HE', 'CZ', 'NH1', 'NH2',
                      'HH11', 'HH12', 'HH21', 'HH22'],
    },

}


def get_atom_type(res, atm):
    try:
        for k in ATOM_DICT[res]:
            if atm in ATOM_DICT[res][k]:
                return k
    except KeyError:
        pass
    return None


TLC = {
    'A': 'ALA',
    'G': 'GLY',
    'S': 'SER',
    'T': 'THR',
    'C': 'CYS',
    'V': 'VAL',
    'L': 'LEU',
    'I': 'ILE',
    'M': 'MET',
    'P': 'PRO',
    'F': 'PHE',
    'Y': 'TYR',
    'W': 'TRP',
    'D': 'ASP',
    'E': 'GLU',
    'N': 'ASN',
    'Q': 'GLN',
    'H': 'HIS',
    'K': 'LYS',
    'R': 'ARG',
}

OLC = {v: k for k, v in TLC.items()}

if __name__ == "__main__":
    c = ['Total', 'H', 'C', 'N']
    sc = ['overall', 'backbone', 'sidechain', 'aromatic', 'sugar', 'base']
    x = check_chemical_shifts('./bmr52130_3.str', '52130')
    out = '52188'
    for l in x['completeness']:
        out = f'{out},{l}'
        for e in x['completeness'][l]:
            out = f'{out},{e}'
            for c1 in c:
                for sc1 in sc:
                    if x["completeness"][l][e][c1][sc1][0] > 0:
                        out = f'{out},{x["completeness"][l][e][c1][sc1][0]},{x["completeness"][l][e][c1][sc1][1]},{round(100.0 * (float(x["completeness"][l][e][c1][sc1][1]) / float(x["completeness"][l][e][c1][sc1][0])), 2)}'
                    else:
                        out = f'{out},{x["completeness"][l][e][c1][sc1][0]},{x["completeness"][l][e][c1][sc1][1]},0.00'

            print(out)
    entries = ["3", "5", "6", "7", "8", "9", "18", "19", "21", "23", "24", "25", "26", "27", "29", "30", "32", "36",
               "39", "40", "41", "42", "44", "45", "46", "48", "49", "50", "53", "55", "58", "60", "62", "63", "64",
               "65", "66", "68", "76", "79", "80", "81", "86", "87", "88", "90", "94", "95", "96", "97", "102", "103",
               "104", "105", "106", "107", "108", "109", "112", "114", "115", "116", "117", "120", "121", "122", "123",
               "126", "127", "131", "132", "133", "136", "139", "140", "141", "142", "144", "146", "148", "149", "150",
               "151", "161", "162", "169", "171", "180", "181", "182", "183", "185", "186", "187", "188", "189", "190",
               "191", "192", "193", "194", "195", "196", "199", "200", "216", "217", "218", "219", "220", "221", "222",
               "223", "224", "225", "226", "227", "228", "229", "233", "234", "235", "236", "237", "238", "239", "240",
               "241", "242", "243", "244", "245", "246", "247", "248", "249", "250", "252", "253", "254", "255", "256",
               "257", "258", "259", "260", "261", "262", "263", "264", "274", "275", "276", "280", "283", "284", "285",
               "286", "287", "288", "289", "290", "291", "292", "293", "294", "295", "298", "300", "301", "302", "306",
               "307", "308", "309", "311", "314", "316", "317", "321", "325", "326", "327", "329", "330", "331", "335",
               "336", "337", "338", "339", "340", "341", "342", "345", "346", "349", "350", "351", "353", "354", "355",
               "356", "357", "358", "359", "360", "361", "362", "363", "364", "367", "369", "371", "372", "373", "374",
               "375", "376", "385", "386", "387", "388", "389", "390", "391", "392", "393", "394", "395", "397", "398",
               "401", "402", "405", "407", "408", "410", "411", "412", "413", "414", "415", "416", "417", "418", "419",
               "420", "421", "422", "423", "426", "428", "429", "430", "432", "433", "434", "435", "436", "437", "438",
               "439", "440", "441", "442", "443", "444", "447", "448", "449", "450", "451", "452", "453", "454", "455",
               "456", "457", "458", "459", "460", "461", "462", "463", "464", "465", "466", "467", "468", "469", "470",
               "471", "472", "473", "474", "475", "476", "477", "478", "479", "480", "481", "485", "486", "488", "489",
               "490", "491", "494", "495", "496", "497", "498", "499", "500", "505", "506", "507", "508", "509", "510",
               "511", "512", "513", "514", "515", "516", "517", "530", "531", "532", "537", "538", "539", "540", "541",
               "543", "544", "545", "546", "547", "548", "554", "555", "556", "557", "558", "559", "560", "561", "562",
               "563", "568", "569", "573", "595", "596", "597", "598", "599", "600", "606", "611", "613", "614", "623",
               "630", "644", "645", "653", "660", "661", "662", "664", "665", "666", "667", "668", "672", "673", "674",
               "676", "677", "682", "695", "699", "700", "704", "705", "706", "707", "736", "749", "750", "753", "754",
               "758", "759", "760", "761", "779", "780", "781", "791", "810", "811", "812", "813", "814", "821", "822",
               "823", "824", "846", "848", "849", "862", "877", "878", "883", "884", "887", "888", "893", "907", "908",
               "909", "910", "911", "912", "913", "914", "915", "916", "918", "919", "920", "922", "923", "932", "934",
               "936", "937", "939", "940", "944", "945", "946", "947", "948", "949", "956", "957", "960", "961", "962",
               "963", "964", "965", "966", "967", "968", "969", "970", "971", "975", "979", "994", "995", "996", "997",
               "998", "999", "1000", "1001", "1002", "1003", "1004", "1005", "1006", "1007", "1008", "1009", "1010",
               "1011", "1012", "1013", "1014", "1015", "1016", "1017", "1018", "1019", "1020", "1021", "1022", "1023",
               "1024", "1025", "1027", "1028", "1029", "1030", "1031", "1032", "1033", "1037", "1039", "1042", "1043",
               "1044", "1045", "1046", "1049", "1053", "1056", "1058", "1061", "1062", "1065", "1066", "1071", "1072",
               "1093", "1095", "1096", "1097", "1098", "1099", "1100", "1101", "1102", "1104", "1105", "1107", "1108",
               "1109", "1110", "1111", "1112", "1113", "1114", "1115", "1116", "1117", "1120", "1127", "1128", "1129",
               "1133", "1134", "1135", "1136", "1137", "1138", "1139", "1140", "1144", "1145", "1147", "1148", "1149",
               "1150", "1151", "1152", "1153", "1154", "1155", "1156", "1158", "1159", "1160", "1161", "1162", "1163",
               "1164", "1165", "1166", "1167", "1168", "1170", "1171", "1174", "1175", "1176", "1177", "1178", "1179",
               "1181", "1182", "1183", "1184", "1192", "1194", "1198", "1200", "1201", "1202", "1203", "1204", "1205",
               "1206", "1207", "1208", "1209", "1210", "1211", "1212", "1213", "1214", "1215", "1216", "1217", "1218",
               "1219", "1220", "1221", "1222", "1225", "1226", "1227", "1228", "1230", "1323", "1324", "1333", "1334",
               "1336", "1337", "1338", "1339", "1343", "1344", "1346", "1347", "1348", "1358", "1360", "1374", "1375",
               "1376", "1377", "1378", "1379", "1381", "1393", "1394", "1396", "1397", "1398", "1399", "1400", "1403",
               "1404", "1413", "1414", "1416", "1439", "1440", "1441", "1442", "1443", "1444", "1445", "1446", "1447",
               "1448", "1449", "1450", "1451", "1452", "1453", "1454", "1455", "1457", "1459", "1461", "1463", "1465",
               "1467", "1469", "1471", "1474", "1475", "1476", "1477", "1478", "1479", "1480", "1481", "1482", "1483",
               "1484", "1485", "1486", "1487", "1494", "1495", "1496", "1497", "1500", "1503", "1515", "1516", "1517",
               "1518", "1519", "1520", "1538", "1539", "1541", "1542", "1543", "1544", "1545", "1546", "1551", "1552",
               "1553", "1562", "1566", "1567", "1571", "1573", "1574", "1575", "1576", "1577", "1578", "1580", "1581",
               "1582", "1583", "1584", "1585", "1586", "1587", "1588", "1589", "1590", "1591", "1592", "1616", "1620",
               "1624", "1632", "1633", "1634", "1635", "1638", "1639", "1640", "1642", "1643", "1646", "1647", "1648",
               "1649", "1650", "1651", "1652", "1653", "1654", "1655", "1656", "1657", "1658", "1660", "1661", "1662",
               "1663", "1664", "1665", "1666", "1669", "1670", "1672", "1673", "1674", "1675", "1676", "1677", "1692",
               "1693", "1696", "1697", "1698", "1699", "1700", "1702", "1704", "1705", "1706", "1707", "1709", "1710",
               "1712", "1713", "1719", "1720", "1722", "1723", "1724", "1725", "1726", "1727", "1728", "1729", "1734",
               "1736", "1739", "1743", "1744", "1745", "1751", "1752", "1756", "1757", "1760", "1761", "1762", "1766",
               "1767", "1769", "1770", "1771", "1772", "1773", "1775", "1777", "1779", "1781", "1783", "1785", "1787",
               "1789", "1791", "1793", "1795", "1797", "1798", "1799", "1800", "1801", "1802", "1803", "1808", "1811",
               "1812", "1813", "1814", "1815", "1816", "1817", "1819", "1820", "1821", "1822", "1825", "1826", "1839",
               "1841", "1843", "1844", "1845", "1846", "1851", "1854", "1855", "1856", "1857", "1858", "1859", "1865",
               "1866", "1867", "1869", "1870", "1871", "1872", "1874", "1875", "1876", "1877", "1878", "1881", "1882",
               "1883", "1884", "1885", "1886", "1887", "1889", "1891", "1892", "1893", "1895", "1911", "1917", "1918",
               "1919", "1920", "1921", "1968", "1975", "1977", "1991", "1995", "1996", "1997", "1998", "1999", "2000",
               "2001", "2002", "2003", "2006", "2012", "2013", "2014", "2024", "2030", "2031", "2038", "2039", "2040",
               "2042", "2043", "2047", "2048", "2049", "2050", "2051", "2052", "2053", "2059", "2060", "2061", "2062",
               "2063", "2065", "2066", "2067", "2068", "2069", "2070", "2071", "2074", "2075", "2076", "2077", "2078",
               "2113", "2114", "2115", "2116", "2117", "2118", "2151", "2169", "2173", "2192", "2193", "2196", "2197",
               "2198", "2199", "2200", "2201", "2202", "2203", "2204", "2205", "2206", "2207", "2208", "2209", "2217",
               "2218", "2219", "2220", "2221", "2222", "2223", "2224", "2225", "2226", "2227", "2228", "2229", "2231",
               "2232", "2233", "2234", "2235", "2236", "2247", "2249", "2261", "2262", "2263", "2265", "2266", "2267",
               "2268", "2269", "2270", "2271", "2272", "2274", "2275", "2276", "2277", "2278", "2279", "2281", "2282",
               "2283", "2284", "2285", "2288", "2324", "2325", "2326", "2327", "2328", "2329", "2331", "2338", "2345",
               "2346", "2347", "2348", "2352", "2353", "2366", "2367", "2368", "2371", "2384", "2395", "2396", "2410",
               "2425", "2426", "2431", "2432", "2433", "2434", "2435", "2436", "2437", "2438", "2439", "2440", "2442",
               "2443", "2446", "2454", "2455", "2457", "2463", "2472", "2473", "2474", "2475", "2476", "2488", "2489",
               "2498", "2506", "2512", "2513", "2527", "2528", "2529", "2539", "2542", "2546", "2547", "2573", "2574",
               "2575", "2580", "2590", "2591", "2592", "2606", "2607", "2608", "2609", "2610", "2611", "2612", "2613",
               "2614", "2615", "2616", "2617", "2618", "2619", "2620", "2621", "2622", "2623", "2624", "2625", "2627",
               "2628", "2630", "2631", "2632", "2638", "2639", "2640", "2641", "2642", "2651", "2652", "2707", "2708",
               "2709", "2710", "2711", "2712", "2718", "2719", "2763", "2764", "2768", "2784", "2785", "2786", "2787",
               "2798", "2799", "2802", "2804", "2852", "2853", "2856", "2858", "2868", "2917", "2928", "2935", "2936",
               "2939", "2940", "2941", "2948", "2950", "2951", "2956", "2957", "2958", "2959", "2960", "2968", "2990",
               "2991", "2999", "3000", "3032", "3036", "3065", "3078", "3079", "3084", "3085", "3118", "3119", "3127",
               "3128", "3322", "3394", "3404", "3427", "3433", "3435", "3436", "3437", "3440", "3441", "3442", "3449",
               "3456", "3466", "3471", "3472", "3485", "3524", "3525", "3548", "4001", "4010", "4011", "4012", "4019",
               "4020", "4022", "4023", "4024", "4027", "4028", "4029", "4030", "4031", "4032", "4033", "4034", "4035",
               "4036", "4037", "4038", "4039", "4040", "4041", "4042", "4043", "4044", "4045", "4046", "4047", "4048",
               "4049", "4050", "4051", "4052", "4053", "4054", "4055", "4056", "4057", "4058", "4059", "4060", "4061",
               "4062", "4063", "4064", "4065", "4066", "4067", "4068", "4069", "4070", "4071", "4072", "4073", "4074",
               "4075", "4076", "4078", "4079", "4080", "4081", "4082", "4083", "4084", "4085", "4086", "4087", "4088",
               "4089", "4090", "4091", "4092", "4093", "4094", "4095", "4097", "4098", "4099", "4100", "4101", "4102",
               "4103", "4104", "4105", "4106", "4107", "4108", "4109", "4110", "4111", "4112", "4113", "4114", "4115",
               "4116", "4117", "4120", "4121", "4122", "4123", "4124", "4125", "4126", "4127", "4128", "4129", "4130",
               "4131", "4132", "4133", "4134", "4135", "4136", "4137", "4138", "4140", "4141", "4142", "4143", "4144",
               "4145", "4146", "4147", "4148", "4149", "4150", "4151", "4152", "4153", "4154", "4155", "4156", "4157",
               "4158", "4159", "4160", "4161", "4162", "4163", "4164", "4165", "4166", "4167", "4168", "4169", "4170",
               "4171", "4172", "4173", "4174", "4175", "4176", "4177", "4179", "4180", "4181", "4182", "4183", "4184",
               "4185", "4186", "4187", "4188", "4189", "4190", "4191", "4192", "4193", "4194", "4195", "4197", "4198",
               "4199", "4200", "4201", "4202", "4203", "4204", "4205", "4206", "4207", "4208", "4209", "4210", "4211",
               "4212", "4213", "4214", "4215", "4216", "4217", "4218", "4219", "4220", "4221", "4222", "4223", "4224",
               "4225", "4226", "4227", "4228", "4229", "4230", "4231", "4232", "4233", "4234", "4235", "4236", "4237",
               "4239", "4240", "4241", "4242", "4243", "4244", "4245", "4246", "4247", "4248", "4249", "4250", "4251",
               "4252", "4253", "4254", "4255", "4256", "4257", "4258", "4259", "4260", "4261", "4262", "4263", "4264",
               "4265", "4266", "4267", "4268", "4269", "4270", "4271", "4272", "4273", "4276", "4278", "4279", "4280",
               "4281", "4282", "4283", "4284", "4285", "4286", "4287", "4288", "4289", "4290", "4291", "4292", "4293",
               "4294", "4295", "4296", "4297", "4298", "4299", "4300", "4301", "4302", "4303", "4304", "4305", "4306",
               "4307", "4308", "4309", "4310", "4311", "4312", "4313", "4314", "4315", "4316", "4317", "4318", "4319",
               "4320", "4321", "4322", "4323", "4324", "4325", "4326", "4327", "4328", "4329", "4330", "4331", "4332",
               "4333", "4334", "4335", "4336", "4337", "4338", "4339", "4340", "4341", "4342", "4343", "4344", "4345",
               "4346", "4347", "4348", "4349", "4350", "4351", "4352", "4353", "4354", "4355", "4357", "4358", "4359",
               "4360", "4361", "4362", "4363", "4364", "4365", "4366", "4367", "4368", "4369", "4370", "4371", "4372",
               "4373", "4374", "4375", "4376", "4377", "4378", "4379", "4380", "4381", "4382", "4383", "4384", "4385",
               "4386", "4387", "4388", "4389", "4390", "4391", "4392", "4393", "4394", "4395", "4396", "4397", "4398",
               "4399", "4400", "4401", "4402", "4403", "4404", "4405", "4406", "4407", "4408", "4409", "4410", "4411",
               "4412", "4413", "4414", "4415", "4416", "4417", "4418", "4419", "4420", "4421", "4422", "4423", "4424",
               "4425", "4426", "4427", "4428", "4429", "4430", "4431", "4432", "4433", "4434", "4435", "4436", "4437",
               "4438", "4439", "4440", "4441", "4442", "4443", "4444", "4445", "4446", "4447", "4448", "4449", "4451",
               "4452", "4453", "4454", "4455", "4456", "4457", "4458", "4459", "4460", "4461", "4462", "4463", "4464",
               "4465", "4466", "4467", "4468", "4469", "4470", "4471", "4472", "4473", "4475", "4476", "4477", "4478",
               "4483", "4486", "4487", "4488", "4490", "4491", "4492", "4493", "4494", "4496", "4497", "4498", "4500",
               "4503", "4506", "4507", "4509", "4510", "4514", "4516", "4519", "4522", "4524", "4526", "4527", "4528",
               "4536", "4540", "4541", "4542", "4547", "4550", "4551", "4552", "4553", "4554", "4555", "4556", "4557",
               "4558", "4559", "4560", "4561", "4562", "4563", "4564", "4565", "4566", "4567", "4568", "4570", "4571",
               "4572", "4573", "4574", "4575", "4576", "4577", "4578", "4579", "4580", "4581", "4583", "4584", "4585",
               "4587", "4588", "4589", "4590", "4591", "4592", "4593", "4595", "4596", "4597", "4598", "4599", "4600",
               "4601", "4602", "4603", "4604", "4607", "4609", "4610", "4612", "4614", "4615", "4616", "4617", "4618",
               "4619", "4620", "4621", "4622", "4623", "4628", "4630", "4631", "4632", "4633", "4634", "4635", "4636",
               "4637", "4638", "4639", "4640", "4641", "4642", "4643", "4644", "4645", "4646", "4647", "4648", "4649",
               "4650", "4651", "4652", "4653", "4654", "4655", "4656", "4661", "4663", "4664", "4666", "4668", "4669",
               "4670", "4671", "4673", "4674", "4675", "4676", "4677", "4678", "4679", "4680", "4681", "4682", "4683",
               "4684", "4685", "4686", "4687", "4688", "4689", "4692", "4694", "4695", "4696", "4697", "4698", "4699",
               "4700", "4701", "4702", "4703", "4704", "4706", "4707", "4708", "4709", "4710", "4711", "4712", "4713",
               "4714", "4715", "4716", "4717", "4718", "4719", "4720", "4721", "4722", "4723", "4724", "4725", "4726",
               "4727", "4728", "4729", "4731", "4732", "4733", "4734", "4735", "4736", "4737", "4738", "4739", "4740",
               "4741", "4742", "4743", "4744", "4745", "4746", "4747", "4748", "4749", "4750", "4751", "4752", "4753",
               "4754", "4755", "4757", "4758", "4759", "4760", "4761", "4762", "4763", "4765", "4766", "4767", "4768",
               "4769", "4770", "4771", "4772", "4773", "4774", "4775", "4777", "4778", "4779", "4780", "4781", "4782",
               "4784", "4785", "4786", "4787", "4788", "4789", "4790", "4791", "4792", "4793", "4794", "4795", "4796",
               "4797", "4798", "4799", "4800", "4801", "4802", "4803", "4804", "4805", "4806", "4807", "4808", "4809",
               "4810", "4811", "4812", "4813", "4814", "4815", "4816", "4817", "4818", "4819", "4820", "4821", "4822",
               "4823", "4824", "4825", "4827", "4828", "4829", "4830", "4831", "4832", "4833", "4834", "4835", "4836",
               "4837", "4838", "4839", "4840", "4841", "4842", "4843", "4844", "4845", "4846", "4847", "4848", "4849",
               "4850", "4851", "4852", "4853", "4854", "4855", "4856", "4857", "4858", "4859", "4860", "4862", "4863",
               "4864", "4865", "4866", "4867", "4868", "4869", "4870", "4871", "4872", "4873", "4874", "4875", "4876",
               "4877", "4878", "4879", "4880", "4881", "4882", "4883", "4884", "4885", "4886", "4887", "4888", "4889",
               "4890", "4891", "4892", "4893", "4894", "4895", "4896", "4897", "4898", "4899", "4900", "4901", "4902",
               "4905", "4906", "4907", "4908", "4909", "4910", "4911", "4912", "4913", "4914", "4915", "4916", "4917",
               "4918", "4919", "4920", "4921", "4922", "4923", "4924", "4925", "4926", "4927", "4928", "4929", "4930",
               "4931", "4932", "4933", "4934", "4935", "4936", "4937", "4938", "4939", "4940", "4941", "4942", "4943",
               "4944", "4945", "4946", "4947", "4948", "4950", "4951", "4952", "4953", "4954", "4955", "4956", "4957",
               "4958", "4959", "4960", "4961", "4963", "4964", "4965", "4966", "4967", "4968", "4969", "4970", "4971",
               "4972", "4973", "4974", "4975", "4976", "4977", "4978", "4979", "4980", "4981", "4982", "4983", "4984",
               "4985", "4986", "4987", "4988", "4989", "4990", "4991", "4992", "4993", "4994", "4995", "4996", "4997",
               "4998", "4999", "5000", "5003", "5004", "5005", "5006", "5007", "5008", "5009", "5010", "5011", "5012",
               "5013", "5014", "5018", "5019", "5020", "5021", "5022", "5023", "5024", "5025", "5026", "5027", "5028",
               "5029", "5030", "5031", "5032", "5033", "5034", "5036", "5037", "5038", "5039", "5040", "5041", "5042",
               "5043", "5044", "5046", "5047", "5048", "5049", "5050", "5051", "5052", "5053", "5054", "5055", "5056",
               "5057", "5058", "5059", "5060", "5061", "5062", "5064", "5065", "5066", "5067", "5068", "5069", "5070",
               "5071", "5072", "5073", "5074", "5075", "5076", "5077", "5078", "5079", "5080", "5081", "5082", "5083",
               "5084", "5085", "5086", "5087", "5088", "5089", "5090", "5091", "5092", "5093", "5094", "5096", "5097",
               "5098", "5099", "5100", "5101", "5102", "5103", "5104", "5105", "5106", "5107", "5108", "5109", "5110",
               "5111", "5112", "5113", "5114", "5115", "5116", "5117", "5119", "5120", "5121", "5122", "5123", "5124",
               "5125", "5126", "5127", "5128", "5129", "5130", "5131", "5132", "5133", "5134", "5135", "5136", "5137",
               "5138", "5139", "5140", "5141", "5142", "5143", "5144", "5145", "5147", "5148", "5149", "5150", "5151",
               "5152", "5153", "5154", "5155", "5156", "5157", "5158", "5159", "5160", "5161", "5162", "5163", "5164",
               "5165", "5166", "5167", "5169", "5170", "5171", "5172", "5173", "5174", "5175", "5176", "5177", "5178",
               "5179", "5180", "5181", "5182", "5183", "5184", "5185", "5186", "5187", "5188", "5189", "5190", "5191",
               "5192", "5193", "5194", "5196", "5197", "5198", "5199", "5200", "5201", "5202", "5203", "5204", "5205",
               "5206", "5207", "5208", "5209", "5210", "5211", "5212", "5213", "5214", "5215", "5216", "5217", "5218",
               "5219", "5220", "5221", "5222", "5223", "5224", "5225", "5226", "5227", "5228", "5231", "5232", "5233",
               "5234", "5235", "5236", "5237", "5238", "5239", "5240", "5241", "5242", "5243", "5244", "5245", "5246",
               "5247", "5248", "5249", "5250", "5251", "5252", "5253", "5254", "5255", "5256", "5257", "5259", "5260",
               "5261", "5262", "5263", "5264", "5265", "5266", "5267", "5268", "5269", "5270", "5271", "5272", "5273",
               "5274", "5275", "5276", "5277", "5278", "5279", "5280", "5281", "5282", "5283", "5284", "5285", "5286",
               "5287", "5288", "5289", "5290", "5291", "5292", "5293", "5294", "5295", "5296", "5297", "5298", "5299",
               "5300", "5301", "5302", "5303", "5304", "5305", "5306", "5307", "5308", "5309", "5310", "5311", "5312",
               "5313", "5314", "5315", "5316", "5317", "5318", "5319", "5320", "5321", "5322", "5323", "5324", "5325",
               "5326", "5327", "5328", "5329", "5330", "5331", "5332", "5333", "5334", "5335", "5337", "5338", "5339",
               "5340", "5341", "5342", "5343", "5344", "5345", "5346", "5347", "5348", "5349", "5350", "5351", "5352",
               "5353", "5354", "5355", "5356", "5357", "5358", "5359", "5360", "5361", "5362", "5363", "5364", "5365",
               "5366", "5368", "5369", "5370", "5371", "5372", "5373", "5374", "5375", "5376", "5377", "5378", "5379",
               "5381", "5382", "5383", "5384", "5385", "5386", "5387", "5388", "5389", "5390", "5391", "5392", "5393",
               "5394", "5395", "5396", "5397", "5398", "5399", "5400", "5401", "5402", "5403", "5404", "5405", "5406",
               "5407", "5408", "5409", "5410", "5411", "5412", "5413", "5414", "5415", "5416", "5417", "5418", "5419",
               "5420", "5421", "5422", "5423", "5424", "5425", "5426", "5427", "5428", "5429", "5430", "5431", "5432",
               "5433", "5434", "5435", "5436", "5437", "5438", "5439", "5440", "5441", "5442", "5443", "5444", "5445",
               "5446", "5447", "5448", "5449", "5450", "5451", "5452", "5453", "5454", "5455", "5456", "5457", "5458",
               "5459", "5461", "5462", "5463", "5464", "5465", "5466", "5467", "5468", "5469", "5470", "5471", "5472",
               "5473", "5474", "5475", "5476", "5477", "5478", "5479", "5480", "5481", "5482", "5483", "5484", "5485",
               "5486", "5487", "5488", "5489", "5490", "5491", "5492", "5493", "5494", "5495", "5496", "5497", "5498",
               "5499", "5500", "5501", "5502", "5503", "5504", "5505", "5506", "5507", "5508", "5509", "5510", "5511",
               "5512", "5513", "5514", "5515", "5516", "5517", "5518", "5519", "5520", "5521", "5522", "5523", "5524",
               "5525", "5527", "5528", "5529", "5530", "5531", "5532", "5534", "5535", "5536", "5537", "5538", "5539",
               "5540", "5541", "5542", "5543", "5544", "5545", "5547", "5548", "5549", "5550", "5551", "5552", "5553",
               "5554", "5555", "5556", "5557", "5558", "5559", "5560", "5561", "5562", "5563", "5564", "5565", "5566",
               "5567", "5568", "5569", "5570", "5571", "5572", "5573", "5574", "5575", "5576", "5577", "5578", "5579",
               "5580", "5581", "5582", "5583", "5584", "5585", "5586", "5587", "5588", "5589", "5590", "5591", "5592",
               "5593", "5594", "5595", "5596", "5597", "5598", "5599", "5600", "5601", "5602", "5603", "5604", "5605",
               "5606", "5607", "5608", "5609", "5610", "5611", "5612", "5613", "5614", "5615", "5616", "5617", "5618",
               "5619", "5620", "5621", "5622", "5623", "5624", "5625", "5626", "5627", "5628", "5629", "5630", "5631",
               "5632", "5633", "5634", "5635", "5636", "5637", "5638", "5639", "5640", "5641", "5642", "5643", "5644",
               "5645", "5646", "5647", "5648", "5649", "5650", "5651", "5652", "5653", "5654", "5655", "5656", "5657",
               "5658", "5659", "5660", "5661", "5662", "5663", "5664", "5665", "5666", "5667", "5668", "5669", "5670",
               "5671", "5672", "5673", "5674", "5675", "5676", "5677", "5678", "5679", "5680", "5681", "5682", "5683",
               "5684", "5685", "5686", "5687", "5688", "5689", "5690", "5691", "5692", "5693", "5694", "5695", "5696",
               "5697", "5698", "5699", "5700", "5701", "5702", "5703", "5704", "5705", "5706", "5707", "5708", "5709",
               "5710", "5711", "5712", "5713", "5714", "5715", "5716", "5717", "5718", "5719", "5720", "5721", "5722",
               "5723", "5724", "5725", "5726", "5727", "5728", "5729", "5730", "5731", "5732", "5733", "5734", "5735",
               "5736", "5737", "5738", "5739", "5740", "5741", "5742", "5743", "5744", "5745", "5746", "5747", "5748",
               "5749", "5750", "5751", "5752", "5753", "5754", "5755", "5756", "5757", "5758", "5759", "5760", "5761",
               "5762", "5763", "5764", "5765", "5766", "5767", "5768", "5769", "5770", "5771", "5772", "5773", "5774",
               "5775", "5776", "5777", "5778", "5779", "5780", "5781", "5782", "5783", "5784", "5785", "5786", "5787",
               "5788", "5789", "5790", "5791", "5792", "5793", "5794", "5795", "5796", "5797", "5798", "5799", "5800",
               "5801", "5802", "5803", "5804", "5805", "5806", "5807", "5808", "5809", "5810", "5811", "5812", "5813",
               "5814", "5815", "5816", "5817", "5818", "5819", "5820", "5821", "5822", "5823", "5824", "5825", "5826",
               "5827", "5828", "5829", "5830", "5832", "5833", "5834", "5835", "5836", "5837", "5838", "5839", "5840",
               "5841", "5842", "5843", "5844", "5845", "5846", "5847", "5848", "5849", "5850", "5851", "5852", "5853",
               "5854", "5856", "5857", "5858", "5859", "5860", "5861", "5862", "5863", "5864", "5865", "5866", "5867",
               "5868", "5869", "5870", "5871", "5872", "5873", "5874", "5875", "5876", "5877", "5878", "5879", "5880",
               "5881", "5882", "5883", "5884", "5885", "5886", "5887", "5888", "5889", "5890", "5891", "5892", "5893",
               "5894", "5895", "5896", "5897", "5898", "5899", "5900", "5901", "5902", "5903", "5904", "5905", "5906",
               "5907", "5908", "5909", "5910", "5911", "5912", "5913", "5914", "5915", "5916", "5917", "5918", "5919",
               "5920", "5921", "5922", "5923", "5924", "5925", "5926", "5927", "5928", "5929", "5930", "5931", "5932",
               "5933", "5934", "5935", "5936", "5937", "5938", "5939", "5940", "5941", "5942", "5943", "5944", "5945",
               "5946", "5947", "5948", "5949", "5950", "5951", "5952", "5953", "5954", "5955", "5956", "5958", "5959",
               "5960", "5961", "5962", "5963", "5964", "5965", "5966", "5967", "5968", "5969", "5970", "5971", "5972",
               "5973", "5974", "5975", "5976", "5977", "5978", "5979", "5980", "5981", "5982", "5983", "5984", "5985",
               "5986", "5987", "5988", "5989", "5991", "5992", "5993", "5994", "5995", "5996", "5997", "5998", "5999",
               "6000", "6001", "6002", "6003", "6004", "6005", "6006", "6007", "6008", "6009", "6010", "6011", "6012",
               "6013", "6014", "6015", "6016", "6017", "6018", "6019", "6020", "6021", "6022", "6023", "6024", "6025",
               "6026", "6027", "6028", "6029", "6030", "6031", "6032", "6033", "6034", "6035", "6037", "6038", "6039",
               "6040", "6041", "6042", "6043", "6044", "6045", "6046", "6047", "6048", "6049", "6050", "6051", "6052",
               "6053", "6054", "6055", "6056", "6057", "6058", "6059", "6060", "6061", "6062", "6063", "6064", "6066",
               "6067", "6068", "6069", "6070", "6071", "6072", "6073", "6074", "6075", "6076", "6077", "6078", "6079",
               "6080", "6081", "6082", "6083", "6084", "6085", "6086", "6087", "6088", "6089", "6090", "6091", "6092",
               "6093", "6094", "6095", "6096", "6097", "6098", "6099", "6100", "6101", "6102", "6103", "6104", "6105",
               "6106", "6107", "6108", "6109", "6110", "6111", "6112", "6113", "6114", "6115", "6116", "6117", "6118",
               "6120", "6121", "6122", "6123", "6125", "6126", "6127", "6128", "6129", "6130", "6131", "6132", "6133",
               "6134", "6135", "6136", "6137", "6138", "6139", "6140", "6141", "6142", "6143", "6144", "6145", "6146",
               "6147", "6148", "6149", "6150", "6151", "6152", "6154", "6155", "6156", "6157", "6158", "6159", "6160",
               "6161", "6162", "6163", "6165", "6166", "6167", "6168", "6169", "6170", "6171", "6172", "6173", "6174",
               "6175", "6176", "6177", "6178", "6179", "6180", "6181", "6182", "6183", "6184", "6185", "6186", "6187",
               "6188", "6189", "6190", "6191", "6192", "6193", "6194", "6195", "6196", "6197", "6198", "6199", "6200",
               "6201", "6202", "6203", "6204", "6205", "6207", "6208", "6209", "6210", "6211", "6212", "6214", "6216",
               "6217", "6218", "6219", "6220", "6221", "6222", "6223", "6224", "6225", "6226", "6227", "6228", "6229",
               "6230", "6231", "6232", "6233", "6234", "6235", "6236", "6237", "6238", "6239", "6240", "6241", "6242",
               "6243", "6244", "6245", "6246", "6247", "6248", "6249", "6250", "6251", "6252", "6253", "6254", "6255",
               "6256", "6257", "6258", "6259", "6260", "6261", "6262", "6263", "6264", "6265", "6266", "6267", "6268",
               "6269", "6270", "6271", "6272", "6273", "6274", "6275", "6276", "6277", "6278", "6279", "6280", "6281",
               "6282", "6283", "6284", "6285", "6286", "6287", "6288", "6289", "6290", "6291", "6292", "6293", "6294",
               "6295", "6296", "6297", "6298", "6299", "6300", "6301", "6302", "6303", "6304", "6305", "6306", "6307",
               "6308", "6309", "6310", "6311", "6312", "6313", "6314", "6315", "6317", "6318", "6319", "6320", "6321",
               "6322", "6324", "6325", "6326", "6327", "6328", "6329", "6330", "6331", "6332", "6333", "6334", "6335",
               "6336", "6337", "6338", "6339", "6340", "6341", "6342", "6343", "6344", "6345", "6346", "6347", "6348",
               "6349", "6350", "6351", "6352", "6353", "6354", "6355", "6356", "6357", "6358", "6359", "6360", "6361",
               "6362", "6363", "6364", "6365", "6366", "6367", "6368", "6369", "6370", "6371", "6372", "6373", "6374",
               "6375", "6376", "6377", "6378", "6379", "6380", "6381", "6382", "6383", "6384", "6385", "6386", "6387",
               "6388", "6389", "6390", "6391", "6392", "6393", "6395", "6396", "6398", "6399", "6400", "6401", "6402",
               "6403", "6404", "6405", "6406", "6407", "6408", "6409", "6410", "6411", "6412", "6413", "6414", "6415",
               "6416", "6418", "6419", "6420", "6421", "6422", "6423", "6424", "6425", "6426", "6427", "6428", "6429",
               "6430", "6431", "6432", "6433", "6434", "6436", "6437", "6438", "6439", "6440", "6441", "6442", "6443",
               "6444", "6445", "6446", "6447", "6448", "6449", "6451", "6452", "6453", "6454", "6455", "6456", "6457",
               "6458", "6459", "6460", "6462", "6463", "6464", "6465", "6466", "6467", "6468", "6469", "6470", "6473",
               "6474", "6475", "6476", "6477", "6478", "6479", "6480", "6481", "6482", "6483", "6484", "6485", "6487",
               "6488", "6489", "6491", "6493", "6494", "6495", "6496", "6497", "6498", "6499", "6500", "6501", "6502",
               "6503", "6504", "6505", "6506", "6507", "6508", "6509", "6510", "6511", "6512", "6513", "6514", "6515",
               "6516", "6517", "6518", "6519", "6520", "6521", "6522", "6524", "6525", "6526", "6527", "6528", "6529",
               "6530", "6531", "6532", "6533", "6535", "6536", "6537", "6538", "6539", "6540", "6541", "6542", "6543",
               "6546", "6547", "6549", "6551", "6552", "6553", "6554", "6555", "6556", "6557", "6558", "6559", "6560",
               "6561", "6562", "6563", "6564", "6565", "6566", "6567", "6568", "6569", "6570", "6571", "6572", "6573",
               "6574", "6575", "6576", "6577", "6578", "6579", "6580", "6581", "6582", "6583", "6584", "6585", "6586",
               "6587", "6588", "6589", "6590", "6591", "6592", "6593", "6594", "6596", "6597", "6598", "6599", "6600",
               "6601", "6602", "6603", "6604", "6605", "6606", "6607", "6608", "6609", "6610", "6611", "6612", "6613",
               "6615", "6616", "6617", "6618", "6619", "6620", "6621", "6622", "6623", "6624", "6625", "6626", "6627",
               "6628", "6629", "6631", "6633", "6634", "6635", "6636", "6637", "6638", "6639", "6640", "6642", "6643",
               "6644", "6645", "6646", "6647", "6648", "6649", "6650", "6651", "6652", "6653", "6654", "6655", "6656",
               "6657", "6658", "6659", "6660", "6661", "6662", "6663", "6664", "6665", "6666", "6667", "6668", "6669",
               "6670", "6671", "6672", "6673", "6674", "6675", "6676", "6677", "6678", "6679", "6680", "6681", "6682",
               "6683", "6685", "6686", "6687", "6688", "6689", "6690", "6691", "6692", "6693", "6695", "6696", "6698",
               "6699", "6700", "6702", "6705", "6707", "6709", "6710", "6711", "6712", "6713", "6714", "6715", "6716",
               "6717", "6718", "6719", "6720", "6721", "6722", "6723", "6724", "6725", "6726", "6727", "6728", "6729",
               "6730", "6731", "6732", "6733", "6734", "6735", "6736", "6737", "6738", "6739", "6740", "6741", "6742",
               "6743", "6744", "6745", "6746", "6747", "6748", "6749", "6750", "6751", "6752", "6753", "6754", "6755",
               "6756", "6757", "6758", "6759", "6760", "6761", "6762", "6763", "6766", "6769", "6771", "6774", "6775",
               "6776", "6777", "6778", "6779", "6780", "6781", "6782", "6783", "6784", "6785", "6786", "6787", "6788",
               "6789", "6790", "6791", "6792", "6793", "6794", "6796", "6797", "6798", "6799", "6800", "6801", "6802",
               "6803", "6804", "6805", "6806", "6807", "6808", "6809", "6810", "6811", "6812", "6814", "6815", "6816",
               "6817", "6818", "6820", "6821", "6822", "6823", "6824", "6825", "6826", "6827", "6828", "6829", "6830",
               "6831", "6832", "6833", "6834", "6835", "6836", "6837", "6838", "6839", "6840", "6841", "6843", "6844",
               "6845", "6846", "6847", "6849", "6850", "6851", "6853", "6854", "6855", "6856", "6857", "6858", "6859",
               "6860", "6861", "6862", "6863", "6864", "6865", "6866", "6867", "6868", "6869", "6870", "6872", "6873",
               "6874", "6875", "6876", "6877", "6878", "6879", "6880", "6881", "6882", "6883", "6884", "6885", "6886",
               "6888", "6890", "6891", "6892", "6893", "6894", "6895", "6896", "6897", "6898", "6899", "6900", "6901",
               "6902", "6903", "6904", "6905", "6906", "6907", "6908", "6909", "6910", "6911", "6912", "6913", "6914",
               "6915", "6916", "6917", "6918", "6919", "6920", "6921", "6922", "6923", "6924", "6925", "6926", "6927",
               "6928", "6929", "6930", "6932", "6933", "6934", "6935", "6936", "6938", "6939", "6940", "6941", "6942",
               "6943", "6944", "6945", "6946", "6947", "6948", "6949", "6950", "6951", "6952", "6953", "6954", "6955",
               "6956", "6957", "6960", "6962", "6963", "6964", "6965", "6966", "6967", "6968", "6969", "6970", "6971",
               "6972", "6973", "6974", "6975", "6976", "6979", "6980", "6981", "6982", "6983", "6984", "6985", "6986",
               "6987", "6988", "6989", "6990", "6991", "6992", "6993", "6995", "6996", "6997", "6998", "6999", "7000",
               "7001", "7002", "7003", "7004", "7005", "7006", "7007", "7008", "7009", "7010", "7011", "7012", "7013",
               "7014", "7015", "7016", "7017", "7018", "7019", "7020", "7021", "7022", "7023", "7024", "7025", "7028",
               "7029", "7030", "7031", "7032", "7033", "7034", "7035", "7036", "7049", "7050", "7051", "7053", "7054",
               "7055", "7056", "7057", "7058", "7059", "7061", "7063", "7064", "7065", "7066", "7067", "7068", "7069",
               "7070", "7071", "7072", "7073", "7074", "7075", "7078", "7079", "7080", "7081", "7082", "7083", "7084",
               "7085", "7086", "7087", "7088", "7089", "7090", "7091", "7092", "7093", "7094", "7095", "7097", "7098",
               "7099", "7101", "7102", "7103", "7104", "7105", "7106", "7107", "7108", "7109", "7110", "7111", "7112",
               "7113", "7114", "7115", "7116", "7117", "7118", "7119", "7120", "7121", "7122", "7123", "7124", "7125",
               "7126", "7127", "7128", "7129", "7130", "7131", "7132", "7133", "7134", "7135", "7136", "7137", "7139",
               "7140", "7141", "7142", "7144", "7147", "7149", "7150", "7151", "7152", "7153", "7154", "7155", "7156",
               "7157", "7158", "7159", "7160", "7161", "7162", "7165", "7166", "7167", "7168", "7169", "7170", "7171",
               "7172", "7173", "7174", "7175", "7176", "7177", "7178", "7179", "7180", "7181", "7182", "7184", "7185",
               "7186", "7187", "7188", "7189", "7190", "7191", "7192", "7193", "7194", "7196", "7199", "7200", "7201",
               "7202", "7203", "7204", "7205", "7206", "7207", "7208", "7209", "7210", "7211", "7212", "7213", "7214",
               "7215", "7216", "7218", "7219", "7220", "7221", "7222", "7223", "7224", "7225", "7226", "7227", "7228",
               "7229", "7230", "7231", "7232", "7233", "7234", "7235", "7236", "7237", "7238", "7239", "7240", "7241",
               "7242", "7243", "7244", "7245", "7246", "7247", "7248", "7249", "7250", "7251", "7252", "7253", "7254",
               "7255", "7256", "7257", "7258", "7259", "7260", "7261", "7262", "7263", "7264", "7266", "7268", "7269",
               "7270", "7271", "7272", "7273", "7274", "7276", "7277", "7278", "7279", "7280", "7281", "7282", "7283",
               "7284", "7285", "7286", "7287", "7288", "7289", "7290", "7291", "7292", "7293", "7294", "7295", "7296",
               "7297", "7298", "7299", "7300", "7301", "7302", "7303", "7304", "7305", "7306", "7307", "7308", "7309",
               "7310", "7311", "7312", "7313", "7314", "7315", "7316", "7317", "7318", "7319", "7320", "7321", "7322",
               "7323", "7324", "7325", "7326", "7327", "7330", "7339", "7340", "7341", "7342", "7349", "7350", "7351",
               "7352", "7354", "7356", "7357", "7358", "7359", "7360", "7361", "7362", "7364", "7365", "7366", "7367",
               "7370", "7371", "7375", "7376", "7377", "7381", "7382", "7383", "7384", "7386", "7387", "7388", "7389",
               "7390", "7391", "7392", "7395", "7396", "7397", "7399", "7400", "7401", "7402", "7403", "7404", "7405",
               "7406", "7407", "7408", "7409", "7410", "7411", "7412", "7413", "7414", "7415", "7416", "7417", "7418",
               "7419", "7420", "7421", "7422", "7423", "7424", "7425", "7426", "7427", "7428", "7429", "7430", "7432",
               "7433", "7434", "7435", "9500", "10001", "10002", "10004", "10005", "10006", "10008", "10009", "10010",
               "10011", "10012", "10013", "10014", "10015", "10016", "10017", "10018", "10019", "10021", "10022",
               "10023", "10024", "10025", "10026", "10027", "10028", "10029", "10030", "10031", "10032", "10033",
               "10034", "10035", "10036", "10037", "10038", "10039", "10040", "10041", "10042", "10043", "10044",
               "10045", "10046", "10047", "10048", "10049", "10050", "10051", "10052", "10053", "10054", "10055",
               "10056", "10057", "10058", "10059", "10060", "10061", "10062", "10063", "10064", "10065", "10066",
               "10067", "10068", "10069", "10070", "10071", "10072", "10073", "10074", "10075", "10076", "10077",
               "10078", "10079", "10080", "10081", "10082", "10083", "10084", "10085", "10086", "10087", "10088",
               "10089", "10090", "10091", "10092", "10093", "10094", "10095", "10096", "10097", "10098", "10099",
               "10100", "10101", "10102", "10103", "10104", "10105", "10106", "10107", "10108", "10109", "10110",
               "10111", "10112", "10113", "10114", "10115", "10116", "10117", "10118", "10119", "10120", "10121",
               "10122", "10123", "10124", "10125", "10126", "10127", "10128", "10129", "10130", "10131", "10132",
               "10133", "10134", "10135", "10136", "10137", "10138", "10139", "10140", "10141", "10142", "10143",
               "10144", "10145", "10146", "10147", "10148", "10149", "10150", "10151", "10152", "10153", "10154",
               "10155", "10156", "10157", "10158", "10159", "10160", "10161", "10162", "10163", "10164", "10165",
               "10166", "10167", "10168", "10169", "10170", "10171", "10172", "10173", "10174", "10175", "10176",
               "10177", "10178", "10179", "10180", "10181", "10182", "10183", "10184", "10185", "10186", "10187",
               "10188", "10189", "10190", "10191", "10192", "10193", "10194", "10195", "10196", "10197", "10198",
               "10199", "10200", "10201", "10202", "10203", "10204", "10205", "10206", "10207", "10208", "10209",
               "10210", "10211", "10212", "10213", "10214", "10215", "10216", "10217", "10218", "10219", "10220",
               "10221", "10222", "10223", "10224", "10225", "10226", "10227", "10228", "10229", "10230", "10231",
               "10232", "10233", "10234", "10235", "10236", "10237", "10238", "10239", "10240", "10241", "10242",
               "10243", "10244", "10245", "10246", "10247", "10248", "10249", "10250", "10251", "10252", "10253",
               "10254", "10255", "10256", "10257", "10258", "10259", "10260", "10261", "10262", "10263", "10264",
               "10265", "10266", "10267", "10268", "10269", "10270", "10271", "10272", "10273", "10274", "10275",
               "10276", "10277", "10278", "10279", "10280", "10281", "10282", "10283", "10284", "10285", "10286",
               "10287", "10288", "10289", "10290", "10291", "10292", "10293", "10294", "10295", "10296", "10297",
               "10298", "10299", "10300", "10301", "10302", "10303", "10304", "10305", "10306", "10307", "10308",
               "10309", "10310", "10311", "10312", "10313", "10314", "10315", "10316", "10317", "10318", "10319",
               "10320", "10321", "10322", "10323", "10324", "10325", "10326", "10327", "10328", "10329", "10330",
               "10331", "10332", "10333", "10334", "10335", "10336", "10337", "10338", "11000", "11001", "11002",
               "11003", "11004", "11005", "11006", "11007", "11008", "11009", "11010", "11011", "11012", "11013",
               "11014", "11016", "11017", "11018", "11019", "11020", "11021", "11022", "11024", "11026", "11027",
               "11028", "11029", "11030", "11031", "11032", "11033", "11034", "11035", "11036", "11037", "11038",
               "11040", "11041", "11042", "11043", "11044", "11045", "11046", "11047", "11048", "11049", "11050",
               "11051", "11052", "11053", "11054", "11055", "11056", "11057", "11058", "11059", "11060", "11061",
               "11062", "11063", "11064", "11065", "11067", "11068", "11069", "11070", "11072", "11073", "11074",
               "11075", "11076", "11077", "11078", "11080", "11081", "11082", "11083", "11084", "11085", "11086",
               "11087", "11088", "11089", "11090", "11091", "11092", "11093", "11094", "11095", "11096", "11097",
               "11098", "11099", "11100", "11101", "11102", "11103", "11104", "11105", "11106", "11107", "11108",
               "11109", "11110", "11111", "11112", "11113", "11114", "11115", "11116", "11117", "11118", "11119",
               "11120", "11121", "11122", "11123", "11124", "11125", "11126", "11127", "11128", "11129", "11130",
               "11131", "11132", "11133", "11134", "11135", "11136", "11137", "11138", "11139", "11140", "11141",
               "11142", "11143", "11144", "11145", "11146", "11147", "11148", "11149", "11150", "11151", "11152",
               "11153", "11154", "11155", "11156", "11157", "11158", "11159", "11160", "11161", "11162", "11163",
               "11164", "11165", "11166", "11167", "11168", "11169", "11170", "11171", "11172", "11173", "11175",
               "11176", "11177", "11178", "11179", "11180", "11181", "11182", "11183", "11184", "11185", "11186",
               "11187", "11188", "11189", "11190", "11191", "11192", "11193", "11194", "11195", "11196", "11197",
               "11198", "11199", "11200", "11201", "11202", "11203", "11204", "11205", "11206", "11207", "11208",
               "11209", "11210", "11211", "11212", "11213", "11214", "11215", "11216", "11217", "11218", "11219",
               "11220", "11221", "11222", "11223", "11224", "11225", "11226", "11227", "11228", "11229", "11230",
               "11231", "11232", "11233", "11234", "11235", "11236", "11237", "11238", "11239", "11240", "11242",
               "11244", "11245", "11246", "11247", "11248", "11249", "11250", "11251", "11252", "11253", "11254",
               "11255", "11256", "11257", "11258", "11259", "11260", "11261", "11262", "11263", "11264", "11265",
               "11266", "11267", "11268", "11269", "11270", "11271", "11272", "11273", "11274", "11275", "11276",
               "11277", "11278", "11279", "11280", "11281", "11282", "11283", "11284", "11285", "11286", "11287",
               "11288", "11289", "11290", "11291", "11292", "11293", "11294", "11295", "11296", "11297", "11298",
               "11299", "11300", "11301", "11302", "11303", "11304", "11305", "11306", "11307", "11308", "11309",
               "11310", "11311", "11312", "11313", "11314", "11315", "11316", "11317", "11318", "11319", "11320",
               "11321", "11322", "11323", "11324", "11325", "11326", "11327", "11328", "11329", "11330", "11331",
               "11332", "11333", "11334", "11335", "11336", "11337", "11338", "11339", "11340", "11341", "11342",
               "11343", "11344", "11345", "11346", "11347", "11348", "11349", "11350", "11351", "11352", "11353",
               "11354", "11355", "11356", "11357", "11358", "11359", "11360", "11361", "11362", "11363", "11364",
               "11365", "11366", "11367", "11368", "11369", "11370", "11371", "11372", "11373", "11374", "11375",
               "11376", "11377", "11378", "11379", "11380", "11381", "11382", "11383", "11384", "11385", "11386",
               "11387", "11388", "11389", "11390", "11391", "11392", "11393", "11394", "11395", "11396", "11397",
               "11398", "11399", "11400", "11401", "11402", "11403", "11404", "11405", "11406", "11407", "11408",
               "11409", "11410", "11411", "11412", "11413", "11414", "11415", "11416", "11417", "11419", "11420",
               "11422", "11423", "11424", "11425", "11426", "11427", "11428", "11429", "11430", "11431", "11434",
               "11435", "11436", "11437", "11438", "11439", "11440", "11441", "11442", "11443", "11450", "11451",
               "11452", "11453", "11454", "11456", "11457", "11458", "11459", "11460", "11461", "11462", "11463",
               "11466", "11467", "11468", "11469", "11470", "11471", "11472", "11473", "11474", "11475", "11476",
               "11477", "11478", "11479", "11480", "11481", "11482", "11483", "11484", "11485", "11486", "11487",
               "11488", "11489", "11490", "11491", "11492", "11493", "11494", "11495", "11496", "11497", "11498",
               "11499", "11500", "11501", "11502", "11504", "11505", "11506", "11507", "11508", "11511", "11512",
               "11513", "11514", "11515", "11516", "11519", "11520", "11521", "11522", "11523", "11524", "11525",
               "11526", "11527", "11528", "11529", "11530", "11531", "11532", "11534", "11535", "11536", "11537",
               "11538", "11539", "11540", "11541", "11542", "11543", "11544", "11545", "11546", "11547", "11548",
               "11549", "11550", "11551", "11552", "11553", "11554", "11555", "11556", "11557", "11558", "11559",
               "11560", "11561", "11562", "11563", "11564", "11565", "11566", "11567", "11568", "11569", "11570",
               "11572", "11573", "11574", "11575", "11578", "11580", "11581", "11582", "11583", "11584", "11585",
               "11587", "11588", "11589", "11590", "11591", "11592", "11593", "11594", "11595", "11596", "11597",
               "11599", "11601", "11604", "11605", "11606", "11607", "11608", "11609", "12000", "12002", "12003",
               "12009", "12010", "12011", "12014", "12016", "12018", "12019", "12020", "12021", "12022", "12024",
               "12026", "12032", "12033", "12035", "12036", "12038", "12039", "12040", "12041", "15000", "15001",
               "15002", "15003", "15007", "15008", "15009", "15010", "15011", "15012", "15013", "15014", "15016",
               "15019", "15020", "15021", "15023", "15025", "15026", "15027", "15028", "15031", "15032", "15033",
               "15034", "15035", "15036", "15037", "15038", "15039", "15040", "15041", "15042", "15043", "15044",
               "15045", "15046", "15047", "15048", "15049", "15050", "15051", "15052", "15053", "15054", "15055",
               "15056", "15057", "15058", "15059", "15060", "15061", "15063", "15064", "15065", "15066", "15067",
               "15068", "15069", "15070", "15071", "15072", "15073", "15074", "15075", "15076", "15077", "15078",
               "15079", "15080", "15081", "15082", "15083", "15084", "15085", "15086", "15087", "15088", "15089",
               "15090", "15091", "15092", "15093", "15094", "15095", "15097", "15098", "15099", "15100", "15101",
               "15102", "15103", "15104", "15105", "15106", "15107", "15108", "15109", "15110", "15111", "15112",
               "15113", "15114", "15115", "15116", "15117", "15118", "15120", "15121", "15122", "15123", "15124",
               "15125", "15126", "15127", "15128", "15129", "15130", "15131", "15132", "15133", "15134", "15135",
               "15136", "15137", "15138", "15139", "15140", "15141", "15142", "15143", "15144", "15145", "15148",
               "15149", "15150", "15152", "15153", "15154", "15156", "15157", "15158", "15159", "15160", "15161",
               "15162", "15163", "15164", "15166", "15167", "15168", "15169", "15170", "15171", "15172", "15173",
               "15174", "15175", "15176", "15177", "15178", "15179", "15180", "15183", "15184", "15185", "15186",
               "15187", "15188", "15189", "15190", "15191", "15192", "15193", "15194", "15195", "15196", "15197",
               "15198", "15199", "15200", "15201", "15202", "15203", "15204", "15206", "15207", "15208", "15209",
               "15210", "15211", "15212", "15213", "15214", "15215", "15217", "15218", "15219", "15221", "15222",
               "15223", "15224", "15225", "15226", "15227", "15228", "15229", "15230", "15231", "15232", "15233",
               "15234", "15235", "15236", "15237", "15238", "15239", "15240", "15241", "15242", "15243", "15244",
               "15245", "15246", "15247", "15248", "15249", "15250", "15252", "15253", "15254", "15255", "15256",
               "15257", "15258", "15259", "15261", "15263", "15264", "15265", "15266", "15267", "15268", "15269",
               "15270", "15271", "15272", "15273", "15274", "15275", "15276", "15277", "15278", "15279", "15280",
               "15281", "15282", "15283", "15284", "15285", "15286", "15287", "15288", "15289", "15290", "15291",
               "15292", "15293", "15294", "15295", "15296", "15298", "15299", "15300", "15301", "15302", "15303",
               "15304", "15305", "15306", "15309", "15312", "15313", "15314", "15315", "15316", "15317", "15318",
               "15319", "15320", "15322", "15323", "15324", "15325", "15326", "15327", "15329", "15330", "15331",
               "15332", "15333", "15334", "15335", "15336", "15337", "15338", "15339", "15340", "15341", "15342",
               "15343", "15344", "15345", "15346", "15347", "15348", "15349", "15350", "15351", "15352", "15353",
               "15354", "15355", "15356", "15357", "15358", "15359", "15360", "15361", "15362", "15363", "15364",
               "15366", "15367", "15368", "15369", "15370", "15371", "15372", "15373", "15374", "15375", "15376",
               "15377", "15378", "15379", "15380", "15381", "15382", "15383", "15384", "15385", "15386", "15388",
               "15389", "15390", "15391", "15392", "15393", "15394", "15395", "15396", "15397", "15398", "15399",
               "15400", "15401", "15402", "15403", "15404", "15405", "15406", "15407", "15408", "15409", "15410",
               "15411", "15412", "15413", "15415", "15417", "15418", "15419", "15420", "15421", "15422", "15423",
               "15424", "15425", "15426", "15427", "15429", "15430", "15431", "15433", "15434", "15435", "15436",
               "15437", "15438", "15439", "15440", "15441", "15442", "15443", "15444", "15445", "15446", "15447",
               "15448", "15449", "15450", "15451", "15452", "15453", "15454", "15455", "15456", "15457", "15458",
               "15459", "15460", "15461", "15462", "15463", "15464", "15465", "15466", "15467", "15468", "15469",
               "15470", "15471", "15473", "15474", "15475", "15476", "15477", "15478", "15479", "15480", "15481",
               "15482", "15483", "15485", "15486", "15487", "15488", "15489", "15490", "15491", "15492", "15493",
               "15494", "15495", "15496", "15497", "15498", "15499", "15500", "15501", "15502", "15503", "15504",
               "15505", "15506", "15507", "15508", "15509", "15510", "15511", "15512", "15513", "15514", "15515",
               "15517", "15518", "15519", "15520", "15521", "15522", "15524", "15525", "15527", "15528", "15529",
               "15530", "15531", "15532", "15533", "15534", "15535", "15536", "15537", "15538", "15539", "15540",
               "15541", "15542", "15543", "15544", "15545", "15546", "15547", "15548", "15549", "15550", "15551",
               "15552", "15553", "15554", "15555", "15557", "15558", "15559", "15560", "15561", "15562", "15563",
               "15564", "15565", "15566", "15567", "15568", "15569", "15570", "15571", "15572", "15573", "15574",
               "15575", "15576", "15577", "15578", "15579", "15580", "15581", "15582", "15583", "15584", "15585",
               "15586", "15587", "15589", "15590", "15591", "15592", "15593", "15594", "15595", "15596", "15597",
               "15598", "15599", "15601", "15602", "15603", "15604", "15605", "15606", "15607", "15608", "15609",
               "15610", "15611", "15612", "15613", "15614", "15615", "15616", "15617", "15618", "15620", "15621",
               "15622", "15623", "15624", "15625", "15626", "15627", "15628", "15629", "15632", "15633", "15634",
               "15635", "15636", "15637", "15638", "15639", "15640", "15641", "15642", "15643", "15644", "15645",
               "15646", "15647", "15648", "15649", "15650", "15651", "15652", "15653", "15654", "15655", "15656",
               "15658", "15659", "15660", "15661", "15662", "15663", "15664", "15665", "15666", "15667", "15669",
               "15670", "15671", "15672", "15673", "15674", "15676", "15677", "15678", "15679", "15680", "15681",
               "15682", "15683", "15684", "15685", "15687", "15688", "15689", "15690", "15691", "15692", "15693",
               "15694", "15695", "15697", "15698", "15700", "15701", "15702", "15703", "15704", "15705", "15706",
               "15707", "15708", "15710", "15711", "15712", "15713", "15714", "15715", "15716", "15717", "15718",
               "15719", "15720", "15721", "15722", "15723", "15724", "15725", "15726", "15727", "15728", "15729",
               "15730", "15731", "15732", "15733", "15734", "15735", "15736", "15737", "15738", "15739", "15740",
               "15741", "15742", "15743", "15744", "15745", "15746", "15747", "15748", "15749", "15750", "15751",
               "15752", "15753", "15755", "15756", "15757", "15758", "15759", "15760", "15761", "15762", "15763",
               "15764", "15765", "15766", "15767", "15768", "15769", "15770", "15773", "15774", "15775", "15776",
               "15777", "15778", "15779", "15780", "15781", "15782", "15783", "15784", "15785", "15786", "15787",
               "15788", "15789", "15790", "15791", "15792", "15793", "15795", "15796", "15797", "15798", "15799",
               "15800", "15801", "15802", "15803", "15804", "15805", "15806", "15807", "15808", "15809", "15810",
               "15811", "15812", "15813", "15814", "15816", "15817", "15818", "15819", "15820", "15821", "15822",
               "15823", "15824", "15825", "15826", "15827", "15828", "15829", "15830", "15831", "15832", "15833",
               "15834", "15835", "15836", "15837", "15839", "15840", "15841", "15843", "15844", "15845", "15846",
               "15847", "15848", "15849", "15850", "15851", "15852", "15853", "15854", "15855", "15856", "15857",
               "15858", "15859", "15860", "15863", "15864", "15865", "15866", "15867", "15868", "15869", "15870",
               "15871", "15873", "15874", "15875", "15876", "15877", "15878", "15879", "15880", "15881", "15882",
               "15883", "15884", "15885", "15887", "15888", "15889", "15890", "15891", "15892", "15893", "15894",
               "15895", "15896", "15897", "15898", "15899", "15900", "15901", "15902", "15903", "15904", "15905",
               "15906", "15907", "15908", "15909", "15911", "15912", "15913", "15914", "15915", "15916", "15917",
               "15918", "15919", "15921", "15922", "15923", "15924", "15925", "15926", "15927", "15928", "15930",
               "15931", "15932", "15933", "15934", "15935", "15936", "15937", "15938", "15939", "15940", "15941",
               "15942", "15943", "15944", "15945", "15946", "15948", "15949", "15950", "15951", "15952", "15953",
               "15954", "15955", "15956", "15957", "15958", "15959", "15960", "15961", "15962", "15963", "15964",
               "15965", "15966", "15967", "15968", "15969", "15970", "15971", "15972", "15973", "15974", "15975",
               "15977", "15978", "15979", "15980", "15981", "15983", "15986", "15987", "15988", "15989", "15990",
               "15991", "15992", "15993", "15994", "15995", "15996", "15997", "15998", "15999", "16000", "16001",
               "16002", "16003", "16005", "16006", "16007", "16008", "16009", "16010", "16011", "16012", "16014",
               "16015", "16016", "16018", "16019", "16020", "16021", "16022", "16023", "16024", "16025", "16026",
               "16027", "16028", "16029", "16030", "16031", "16032", "16033", "16034", "16035", "16037", "16038",
               "16039", "16040", "16041", "16042", "16043", "16044", "16045", "16046", "16047", "16048", "16049",
               "16050", "16051", "16052", "16053", "16054", "16055", "16056", "16057", "16058", "16059", "16060",
               "16061", "16062", "16063", "16064", "16065", "16066", "16067", "16068", "16069", "16070", "16071",
               "16072", "16073", "16074", "16075", "16076", "16077", "16078", "16079", "16080", "16082", "16083",
               "16084", "16085", "16087", "16088", "16089", "16090", "16091", "16093", "16094", "16096", "16097",
               "16098", "16099", "16100", "16101", "16102", "16103", "16104", "16105", "16107", "16108", "16109",
               "16110", "16111", "16112", "16113", "16114", "16115", "16116", "16117", "16118", "16119", "16120",
               "16121", "16122", "16123", "16124", "16126", "16127", "16129", "16130", "16131", "16133", "16134",
               "16135", "16136", "16137", "16138", "16139", "16140", "16141", "16142", "16143", "16144", "16145",
               "16146", "16147", "16148", "16149", "16150", "16151", "16152", "16153", "16154", "16155", "16156",
               "16157", "16158", "16159", "16160", "16161", "16162", "16163", "16164", "16165", "16166", "16167",
               "16168", "16169", "16170", "16171", "16172", "16173", "16174", "16175", "16176", "16177", "16178",
               "16179", "16180", "16181", "16182", "16183", "16184", "16185", "16186", "16187", "16188", "16189",
               "16190", "16191", "16192", "16193", "16194", "16195", "16196", "16197", "16198", "16199", "16200",
               "16201", "16202", "16203", "16204", "16205", "16206", "16207", "16208", "16209", "16210", "16211",
               "16212", "16213", "16214", "16215", "16216", "16217", "16218", "16219", "16221", "16222", "16223",
               "16225", "16226", "16228", "16229", "16230", "16231", "16233", "16234", "16235", "16236", "16237",
               "16238", "16239", "16240", "16241", "16243", "16244", "16245", "16246", "16247", "16248", "16249",
               "16250", "16251", "16252", "16253", "16254", "16255", "16256", "16257", "16258", "16259", "16260",
               "16261", "16262", "16263", "16264", "16265", "16266", "16267", "16268", "16269", "16270", "16271",
               "16272", "16273", "16274", "16275", "16276", "16277", "16278", "16279", "16280", "16281", "16282",
               "16283", "16284", "16285", "16286", "16287", "16288", "16289", "16290", "16291", "16292", "16293",
               "16294", "16295", "16296", "16297", "16298", "16299", "16300", "16301", "16302", "16303", "16304",
               "16305", "16306", "16307", "16309", "16310", "16311", "16312", "16313", "16314", "16315", "16316",
               "16317", "16318", "16319", "16320", "16321", "16322", "16323", "16325", "16327", "16328", "16329",
               "16330", "16332", "16333", "16334", "16335", "16336", "16337", "16338", "16339", "16340", "16341",
               "16342", "16343", "16344", "16345", "16347", "16348", "16349", "16350", "16352", "16353", "16354",
               "16355", "16356", "16357", "16359", "16360", "16361", "16362", "16363", "16364", "16365", "16366",
               "16367", "16368", "16370", "16371", "16372", "16373", "16374", "16375", "16376", "16377", "16378",
               "16379", "16380", "16381", "16382", "16384", "16385", "16386", "16387", "16388", "16389", "16390",
               "16391", "16392", "16393", "16394", "16395", "16396", "16397", "16398", "16399", "16400", "16401",
               "16402", "16403", "16404", "16405", "16406", "16407", "16408", "16409", "16411", "16412", "16413",
               "16414", "16415", "16416", "16417", "16418", "16419", "16420", "16421", "16422", "16423", "16424",
               "16425", "16426", "16428", "16430", "16431", "16435", "16436", "16437", "16438", "16439", "16440",
               "16441", "16442", "16443", "16444", "16445", "16446", "16447", "16448", "16449", "16450", "16451",
               "16452", "16453", "16454", "16455", "16456", "16457", "16458", "16459", "16460", "16461", "16462",
               "16463", "16464", "16465", "16466", "16467", "16468", "16469", "16470", "16471", "16472", "16473",
               "16475", "16476", "16477", "16478", "16479", "16480", "16481", "16482", "16483", "16484", "16485",
               "16486", "16487", "16489", "16490", "16491", "16492", "16493", "16494", "16495", "16496", "16497",
               "16498", "16499", "16500", "16501", "16502", "16503", "16504", "16505", "16506", "16507", "16508",
               "16509", "16510", "16511", "16512", "16514", "16515", "16516", "16517", "16518", "16519", "16520",
               "16521", "16522", "16523", "16524", "16525", "16526", "16527", "16528", "16529", "16530", "16531",
               "16532", "16533", "16534", "16536", "16537", "16538", "16539", "16540", "16541", "16542", "16543",
               "16544", "16545", "16546", "16547", "16548", "16549", "16550", "16551", "16555", "16556", "16557",
               "16558", "16559", "16560", "16561", "16562", "16563", "16564", "16565", "16566", "16567", "16568",
               "16569", "16570", "16571", "16572", "16574", "16575", "16576", "16577", "16578", "16579", "16580",
               "16581", "16582", "16583", "16584", "16585", "16586", "16587", "16588", "16589", "16590", "16591",
               "16592", "16593", "16594", "16595", "16596", "16597", "16598", "16599", "16600", "16601", "16602",
               "16603", "16604", "16605", "16606", "16607", "16608", "16609", "16610", "16611", "16612", "16613",
               "16614", "16615", "16616", "16617", "16618", "16619", "16620", "16621", "16622", "16623", "16624",
               "16626", "16627", "16628", "16629", "16630", "16632", "16633", "16634", "16635", "16636", "16637",
               "16638", "16639", "16640", "16641", "16642", "16643", "16646", "16647", "16648", "16649", "16652",
               "16653", "16654", "16655", "16656", "16657", "16658", "16659", "16660", "16661", "16662", "16663",
               "16664", "16665", "16666", "16667", "16668", "16669", "16670", "16671", "16672", "16673", "16674",
               "16675", "16676", "16677", "16678", "16679", "16680", "16681", "16682", "16683", "16684", "16685",
               "16686", "16687", "16688", "16689", "16690", "16691", "16692", "16693", "16694", "16695", "16696",
               "16697", "16698", "16699", "16700", "16701", "16703", "16704", "16705", "16706", "16707", "16708",
               "16709", "16710", "16711", "16712", "16713", "16714", "16715", "16716", "16717", "16718", "16719",
               "16720", "16721", "16722", "16723", "16726", "16727", "16729", "16731", "16732", "16733", "16734",
               "16735", "16736", "16737", "16738", "16739", "16740", "16741", "16742", "16743", "16744", "16745",
               "16746", "16747", "16748", "16749", "16750", "16751", "16752", "16753", "16754", "16755", "16756",
               "16757", "16758", "16759", "16760", "16761", "16762", "16763", "16764", "16766", "16767", "16768",
               "16769", "16770", "16771", "16772", "16773", "16774", "16775", "16776", "16777", "16778", "16779",
               "16780", "16781", "16782", "16783", "16784", "16785", "16786", "16787", "16788", "16789", "16790",
               "16791", "16792", "16794", "16795", "16796", "16797", "16798", "16799", "16800", "16801", "16802",
               "16803", "16804", "16805", "16806", "16807", "16808", "16809", "16810", "16811", "16812", "16813",
               "16814", "16815", "16816", "16817", "16818", "16819", "16820", "16821", "16822", "16824", "16831",
               "16832", "16833", "16834", "16835", "16837", "16838", "16839", "16840", "16841", "16842", "16843",
               "16845", "16846", "16847", "16848", "16849", "16850", "16851", "16852", "16853", "16854", "16856",
               "16857", "16858", "16859", "16860", "16861", "16862", "16863", "16864", "16865", "16866", "16867",
               "16868", "16869", "16870", "16871", "16872", "16873", "16874", "16876", "16877", "16878", "16879",
               "16880", "16881", "16882", "16883", "16884", "16885", "16886", "16887", "16888", "16889", "16890",
               "16891", "16892", "16893", "16894", "16895", "16897", "16898", "16899", "16900", "16901", "16902",
               "16904", "16905", "16907", "16908", "16909", "16910", "16911", "16912", "16913", "16914", "16915",
               "16916", "16917", "16918", "16919", "16920", "16921", "16922", "16923", "16925", "16926", "16927",
               "16928", "16929", "16930", "16931", "16932", "16933", "16934", "16935", "16936", "16937", "16939",
               "16940", "16941", "16942", "16943", "16944", "16945", "16946", "16947", "16948", "16949", "16950",
               "16951", "16952", "16953", "16954", "16955", "16956", "16957", "16958", "16959", "16960", "16961",
               "16962", "16963", "16964", "16965", "16966", "16967", "16968", "16969", "16970", "16971", "16977",
               "16978", "16979", "16980", "16981", "16982", "16983", "16984", "16986", "16988", "16989", "16991",
               "16992", "16994", "16995", "16996", "16997", "16998", "16999", "17000", "17001", "17002", "17003",
               "17005", "17006", "17007", "17008", "17009", "17010", "17011", "17012", "17013", "17014", "17015",
               "17016", "17017", "17018", "17019", "17020", "17021", "17022", "17023", "17024", "17025", "17026",
               "17027", "17028", "17029", "17030", "17031", "17032", "17033", "17034", "17035", "17036", "17038",
               "17039", "17040", "17041", "17042", "17043", "17044", "17045", "17046", "17047", "17048", "17050",
               "17052", "17053", "17054", "17055", "17056", "17057", "17058", "17059", "17060", "17061", "17062",
               "17063", "17064", "17065", "17066", "17067", "17068", "17069", "17070", "17071", "17072", "17073",
               "17074", "17075", "17076", "17077", "17078", "17079", "17080", "17081", "17082", "17083", "17084",
               "17085", "17086", "17087", "17088", "17089", "17090", "17091", "17092", "17093", "17094", "17095",
               "17096", "17097", "17098", "17099", "17100", "17101", "17102", "17103", "17104", "17105", "17106",
               "17107", "17108", "17109", "17110", "17111", "17112", "17113", "17114", "17115", "17116", "17118",
               "17119", "17120", "17121", "17122", "17123", "17124", "17125", "17126", "17127", "17128", "17129",
               "17130", "17131", "17132", "17133", "17134", "17135", "17136", "17137", "17138", "17139", "17141",
               "17143", "17144", "17145", "17146", "17147", "17148", "17149", "17150", "17151", "17152", "17153",
               "17154", "17155", "17156", "17157", "17158", "17159", "17160", "17161", "17162", "17163", "17165",
               "17166", "17169", "17170", "17171", "17172", "17173", "17174", "17175", "17176", "17177", "17178",
               "17179", "17180", "17181", "17182", "17183", "17184", "17185", "17186", "17187", "17188", "17189",
               "17190", "17191", "17192", "17193", "17194", "17195", "17196", "17199", "17200", "17201", "17202",
               "17203", "17204", "17205", "17206", "17207", "17208", "17209", "17210", "17211", "17212", "17213",
               "17214", "17215", "17216", "17217", "17218", "17219", "17220", "17221", "17222", "17223", "17224",
               "17225", "17226", "17227", "17228", "17229", "17231", "17232", "17233", "17234", "17235", "17236",
               "17237", "17238", "17239", "17240", "17241", "17242", "17243", "17244", "17245", "17246", "17247",
               "17249", "17250", "17251", "17252", "17253", "17254", "17255", "17256", "17257", "17258", "17260",
               "17261", "17262", "17263", "17264", "17265", "17266", "17267", "17268", "17270", "17271", "17272",
               "17273", "17274", "17275", "17276", "17277", "17278", "17279", "17280", "17281", "17282", "17283",
               "17284", "17285", "17286", "17287", "17288", "17289", "17290", "17291", "17292", "17293", "17294",
               "17295", "17296", "17297", "17298", "17299", "17300", "17301", "17302", "17303", "17304", "17305",
               "17306", "17307", "17308", "17309", "17310", "17311", "17312", "17313", "17314", "17315", "17316",
               "17318", "17319", "17320", "17321", "17322", "17323", "17324", "17325", "17326", "17327", "17328",
               "17329", "17330", "17331", "17332", "17333", "17334", "17338", "17339", "17340", "17341", "17343",
               "17344", "17345", "17346", "17347", "17350", "17351", "17352", "17353", "17354", "17355", "17356",
               "17357", "17358", "17359", "17360", "17361", "17362", "17363", "17364", "17365", "17367", "17368",
               "17369", "17370", "17371", "17373", "17374", "17375", "17376", "17377", "17378", "17379", "17380",
               "17381", "17382", "17383", "17384", "17385", "17386", "17387", "17388", "17389", "17390", "17391",
               "17392", "17393", "17394", "17395", "17396", "17397", "17398", "17399", "17400", "17401", "17402",
               "17403", "17404", "17405", "17406", "17407", "17408", "17409", "17410", "17411", "17412", "17413",
               "17414", "17415", "17416", "17417", "17418", "17419", "17420", "17421", "17422", "17423", "17424",
               "17425", "17426", "17427", "17428", "17429", "17430", "17431", "17432", "17433", "17434", "17435",
               "17436", "17437", "17438", "17439", "17440", "17441", "17442", "17443", "17444", "17445", "17446",
               "17447", "17448", "17449", "17450", "17451", "17452", "17453", "17454", "17455", "17456", "17457",
               "17458", "17459", "17460", "17461", "17462", "17463", "17464", "17465", "17466", "17467", "17468",
               "17469", "17470", "17471", "17472", "17473", "17474", "17475", "17476", "17477", "17478", "17479",
               "17481", "17482", "17483", "17484", "17485", "17487", "17488", "17489", "17490", "17491", "17492",
               "17493", "17494", "17495", "17496", "17497", "17498", "17499", "17500", "17501", "17502", "17503",
               "17504", "17505", "17506", "17507", "17508", "17509", "17510", "17511", "17512", "17513", "17514",
               "17515", "17517", "17518", "17519", "17520", "17521", "17523", "17524", "17525", "17526", "17527",
               "17528", "17529", "17530", "17531", "17532", "17533", "17534", "17535", "17536", "17537", "17538",
               "17539", "17540", "17541", "17542", "17543", "17544", "17545", "17546", "17547", "17549", "17550",
               "17551", "17552", "17553", "17554", "17555", "17556", "17557", "17558", "17559", "17560", "17561",
               "17562", "17563", "17564", "17565", "17566", "17567", "17568", "17569", "17570", "17571", "17572",
               "17573", "17574", "17575", "17576", "17577", "17578", "17579", "17580", "17581", "17582", "17583",
               "17584", "17585", "17586", "17588", "17589", "17590", "17591", "17592", "17593", "17594", "17595",
               "17596", "17597", "17598", "17599", "17600", "17601", "17602", "17603", "17604", "17606", "17607",
               "17608", "17609", "17610", "17611", "17612", "17613", "17614", "17615", "17616", "17617", "17618",
               "17619", "17620", "17621", "17622", "17623", "17625", "17626", "17627", "17628", "17629", "17630",
               "17631", "17632", "17633", "17634", "17635", "17636", "17637", "17638", "17639", "17640", "17641",
               "17642", "17643", "17644", "17645", "17646", "17647", "17648", "17649", "17650", "17651", "17652",
               "17653", "17654", "17655", "17656", "17657", "17658", "17659", "17660", "17661", "17662", "17663",
               "17664", "17665", "17667", "17668", "17669", "17670", "17671", "17672", "17673", "17674", "17675",
               "17678", "17679", "17680", "17681", "17682", "17683", "17685", "17686", "17687", "17688", "17689",
               "17690", "17691", "17692", "17693", "17694", "17695", "17697", "17698", "17699", "17700", "17701",
               "17702", "17703", "17704", "17705", "17706", "17707", "17708", "17709", "17710", "17711", "17713",
               "17714", "17715", "17716", "17717", "17718", "17719", "17720", "17721", "17722", "17723", "17724",
               "17725", "17726", "17727", "17728", "17729", "17731", "17732", "17733", "17734", "17735", "17736",
               "17737", "17738", "17739", "17740", "17741", "17742", "17743", "17744", "17745", "17746", "17747",
               "17748", "17749", "17750", "17751", "17752", "17753", "17754", "17755", "17756", "17757", "17758",
               "17759", "17760", "17761", "17762", "17763", "17764", "17765", "17766", "17767", "17768", "17769",
               "17770", "17771", "17772", "17773", "17774", "17775", "17776", "17777", "17779", "17780", "17781",
               "17782", "17783", "17784", "17785", "17786", "17787", "17788", "17789", "17790", "17791", "17792",
               "17793", "17794", "17795", "17796", "17797", "17798", "17799", "17803", "17804", "17805", "17806",
               "17807", "17808", "17809", "17810", "17811", "17812", "17813", "17814", "17815", "17816", "17817",
               "17818", "17819", "17820", "17821", "17822", "17823", "17824", "17825", "17826", "17827", "17828",
               "17829", "17830", "17831", "17832", "17833", "17834", "17835", "17836", "17837", "17838", "17839",
               "17840", "17841", "17842", "17843", "17844", "17845", "17846", "17847", "17848", "17849", "17850",
               "17851", "17852", "17853", "17855", "17856", "17857", "17858", "17859", "17860", "17861", "17862",
               "17863", "17864", "17865", "17866", "17867", "17868", "17869", "17870", "17871", "17872", "17873",
               "17874", "17875", "17876", "17877", "17878", "17879", "17880", "17881", "17882", "17883", "17884",
               "17885", "17886", "17887", "17888", "17889", "17890", "17891", "17892", "17893", "17894", "17895",
               "17896", "17897", "17898", "17899", "17900", "17901", "17902", "17903", "17904", "17905", "17906",
               "17907", "17908", "17909", "17910", "17911", "17912", "17913", "17914", "17915", "17916", "17917",
               "17918", "17919", "17920", "17921", "17922", "17923", "17924", "17925", "17926", "17927", "17928",
               "17929", "17930", "17931", "17932", "17934", "17935", "17936", "17937", "17938", "17939", "17940",
               "17941", "17942", "17943", "17944", "17945", "17946", "17947", "17948", "17949", "17950", "17952",
               "17953", "17955", "17956", "17957", "17958", "17959", "17960", "17961", "17962", "17963", "17964",
               "17965", "17966", "17967", "17968", "17969", "17970", "17971", "17972", "17973", "17974", "17975",
               "17976", "17977", "17978", "17979", "17980", "17981", "17982", "17983", "17985", "17986", "17987",
               "17988", "17989", "17990", "17991", "17992", "17993", "17994", "17995", "17996", "17997", "17998",
               "17999", "18000", "18001", "18002", "18003", "18004", "18005", "18006", "18007", "18008", "18009",
               "18010", "18011", "18012", "18013", "18014", "18015", "18016", "18017", "18018", "18019", "18020",
               "18021", "18022", "18023", "18024", "18025", "18026", "18027", "18028", "18029", "18031", "18032",
               "18034", "18035", "18036", "18037", "18039", "18040", "18041", "18042", "18043", "18044", "18045",
               "18046", "18047", "18048", "18049", "18050", "18051", "18052", "18053", "18054", "18055", "18056",
               "18057", "18058", "18059", "18060", "18061", "18062", "18063", "18064", "18065", "18066", "18067",
               "18068", "18069", "18070", "18071", "18072", "18073", "18074", "18075", "18076", "18077", "18078",
               "18079", "18080", "18081", "18082", "18083", "18084", "18085", "18086", "18087", "18088", "18089",
               "18090", "18091", "18092", "18093", "18094", "18095", "18096", "18097", "18098", "18099", "18101",
               "18102", "18103", "18104", "18105", "18106", "18107", "18108", "18109", "18110", "18111", "18112",
               "18113", "18114", "18115", "18116", "18118", "18119", "18121", "18122", "18123", "18124", "18125",
               "18126", "18127", "18128", "18129", "18130", "18131", "18132", "18133", "18134", "18135", "18137",
               "18138", "18141", "18142", "18145", "18146", "18147", "18148", "18151", "18152", "18153", "18154",
               "18155", "18156", "18157", "18158", "18159", "18160", "18161", "18162", "18163", "18164", "18165",
               "18166", "18167", "18169", "18170", "18171", "18175", "18176", "18177", "18178", "18179", "18180",
               "18181", "18182", "18183", "18184", "18185", "18186", "18187", "18188", "18189", "18190", "18191",
               "18192", "18193", "18194", "18195", "18196", "18197", "18198", "18199", "18200", "18201", "18202",
               "18203", "18204", "18205", "18206", "18207", "18208", "18209", "18210", "18211", "18214", "18215",
               "18216", "18217", "18218", "18219", "18220", "18221", "18222", "18223", "18224", "18225", "18226",
               "18227", "18228", "18229", "18230", "18231", "18232", "18234", "18235", "18236", "18237", "18238",
               "18239", "18240", "18242", "18243", "18244", "18245", "18246", "18248", "18249", "18250", "18251",
               "18252", "18253", "18254", "18255", "18256", "18257", "18260", "18261", "18262", "18263", "18265",
               "18266", "18267", "18268", "18269", "18276", "18277", "18278", "18279", "18280", "18281", "18282",
               "18283", "18284", "18285", "18286", "18287", "18288", "18289", "18290", "18291", "18292", "18293",
               "18294", "18295", "18296", "18297", "18298", "18299", "18300", "18301", "18302", "18303", "18304",
               "18305", "18306", "18307", "18308", "18309", "18310", "18312", "18313", "18314", "18315", "18316",
               "18317", "18318", "18319", "18320", "18321", "18322", "18323", "18324", "18325", "18326", "18327",
               "18328", "18329", "18331", "18332", "18333", "18334", "18335", "18336", "18337", "18338", "18339",
               "18340", "18341", "18342", "18344", "18345", "18346", "18347", "18348", "18349", "18350", "18351",
               "18352", "18353", "18354", "18355", "18356", "18357", "18358", "18359", "18360", "18361", "18362",
               "18363", "18364", "18365", "18366", "18367", "18368", "18369", "18370", "18371", "18372", "18373",
               "18374", "18375", "18376", "18377", "18378", "18379", "18380", "18381", "18385", "18386", "18387",
               "18388", "18389", "18390", "18391", "18392", "18393", "18394", "18395", "18396", "18397", "18398",
               "18399", "18400", "18403", "18404", "18405", "18406", "18407", "18408", "18409", "18410", "18411",
               "18412", "18413", "18414", "18415", "18416", "18417", "18418", "18419", "18420", "18421", "18422",
               "18423", "18424", "18425", "18426", "18427", "18428", "18429", "18430", "18431", "18432", "18433",
               "18434", "18435", "18437", "18438", "18439", "18440", "18441", "18442", "18443", "18445", "18446",
               "18447", "18448", "18449", "18452", "18453", "18454", "18455", "18458", "18459", "18461", "18462",
               "18463", "18464", "18465", "18466", "18467", "18468", "18469", "18470", "18471", "18472", "18473",
               "18474", "18475", "18477", "18478", "18479", "18480", "18481", "18484", "18485", "18486", "18487",
               "18488", "18489", "18490", "18491", "18492", "18493", "18494", "18495", "18496", "18497", "18498",
               "18499", "18500", "18501", "18502", "18503", "18504", "18505", "18506", "18507", "18508", "18509",
               "18510", "18511", "18513", "18514", "18515", "18516", "18517", "18518", "18519", "18520", "18521",
               "18522", "18523", "18524", "18526", "18527", "18528", "18529", "18530", "18531", "18532", "18533",
               "18534", "18535", "18536", "18537", "18538", "18539", "18540", "18541", "18542", "18543", "18544",
               "18545", "18546", "18547", "18548", "18549", "18550", "18551", "18552", "18553", "18555", "18556",
               "18557", "18558", "18559", "18560", "18561", "18562", "18563", "18565", "18566", "18567", "18568",
               "18569", "18570", "18571", "18572", "18573", "18574", "18575", "18576", "18577", "18578", "18579",
               "18580", "18581", "18582", "18583", "18584", "18585", "18586", "18587", "18588", "18589", "18590",
               "18591", "18592", "18593", "18595", "18596", "18598", "18599", "18600", "18601", "18602", "18603",
               "18604", "18605", "18606", "18607", "18608", "18609", "18610", "18611", "18612", "18613", "18614",
               "18615", "18616", "18617", "18618", "18619", "18620", "18621", "18622", "18623", "18624", "18625",
               "18626", "18627", "18628", "18629", "18630", "18631", "18632", "18633", "18634", "18635", "18636",
               "18637", "18638", "18639", "18640", "18641", "18642", "18643", "18644", "18645", "18647", "18648",
               "18649", "18650", "18651", "18652", "18653", "18654", "18655", "18656", "18657", "18658", "18659",
               "18661", "18662", "18663", "18664", "18665", "18666", "18667", "18668", "18669", "18670", "18671",
               "18672", "18673", "18676", "18677", "18678", "18679", "18680", "18681", "18682", "18683", "18684",
               "18685", "18686", "18687", "18688", "18689", "18690", "18691", "18692", "18693", "18694", "18695",
               "18696", "18697", "18698", "18699", "18700", "18701", "18702", "18703", "18704", "18705", "18706",
               "18707", "18708", "18709", "18710", "18711", "18712", "18713", "18714", "18715", "18716", "18717",
               "18718", "18719", "18720", "18721", "18722", "18723", "18724", "18725", "18726", "18728", "18729",
               "18730", "18731", "18732", "18733", "18734", "18735", "18736", "18737", "18738", "18740", "18741",
               "18748", "18749", "18750", "18753", "18754", "18755", "18756", "18757", "18758", "18759", "18760",
               "18761", "18762", "18763", "18764", "18765", "18766", "18767", "18768", "18769", "18770", "18771",
               "18772", "18773", "18774", "18775", "18776", "18777", "18778", "18779", "18780", "18781", "18782",
               "18783", "18784", "18785", "18786", "18787", "18788", "18789", "18791", "18792", "18793", "18794",
               "18795", "18796", "18797", "18798", "18799", "18800", "18801", "18802", "18803", "18804", "18805",
               "18806", "18807", "18808", "18809", "18811", "18812", "18813", "18814", "18815", "18816", "18817",
               "18818", "18819", "18820", "18821", "18822", "18823", "18824", "18825", "18826", "18827", "18828",
               "18829", "18830", "18831", "18832", "18833", "18834", "18835", "18836", "18837", "18838", "18839",
               "18840", "18841", "18842", "18843", "18844", "18845", "18846", "18847", "18848", "18849", "18850",
               "18851", "18852", "18853", "18854", "18855", "18856", "18857", "18858", "18859", "18860", "18861",
               "18862", "18863", "18864", "18865", "18867", "18868", "18869", "18870", "18871", "18872", "18873",
               "18874", "18875", "18876", "18877", "18878", "18879", "18880", "18881", "18882", "18883", "18884",
               "18885", "18887", "18888", "18889", "18890", "18891", "18892", "18893", "18894", "18895", "18896",
               "18897", "18898", "18899", "18900", "18901", "18902", "18903", "18904", "18905", "18906", "18907",
               "18908", "18909", "18910", "18911", "18912", "18913", "18914", "18915", "18916", "18917", "18918",
               "18919", "18920", "18921", "18922", "18923", "18924", "18925", "18926", "18927", "18928", "18929",
               "18930", "18931", "18933", "18934", "18935", "18937", "18938", "18939", "18940", "18941", "18942",
               "18943", "18944", "18945", "18946", "18947", "18948", "18949", "18950", "18951", "18952", "18953",
               "18954", "18955", "18956", "18958", "18959", "18961", "18963", "18964", "18965", "18966", "18967",
               "18968", "18969", "18970", "18971", "18972", "18973", "18974", "18975", "18976", "18977", "18978",
               "18979", "18980", "18981", "18982", "18984", "18985", "18986", "18987", "18988", "18989", "18990",
               "18991", "18992", "18993", "18994", "18995", "18998", "18999", "19000", "19001", "19002", "19003",
               "19004", "19005", "19006", "19007", "19008", "19009", "19010", "19011", "19012", "19013", "19014",
               "19015", "19017", "19018", "19023", "19024", "19025", "19026", "19027", "19028", "19029", "19030",
               "19031", "19032", "19033", "19034", "19035", "19036", "19037", "19038", "19039", "19040", "19041",
               "19042", "19043", "19044", "19045", "19046", "19047", "19048", "19049", "19050", "19051", "19052",
               "19053", "19054", "19055", "19056", "19057", "19058", "19059", "19060", "19061", "19062", "19063",
               "19064", "19065", "19066", "19067", "19068", "19069", "19070", "19071", "19072", "19074", "19075",
               "19076", "19077", "19078", "19079", "19080", "19081", "19082", "19083", "19084", "19085", "19086",
               "19087", "19088", "19089", "19090", "19091", "19092", "19093", "19094", "19095", "19098", "19099",
               "19101", "19102", "19103", "19104", "19105", "19106", "19107", "19108", "19109", "19110", "19111",
               "19112", "19113", "19114", "19115", "19116", "19117", "19118", "19119", "19120", "19121", "19122",
               "19123", "19124", "19125", "19126", "19127", "19128", "19129", "19130", "19131", "19132", "19133",
               "19134", "19135", "19136", "19137", "19138", "19139", "19140", "19141", "19142", "19143", "19144",
               "19145", "19146", "19147", "19148", "19149", "19150", "19151", "19152", "19153", "19154", "19155",
               "19156", "19157", "19158", "19159", "19160", "19161", "19162", "19163", "19164", "19165", "19167",
               "19168", "19169", "19170", "19171", "19172", "19173", "19174", "19175", "19176", "19177", "19178",
               "19179", "19180", "19181", "19182", "19183", "19184", "19185", "19186", "19187", "19188", "19189",
               "19190", "19191", "19192", "19193", "19194", "19195", "19196", "19197", "19198", "19200", "19201",
               "19202", "19203", "19204", "19205", "19206", "19207", "19208", "19209", "19210", "19211", "19212",
               "19213", "19214", "19215", "19216", "19217", "19218", "19219", "19220", "19221", "19222", "19223",
               "19224", "19225", "19226", "19227", "19228", "19229", "19230", "19231", "19232", "19233", "19234",
               "19235", "19236", "19237", "19238", "19239", "19240", "19241", "19242", "19243", "19244", "19245",
               "19246", "19247", "19248", "19249", "19250", "19251", "19252", "19253", "19254", "19255", "19256",
               "19257", "19258", "19259", "19260", "19261", "19262", "19263", "19264", "19266", "19267", "19268",
               "19269", "19270", "19271", "19272", "19273", "19274", "19275", "19276", "19277", "19278", "19279",
               "19280", "19281", "19282", "19283", "19284", "19285", "19286", "19287", "19288", "19290", "19291",
               "19292", "19293", "19294", "19295", "19296", "19297", "19298", "19299", "19300", "19301", "19302",
               "19303", "19304", "19305", "19306", "19307", "19308", "19309", "19311", "19312", "19314", "19315",
               "19317", "19318", "19319", "19320", "19321", "19322", "19323", "19324", "19325", "19327", "19328",
               "19329", "19330", "19331", "19332", "19333", "19334", "19335", "19336", "19337", "19338", "19339",
               "19340", "19342", "19344", "19345", "19346", "19347", "19348", "19349", "19350", "19351", "19353",
               "19354", "19355", "19356", "19357", "19358", "19361", "19362", "19363", "19364", "19365", "19366",
               "19367", "19368", "19369", "19370", "19371", "19372", "19374", "19375", "19376", "19377", "19379",
               "19380", "19381", "19382", "19383", "19384", "19385", "19386", "19387", "19388", "19389", "19390",
               "19391", "19392", "19393", "19394", "19395", "19396", "19397", "19398", "19399", "19400", "19402",
               "19403", "19404", "19406", "19407", "19408", "19409", "19410", "19411", "19412", "19413", "19414",
               "19415", "19416", "19417", "19418", "19419", "19421", "19422", "19423", "19424", "19425", "19426",
               "19427", "19428", "19429", "19430", "19431", "19432", "19433", "19435", "19436", "19437", "19438",
               "19439", "19440", "19441", "19442", "19443", "19444", "19446", "19447", "19448", "19449", "19450",
               "19451", "19452", "19453", "19454", "19455", "19456", "19457", "19458", "19459", "19460", "19461",
               "19462", "19463", "19464", "19466", "19467", "19468", "19469", "19470", "19471", "19472", "19473",
               "19474", "19476", "19477", "19478", "19479", "19480", "19481", "19482", "19483", "19484", "19485",
               "19486", "19487", "19488", "19489", "19490", "19491", "19493", "19494", "19495", "19496", "19497",
               "19498", "19499", "19500", "19501", "19502", "19503", "19504", "19505", "19506", "19507", "19508",
               "19510", "19511", "19512", "19513", "19514", "19515", "19516", "19517", "19518", "19519", "19520",
               "19521", "19522", "19523", "19524", "19525", "19526", "19527", "19528", "19530", "19531", "19532",
               "19533", "19534", "19535", "19536", "19538", "19539", "19540", "19541", "19542", "19544", "19545",
               "19546", "19547", "19548", "19549", "19550", "19551", "19552", "19553", "19554", "19555", "19556",
               "19557", "19558", "19559", "19560", "19561", "19562", "19563", "19564", "19565", "19566", "19567",
               "19568", "19569", "19570", "19571", "19572", "19573", "19575", "19576", "19577", "19578", "19579",
               "19580", "19581", "19582", "19583", "19584", "19585", "19586", "19587", "19589", "19590", "19591",
               "19592", "19593", "19594", "19595", "19596", "19597", "19598", "19599", "19600", "19601", "19602",
               "19603", "19604", "19605", "19606", "19607", "19608", "19609", "19610", "19611", "19613", "19614",
               "19615", "19616", "19617", "19618", "19619", "19620", "19621", "19622", "19623", "19624", "19625",
               "19626", "19627", "19628", "19629", "19632", "19633", "19634", "19635", "19637", "19638", "19641",
               "19642", "19643", "19644", "19645", "19646", "19648", "19649", "19650", "19651", "19653", "19654",
               "19655", "19656", "19657", "19658", "19659", "19660", "19661", "19662", "19663", "19664", "19666",
               "19667", "19668", "19669", "19670", "19671", "19672", "19673", "19674", "19675", "19677", "19678",
               "19679", "19680", "19681", "19682", "19683", "19684", "19685", "19686", "19687", "19688", "19689",
               "19692", "19693", "19694", "19695", "19696", "19697", "19698", "19699", "19700", "19701", "19702",
               "19703", "19707", "19708", "19709", "19710", "19711", "19712", "19713", "19714", "19715", "19716",
               "19717", "19718", "19719", "19720", "19721", "19722", "19723", "19724", "19725", "19726", "19727",
               "19728", "19729", "19730", "19731", "19732", "19733", "19734", "19735", "19736", "19737", "19738",
               "19739", "19740", "19741", "19742", "19743", "19744", "19745", "19746", "19747", "19748", "19749",
               "19750", "19751", "19752", "19753", "19754", "19755", "19757", "19758", "19759", "19760", "19763",
               "19764", "19765", "19766", "19767", "19768", "19769", "19770", "19771", "19773", "19774", "19775",
               "19776", "19777", "19778", "19779", "19782", "19783", "19784", "19787", "19788", "19789", "19791",
               "19792", "19796", "19797", "19798", "19799", "19800", "19801", "19803", "19805", "19806", "19807",
               "19808", "19809", "19810", "19811", "19813", "19814", "19815", "19816", "19817", "19818", "19819",
               "19821", "19822", "19823", "19824", "19825", "19826", "19828", "19829", "19830", "19831", "19832",
               "19833", "19834", "19835", "19836", "19837", "19838", "19839", "19840", "19841", "19842", "19843",
               "19844", "19845", "19846", "19847", "19848", "19849", "19850", "19851", "19852", "19853", "19854",
               "19855", "19856", "19857", "19858", "19859", "19860", "19861", "19862", "19863", "19864", "19865",
               "19866", "19867", "19869", "19870", "19872", "19873", "19874", "19875", "19876", "19878", "19879",
               "19880", "19881", "19882", "19883", "19884", "19885", "19886", "19887", "19888", "19889", "19890",
               "19891", "19892", "19893", "19901", "19902", "19904", "19905", "19906", "19907", "19908", "19910",
               "19911", "19912", "19913", "19914", "19915", "19916", "19917", "19921", "19922", "19923", "19925",
               "19926", "19927", "19928", "19929", "19930", "19931", "19932", "19933", "19934", "19935", "19936",
               "19937", "19938", "19939", "19940", "19941", "19942", "19943", "19945", "19946", "19947", "19948",
               "19949", "19951", "19952", "19953", "19954", "19955", "19957", "19958", "19959", "19960", "19961",
               "19962", "19963", "19966", "19967", "19968", "19970", "19971", "19972", "19973", "19974", "19976",
               "19977", "19978", "19979", "19980", "19981", "19982", "19984", "19985", "19986", "19987", "19988",
               "19989", "19990", "19991", "19992", "19993", "19994", "19995", "19996", "19998", "19999", "20001",
               "20002", "20003", "20004", "20005", "20007", "20008", "20009", "20010", "20011", "20012", "20013",
               "20014", "20015", "20016", "20017", "20018", "20019", "20020", "20022", "20023", "20024", "20025",
               "20026", "20027", "20028", "20029", "20030", "20031", "20032", "20033", "20034", "20036", "20037",
               "20038", "20039", "20040", "20041", "20042", "20043", "20044", "20045", "20046", "20047", "20048",
               "20049", "20050", "20051", "20052", "20053", "20054", "20055", "20056", "20057", "20058", "20059",
               "20060", "20061", "20062", "20063", "20064", "20065", "20066", "20067", "20068", "20069", "20070",
               "20071", "20072", "20073", "20074", "20075", "20076", "20078", "20079", "20080", "20081", "20082",
               "20084", "20085", "20086", "20087", "20088", "20089", "20090", "20091", "20092", "20093", "20094",
               "20095", "20098", "20101", "20102", "20103", "20104", "20105", "20107", "20108", "20109", "20110",
               "20111", "20112", "20113", "20114", "20115", "20116", "20117", "20118", "20119", "20121", "20122",
               "20123", "20124", "20125", "20126", "20127", "20128", "21000", "21001", "21002", "21006", "21007",
               "21008", "21009", "21010", "21011", "21013", "21014", "21015", "21018", "21019", "21022", "21023",
               "21024", "21025", "21026", "21027", "21028", "21031", "21032", "21033", "21034", "21035", "21036",
               "21037", "21038", "21039", "21040", "21041", "21042", "21045", "21046", "21047", "21053", "21054",
               "21056", "21057", "21058", "21059", "21060", "21062", "21064", "21065", "21066", "21067", "21068",
               "21069", "21070", "21071", "21072", "21073", "21074", "21075", "21076", "21077", "21078", "21079",
               "21082", "21084", "21085", "21086", "21097", "21098", "21099", "21100", "21101", "21102", "25000",
               "25001", "25002", "25003", "25004", "25005", "25007", "25008", "25009", "25010", "25011", "25012",
               "25013", "25014", "25015", "25016", "25018", "25019", "25020", "25021", "25022", "25023", "25024",
               "25025", "25026", "25027", "25028", "25029", "25030", "25031", "25032", "25033", "25034", "25035",
               "25036", "25037", "25038", "25039", "25040", "25041", "25042", "25043", "25044", "25046", "25047",
               "25048", "25049", "25050", "25052", "25059", "25060", "25061", "25062", "25063", "25064", "25065",
               "25066", "25067", "25068", "25069", "25070", "25071", "25072", "25076", "25077", "25078", "25079",
               "25080", "25081", "25082", "25083", "25084", "25085", "25086", "25087", "25088", "25092", "25093",
               "25094", "25096", "25098", "25099", "25100", "25101", "25102", "25104", "25105", "25106", "25107",
               "25108", "25109", "25110", "25111", "25112", "25113", "25114", "25115", "25116", "25117", "25118",
               "25119", "25120", "25121", "25122", "25123", "25124", "25125", "25127", "25128", "25129", "25130",
               "25131", "25132", "25133", "25134", "25135", "25136", "25137", "25138", "25139", "25140", "25141",
               "25142", "25143", "25145", "25146", "25147", "25148", "25149", "25150", "25151", "25152", "25153",
               "25154", "25155", "25156", "25157", "25158", "25159", "25160", "25161", "25162", "25163", "25164",
               "25165", "25166", "25167", "25168", "25169", "25170", "25171", "25172", "25173", "25174", "25175",
               "25176", "25177", "25178", "25179", "25180", "25181", "25182", "25183", "25184", "25185", "25186",
               "25187", "25188", "25189", "25190", "25191", "25192", "25193", "25194", "25195", "25196", "25197",
               "25199", "25200", "25201", "25202", "25203", "25204", "25205", "25206", "25207", "25208", "25209",
               "25210", "25211", "25212", "25213", "25214", "25215", "25217", "25218", "25219", "25220", "25221",
               "25222", "25223", "25224", "25225", "25226", "25227", "25228", "25229", "25230", "25231", "25232",
               "25233", "25234", "25235", "25236", "25237", "25238", "25239", "25240", "25241", "25242", "25243",
               "25244", "25247", "25248", "25249", "25250", "25251", "25253", "25254", "25255", "25257", "25259",
               "25260", "25261", "25263", "25264", "25265", "25266", "25267", "25268", "25269", "25270", "25271",
               "25272", "25273", "25274", "25275", "25276", "25277", "25278", "25279", "25280", "25281", "25282",
               "25283", "25284", "25285", "25286", "25287", "25288", "25289", "25290", "25291", "25292", "25293",
               "25294", "25295", "25296", "25297", "25298", "25299", "25300", "25301", "25302", "25303", "25304",
               "25305", "25306", "25307", "25308", "25309", "25310", "25311", "25312", "25313", "25314", "25315",
               "25316", "25317", "25318", "25319", "25320", "25321", "25322", "25323", "25324", "25325", "25326",
               "25327", "25328", "25329", "25330", "25331", "25332", "25333", "25334", "25335", "25336", "25337",
               "25338", "25339", "25340", "25341", "25342", "25343", "25344", "25345", "25346", "25347", "25348",
               "25349", "25350", "25351", "25352", "25353", "25354", "25355", "25356", "25357", "25358", "25359",
               "25360", "25361", "25362", "25363", "25364", "25365", "25366", "25367", "25368", "25369", "25370",
               "25371", "25372", "25375", "25376", "25377", "25378", "25379", "25380", "25381", "25382", "25383",
               "25384", "25385", "25386", "25387", "25389", "25390", "25391", "25393", "25394", "25395", "25396",
               "25397", "25398", "25399", "25400", "25401", "25402", "25403", "25404", "25405", "25406", "25407",
               "25408", "25409", "25410", "25411", "25413", "25414", "25415", "25416", "25417", "25418", "25419",
               "25420", "25421", "25422", "25423", "25424", "25425", "25426", "25427", "25428", "25429", "25430",
               "25432", "25433", "25434", "25435", "25436", "25437", "25439", "25440", "25441", "25442", "25443",
               "25445", "25446", "25447", "25448", "25449", "25451", "25452", "25453", "25454", "25455", "25456",
               "25457", "25458", "25459", "25460", "25461", "25462", "25463", "25464", "25465", "25466", "25467",
               "25468", "25469", "25470", "25471", "25472", "25473", "25474", "25475", "25476", "25477", "25478",
               "25479", "25481", "25482", "25483", "25484", "25485", "25486", "25487", "25488", "25489", "25490",
               "25491", "25492", "25493", "25494", "25495", "25496", "25497", "25498", "25499", "25500", "25501",
               "25502", "25504", "25505", "25506", "25507", "25508", "25509", "25510", "25511", "25512", "25513",
               "25514", "25515", "25516", "25517", "25518", "25519", "25520", "25521", "25522", "25523", "25524",
               "25525", "25526", "25527", "25528", "25529", "25530", "25531", "25532", "25533", "25534", "25535",
               "25536", "25537", "25538", "25539", "25540", "25541", "25542", "25543", "25544", "25545", "25546",
               "25548", "25549", "25551", "25552", "25553", "25554", "25555", "25556", "25557", "25558", "25559",
               "25560", "25561", "25562", "25564", "25565", "25566", "25567", "25568", "25569", "25570", "25571",
               "25572", "25573", "25574", "25575", "25576", "25577", "25582", "25584", "25585", "25586", "25587",
               "25588", "25590", "25591", "25592", "25593", "25595", "25596", "25597", "25598", "25599", "25600",
               "25601", "25602", "25603", "25604", "25605", "25606", "25607", "25609", "25610", "25611", "25612",
               "25613", "25614", "25615", "25616", "25617", "25618", "25619", "25620", "25621", "25622", "25623",
               "25624", "25625", "25626", "25627", "25628", "25629", "25630", "25631", "25632", "25634", "25636",
               "25638", "25639", "25640", "25642", "25643", "25645", "25646", "25647", "25648", "25649", "25650",
               "25651", "25652", "25653", "25654", "25655", "25656", "25657", "25658", "25659", "25660", "25661",
               "25662", "25663", "25664", "25665", "25666", "25667", "25668", "25669", "25670", "25671", "25672",
               "25673", "25674", "25675", "25676", "25677", "25678", "25679", "25680", "25681", "25682", "25683",
               "25684", "25685", "25686", "25688", "25689", "25690", "25691", "25692", "25693", "25694", "25696",
               "25697", "25698", "25699", "25700", "25701", "25702", "25703", "25704", "25705", "25706", "25707",
               "25708", "25709", "25710", "25711", "25712", "25713", "25714", "25715", "25716", "25717", "25718",
               "25719", "25720", "25721", "25722", "25723", "25724", "25725", "25726", "25727", "25728", "25729",
               "25730", "25731", "25732", "25734", "25735", "25736", "25737", "25738", "25739", "25740", "25741",
               "25742", "25743", "25744", "25745", "25746", "25748", "25749", "25750", "25751", "25752", "25753",
               "25754", "25755", "25756", "25757", "25758", "25759", "25760", "25761", "25762", "25763", "25764",
               "25765", "25766", "25767", "25768", "25769", "25770", "25771", "25772", "25773", "25774", "25775",
               "25776", "25777", "25778", "25779", "25780", "25781", "25782", "25783", "25784", "25785", "25786",
               "25787", "25788", "25789", "25790", "25791", "25792", "25793", "25794", "25796", "25797", "25798",
               "25799", "25800", "25801", "25802", "25803", "25804", "25805", "25806", "25807", "25808", "25809",
               "25810", "25811", "25813", "25814", "25816", "25817", "25818", "25819", "25820", "25821", "25825",
               "25826", "25827", "25828", "25829", "25830", "25831", "25832", "25833", "25834", "25835", "25836",
               "25837", "25838", "25839", "25840", "25844", "25845", "25846", "25847", "25848", "25849", "25850",
               "25851", "25852", "25853", "25854", "25855", "25856", "25857", "25858", "25859", "25860", "25861",
               "25863", "25864", "25865", "25866", "25867", "25868", "25869", "25870", "25871", "25872", "25873",
               "25874", "25875", "25877", "25878", "25879", "25880", "25881", "25882", "25883", "25885", "25886",
               "25887", "25888", "25889", "25890", "25891", "25892", "25893", "25894", "25895", "25896", "25897",
               "25898", "25899", "25900", "25901", "25902", "25903", "25904", "25905", "25906", "25907", "25908",
               "25909", "25910", "25911", "25912", "25913", "25914", "25915", "25916", "25917", "25918", "25919",
               "25920", "25921", "25922", "25923", "25924", "25925", "25926", "25927", "25928", "25929", "25930",
               "25931", "25932", "25933", "25934", "25935", "25939", "25941", "25942", "25943", "25944", "25945",
               "25946", "25947", "25948", "25949", "25950", "25951", "25952", "25953", "25954", "25955", "25956",
               "25957", "25958", "25959", "25961", "25962", "25963", "25964", "25965", "25966", "25967", "25968",
               "25969", "25970", "25971", "25972", "25973", "25974", "25975", "25976", "25979", "25982", "25983",
               "25984", "25985", "25986", "25987", "25988", "25989", "25990", "25991", "25992", "25993", "25994",
               "25995", "25996", "25997", "25998", "25999", "26000", "26001", "26002", "26003", "26004", "26005",
               "26006", "26007", "26008", "26009", "26010", "26011", "26012", "26016", "26021", "26022", "26024",
               "26025", "26026", "26027", "26028", "26029", "26030", "26031", "26032", "26033", "26034", "26035",
               "26036", "26037", "26038", "26040", "26041", "26042", "26043", "26044", "26045", "26046", "26047",
               "26048", "26049", "26050", "26051", "26052", "26053", "26054", "26055", "26056", "26057", "26058",
               "26059", "26060", "26061", "26062", "26063", "26064", "26065", "26066", "26068", "26069", "26300",
               "26301", "26307", "26308", "26309", "26312", "26314", "26315", "26316", "26317", "26318", "26319",
               "26320", "26321", "26322", "26323", "26328", "26329", "26331", "26333", "26334", "26335", "26336",
               "26503", "26504", "26505", "26506", "26507", "26508", "26509", "26510", "26511", "26512", "26513",
               "26514", "26515", "26516", "26517", "26518", "26519", "26520", "26525", "26526", "26527", "26528",
               "26529", "26530", "26531", "26532", "26533", "26534", "26535", "26536", "26537", "26538", "26541",
               "26543", "26544", "26545", "26546", "26547", "26548", "26549", "26550", "26551", "26552", "26553",
               "26554", "26557", "26558", "26562", "26564", "26565", "26568", "26569", "26570", "26571", "26573",
               "26574", "26575", "26576", "26577", "26578", "26579", "26580", "26582", "26583", "26584", "26585",
               "26586", "26587", "26588", "26589", "26590", "26591", "26592", "26593", "26594", "26595", "26597",
               "26598", "26599", "26600", "26601", "26602", "26603", "26604", "26605", "26607", "26608", "26609",
               "26610", "26611", "26612", "26613", "26614", "26615", "26616", "26617", "26618", "26619", "26620",
               "26623", "26624", "26625", "26626", "26627", "26628", "26629", "26630", "26631", "26632", "26634",
               "26635", "26636", "26637", "26638", "26639", "26640", "26641", "26642", "26643", "26644", "26645",
               "26646", "26647", "26648", "26649", "26650", "26651", "26652", "26653", "26654", "26655", "26656",
               "26657", "26658", "26659", "26660", "26661", "26662", "26663", "26664", "26667", "26670", "26672",
               "26674", "26675", "26676", "26677", "26678", "26679", "26680", "26681", "26682", "26683", "26684",
               "26685", "26686", "26687", "26688", "26689", "26690", "26691", "26692", "26693", "26694", "26695",
               "26697", "26698", "26700", "26701", "26702", "26703", "26705", "26706", "26708", "26709", "26710",
               "26711", "26712", "26713", "26714", "26715", "26716", "26717", "26718", "26719", "26720", "26721",
               "26722", "26723", "26724", "26725", "26726", "26727", "26728", "26730", "26731", "26732", "26733",
               "26734", "26735", "26738", "26739", "26740", "26741", "26742", "26743", "26744", "26745", "26746",
               "26748", "26750", "26751", "26752", "26753", "26754", "26755", "26756", "26757", "26758", "26759",
               "26760", "26761", "26762", "26763", "26764", "26765", "26766", "26767", "26768", "26769", "26770",
               "26771", "26772", "26773", "26774", "26775", "26776", "26777", "26778", "26779", "26780", "26782",
               "26783", "26784", "26785", "26786", "26787", "26788", "26789", "26790", "26791", "26792", "26793",
               "26794", "26795", "26796", "26797", "26798", "26799", "26800", "26801", "26802", "26803", "26804",
               "26805", "26806", "26807", "26808", "26809", "26810", "26811", "26812", "26813", "26814", "26815",
               "26816", "26818", "26819", "26820", "26821", "26822", "26823", "26824", "26825", "26826", "26827",
               "26828", "26829", "26830", "26831", "26832", "26833", "26834", "26835", "26836", "26837", "26838",
               "26839", "26840", "26841", "26842", "26843", "26844", "26845", "26846", "26847", "26848", "26849",
               "26850", "26851", "26852", "26853", "26854", "26855", "26857", "26858", "26860", "26861", "26862",
               "26863", "26864", "26866", "26867", "26870", "26871", "26873", "26874", "26875", "26876", "26877",
               "26878", "26879", "26880", "26881", "26882", "26883", "26884", "26885", "26886", "26887", "26888",
               "26889", "26890", "26891", "26892", "26893", "26894", "26896", "26897", "26898", "26900", "26901",
               "26902", "26903", "26904", "26905", "26906", "26907", "26908", "26909", "26910", "26912", "26913",
               "26914", "26915", "26916", "26919", "26920", "26921", "26922", "26923", "26924", "26925", "26926",
               "26927", "26928", "26929", "26931", "26932", "26933", "26934", "26935", "26936", "26937", "26938",
               "26939", "26940", "26942", "26943", "26944", "26945", "26946", "26947", "26948", "26949", "26950",
               "26951", "26952", "26953", "26954", "26955", "26956", "26957", "26958", "26959", "26960", "26961",
               "26962", "26963", "26964", "26965", "26966", "26967", "26968", "26969", "26970", "26971", "26972",
               "26973", "26974", "26975", "26976", "26977", "26978", "26979", "26980", "26981", "26982", "26983",
               "26984", "26985", "26986", "26987", "26989", "26990", "26991", "26992", "26993", "26994", "26995",
               "26996", "26997", "26998", "26999", "27000", "27001", "27002", "27003", "27004", "27005", "27007",
               "27008", "27010", "27011", "27012", "27013", "27014", "27015", "27016", "27017", "27019", "27020",
               "27021", "27022", "27023", "27024", "27025", "27026", "27027", "27028", "27029", "27030", "27031",
               "27032", "27033", "27034", "27035", "27036", "27037", "27038", "27039", "27040", "27041", "27042",
               "27043", "27044", "27045", "27046", "27047", "27048", "27049", "27050", "27051", "27052", "27053",
               "27055", "27056", "27057", "27058", "27059", "27060", "27061", "27062", "27063", "27064", "27065",
               "27066", "27067", "27069", "27070", "27071", "27072", "27073", "27074", "27075", "27076", "27077",
               "27078", "27079", "27080", "27081", "27082", "27083", "27084", "27085", "27086", "27087", "27088",
               "27089", "27090", "27091", "27092", "27093", "27094", "27095", "27096", "27097", "27098", "27099",
               "27102", "27103", "27104", "27106", "27107", "27108", "27109", "27110", "27111", "27112", "27113",
               "27119", "27120", "27121", "27122", "27123", "27124", "27125", "27126", "27127", "27128", "27129",
               "27130", "27131", "27133", "27134", "27135", "27136", "27137", "27138", "27139", "27140", "27141",
               "27142", "27143", "27144", "27145", "27146", "27147", "27148", "27151", "27152", "27153", "27154",
               "27155", "27156", "27157", "27158", "27159", "27160", "27161", "27162", "27163", "27164", "27165",
               "27166", "27167", "27168", "27169", "27170", "27171", "27172", "27173", "27174", "27175", "27176",
               "27177", "27179", "27180", "27183", "27184", "27185", "27186", "27187", "27188", "27189", "27190",
               "27191", "27192", "27193", "27194", "27195", "27196", "27197", "27198", "27199", "27201", "27202",
               "27203", "27204", "27205", "27206", "27207", "27208", "27209", "27210", "27211", "27212", "27213",
               "27214", "27215", "27216", "27217", "27218", "27219", "27220", "27221", "27224", "27225", "27226",
               "27227", "27228", "27229", "27230", "27231", "27232", "27233", "27234", "27235", "27236", "27237",
               "27238", "27239", "27240", "27241", "27242", "27243", "27244", "27245", "27246", "27247", "27248",
               "27249", "27250", "27251", "27252", "27254", "27255", "27256", "27257", "27258", "27259", "27260",
               "27261", "27262", "27263", "27264", "27265", "27266", "27267", "27268", "27269", "27270", "27271",
               "27272", "27273", "27274", "27275", "27276", "27277", "27278", "27279", "27280", "27281", "27282",
               "27284", "27285", "27286", "27287", "27288", "27289", "27290", "27292", "27293", "27294", "27295",
               "27296", "27297", "27298", "27299", "27300", "27301", "27302", "27303", "27304", "27305", "27306",
               "27307", "27312", "27313", "27314", "27315", "27316", "27317", "27320", "27321", "27322", "27323",
               "27324", "27325", "27326", "27327", "27328", "27330", "27331", "27332", "27333", "27334", "27336",
               "27337", "27338", "27339", "27341", "27342", "27343", "27344", "27345", "27346", "27347", "27348",
               "27349", "27351", "27352", "27353", "27354", "27355", "27356", "27357", "27358", "27359", "27360",
               "27361", "27362", "27363", "27364", "27365", "27366", "27367", "27368", "27369", "27370", "27371",
               "27372", "27373", "27374", "27375", "27376", "27377", "27378", "27379", "27380", "27381", "27382",
               "27383", "27384", "27385", "27386", "27387", "27388", "27389", "27390", "27391", "27392", "27393",
               "27394", "27395", "27396", "27397", "27398", "27399", "27400", "27401", "27402", "27403", "27404",
               "27405", "27406", "27407", "27408", "27409", "27410", "27411", "27412", "27413", "27414", "27415",
               "27416", "27417", "27418", "27419", "27420", "27421", "27422", "27424", "27425", "27426", "27427",
               "27428", "27429", "27430", "27431", "27432", "27433", "27434", "27435", "27436", "27437", "27438",
               "27439", "27440", "27441", "27442", "27443", "27445", "27446", "27447", "27448", "27449", "27450",
               "27451", "27452", "27453", "27455", "27457", "27458", "27459", "27460", "27461", "27462", "27463",
               "27464", "27465", "27466", "27467", "27468", "27469", "27470", "27471", "27472", "27473", "27475",
               "27476", "27477", "27478", "27479", "27480", "27481", "27482", "27483", "27484", "27485", "27488",
               "27489", "27490", "27491", "27492", "27493", "27494", "27495", "27496", "27497", "27498", "27499",
               "27500", "27501", "27502", "27503", "27506", "27508", "27509", "27510", "27511", "27512", "27513",
               "27514", "27515", "27516", "27517", "27518", "27519", "27520", "27521", "27522", "27523", "27524",
               "27525", "27526", "27527", "27528", "27529", "27530", "27531", "27532", "27533", "27534", "27536",
               "27537", "27538", "27539", "27540", "27541", "27542", "27543", "27544", "27545", "27546", "27547",
               "27548", "27549", "27550", "27552", "27553", "27554", "27555", "27556", "27557", "27558", "27559",
               "27560", "27561", "27562", "27563", "27564", "27565", "27566", "27567", "27569", "27570", "27571",
               "27572", "27573", "27574", "27575", "27576", "27577", "27578", "27579", "27580", "27581", "27582",
               "27583", "27584", "27585", "27586", "27587", "27588", "27589", "27590", "27591", "27592", "27593",
               "27594", "27595", "27596", "27597", "27598", "27599", "27600", "27601", "27602", "27603", "27605",
               "27606", "27607", "27608", "27609", "27610", "27611", "27612", "27613", "27614", "27616", "27617",
               "27618", "27619", "27620", "27621", "27622", "27623", "27624", "27625", "27626", "27627", "27628",
               "27629", "27630", "27631", "27632", "27633", "27634", "27635", "27636", "27637", "27638", "27639",
               "27640", "27641", "27642", "27643", "27645", "27646", "27647", "27648", "27649", "27651", "27652",
               "27653", "27654", "27655", "27656", "27659", "27661", "27662", "27663", "27664", "27666", "27667",
               "27668", "27669", "27670", "27671", "27672", "27673", "27674", "27676", "27677", "27678", "27679",
               "27680", "27681", "27682", "27683", "27684", "27687", "27688", "27689", "27690", "27691", "27692",
               "27693", "27694", "27695", "27696", "27697", "27698", "27699", "27700", "27701", "27702", "27703",
               "27704", "27705", "27706", "27707", "27708", "27709", "27711", "27712", "27713", "27714", "27715",
               "27716", "27717", "27718", "27719", "27720", "27721", "27722", "27724", "27725", "27726", "27727",
               "27728", "27729", "27730", "27731", "27732", "27733", "27734", "27737", "27738", "27739", "27740",
               "27741", "27742", "27743", "27744", "27745", "27746", "27747", "27750", "27751", "27752", "27753",
               "27754", "27755", "27756", "27757", "27759", "27761", "27762", "27763", "27764", "27767", "27768",
               "27769", "27770", "27771", "27772", "27774", "27775", "27776", "27777", "27778", "27779", "27780",
               "27781", "27782", "27783", "27784", "27785", "27786", "27787", "27788", "27789", "27790", "27791",
               "27792", "27793", "27794", "27795", "27796", "27797", "27798", "27799", "27800", "27801", "27802",
               "27803", "27804", "27805", "27806", "27807", "27808", "27809", "27810", "27811", "27812", "27813",
               "27814", "27815", "27816", "27817", "27818", "27819", "27820", "27821", "27825", "27826", "27827",
               "27828", "27829", "27830", "27831", "27832", "27833", "27834", "27835", "27836", "27837", "27838",
               "27839", "27840", "27841", "27842", "27843", "27844", "27845", "27847", "27848", "27850", "27851",
               "27852", "27853", "27854", "27856", "27857", "27858", "27859", "27860", "27861", "27862", "27863",
               "27864", "27865", "27866", "27867", "27868", "27869", "27870", "27872", "27873", "27874", "27875",
               "27876", "27877", "27878", "27879", "27880", "27881", "27882", "27883", "27884", "27885", "27886",
               "27887", "27888", "27890", "27893", "27894", "27895", "27900", "27901", "27902", "27903", "27904",
               "27905", "27911", "27912", "27913", "27914", "27915", "27916", "27917", "27920", "27922", "27923",
               "27924", "27925", "27926", "27927", "27928", "27929", "27930", "27931", "27932", "27934", "27935",
               "27936", "27937", "27938", "27939", "27943", "27944", "27945", "27946", "27947", "27948", "27949",
               "27950", "27951", "27952", "27953", "27954", "27955", "27956", "27957", "27958", "27959", "27960",
               "27961", "27962", "27963", "27964", "27965", "27967", "27968", "27969", "27970", "27971", "27972",
               "27973", "27974", "27975", "27976", "27977", "27978", "27979", "27980", "27981", "27982", "27983",
               "27984", "27985", "27986", "27987", "27988", "27989", "27990", "27991", "27992", "27993", "27994",
               "27995", "27996", "27997", "27998", "27999", "28000", "28001", "28002", "28003", "28005", "28006",
               "28008", "28009", "28011", "28012", "28014", "28015", "28016", "28017", "28019", "28020", "28021",
               "28022", "28025", "28026", "28027", "28028", "28029", "28030", "28031", "28032", "28033", "28034",
               "28035", "28036", "28037", "28038", "28039", "28040", "28041", "28042", "28043", "28045", "28046",
               "28047", "28048", "28049", "28050", "28051", "28052", "28053", "28054", "28055", "28056", "28057",
               "28058", "28059", "28060", "28061", "28062", "28063", "28064", "28065", "28066", "28069", "28070",
               "28071", "28072", "28073", "28074", "28075", "28076", "28077", "28079", "28080", "28081", "28082",
               "28084", "28085", "28086", "28087", "28088", "28089", "28090", "28091", "28092", "28094", "28095",
               "28096", "28097", "28098", "28099", "28100", "28101", "28102", "28103", "28104", "28105", "28106",
               "28107", "28108", "28109", "28110", "28111", "28112", "28113", "28114", "28115", "28116", "28117",
               "28118", "28121", "28122", "28123", "28124", "28125", "28126", "28127", "28128", "28129", "28130",
               "28131", "28132", "28133", "28134", "28135", "28136", "28137", "28138", "28139", "30000", "30002",
               "30003", "30004", "30005", "30006", "30007", "30008", "30009", "30010", "30011", "30012", "30013",
               "30015", "30016", "30017", "30019", "30020", "30021", "30022", "30023", "30024", "30025", "30026",
               "30027", "30028", "30029", "30030", "30031", "30032", "30033", "30034", "30035", "30037", "30038",
               "30039", "30040", "30042", "30043", "30044", "30045", "30046", "30047", "30048", "30049", "30050",
               "30051", "30052", "30053", "30054", "30055", "30056", "30057", "30058", "30059", "30060", "30061",
               "30062", "30063", "30064", "30065", "30066", "30067", "30068", "30069", "30070", "30071", "30072",
               "30073", "30074", "30075", "30076", "30077", "30078", "30079", "30080", "30081", "30082", "30083",
               "30084", "30085", "30086", "30087", "30088", "30089", "30090", "30091", "30092", "30093", "30094",
               "30097", "30098", "30099", "30100", "30101", "30102", "30105", "30106", "30107", "30108", "30109",
               "30110", "30111", "30112", "30113", "30114", "30115", "30116", "30117", "30118", "30119", "30120",
               "30121", "30123", "30124", "30126", "30127", "30128", "30129", "30130", "30131", "30132", "30133",
               "30134", "30135", "30136", "30137", "30138", "30139", "30140", "30141", "30142", "30143", "30144",
               "30145", "30146", "30147", "30148", "30149", "30150", "30151", "30152", "30153", "30154", "30155",
               "30156", "30157", "30158", "30159", "30160", "30161", "30162", "30163", "30164", "30165", "30170",
               "30171", "30172", "30176", "30177", "30178", "30179", "30180", "30181", "30184", "30185", "30186",
               "30188", "30189", "30190", "30191", "30192", "30193", "30194", "30195", "30196", "30197", "30198",
               "30199", "30200", "30201", "30202", "30203", "30204", "30205", "30206", "30207", "30208", "30209",
               "30210", "30211", "30212", "30213", "30214", "30215", "30216", "30217", "30218", "30219", "30220",
               "30221", "30222", "30223", "30224", "30225", "30226", "30227", "30229", "30230", "30231", "30232",
               "30233", "30234", "30235", "30236", "30237", "30238", "30239", "30240", "30241", "30242", "30243",
               "30244", "30245", "30246", "30247", "30248", "30249", "30250", "30251", "30252", "30253", "30254",
               "30255", "30256", "30257", "30258", "30259", "30260", "30261", "30262", "30263", "30265", "30267",
               "30268", "30270", "30271", "30273", "30274", "30275", "30276", "30281", "30282", "30283", "30284",
               "30285", "30286", "30287", "30288", "30289", "30290", "30291", "30292", "30293", "30295", "30296",
               "30297", "30298", "30300", "30301", "30303", "30304", "30305", "30306", "30307", "30308", "30309",
               "30310", "30311", "30312", "30313", "30314", "30315", "30316", "30317", "30318", "30319", "30320",
               "30321", "30322", "30323", "30324", "30325", "30326", "30327", "30328", "30329", "30330", "30331",
               "30332", "30333", "30334", "30335", "30336", "30337", "30338", "30339", "30340", "30341", "30342",
               "30343", "30344", "30345", "30346", "30347", "30348", "30349", "30350", "30351", "30352", "30353",
               "30354", "30355", "30356", "30357", "30358", "30359", "30360", "30361", "30362", "30363", "30364",
               "30365", "30366", "30367", "30368", "30370", "30371", "30372", "30373", "30374", "30375", "30376",
               "30377", "30378", "30379", "30380", "30381", "30382", "30383", "30385", "30386", "30388", "30389",
               "30390", "30391", "30394", "30395", "30396", "30397", "30398", "30399", "30400", "30401", "30402",
               "30403", "30404", "30405", "30406", "30407", "30408", "30409", "30410", "30411", "30412", "30413",
               "30414", "30415", "30416", "30417", "30418", "30419", "30420", "30421", "30422", "30423", "30424",
               "30425", "30426", "30427", "30428", "30429", "30430", "30431", "30432", "30433", "30434", "30435",
               "30436", "30437", "30438", "30439", "30440", "30441", "30442", "30443", "30444", "30445", "30446",
               "30447", "30448", "30449", "30450", "30451", "30452", "30453", "30454", "30456", "30457", "30458",
               "30459", "30460", "30461", "30462", "30463", "30464", "30465", "30466", "30467", "30468", "30469",
               "30470", "30471", "30472", "30473", "30474", "30475", "30476", "30477", "30478", "30479", "30480",
               "30481", "30482", "30483", "30484", "30485", "30486", "30487", "30488", "30489", "30490", "30491",
               "30492", "30493", "30494", "30495", "30496", "30497", "30498", "30499", "30500", "30501", "30502",
               "30503", "30504", "30505", "30506", "30507", "30508", "30509", "30510", "30511", "30512", "30513",
               "30514", "30515", "30516", "30517", "30518", "30519", "30520", "30521", "30522", "30523", "30524",
               "30525", "30527", "30528", "30529", "30530", "30531", "30532", "30533", "30534", "30535", "30536",
               "30537", "30538", "30543", "30544", "30545", "30546", "30547", "30548", "30550", "30551", "30552",
               "30553", "30554", "30555", "30556", "30557", "30558", "30559", "30560", "30561", "30562", "30565",
               "30566", "30567", "30568", "30569", "30570", "30571", "30572", "30573", "30574", "30575", "30576",
               "30577", "30579", "30580", "30583", "30584", "30585", "30586", "30587", "30588", "30590", "30591",
               "30592", "30593", "30594", "30595", "30596", "30597", "30598", "30599", "30600", "30601", "30602",
               "30603", "30604", "30605", "30606", "30607", "30608", "30609", "30610", "30611", "30613", "30614",
               "30615", "30616", "30617", "30618", "30619", "30620", "30621", "30622", "30623", "30625", "30626",
               "30627", "30628", "30629", "30630", "30631", "30632", "30633", "30634", "30635", "30636", "30637",
               "30638", "30639", "30640", "30641", "30642", "30643", "30644", "30645", "30646", "30647", "30648",
               "30650", "30651", "30652", "30653", "30654", "30655", "30656", "30657", "30658", "30659", "30660",
               "30661", "30662", "30663", "30664", "30665", "30666", "30667", "30668", "30669", "30670", "30671",
               "30672", "30673", "30674", "30675", "30676", "30677", "30678", "30679", "30680", "30681", "30682",
               "30683", "30684", "30686", "30687", "30688", "30690", "30691", "30692", "30693", "30694", "30695",
               "30696", "30697", "30698", "30699", "30700", "30701", "30702", "30703", "30704", "30705", "30706",
               "30707", "30708", "30710", "30711", "30712", "30713", "30714", "30716", "30717", "30718", "30719",
               "30720", "30721", "30722", "30723", "30724", "30725", "30726", "30727", "30728", "30729", "30730",
               "30731", "30734", "30735", "30736", "30737", "30738", "30739", "30740", "30741", "30744", "30745",
               "30746", "30747", "30748", "30749", "30750", "30752", "30753", "30754", "30755", "30756", "30757",
               "30758", "30759", "30760", "30761", "30763", "30765", "30767", "30768", "30769", "30770", "30771",
               "30772", "30773", "30774", "30775", "30776", "30777", "30778", "30779", "30782", "30783", "30784",
               "30785", "30786", "30787", "30788", "30789", "30790", "30791", "30792", "30793", "30795", "30796",
               "30797", "30798", "30799", "30800", "30801", "30802", "30803", "30804", "30805", "30807", "30808",
               "30809", "30810", "30811", "30812", "30813", "30814", "30816", "30817", "30818", "30819", "30820",
               "30821", "30822", "30823", "30824", "30825", "30826", "30827", "30828", "30829", "30830", "30832",
               "30833", "30834", "30838", "30839", "30840", "30842", "30843", "30844", "30845", "30846", "30847",
               "30848", "30849", "30850", "30851", "30852", "30853", "30854", "30855", "30856", "30857", "30858",
               "30861", "30863", "30864", "30866", "30867", "30868", "30869", "30870", "30871", "30872", "30873",
               "30874", "30875", "30876", "30877", "30878", "30879", "30880", "30881", "30882", "30883", "30884",
               "30885", "30886", "30887", "30888", "30889", "30890", "30891", "30892", "30893", "30894", "30895",
               "30896", "30899", "30900", "30901", "30902", "30903", "30904", "30905", "30906", "30907", "30908",
               "30909", "30910", "30911", "30912", "30913", "30915", "30916", "30917", "30918", "30919", "30920",
               "30921", "30922", "30923", "30924", "30925", "30926", "30927", "30928", "30929", "30930", "30931",
               "30932", "30933", "30934", "30935", "30936", "30937", "30938", "30939", "30940", "30941", "30942",
               "30943", "30944", "30945", "30946", "30947", "30948", "30949", "30950", "30951", "30952", "30953",
               "30954", "30955", "30956", "30957", "30958", "30959", "30960", "30961", "30962", "30963", "30965",
               "30966", "30967", "30968", "30969", "30970", "30971", "30972", "30973", "30974", "30975", "30976",
               "30977", "30978", "30979", "30980", "30981", "30982", "30983", "30984", "30985", "30986", "30987",
               "30988", "30989", "30990", "30991", "30992", "30993", "30994", "30995", "30997", "30998", "30999",
               "31000", "31001", "31002", "31003", "31004", "31005", "31006", "31007", "31008", "31009", "31010",
               "31011", "31012", "31013", "31014", "31015", "31016", "31017", "31018", "31019", "31020", "31021",
               "31022", "31023", "31024", "31025", "31026", "31027", "31028", "31029", "31030", "31031", "31032",
               "31033", "31034", "31035", "31036", "31037", "31038", "31039", "31040", "31043", "31045", "31047",
               "31048", "31049", "31050", "31051", "31052", "31053", "31054", "31055", "31056", "31057", "31058",
               "31059", "31060", "31061", "31062", "31063", "31064", "31065", "31066", "31068", "31069", "31070",
               "31071", "31072", "31073", "31074", "31075", "31076", "31077", "31078", "31079", "31085", "31086",
               "31089", "31090", "31091", "31094", "31095", "31096", "31097", "31099", "31100", "31101", "31102",
               "31103", "31104", "31106", "31107", "31111", "31112", "31113", "31124", "31126", "31134", "34000",
               "34001", "34002", "34003", "34004", "34005", "34006", "34007", "34008", "34009", "34010", "34011",
               "34012", "34013", "34014", "34015", "34016", "34017", "34018", "34019", "34022", "34024", "34025",
               "34026", "34027", "34028", "34029", "34030", "34031", "34032", "34033", "34034", "34035", "34036",
               "34037", "34038", "34039", "34040", "34041", "34042", "34043", "34044", "34045", "34046", "34047",
               "34048", "34049", "34050", "34051", "34052", "34053", "34054", "34055", "34056", "34057", "34058",
               "34059", "34060", "34061", "34062", "34063", "34064", "34065", "34066", "34067", "34068", "34069",
               "34070", "34071", "34072", "34073", "34074", "34075", "34076", "34077", "34078", "34079", "34080",
               "34081", "34082", "34083", "34084", "34085", "34086", "34087", "34088", "34089", "34090", "34091",
               "34092", "34093", "34094", "34095", "34098", "34099", "34100", "34101", "34102", "34103", "34104",
               "34105", "34106", "34108", "34109", "34110", "34111", "34112", "34113", "34114", "34115", "34116",
               "34118", "34119", "34120", "34121", "34122", "34123", "34124", "34125", "34126", "34127", "34129",
               "34130", "34131", "34132", "34133", "34134", "34135", "34136", "34137", "34138", "34139", "34140",
               "34141", "34142", "34143", "34144", "34145", "34146", "34149", "34151", "34152", "34153", "34154",
               "34155", "34157", "34158", "34159", "34160", "34161", "34162", "34163", "34164", "34165", "34166",
               "34167", "34168", "34169", "34170", "34171", "34172", "34173", "34174", "34175", "34176", "34178",
               "34179", "34180", "34181", "34182", "34183", "34184", "34185", "34186", "34187", "34188", "34189",
               "34190", "34191", "34193", "34194", "34195", "34196", "34198", "34199", "34200", "34201", "34202",
               "34203", "34206", "34207", "34208", "34209", "34210", "34211", "34212", "34213", "34214", "34215",
               "34216", "34217", "34218", "34219", "34220", "34221", "34222", "34223", "34224", "34225", "34226",
               "34227", "34228", "34229", "34231", "34232", "34233", "34234", "34235", "34236", "34237", "34238",
               "34239", "34240", "34243", "34244", "34245", "34246", "34247", "34248", "34249", "34250", "34251",
               "34252", "34253", "34254", "34255", "34256", "34257", "34258", "34259", "34260", "34261", "34262",
               "34263", "34264", "34265", "34266", "34267", "34268", "34269", "34270", "34271", "34272", "34273",
               "34274", "34276", "34277", "34278", "34279", "34280", "34281", "34282", "34283", "34284", "34285",
               "34286", "34287", "34288", "34289", "34290", "34291", "34292", "34293", "34294", "34295", "34296",
               "34297", "34298", "34299", "34300", "34301", "34302", "34303", "34304", "34305", "34306", "34307",
               "34308", "34309", "34311", "34312", "34313", "34314", "34315", "34316", "34317", "34318", "34319",
               "34320", "34321", "34322", "34323", "34324", "34325", "34326", "34328", "34329", "34330", "34331",
               "34332", "34333", "34334", "34335", "34336", "34338", "34339", "34340", "34341", "34342", "34343",
               "34344", "34345", "34346", "34347", "34348", "34349", "34350", "34351", "34352", "34353", "34354",
               "34355", "34356", "34357", "34358", "34359", "34360", "34361", "34362", "34363", "34365", "34366",
               "34367", "34369", "34370", "34371", "34372", "34373", "34374", "34376", "34378", "34379", "34380",
               "34381", "34383", "34384", "34385", "34386", "34387", "34388", "34389", "34390", "34391", "34392",
               "34393", "34394", "34395", "34396", "34397", "34398", "34399", "34400", "34401", "34402", "34403",
               "34404", "34405", "34406", "34407", "34408", "34409", "34410", "34411", "34412", "34413", "34414",
               "34416", "34417", "34418", "34419", "34420", "34421", "34422", "34423", "34424", "34425", "34427",
               "34428", "34429", "34430", "34431", "34432", "34433", "34434", "34435", "34436", "34437", "34438",
               "34439", "34440", "34441", "34442", "34443", "34444", "34445", "34446", "34447", "34448", "34449",
               "34450", "34451", "34452", "34453", "34454", "34455", "34456", "34457", "34459", "34460", "34461",
               "34462", "34463", "34464", "34465", "34466", "34467", "34468", "34469", "34470", "34471", "34472",
               "34473", "34474", "34475", "34476", "34477", "34478", "34479", "34480", "34481", "34482", "34483",
               "34484", "34485", "34486", "34487", "34488", "34489", "34490", "34491", "34492", "34493", "34494",
               "34495", "34496", "34497", "34498", "34499", "34500", "34502", "34503", "34504", "34506", "34507",
               "34508", "34509", "34510", "34511", "34512", "34513", "34514", "34515", "34516", "34517", "34518",
               "34519", "34520", "34521", "34522", "34523", "34524", "34525", "34526", "34527", "34528", "34529",
               "34531", "34532", "34533", "34534", "34535", "34536", "34537", "34538", "34539", "34540", "34541",
               "34542", "34543", "34544", "34545", "34546", "34547", "34548", "34549", "34550", "34551", "34552",
               "34553", "34554", "34555", "34556", "34557", "34558", "34559", "34560", "34564", "34565", "34566",
               "34567", "34568", "34569", "34570", "34571", "34572", "34573", "34574", "34575", "34576", "34577",
               "34578", "34579", "34580", "34581", "34582", "34583", "34584", "34585", "34586", "34587", "34588",
               "34589", "34590", "34591", "34592", "34593", "34594", "34595", "34596", "34598", "34599", "34600",
               "34602", "34603", "34604", "34605", "34606", "34607", "34608", "34609", "34610", "34611", "34612",
               "34613", "34614", "34615", "34616", "34617", "34618", "34619", "34620", "34621", "34622", "34623",
               "34624", "34625", "34626", "34627", "34628", "34629", "34630", "34631", "34632", "34633", "34634",
               "34635", "34636", "34637", "34638", "34639", "34640", "34641", "34642", "34643", "34645", "34646",
               "34647", "34648", "34649", "34650", "34651", "34652", "34653", "34654", "34655", "34656", "34657",
               "34660", "34661", "34662", "34663", "34664", "34665", "34666", "34667", "34668", "34669", "34670",
               "34671", "34673", "34674", "34675", "34676", "34677", "34678", "34679", "34680", "34681", "34682",
               "34683", "34685", "34686", "34688", "34689", "34691", "34692", "34693", "34694", "34695", "34696",
               "34697", "34698", "34699", "34700", "34701", "34702", "34703", "34704", "34705", "34706", "34707",
               "34708", "34709", "34710", "34711", "34712", "34713", "34714", "34715", "34716", "34717", "34718",
               "34719", "34720", "34721", "34722", "34723", "34724", "34725", "34726", "34727", "34728", "34729",
               "34730", "34731", "34732", "34734", "34735", "34736", "34737", "34738", "34739", "34740", "34741",
               "34742", "34743", "34744", "34745", "34746", "34747", "34748", "34749", "34750", "34751", "34752",
               "34753", "34754", "34755", "34756", "34757", "34758", "34759", "34760", "34762", "34763", "34764",
               "34765", "34766", "34767", "34768", "34769", "34770", "34771", "34772", "34773", "34774", "34775",
               "34776", "34777", "34778", "34779", "34780", "34781", "34782", "34783", "34784", "34785", "34786",
               "34787", "34790", "34791", "34792", "34793", "34794", "34795", "34796", "34797", "34798", "34799",
               "34800", "34802", "34803", "34804", "34805", "34806", "34807", "34808", "34812", "34816", "34817",
               "34818", "34826", "34831", "34832", "34833", "34834", "34835", "34836", "34837", "34838", "34843",
               "34844", "34845", "34846", "34848", "34849", "34850", "34863", "34864", "34865", "34866", "34876",
               "34877", "34878", "34879", "34880", "34881", "34882", "34883", "34885", "34887", "34895", "36000",
               "36001", "36005", "36006", "36007", "36008", "36009", "36011", "36012", "36013", "36014", "36016",
               "36017", "36018", "36019", "36020", "36021", "36022", "36023", "36024", "36025", "36026", "36027",
               "36034", "36037", "36038", "36040", "36041", "36044", "36045", "36047", "36048", "36049", "36050",
               "36051", "36052", "36053", "36054", "36055", "36056", "36057", "36058", "36059", "36060", "36061",
               "36062", "36064", "36066", "36068", "36069", "36070", "36071", "36072", "36073", "36074", "36075",
               "36077", "36078", "36079", "36080", "36081", "36082", "36083", "36084", "36085", "36086", "36087",
               "36096", "36097", "36098", "36101", "36103", "36105", "36106", "36107", "36109", "36110", "36111",
               "36112", "36114", "36115", "36116", "36117", "36119", "36124", "36125", "36126", "36129", "36133",
               "36136", "36142", "36143", "36144", "36146", "36147", "36148", "36149", "36150", "36157", "36158",
               "36159", "36160", "36161", "36162", "36163", "36165", "36166", "36167", "36168", "36170", "36171",
               "36172", "36174", "36175", "36176", "36177", "36179", "36180", "36181", "36185", "36186", "36187",
               "36207", "36220", "36221", "36228", "36243", "36244", "36258", "36263", "36264", "36265", "36266",
               "36273", "36284", "36288", "36293", "36294", "36299", "36309", "36320", "36321", "36327", "36328",
               "36329", "36330", "36331", "36332", "36333", "36334", "36335", "36336", "36337", "36338", "36339",
               "36341", "36342", "36343", "36345", "36368", "36378", "36379", "36380", "36385", "36404", "36405",
               "36411", "36417", "36419", "36422", "36426", "36427", "36430", "36431", "36445", "36447", "36471",
               "36473", "36489", "36515", "36516", "36589", "50000", "50001", "50002", "50003", "50004", "50007",
               "50008", "50009", "50010", "50011", "50012", "50013", "50014", "50015", "50017", "50018", "50019",
               "50020", "50021", "50022", "50023", "50024", "50025", "50026", "50027", "50028", "50029", "50030",
               "50031", "50032", "50033", "50034", "50035", "50036", "50037", "50038", "50039", "50040", "50041",
               "50042", "50043", "50044", "50045", "50046", "50047", "50048", "50049", "50050", "50051", "50052",
               "50053", "50054", "50055", "50056", "50057", "50058", "50059", "50060", "50061", "50062", "50063",
               "50064", "50065", "50066", "50067", "50068", "50069", "50070", "50071", "50072", "50073", "50074",
               "50075", "50077", "50078", "50079", "50080", "50082", "50084", "50085", "50086", "50088", "50089",
               "50090", "50091", "50093", "50094", "50095", "50096", "50097", "50098", "50099", "50100", "50101",
               "50102", "50103", "50104", "50105", "50106", "50107", "50108", "50109", "50110", "50113", "50114",
               "50115", "50117", "50118", "50119", "50120", "50121", "50122", "50123", "50124", "50125", "50126",
               "50127", "50128", "50129", "50130", "50131", "50132", "50133", "50134", "50135", "50138", "50139",
               "50141", "50142", "50143", "50145", "50146", "50147", "50148", "50149", "50151", "50152", "50153",
               "50154", "50156", "50157", "50158", "50159", "50160", "50161", "50162", "50163", "50164", "50165",
               "50166", "50167", "50168", "50169", "50172", "50173", "50174", "50175", "50176", "50177", "50178",
               "50179", "50180", "50181", "50182", "50183", "50185", "50186", "50187", "50188", "50189", "50190",
               "50191", "50192", "50193", "50195", "50196", "50198", "50199", "50201", "50202", "50203", "50204",
               "50205", "50206", "50207", "50208", "50210", "50211", "50212", "50213", "50214", "50215", "50216",
               "50217", "50218", "50219", "50220", "50221", "50228", "50229", "50231", "50233", "50234", "50235",
               "50236", "50237", "50238", "50239", "50240", "50241", "50242", "50243", "50244", "50245", "50246",
               "50247", "50248", "50249", "50250", "50251", "50253", "50254", "50255", "50256", "50257", "50258",
               "50259", "50260", "50261", "50262", "50263", "50264", "50265", "50267", "50268", "50269", "50270",
               "50271", "50273", "50274", "50275", "50276", "50277", "50278", "50283", "50284", "50285", "50297",
               "50298", "50299", "50300", "50301", "50302", "50304", "50305", "50308", "50309", "50310", "50312",
               "50313", "50314", "50315", "50316", "50317", "50318", "50319", "50320", "50321", "50322", "50323",
               "50324", "50325", "50326", "50328", "50329", "50330", "50331", "50332", "50333", "50334", "50335",
               "50336", "50337", "50338", "50339", "50340", "50341", "50342", "50343", "50344", "50345", "50346",
               "50347", "50348", "50349", "50350", "50351", "50352", "50353", "50355", "50356", "50357", "50360",
               "50361", "50362", "50366", "50368", "50369", "50370", "50371", "50372", "50373", "50374", "50375",
               "50376", "50377", "50378", "50379", "50380", "50381", "50382", "50383", "50384", "50385", "50386",
               "50387", "50388", "50389", "50391", "50392", "50393", "50394", "50395", "50396", "50397", "50398",
               "50399", "50400", "50401", "50402", "50403", "50404", "50405", "50407", "50408", "50409", "50410",
               "50411", "50412", "50414", "50415", "50416", "50417", "50418", "50421", "50422", "50423", "50424",
               "50425", "50426", "50427", "50428", "50430", "50431", "50432", "50433", "50434", "50437", "50438",
               "50443", "50444", "50445", "50446", "50447", "50448", "50449", "50450", "50451", "50452", "50453",
               "50454", "50455", "50456", "50457", "50458", "50459", "50460", "50461", "50463", "50464", "50465",
               "50466", "50467", "50468", "50469", "50470", "50471", "50472", "50473", "50474", "50476", "50477",
               "50478", "50479", "50480", "50481", "50482", "50484", "50485", "50486", "50487", "50488", "50489",
               "50490", "50491", "50492", "50493", "50494", "50495", "50496", "50497", "50507", "50509", "50512",
               "50513", "50515", "50516", "50517", "50518", "50519", "50520", "50521", "50522", "50523", "50524",
               "50525", "50526", "50527", "50528", "50529", "50530", "50531", "50532", "50533", "50534", "50535",
               "50536", "50537", "50538", "50539", "50540", "50541", "50542", "50543", "50544", "50545", "50546",
               "50548", "50549", "50550", "50552", "50553", "50554", "50555", "50556", "50557", "50558", "50559",
               "50560", "50561", "50563", "50564", "50565", "50566", "50567", "50568", "50569", "50570", "50571",
               "50572", "50574", "50575", "50576", "50584", "50585", "50591", "50592", "50593", "50594", "50595",
               "50596", "50598", "50599", "50600", "50601", "50602", "50603", "50604", "50606", "50607", "50608",
               "50609", "50610", "50611", "50612", "50613", "50614", "50615", "50616", "50617", "50618", "50619",
               "50620", "50621", "50622", "50625", "50626", "50627", "50628", "50629", "50630", "50631", "50632",
               "50633", "50634", "50635", "50636", "50637", "50638", "50639", "50640", "50641", "50642", "50643",
               "50644", "50645", "50646", "50647", "50648", "50649", "50650", "50651", "50652", "50653", "50654",
               "50657", "50658", "50659", "50660", "50661", "50662", "50663", "50664", "50665", "50666", "50667",
               "50668", "50669", "50670", "50671", "50672", "50673", "50674", "50676", "50678", "50679", "50680",
               "50681", "50682", "50683", "50686", "50687", "50688", "50689", "50691", "50694", "50695", "50696",
               "50697", "50698", "50699", "50700", "50701", "50702", "50707", "50708", "50709", "50710", "50711",
               "50712", "50713", "50714", "50715", "50718", "50719", "50720", "50721", "50722", "50723", "50725",
               "50730", "50731", "50732", "50733", "50734", "50735", "50736", "50737", "50738", "50739", "50740",
               "50741", "50742", "50743", "50744", "50745", "50746", "50747", "50749", "50750", "50751", "50752",
               "50754", "50755", "50756", "50757", "50758", "50759", "50760", "50761", "50762", "50763", "50765",
               "50766", "50767", "50771", "50772", "50773", "50774", "50776", "50777", "50778", "50779", "50780",
               "50781", "50782", "50783", "50784", "50785", "50786", "50787", "50789", "50790", "50793", "50794",
               "50795", "50796", "50797", "50798", "50799", "50801", "50802", "50803", "50804", "50805", "50806",
               "50807", "50808", "50809", "50811", "50812", "50813", "50814", "50815", "50816", "50817", "50818",
               "50819", "50820", "50821", "50824", "50825", "50826", "50827", "50828", "50829", "50830", "50831",
               "50835", "50836", "50837", "50839", "50842", "50845", "50846", "50847", "50848", "50849", "50850",
               "50851", "50852", "50853", "50855", "50858", "50859", "50860", "50861", "50862", "50863", "50864",
               "50865", "50866", "50867", "50868", "50869", "50870", "50871", "50886", "50887", "50888", "50889",
               "50890", "50891", "50892", "50895", "50896", "50897", "50898", "50899", "50900", "50901", "50902",
               "50903", "50904", "50905", "50906", "50907", "50908", "50909", "50910", "50912", "50913", "50914",
               "50915", "50916", "50917", "50918", "50919", "50920", "50921", "50922", "50923", "50924", "50925",
               "50926", "50927", "50928", "50929", "50930", "50931", "50932", "50933", "50934", "50936", "50938",
               "50940", "50941", "50942", "50944", "50945", "50946", "50948", "50949", "50950", "50951", "50952",
               "50953", "50955", "50956", "50957", "50959", "50960", "50961", "50962", "50963", "50964", "50966",
               "50968", "50969", "50970", "50971", "50972", "50974", "50976", "50977", "50978", "50979", "50980",
               "50981", "50982", "50983", "50985", "50986", "50987", "50988", "50989", "50993", "50994", "50996",
               "50997", "50998", "51003", "51004", "51006", "51007", "51008", "51009", "51010", "51011", "51012",
               "51013", "51014", "51015", "51016", "51018", "51019", "51020", "51021", "51022", "51023", "51024",
               "51025", "51026", "51027", "51028", "51029", "51030", "51031", "51032", "51033", "51034", "51035",
               "51036", "51037", "51038", "51039", "51040", "51041", "51042", "51043", "51044", "51045", "51046",
               "51047", "51048", "51050", "51051", "51052", "51054", "51055", "51056", "51057", "51058", "51059",
               "51061", "51062", "51063", "51064", "51065", "51066", "51067", "51068", "51069", "51070", "51071",
               "51072", "51073", "51074", "51075", "51076", "51077", "51078", "51079", "51080", "51081", "51083",
               "51084", "51085", "51086", "51087", "51090", "51091", "51092", "51093", "51094", "51096", "51097",
               "51098", "51099", "51100", "51101", "51102", "51103", "51104", "51105", "51106", "51107", "51108",
               "51109", "51110", "51111", "51112", "51113", "51114", "51115", "51116", "51117", "51118", "51119",
               "51120", "51121", "51122", "51123", "51124", "51125", "51126", "51127", "51128", "51129", "51130",
               "51131", "51132", "51133", "51134", "51135", "51137", "51138", "51139", "51140", "51141", "51142",
               "51143", "51144", "51145", "51146", "51147", "51148", "51149", "51152", "51155", "51156", "51158",
               "51159", "51160", "51161", "51162", "51163", "51164", "51165", "51166", "51167", "51169", "51170",
               "51171", "51172", "51173", "51174", "51175", "51176", "51177", "51178", "51179", "51180", "51181",
               "51182", "51183", "51184", "51186", "51187", "51188", "51192", "51193", "51194", "51195", "51198",
               "51200", "51201", "51202", "51203", "51204", "51205", "51206", "51207", "51208", "51209", "51210",
               "51211", "51222", "51223", "51224", "51226", "51227", "51228", "51230", "51231", "51232", "51233",
               "51234", "51235", "51236", "51237", "51238", "51239", "51240", "51241", "51242", "51243", "51244",
               "51245", "51246", "51247", "51248", "51249", "51253", "51254", "51255", "51256", "51257", "51258",
               "51259", "51260", "51262", "51264", "51265", "51268", "51269", "51270", "51271", "51272", "51273",
               "51274", "51275", "51276", "51277", "51278", "51279", "51280", "51281", "51282", "51283", "51284",
               "51285", "51286", "51287", "51288", "51289", "51290", "51291", "51292", "51293", "51294", "51296",
               "51297", "51299", "51301", "51303", "51304", "51305", "51306", "51307", "51308", "51309", "51310",
               "51311", "51312", "51313", "51314", "51315", "51316", "51317", "51318", "51321", "51322", "51323",
               "51324", "51325", "51327", "51331", "51332", "51333", "51334", "51335", "51336", "51337", "51338",
               "51339", "51340", "51341", "51342", "51343", "51344", "51345", "51347", "51348", "51349", "51350",
               "51351", "51352", "51353", "51354", "51355", "51356", "51358", "51359", "51360", "51362", "51363",
               "51364", "51365", "51366", "51367", "51368", "51369", "51370", "51371", "51372", "51373", "51374",
               "51375", "51376", "51377", "51378", "51379", "51380", "51381", "51382", "51384", "51388", "51389",
               "51390", "51391", "51392", "51393", "51394", "51395", "51396", "51397", "51398", "51399", "51400",
               "51401", "51402", "51403", "51404", "51405", "51406", "51407", "51408", "51410", "51411", "51412",
               "51413", "51414", "51415", "51416", "51417", "51418", "51419", "51420", "51421", "51422", "51423",
               "51424", "51427", "51428", "51429", "51430", "51431", "51432", "51433", "51434", "51435", "51436",
               "51437", "51438", "51439", "51440", "51441", "51442", "51443", "51444", "51445", "51447", "51448",
               "51449", "51450", "51451", "51452", "51453", "51454", "51455", "51456", "51457", "51458", "51459",
               "51460", "51461", "51462", "51463", "51464", "51465", "51466", "51467", "51468", "51469", "51470",
               "51471", "51472", "51473", "51474", "51475", "51476", "51478", "51479", "51480", "51481", "51483",
               "51485", "51486", "51487", "51488", "51489", "51490", "51491", "51492", "51494", "51495", "51497",
               "51500", "51502", "51503", "51504", "51505", "51506", "51508", "51509", "51510", "51512", "51513",
               "51514", "51515", "51516", "51517", "51518", "51519", "51520", "51521", "51522", "51523", "51524",
               "51525", "51526", "51527", "51528", "51529", "51530", "51531", "51532", "51533", "51534", "51535",
               "51536", "51537", "51538", "51539", "51540", "51541", "51546", "51547", "51548", "51550", "51565",
               "51566", "51567", "51568", "51569", "51570", "51571", "51572", "51573", "51574", "51575", "51576",
               "51577", "51578", "51579", "51580", "51581", "51582", "51583", "51584", "51591", "51592", "51593",
               "51594", "51595", "51596", "51597", "51598", "51599", "51600", "51601", "51603", "51604", "51605",
               "51606", "51607", "51608", "51609", "51610", "51611", "51612", "51614", "51616", "51617", "51618",
               "51619", "51620", "51621", "51622", "51623", "51624", "51625", "51627", "51628", "51629", "51630",
               "51631", "51632", "51633", "51634", "51637", "51638", "51639", "51642", "51644", "51645", "51646",
               "51647", "51648", "51649", "51650", "51651", "51652", "51653", "51654", "51655", "51656", "51657",
               "51658", "51659", "51661", "51662", "51663", "51672", "51673", "51674", "51675", "51676", "51677",
               "51678", "51681", "51682", "51683", "51687", "51688", "51689", "51690", "51691", "51692", "51695",
               "51696", "51697", "51698", "51700", "51701", "51702", "51703", "51704", "51705", "51706", "51707",
               "51708", "51709", "51710", "51711", "51713", "51714", "51715", "51716", "51717", "51719", "51720",
               "51721", "51723", "51724", "51725", "51726", "51728", "51729", "51730", "51731", "51735", "51736",
               "51737", "51738", "51739", "51741", "51742", "51743", "51744", "51745", "51746", "51747", "51748",
               "51749", "51753", "51754", "51755", "51756", "51758", "51763", "51766", "51768", "51769", "51770",
               "51771", "51772", "51773", "51774", "51775", "51776", "51777", "51779", "51780", "51781", "51782",
               "51783", "51784", "51785", "51786", "51787", "51788", "51790", "51791", "51792", "51793", "51794",
               "51795", "51796", "51797", "51798", "51800", "51802", "51808", "51809", "51810", "51811", "51812",
               "51813", "51814", "51815", "51819", "51820", "51821", "51822", "51823", "51825", "51827", "51828",
               "51831", "51832", "51833", "51834", "51837", "51840", "51841", "51842", "51843", "51844", "51845",
               "51846", "51847", "51848", "51849", "51850", "51851", "51852", "51853", "51855", "51857", "51858",
               "51859", "51860", "51869", "51870", "51871", "51877", "51878", "51881", "51882", "51885", "51886",
               "51887", "51888", "51890", "51891", "51892", "51893", "51894", "51896", "51899", "51900", "51901",
               "51902", "51905", "51906", "51907", "51908", "51909", "51910", "51911", "51912", "51913", "51914",
               "51915", "51916", "51917", "51919", "51920", "51921", "51922", "51923", "51924", "51925", "51926",
               "51927", "51929", "51930", "51931", "51932", "51933", "51934", "51935", "51937", "51938", "51939",
               "51940", "51941", "51942", "51943", "51950", "51951", "51952", "51953", "51955", "51956", "51957",
               "51959", "51960", "51962", "51963", "51964", "51973", "51974", "51975", "51978", "51979", "51981",
               "51982", "51983", "51984", "51991", "51992", "51993", "51994", "51996", "51997", "51999", "52001",
               "52002", "52003", "52004", "52006", "52007", "52008", "52009", "52010", "52012", "52013", "52014",
               "52016", "52017", "52018", "52021", "52023", "52024", "52028", "52030", "52031", "52032", "52033",
               "52034", "52035", "52040", "52041", "52042", "52043", "52045", "52049", "52051", "52052", "52053",
               "52059", "52060", "52061", "52062", "52064", "52065", "52071", "52075", "52077", "52078", "52079",
               "52080", "52081", "52082", "52087", "52090", "52091", "52092", "52098", "52099", "52100", "52101",
               "52102", "52104", "52105", "52108", "52109", "52115", "52121", "52122", "52123", "52124", "52125",
               "52126", "52127", "52130", "52137", "52138", "52139", "52140", "52141", "52142", "52143", "52145",
               "52146", "52147", "52150", "52152", "52153", "52154", "52155", "52158", "52160", "52162", "52169",
               "52170", "52171", "52174", "52175", "52179", "52187", "52188", "52196", "52197", "52207", "52209",
               "52210", "52211", "52212", "52213", "52214", "52217", "52235", "52237", "52244", "52245", "52251",
               "52253", "52254", "52255", "52256", "52257", "52258", "52266", "52273", "52284", "52289", "52293",
               "52294", "52295", "52296", "52297", "52298", "52300", "52301", "52302", "52305", "52306", "52307",
               "52310", "52339", "52355"]
    # f = open('sc_stat.csv', 'w')
    # for ent in entries:
    #     str_file = f'/reboxitory/2024/03/BMRB/macromolecules/bmr{ent}/bmr{ent}_3.str'
    #     if os.path.isfile(str_file):
    #         x = check_chemical_shifts(str_file, ent)
    #         out = f'{ent}'
    #         for l in x['completeness']:
    #             out = f'{ent},{l}'
    #             for e in x['completeness'][l]:
    #                 out = f'{out},{e}'
    #                 for c1 in c:
    #                     for sc1 in sc:
    #                         if x["completeness"][l][e][c1][sc1][0] > 0:
    #                             out = f'{out},{x["completeness"][l][e][c1][sc1][0]},{x["completeness"][l][e][c1][sc1][1]},{round(100.0 * (float(x["completeness"][l][e][c1][sc1][1]) / float(x["completeness"][l][e][c1][sc1][0])), 2)}'
    #                         else:
    #                             out = f'{out},{x["completeness"][l][e][c1][sc1][0]},{x["completeness"][l][e][c1][sc1][1]},0.00'
    #
    #                 f.write(f'{out}\n')
    #     else:
    #         print(f'file not found {str_file}')
