import pynmrstar
import os
import json
import logging
from urllib.request import urlopen, Request
import numpy
from typing import Union, List, Optional
import Levenshtein as pylev
import csv
from multiprocessing import Pool
def get_entry_list():
    url = Request('http://api.bmrb.io/v2/list_entries?database=macromolecules')
    url.add_header('Application', 'PyBMRB')
    r = urlopen(url)
    entry_list = json.loads(r.read())
    sf={}
    lp={}
    f=open('bmrb_seq.csv','w')
    for bid in entry_list:
        print (bid)
        str_file = f'/reboxitory/2024/03/BMRB/macromolecules/bmr{bid}/bmr{bid}_3.str'
        if os.path.isfile(str_file):
            ent = pynmrstar.Entry.from_file(str_file)
            sq=ent.get_tag('_Entity.Polymer_seq_one_letter_code')
            for i in range(len(sq)):
                sq1=sq[i].replace("\\n","")
                f.write(f'{bid},{i+1},{sq1}\n')
    f.close()

            # for saveframe in ent:
            #     if saveframe.category not in sf:
            #         sf[saveframe.category]={}
            #     if bid not in sf[saveframe.category]:
            #         sf[saveframe.category][bid]=0
            #     sf[saveframe.category][bid]+=1
            #     for loop in saveframe:
            #         if loop.category not in lp:
            #             lp[loop.category]={}
            #         if bid not in lp[loop.category]:
            #             lp[loop.category][bid]=[]
            #         lp[loop.category][bid].append(len(loop.data))
    # f=open('sfinfo.txt','w')
    # for s in sf:
    #     print (s,sf[s])
    #     f.write(f'{s},{sf[s]}\n')
    # f.close()
    # f=open('lpinfo.txt','w')
    # for l in lp:
    #     print (l,len(lp[l]))
    #     f.write(f'{l},{lp[l]}\n')
    # f.close()








def read_disprot_annotations(disprot_file):
    f=open(disprot_file,'r')
    d=json.load(f)
    disp={}
    for i in d['data']:
        if i['disprot_id'] not in disp:
            disp[i['disprot_id']]={}
        #print (i['sequence'])
        #print(i['disprot_id'])
        for j in i['regions']:
            sq=i['sequence'][j['start']-1:j['end']]
            seg=(j['start'],j['end'])
            if seg not in disp[i['disprot_id']]:
                disp[i['disprot_id']][seg]=sq
    return disp

def sequence_match(sq1,sq2):
    n1 = len(sq1)
    n2 = len(sq2)
    if n1>n2:
        long_seq = sq1
        short_seq = sq2

    else:
        long_seq = sq2
        short_seq = sq1
    n=len(short_seq)
    r=0.0
    for i in range(len(long_seq)):
        s1=long_seq[i:i+n]
        s2=short_seq
        if len(s1)>=n:
            if pylev.ratio(s1,s2)>r:
                r=pylev.ratio(s1,s2)
    return r

def read_star(star_file,bid):
    ent = pynmrstar.Entry.from_file(star_file)
    #het_noe = ent.get_loops_by_category('_Heteronucl_NOE')
    het_noe_sf = ent.get_saveframes_by_category('heteronucl_NOEs')
    T1_sf = ent.get_saveframes_by_category('heteronucl_T1_relaxation')
    T2_sf = ent.get_saveframes_by_category('heteronucl_T2_relaxation')
    order_pm_sf = ent.get_saveframes_by_category('order_parameters')
    if len(het_noe_sf):
        het_noe = "yes"
    else:
        het_noe = "no"
    if len(T1_sf):
        T1 = "yes"
    else:
        T1 = "no"
    if len(T2_sf):
        T2 = "yes"
    else:
        T2 = "no"
    if len(order_pm_sf):
        op = "yes"
    else:
        op = "no"
    if "yes" in [het_noe,T1,T2,op]:
        print (bid,het_noe,T1,T2,op)
    return [het_noe,T1,T2,op]

def mathc_bmrb(dat):
    sq=dat[0]
    bmrb=dat[1]
    did=dat[2]
    r=0.75
    match_id = None
    for bid in bmrb:
        bq = bmrb[bid]
        if len(sq) > 8 and len(bq) > 8:
            if  sequence_match(sq, bq) > r:
                r = sequence_match(sq, bq)
                match_id = bid
    return did,match_id,r

def mathc_bmrb2(dat):
    sq=dat[0]
    bmrb=dat[1]
    did=dat[2]
    r=0.75
    match_id = None
    for bid in bmrb:
        bq = bmrb[bid]
        if len(sq) > 8 and len(bq) > 8:
            if  sequence_match(sq, bq) > r:
                r = sequence_match(sq, bq)
                match_id = bid
    return did,match_id,r



def find_seq_match(bmrb):
    disp = read_disprot_annotations('DisProt release_2023_12 with_ambiguous_evidences.json')
    with open(bmrb, mode='r') as file:
        csvFile = csv.reader(file)
        bmrb={}
        for lines in csvFile:
            if len(lines[2])>2:
                k=f'{lines[0]}-{lines[1]}'
                if k not in bmrb:
                    bmrb[k]=lines[2]
                else:
                    print (lines)
    d=[]
    for dpid in disp:
        for k in disp[dpid]:
            sq = disp[dpid][k]
            d.append((sq,bmrb,(dpid,k)))
    with Pool() as pool:
        result = pool.map(mathc_bmrb,d[0:8])
    f=open('bmrb_disprot_mapping.csv','w')
    for res in result:
        f.write(f'{res[0][0]},{res[0][1][0]}-{res[0][1][1]},{res[1]},{res[2]}\n')

    f.close()

def find_seq_match2(bmrb):
    disp = read_disprot_annotations('DisProt release_2023_12 with_ambiguous_evidences.json')
    with open(bmrb, mode='r') as file:
        csvFile = csv.reader(file)
        bmrb={}
        for lines in csvFile:
            if len(lines[2])>2:
                k=f'{lines[0]}-{lines[1]}'
                if k not in bmrb:
                    bmrb[k]=lines[2]
                else:
                    print (lines)
    disp_data = {}
    for k in disp:
        for seg in disp[k]:
            disp_data[f'{k}-{seg[0]}-{seg[1]}']= disp[k][seg]
    d=[]
    for bmid in bmrb:
        d.append((bmrb[bmid],disp_data,bmid))
    with Pool() as pool:
        result = pool.map(mathc_bmrb2,d[0:8])
    f=open('bmrb_disprot_mapping2.csv','w')
    for res in result:
        f.write(f'{res[0]},{res[1]},{res[2]}\n')

    f.close()

def read_disprot_mapping(csvfile,csvfile2,coff=0.8):
    f=open(f'disprot_bmrb_mapping_{coff}.csv','w')
    with open(csvfile, mode='r') as file:
        csvFile = csv.reader(file)
        bmrb={}
        c=0
        p=0
        did={}
        bid=[]
        ddd=[]
        disp_id1=[]
        disp_id0 = []
        for lines in csvFile:
            c+=1
            disp_id0.append(lines[0])
            if lines[0] not in did:
                did[lines[0]]=[]
            did[lines[0]].append(lines[1:])
            if float(lines[3])>=coff:
                ddd.append(f'{lines[0]}-{lines[1]},{lines[2]},{round(float(lines[3]),2)}')
                disp_id1.append(lines[0])
                p+=1
                bid.append(lines[2].split("-")[0])
        p2=0
        for k in did:
            flag='no'
            for k1 in did[k]:
                if float(k1[2])>=coff:
                    flag='yes'
            if flag == 'yes':
                p2+=1
    bid2=[]
    disp_id2=[]
    with open(csvfile2, mode='r') as file:
        csvFile = csv.reader(file)
        for lines in csvFile:
            if float(lines[2]) >=coff:
                disp_id2.append(lines[1].split("-")[0])
                ddd.append(f'{lines[1]},{lines[0]},{round(float(lines[2]),2)}')
                bid2.append(lines[0].split("-")[0])
    print (c)
    print (len(disp_id0),len(disp_id1),len(disp_id2))
    print (len(set(disp_id0)))
    print (len(set(disp_id1)),len(set(disp_id2)))
    print (len(set(disp_id1+disp_id2)))
    b=bid+bid2
    dd=list(set(ddd))
    print (len(dd))
    for l in dd:
        f.write(f'{l}\n')
    f.close()
    print ('er',len(set(b)))
    return set(b)

def all_stat(lpfile):
    catch_list = ['_Heteronucl_T1_experiment','_Heteronucl_T2_experiment','_Heteronucl_NOE_experiment']
    disprot_bmrb_ids=read_disprot_mapping('bmrb_disprot_mapping.csv','bmrb_disprot_mapping2.csv',0.8)
    f=open(lpfile,'r')
    lp_info = json.load(f)
    data_sets={}
    for k in lp_info:
        if k in catch_list:
            data_sets[k]=list(lp_info[k].keys())
    s1=set(data_sets['_Heteronucl_T1_experiment'])
    s2=set(data_sets['_Heteronucl_T2_experiment'])
    s3 =set(data_sets['_Heteronucl_NOE_experiment'])
    s4= set(disprot_bmrb_ids)
    s1s2=s1.intersection(s2)
    s1s2s3=s1s2.intersection(s3)
    s1s2s3s4=s1s2s3.intersection(s4)
    s1s4 = s1.intersection(s4)
    s4s1 = s4.intersection(s1)
    print (s1s2s3)
    print (len(s1),len(s2),len(s3),len(s1s2s3),len(s1s4))

ALL_3 = ['17881', '17947', '5080', '17047', '18758', '19361', '27722', '27655', '27916', '5841', '25013', '25113', '51103', '28132', '4365', '15065', '19160', '15437', '50233', '50745', '5519', '15795', '7432', '4364', '50819', '18362', '4762', '7056', '51860', '26791', '25636', '51423', '18388', '16234', '19388', '18773', '27448', '18389', '16737', '17266', '16218', '50734', '5153', '25035', '15445', '16426', '50438', '17012', '5991', '19188', '51416', '51415', '26742', '16876', '19164', '16657', '51087', '17246', '17308', '18361', '18424', '17783', '18380', '34544', '17226', '51413', '5272', '25525', '51933', '7415', '18617', '4390', '27887', '51414', '51931', '19153', '51230', '16392', '7088', '25390', '27888', '5505', '5548', '7190', '51144', '4970', '17282', '51119', '19127', '34496', '51930', '5330', '19189', '30523', '50495', '51236', '50283', '27321', '51235', '27011', '18230', '18191', '51224', '25389', '15254', '15655', '18257', '26823', '25523', '4689', '25025', '15541', '5518', '15793', '27929', '18864', '5687', '51932', '7414', '18545', '17983', '15064', '15230', '5550', '6470', '15364', '50238', '19072', '36171', '51065', '51934', '27646', '16482', '16659', '17982', '50332', '51066', '34262', '18903', '16845', '25014', '50120', '25034', '5569', '15703', '4366', '15521', '16307', '18192', '19335', '51591', '51935', '19356', '7288', '5707', '6332', '52077', '15097', '17981', '25015', '6495', '27721', '4245', '50410', '11080', '6881', '5520', '15067', '34497', '51223', '50285', '25121', '18087', '6092', '26779', '18971', '6577', '15014', '16907', '18772', '25519', '5549', '27890', '19284', '50212', '6060', '50001', '5154', '50119', '51530', '27194', '16485', '17865', '15728', '26741', '51174', '30834', '17010', '17306', '18461', '27656', '34545', '50333', '17013', '6243', '16306', '50494', '17046', '51117', '18260', '17610', '27179', '50284', '50482', '17857', '16483', '51306', '26724', '5079', '5995', '5839', '4763', '5996', '5746', '19993', '5331', '6474', '4870', '18359', '6758', '51234', '18477', '18306', '17018', '50243', '25852', '51305', '6182', '26506', '52082', '15975', '27063', '5521', '51417', '18231', '16360', '27594', '50020', '5762', '27447', '51418', '51529', '16033', '18360', '17041', '15562', '26507', '17069', '15066', '27245', '15255', '50234', '50115', '15578', '4697', '6880', '18092', '25255', '51964', '15923', '27915', '6494', '6838', '18351', '26511', '18422', '16034', '5720', '28058', '15451', '51307', '17080', '18423', '26723', '26513', '16480']


def get_t1_t2_values(bid):
    str_file = f'/reboxitory/2024/03/BMRB/macromolecules/bmr{bid}/bmr{bid}_3.str'
    ent = pynmrstar.Entry.from_file(str_file)
    #t1_loops = ent.get_loops_by_category('_T1')
    #t2_loops = ent.get_loops_by_category('_T2')
    t1={}
    t2={}
    noe={}
    d={}
    units=[]
    for saveframe in ent:
        if saveframe.category == 'heteronucl_T1_relaxation':
            tag_dict = saveframe.tag_dict
            field_strength = tag_dict['spectrometer_frequency_1h']
            #print (bid,tag_dict)
            unit = tag_dict['t1_val_units']
            if unit not in units:
                units.append(unit)
            if field_strength not in d:
                d[field_strength]={}
            for t1_loop in saveframe:
                if t1_loop.category == '_T1':
                    t1[field_strength]={}
                    t1_columns = t1_loop.get_tag_names()
                    entity_assembly_idx = t1_columns.index('_T1.Entity_assembly_ID')
                    entity_idx = t1_columns.index('_T1.Entity_ID')
                    comp_index_idx = t1_columns.index('_T1.Comp_index_ID')
                    comp_idx = t1_columns.index('_T1.Comp_ID')
                    atom_idx = t1_columns.index('_T1.Atom_ID')
                    val_idx = t1_columns.index('_T1.Val')
                    val_err_idx = t1_columns.index('_T1.Val_err')
                    list_idx = t1_columns.index('_T1.Heteronucl_T1_list_ID')
                    for data in t1_loop.data:
                        if (data[entity_assembly_idx],data[entity_idx],data[comp_index_idx],data[comp_idx]) not in d[field_strength]:
                            d[field_strength][(data[entity_assembly_idx],data[entity_idx],data[comp_index_idx],data[comp_idx])]=[None,None,None]
                        try:
                            # if unit == 's-1' and float(data[val_idx])!=0:
                            #     #d[field_strength][(data[entity_assembly_idx],data[entity_idx],data[comp_index_idx],data[comp_idx])][0]=float(data[val_idx])
                            #     pass
                            # elif unit == 'ms':
                            #     d[field_strength][(
                            #     data[entity_assembly_idx], data[entity_idx], data[comp_index_idx], data[comp_idx])][
                            #         0] = float(data[val_idx])/1000.00
                            # elif unit == 's' and float(data[val_idx])>10:
                            #     d[field_strength][(
                            #         data[entity_assembly_idx], data[entity_idx], data[comp_index_idx], data[comp_idx])][
                            #         0] = float(data[val_idx])/1000.00
                            # else:
                            d[field_strength][(
                            data[entity_assembly_idx], data[entity_idx], data[comp_index_idx], data[comp_idx])][
                                0] = (float(data[val_idx]),unit)
                        except ValueError:
                            pass
                        t1[field_strength][(data[entity_assembly_idx],data[entity_idx],data[comp_index_idx],data[comp_idx],data[atom_idx])] = (data[val_idx],data[val_err_idx])
        if saveframe.category == 'heteronucl_T2_relaxation':
            tag_dict = saveframe.tag_dict
            field_strength = tag_dict['spectrometer_frequency_1h']
            unit = tag_dict['t2_val_units']
            if unit not in units:
                units.append(unit)
            if field_strength not in d:
                d[field_strength]={}
            for t2_loop in saveframe:
                if t2_loop.category == '_T2':
                    t2[field_strength]={}
                    t2_columns = t2_loop.get_tag_names()
                    entity_assembly_idx = t2_columns.index('_T2.Entity_assembly_ID')
                    entity_idx = t2_columns.index('_T2.Entity_ID')
                    comp_index_idx = t2_columns.index('_T2.Comp_index_ID')
                    comp_idx = t2_columns.index('_T2.Comp_ID')
                    atom_idx = t2_columns.index('_T2.Atom_ID')
                    val_idx = t2_columns.index('_T2.T2_val')
                    val_err_idx = t2_columns.index('_T2.T2_val_err')
                    list_idx = t2_columns.index('_T2.Heteronucl_T2_list_ID')
                    for data in t2_loop.data:
                        if (data[entity_assembly_idx],data[entity_idx],data[comp_index_idx],data[comp_idx]) not in d[field_strength]:
                            d[field_strength][(data[entity_assembly_idx],data[entity_idx],data[comp_index_idx],data[comp_idx])]=[None,None,None]
                        try:
                            # if unit == 's-1' and float(data[val_idx])!=0:
                            #     pass
                            #     # d[field_strength][
                            #     # (data[entity_assembly_idx], data[entity_idx], data[comp_index_idx], data[comp_idx])][
                            #     # 1] = float(data[val_idx])
                            # elif unit == 'ms':
                            #     d[field_strength][
                            #         (
                            #         data[entity_assembly_idx], data[entity_idx], data[comp_index_idx], data[comp_idx])][
                            #         1] = float(data[val_idx])/1000.00
                            # elif unit == 's' and float(data[val_idx]) > 10.0:
                            #     d[field_strength][
                            #         (
                            #             data[entity_assembly_idx], data[entity_idx], data[comp_index_idx],
                            #             data[comp_idx])][
                            #         1] = float(data[val_idx])/1000.00
                            # else:
                            d[field_strength][
                                (
                                data[entity_assembly_idx], data[entity_idx], data[comp_index_idx], data[comp_idx])][
                                1] = (float(data[val_idx]),unit)
                        except ValueError:
                            #print (data)
                            #print ('Errir',bid,data[entity_assembly_idx], data[entity_idx], data[comp_index_idx], data[comp_idx])
                            pass
                        t2[field_strength][(data[entity_assembly_idx],data[entity_idx],data[comp_index_idx],data[comp_idx],data[atom_idx])] = (data[val_idx],data[val_err_idx])
        if saveframe.category == 'heteronucl_NOEs':
            tag_dict = saveframe.tag_dict
            field_strength = tag_dict['spectrometer_frequency_1h']
            if field_strength not in d:
                d[field_strength]={}
            for noe_loop in saveframe:
                if noe_loop.category == '_Heteronucl_NOE':
                    noe[field_strength] = {}
                    noe_columns = noe_loop.get_tag_names()
                    entity_assembly_idx1 = noe_columns.index('_Heteronucl_NOE.Entity_assembly_ID_1')
                    entity_idx1 = noe_columns.index('_Heteronucl_NOE.Entity_ID_1')
                    comp_index_idx1 = noe_columns.index('_Heteronucl_NOE.Comp_index_ID_1')
                    comp_idx1 = noe_columns.index('_Heteronucl_NOE.Comp_ID_1')
                    atom_idx1 = noe_columns.index('_Heteronucl_NOE.Atom_ID_1')
                    entity_assembly_idx2 = noe_columns.index('_Heteronucl_NOE.Entity_assembly_ID_2')
                    entity_idx2 = noe_columns.index('_Heteronucl_NOE.Entity_ID_2')
                    comp_index_idx2 = noe_columns.index('_Heteronucl_NOE.Comp_index_ID_2')
                    comp_idx2 = noe_columns.index('_Heteronucl_NOE.Comp_ID_2')
                    atom_idx2 = noe_columns.index('_Heteronucl_NOE.Atom_ID_2')
                    val_idx = noe_columns.index('_Heteronucl_NOE.Val')
                    val_err_idx = noe_columns.index('_Heteronucl_NOE.Val_err')
                    list_idx = noe_columns.index('_Heteronucl_NOE.Heteronucl_NOE_list_ID')
                    for data in noe_loop.data:
                        if (data[entity_assembly_idx1], data[entity_idx1], data[comp_index_idx1], data[comp_idx1]) not in d[
                            field_strength]:
                            d[field_strength][
                                (data[entity_assembly_idx1], data[entity_idx1], data[comp_index_idx1], data[comp_idx1])] = [
                                None, None, None]
                        d[field_strength][
                            (data[entity_assembly_idx1], data[entity_idx1], data[comp_index_idx1], data[comp_idx1])][
                            2] = float(data[val_idx])
                        noe[field_strength][(data[entity_assembly_idx1],data[entity_idx1],data[comp_index_idx1],data[comp_idx1],data[atom_idx1]),
                        (data[entity_assembly_idx2],data[entity_idx2],data[comp_index_idx2],data[comp_idx2],data[atom_idx2])]=(data[val_idx],data[val_err_idx])
    return t1,t2,noe,d,units
if __name__ == "__main__":
    # 16737['ms-1']
    # 4390['Hz']
    # 5505['ns', 's-1']

    c=0
    u=[]
    for bid in ALL_3:
        if True:#bid not in ['16737','4390','5505']:
            t1,t2,noe,d,units=get_t1_t2_values(bid)
            for fs in d:
                for atm in d[fs]:
                    try:
                        print (f'{bid}-{fs}-{"-".join(list(atm))}:{d[fs][atm][0][1]}:{d[fs][atm][1][1]},{d[fs][atm][0][0]},{d[fs][atm][1][0]},{d[fs][atm][2]}')
                    except TypeError:
                        pass
    #print (set(u))


    #all_stat('lpinfo.json')
    #read_disprot_mapping('bmrb_disprot_mapping.csv','bmrb_disprot_mapping2.csv')
    #find_seq_match2('bmrb_seq.csv')
    #print (sequence_match('ac','abcdefghijklmnopqrstuvwxyz'))
    #read_disprot_annotations('DisProt release_2023_12 with_ambiguous_evidences.json')
    # c2=0
    # for ent in entries:
    #
    #     str_file = f'/reboxitory/2024/03/BMRB/macromolecules/bmr{ent}/bmr{ent}_3.str'
    #     if os.path.isfile(str_file):
    #         x=read_star(str_file,ent)
    #         if 'yes' in x:
    #             c1+=1
    #         if x[0] == "yes" and  x[1] == "yes" and  x[2] == "yes":
    #             c2+=1
    # print (c1,c2)
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
