from mmcif.io.PdbxReader import PdbxReader
import plotly.express as px

def read_validation_cif(fname):

    cif_data = []
    ifh = open(fname, 'r')
    pRd = PdbxReader(ifh)
    pRd.read(cif_data)
    ifh.close()
    c0 = cif_data[0]
    status = c0.getObj('pdbx_vrpt_software')
    prog_stat={}
    for row in status.getRowList():
        prog_stat[row[1]]=row[3]
    try:
        if prog_stat['restraintsanalysis'] == 'Y':
            data = {}
            summary=c0.getObj('pdbx_vrpt_summary')
            entry_id = summary.getValue('entry_id')
            deposition_date=summary.getValue('PDB_deposition_date')
            rest_summary = c0.getObj('pdbx_vrpt_restraint_summary')
            dist_violation = c0.getObj('pdbx_vrpt_distance_violation_summary')
            residual_dist_violation = c0.getObj('pdbx_vrpt_residual_distance_violations')
            ang_violation=c0.getObj('pdbx_vrpt_dihedralangle_violations_summary')
            residual_ang_violation = c0.getObj('pdbx_vrpt_residual_angle_violations')
            data['id'] = entry_id
            data['date'] = deposition_date
            data['rest_summary']={}
            for row in rest_summary.getRowList():
                data['rest_summary'][row[1]]= float(row[2])
            data['dist_violation']={}
            col_names = dist_violation.getAttributeList()
            for row in dist_violation.getRowList():
                if row[1] not in data['dist_violation']:
                    data['dist_violation'][row[1]]={}
                if row[2] not in data['dist_violation'][row[1]]:
                    data['dist_violation'][row[1]][row[2]]={}
                for tag in col_names[3:]:
                    data['dist_violation'][row[1]][row[2]][tag]=float(row[col_names.index(tag)])
            data['residual_dist_violation']={}
            data['residual_dist_violation']['max_violation']={}
            data['residual_dist_violation']['violation_per_model']={}
            try:
                for row in residual_dist_violation.getRowList():
                    try:
                        data['residual_dist_violation']['max_violation'][row[1]]=float(row[2])
                    except ValueError:
                        data['residual_dist_violation']['max_violation'][row[1]]=0.0
                    try:
                        data['residual_dist_violation']['violation_per_model'][row[1]]=float(row[3])
                    except ValueError:
                        data['residual_dist_violation']['violation_per_model'][row[1]]=0.0
            except AttributeError:
                pass
            data['residual_ang_violation'] = {}
            data['residual_ang_violation']['max_violation'] = {}
            data['residual_ang_violation']['violation_per_model'] = {}
            try:
                for row in residual_ang_violation.getRowList():
                    try:
                        data['residual_ang_violation']['max_violation'][row[1]] = float(row[2])
                    except ValueError:
                        data['residual_ang_violation']['max_violation'][row[1]] = 0.0
                    try:
                        data['residual_ang_violation']['violation_per_model'][row[1]] = float(row[3])
                    except ValueError:
                        data['residual_ang_violation']['violation_per_model'][row[1]] = 0.0
            except AttributeError:
                pass
            data['ang_violation'] = {}
            try:
                col_names = ang_violation.getAttributeList()
                for row in ang_violation.getRowList():
                    if row[1] not in data['ang_violation']:
                        data['ang_violation'][row[1]] = {}
                    for tag in col_names[2:]:
                        data['ang_violation'][row[1]][tag] = float(row[col_names.index(tag)])
            except AttributeError:
                data['ang_violation'] = None


        else:
            print (f'{fname} has no restraints')
            data = None
    except KeyError:
        data = None
    return data

def collect_statistics(flist):
    f=open(flist,'r').read().split("\n")[:-1]
    total_rest = []
    date = []
    id = []
    rest_per_res=[]
    consistently_violated_count=[]
    consistently_violated_percent=[]
    id_dist=[]
    max_dist_viol=[]
    viol_per_model_dist=[]
    id_ang = []
    max_ang_viol = []
    viol_per_model_ang = []
    dist_tag=[]
    ang_tag=[]
    for fname in f:
        #print (f'/Users/kumaranbaskaran/test/testing2/{fname}')
        rest_data=read_validation_cif(f'/Users/kumaranbaskaran/test/testing2/{fname}')
        if rest_data != None:
            id.append(rest_data['id'])
            date.append(rest_data['date'])
            total_rest.append(rest_data['rest_summary']['Total distance restraints'])
            rest_per_res.append(rest_data['rest_summary']['Number of restraints per residue'])
            consistently_violated_count.append(rest_data['dist_violation']['Total']['all']['consistently_violated_count'])
            consistently_violated_percent.append(rest_data['dist_violation']['Total']['all']['consistently_violated_percent_total'])
            for k in rest_data['residual_dist_violation']['max_violation']:
                id_dist.append(rest_data['id'])
                max_dist_viol.append(rest_data['residual_dist_violation']['max_violation'][k])
                viol_per_model_dist.append(rest_data['residual_dist_violation']['violation_per_model'][k])
                dist_tag.append(k)
            for k in rest_data['residual_ang_violation']['max_violation']:
                id_ang.append(rest_data['id'])
                max_ang_viol.append(rest_data['residual_ang_violation']['max_violation'][k])
                viol_per_model_ang.append(rest_data['residual_ang_violation']['violation_per_model'][k])
                ang_tag.append(k)
        else:
            print (f'{fname} Failed')
    f=open('/Users/kumaranbaskaran/rest_analysis/vali_out.txt','w')
    for i in range(len(id)):
        f.write(f'{id[i]},{rest_per_res[i]},{consistently_violated_percent[i]},{consistently_violated_count[i]},{max_dist_viol[i]},{max_ang_viol[i]}\n')
    f.close()

    fig = px.bar(x=id,y=rest_per_res,labels={'x':'PDB ID','y':'No. of restraints per residue'})
    fig.write_html('/Users/kumaranbaskaran/rest_analysis/rest_per_residue.html')
    fig.show()
    fig2= px.histogram(rest_per_res,labels={'value': 'No. of restraints per residue'})
    fig2.write_html('/Users/kumaranbaskaran/rest_analysis/rest_per_residue_hist.html')
    fig2.show()
    fig3= px.histogram(x=id,y=consistently_violated_percent,labels={'x':'PDB ID','y':'Consistently Violated (%)'})
    fig3.write_html('/Users/kumaranbaskaran/rest_analysis/violations.html')
    fig3.show()
    fig4=px.scatter(x=rest_per_res,y=consistently_violated_percent,hover_name=id,labels={'x':'No. of restraints per residue','y':'Consistently violated'})
    fig4.write_html('/Users/kumaranbaskaran/rest_analysis/rest_vs_violated.html')
    fig4.show()
    fig5=px.scatter(x=max_dist_viol,y=viol_per_model_dist,hover_name=id_dist,color=dist_tag,
                    labels={'x':'Max distance violation','y':'Violations per model'})
    fig5.write_html('/Users/kumaranbaskaran/rest_analysis/residual_dist_violation.html')
    fig5 = px.scatter(x=max_ang_viol, y=viol_per_model_ang, hover_name=id_ang, color=ang_tag,
                      labels={'x': 'Max angle violation', 'y': 'Violations per model'})
    fig5.write_html('/Users/kumaranbaskaran/rest_analysis/residual_ang_violation.html')

if __name__ == "__main__":
    #collect_statistics('/Users/kumaranbaskaran/test/testing2/cif_list.txt')
    read_validation_cif('/Users/kumaranbaskaran/test/testing2/2lg5_val-data.cif')
