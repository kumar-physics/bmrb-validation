import plotly.express as px
import csv

def plot_data(fname):
    tot=[]
    x=[]
    y=[]
    l=[]
    with open(fname, mode='r') as file:
        csvFile = csv.reader(file)
        for lines in csvFile:
            x.append(float(lines[5]))
            y.append(float(lines[8]))
            l.append(lines[0])
            tot.append(float(lines[5]))
            if float(lines[5])>100:
                print (lines)
    fig = px.histogram(tot)
    fig.show()
    # fig2 = px.scatter(x=x,y=y,hover_name=l)
    # fig2.show()


def plot_relx_data(fname):
    t1=[]
    t2=[]
    noe=[]
    tag=[]
    r=[]
    u=[]
    s=[]
    t1_modified = []
    t2_modified = []
    with open(fname, mode='r') as file:
        csvFile = csv.reader(file)
        for lines in csvFile:
            if lines[1] != 'None' and lines[2] != 'None' and lines[3] != 'None':
                if float(lines[1])>-10000:
                    tag.append(lines[0])
                    u1=lines[0].split(":")[1]
                    u2=lines[0].split(":")[2]
                    t1.append(float(lines[1]))
                    t2.append(float(lines[2]))
                    if u1 == 's-1':
                        if float(lines[1])!=0:
                            t1_modified.append(1.0/float(lines[1]))
                        else:
                            t1_modified.append(float(lines[1]))
                    elif u1 == 'ms':
                        t1_modified.append(float(lines[1])/1000.00)
                    elif u1 == 'Hz':
                        t1_modified.append(1.0 / float(lines[1]))
                    elif u1 == 'ms-1':
                        t1_modified.append(float(lines[1]) / 1000.00)
                    else:
                        t1_modified.append(float(lines[1]))
                    if u2 == 's-1':
                        if float(lines[2]) !=0:
                            t2_modified.append(1.0/float(lines[2]))
                        else:
                            t2_modified.append(float(lines[2]))
                    elif u2 == 'ms':
                        t2_modified.append(float(lines[2])/1000.00)
                    elif u1 == 'Hz':
                        t2_modified.append(1.0 / float(lines[2]))
                    elif u1 == 'ms-1':
                        t2_modified.append(float(lines[2]) / 1000.00)
                    else:
                        t2_modified.append(float(lines[2]))

                    noe.append(float(lines[3]))
                    s.append(u1)
                    print (u1,u2)
                    if u1==u2:
                        u.append('Same')
                    else:
                        u.append('Different')
                    if float(lines[1]) > 0:
                        if float(lines[2])/float(lines[1]) < 2.0:
                            r.append(float(lines[2])/float(lines[1]))
    print (len(t1))
    fig = px.scatter(x=t1,y=t2,hover_name=tag,color=u,symbol=s,labels={'x':'T1[s,s-1,ms,ms-1,Hz]','y':'T2[s,s-1,ms,ms-1,Hz]'})
    fig.show()
    fig.write_html('dirty_data.html')
    fig2 = px.scatter(x=t1_modified, y=t2_modified, hover_name=tag, color=u, symbol=s,labels={'x':'T1[s]','y':'T2[s]'})
    fig2.show()
    fig2.write_html('dirty_data2.html')
    fig3 = px.scatter_3d(x=t1_modified, y=t2_modified, z=noe,hover_name=tag, color=u, symbol=s,
                      labels={'x': 'T1[s]', 'y': 'T2[s]','z':'NOE[?]'})

    fig3.write_html('dirty_data3.html')
    # fig2 = px.histogram(x=r)
    # fig2.show()
if __name__ == "__main__":
    #plot_data('./sc_stat.csv')
    plot_relx_data('./t1t2noe.csv')