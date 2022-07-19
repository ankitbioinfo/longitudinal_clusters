
def read_db():

    f=open('sort_3_db_L_R_high_confident.dat','r')
    f=open('comb.dat','r')
    contdb=f.readlines()
    totalLRpairs={}
    ligand_one={}
    receptor_one={}
    combined_one={}

    for j in range(len(contdb)):
        l=contdb[j][0:-1].split()
        flag1=0
        flag2=0
        '''
        for i in range(n):
            t=cont[i][0:-1].split('\t')
            if (t[0].upper()==l[0].upper()):
                flag1=1
            if (t[0].upper()==l[1].upper()):
                flag2=1

        if flag1+flag2==2:
        #if (flag1+flag2)>0:
        '''
        totalLRpairs[l[0].capitalize()+'--'+l[1].capitalize()]=1
        #    ligand_one[l[0]]=1
        #    receptor_one[l[1]]=1
        #    combined_one[l[0]]=1
        #    combined_one[l[1]]=1

            #totalLRpairs[l[1]+'--'+l[0]]=1

    fw=open('temp.dat','w')
    for key in totalLRpairs:
        fw.write(key+'\n')

    return totalLRpairs





def find_pairs(cc,nc,totalLRpairs):
    Found1=[]
    Found2=[]
    for i in range(len(cc)):
        for j in range(len(nc)):
            name1=cc[i]+'--'+nc[j]
            name2=nc[j]+'--'+cc[i]
            if name1 in totalLRpairs:
                Found1.append(name1)
            if name2 in totalLRpairs:
                Found2.append(name2)
    return Found1,Found2

def main():
    totalLRpairs=read_db()
    f=open("PC_loading_prediction_2.dat")
    cont=f.readlines()

    pos=[]
    pos.append(0)
    for i in range(len(cont)):
        if cont[i][0]=='\n':
            pos.append(i)

    pos.append(len(cont))


    listofallLR={}

    for i in range(1,len(pos)):
        matter=cont[pos[i-1]+1:pos[i]+1]
        #if len(matter)>5:
        print('\n')
        cc=[]
        nc=[]
        for j in range(len(matter)):
            name=matter[j][0:3]
            #print(name)
            l=matter[j][0:-1].split(',')
            if name=='CC,':
                for k in range(1,len(l)):
                    cc.append(l[k])
            if name=='NC,':
                for k in range(1,len(l)):
                    nc.append(l[k])
            if name=='CC-':
                start=matter[j][0:-1]


        if len(nc)>0:
            title=matter[0][0:-1]
            if len(title)<30:
                print(title)
            print(start)
            found1,found2=find_pairs(cc,nc,totalLRpairs)
            print(found1)
            print(found2)

            for k in range(len(found1)):
                l=found1[k].split('--')
                listofallLR[l[0]]=1
                listofallLR[l[1]]=1

            for k in range(len(found2)):
                l=found2[k].split('--')
                listofallLR[l[0]]=1
                listofallLR[l[1]]=1

    fw=open('listofallLR.dat','w')
    for key in listofallLR:
        fw.write(key+'\n')

    print(sorted(list(listofallLR.keys())))

main()
