from numpy import *
import pyn_fort as f

class cl_net:

    def __init__(self,file_net=None,file_keys=None):

        self.file_net=file_net
        self.file_keys=file_keys

        self.num_nodes=0
        self.weight_total=0
        self.k_total=0
        self.k_max=0

        if self.file_net!=None:
            self.read_net(self.file_net)
        if self.file_keys!=None:
            self.read_keys(self.file_keys)




    def read_net(self,name_file):

        ff=open(name_file,'r')
        line=ff.readline()
        self.num_nodes=int(line.split()[0])
        self.k_max=int(line.split()[1])
        self.k_total=int(line.split()[2])
        self.k_out=[]
        self.k_in=[]
        self.weight=[]
        self.weight_total=[]
        self.links_out=[]

        self.T_ind=zeros(shape=(self.k_total),dtype=int)
        self.T_start=zeros(shape=(self.num_nodes+1),dtype=int)
        self.T_weight=zeros(shape=(self.k_total),dtype=int)
        
        data=ff.read()
        ff.close()
        data2=[]
        for aa in data.split():
            data2.append(int(aa))


        jj=-1
        sumk=0
        for ii in range(self.num_nodes):
            jj+=1
            node_ind=data2[jj]
            k_out=data2[jj+1]
            weight=data2[jj+2]
            self.k_out.append(k_out)
            self.weight.append(weight)
            self.T_start[ii]=sumk
            self.links_out.append({})
            jj=jj+2
            for kk in range(k_out):
                jj+=1
                neigh=data2[jj]
                jj+=1
                flux=data2[jj]
                self.T_ind[sumk]=neigh
                self.T_weight[sumk]=flux
                sumk+=1
                self.links_out[ii][neigh-1]=flux
        self.T_start[self.num_nodes]=sumk


        print '# Number of nodes:',self.num_nodes
        print '# Number of links:',self.k_total
        print '# Max. conectivity:',self.k_max



    def read_keys(self,name_file):

        ff=open(name_file,'r')
        self.key=[]
        for ii in range(self.num_nodes):
            line=ff.readline()
            mss=line.split()[1]+' |'
            for jj in range(2,6):
                mss=mss+' '+line.split()[jj]
            mss=mss+' |' 
            for jj in range(6,9):
                mss=mss+' '+line.split()[jj]
            mss=mss+' |'
            for jj in range(9,12):
                mss=mss+' '+line.split()[jj]
            mss=mss+' |'
            for jj in range(12,15):
                mss=mss+' '+line.split()[jj]
            mss=mss+' |'
            for jj in range(15,18):
                mss=mss+' '+line.split()[jj]
            self.key.append(mss)

    def symmetrize(self):

        temp=cl_net()
        temp.key=self.key
        temp.file_net=self.file_net
        temp.file_keys=self.file_keys
        temp.num_nodes=self.num_nodes
        aux={}
        for ii in range(self.num_nodes):
            for jj in self.links_out[ii].keys():
                aux[(ii,jj)]=0
                aux[(jj,ii)]=0
        temp.k_total=len(aux)
        del(aux)        
        pfff=f.aux_funcs.symmetrize_net(temp.k_total,self.T_ind,self.T_weight,self.T_start,self.num_nodes,self.k_total)
        temp.k_max=pfff[0]
        temp.T_weight=pfff[1]
        temp.T_ind=pfff[2]
        temp.T_start=pfff[3]
        temp.weight=pfff[4]
        temp.weight_total=sum(temp.weight)
        temp.k_out=zeros(temp.num_nodes)
        temp.links_out=[]
        for ii in range(temp.num_nodes):
            temp.links_out.append({})
            for jj in range(temp.T_start[ii],temp.T_start[ii+1]):
                neigh=temp.T_ind[jj]
                temp.links_out[ii][neigh-1]=temp.T_weight[jj]
                temp.k_out[ii]+=1

        temp.k_in=[]

        print '# Number of nodes:',temp.num_nodes
        print '# Number of links:',temp.k_total
        print '# Max. conectivity:',temp.k_max


        return temp

    def gradient_clusters(self):

        self.num_clusters,self.node_belongs2=f.aux_funcs.grad(self.weight,self.T_ind,self.T_weight,self.T_start,self.num_nodes,self.k_total)

        print '# Number of clusters: ',self.num_clusters
        Clust={}
        Aux={}
        for ii in range(len(self.node_belongs2)):
            try:
                Clust[self.node_belongs2[ii]].append(ii)
            except:
                Clust[self.node_belongs2[ii]]=[]
                Clust[self.node_belongs2[ii]].append(ii)
        a=0
        repre=[]
        for ii in Clust.keys():
            Aux[a]=Clust[ii]
            repre.append(int(ii))
            a+=1
        del Clust
        self.representants=repre
        self.clust_info=Aux
        Aux1={}
        for ii in range(len(self.representants)):
            Aux1[self.representants[ii]]=ii
        for ii in range(len(self.node_belongs2)):
            self.node_belongs2[ii]=Aux1[self.node_belongs2[ii]]
        self.cluster_weight=zeros(self.num_clusters)
        for ii in range(self.num_clusters):
            for jj in self.clust_info[ii]:
                self.cluster_weight[ii]+=self.weight[jj]
