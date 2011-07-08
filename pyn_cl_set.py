from pyn_cl_unit import *
from pyn_cl_coors import *
from pyn_cl_net import *
from os import system
from os import path
from numpy import *
import pylab
import pyn_fort as f
#import pyn_fort_gnm as f_gnm
import pickle as pic

class cl_set:
#--------------------------------------------------------------
    def __init__(self,input_file=None,input_selection=None,download=None):

        # Instantation options:
        self.file=input_file
        self.select=input_selection
        self.file_hbonds=''
        self.file_mss=''
        self.file_shell=''

        # > Topological properties
        self.num_atoms=0
        self.name=''
        self.index=0
        self.pdb_index=0
        self.atom=[]       #list of objects unit
        self.list_atoms=[]
        self.acceptors=[]
        self.donors=[]
        self.dimensionality=0
        self.ss_pdb=[]
        self.chains=[]
        self.water_model=''

        # > Coordinates
        self.recent_frame=0
        self.coors=[]      # list of objects coor
        self.dist_matrix=[]

        ##################################

        # BUILDING THE TOPOLOGY OF THE SET FROM A FILE

        ### Download option:

        if download:
            if not path.exists(download):
                temp='wget -nv http://www.rcsb.org/pdb/files/'+download
                system(temp)

            input_file=download
            self.file=input_file

        ### Reading the file and attaching the atoms to the set

        if self.file:

            if self.file.endswith('pdb'):
                self.read_pdb(self.file)

            elif self.file.endswith('gro'):
                self.read_gro(self.file)

            ### Setting up other atoms' attributes

            ### Setting up the chains
            for aa in self.atom[:]:
                if aa.type_pdb in ['ATOM']:
                    if aa.chain not in self.chains:
                        self.chains.append(aa.chain)


               ## Auxiliary dictionary
            
            before=-99
            ii=-1
            for aa in self.atom:
                if aa.resid_pdb_index != before :
                    before=aa.resid_pdb_index
                    ii+=1
                aa.resid_index=ii


            aux={}

            for aa in self.atom[:]:
                ii=aa.resid_index
                try: 
                    aux[ii][aa.name]=aa.index
                except:
                    aux[ii]={}
                    aux[ii][aa.name]=aa.index

            for aa in self.atom[:]:
                if aa.name in ['OW']:
                    aa.acceptor=True
                    aa.polar_class='acceptor'
                    aa.polarizability=True
                    aa.covalent_bond.append(aux[aa.resid_index]['HW1'])
                    aa.covalent_bond.append(aux[aa.resid_index]['HW2'])
                if aa.name in ['O'] and aa.resid_name in ['HOH','SOL','HO4']:
                    aa.acceptor=True
                    aa.polar_class='acceptor'
                    aa.polarizability=True
                if aa.name in ['HW1','HW2']:
                    aa.donor=True
                    aa.polar_class='donor'
                    aa.polarizability=True
                    aa.covalent_bond.append(aux[aa.resid_index]['OW'])
        
            ### Setting up the residues

            self.residue=[]
            for aa in aux.keys():
                temp_residue=cl_residue()
                temp_residue.index=aa
                temp_residue.list_atoms=aux[aa].values()
                ii=temp_residue.list_atoms[0]
                temp_residue.pdb_index=self.atom[ii].resid_pdb_index
                temp_residue.name=self.atom[ii].resid_name
                self.residue.append(temp_residue)


            """ Check Roman's way to implement residues """

            ### Setting up the waters

            self.water=[]
            for aa in self.residue[:]:
                if aa.name in ['HOH','SOL','HO4']:
                    temp_water=cl_water()
                    temp_water.name=aa.name
                    temp_water.index=aa.index
                    temp_water.pdb_index=aa.pdb_index
                    temp_water.list_atoms=aa.list_atoms
                    if self.water_model=='':
                        self.water_model='tip'+str(len(aa.list_atoms))+'p'
                    if 'HW1' in aux[aa.index].keys():
                        temp_water.H1=aux[aa.index]['HW1']
                        temp_water.H2=aux[aa.index]['HW2']
                        temp_water.O=aux[aa.index]['OW']
                    if 'O' in aux[aa.index].keys():
                        temp_water.O=aux[aa.index]['O']
                    self.water.append(temp_water)

            ### Setting up global attributes
            ii=self.file[::-1].find('.')
            self.name=self.file[:-ii]
            self.num_atoms=len(self.atom)
            self.dimensionality=self.num_atoms*3
            for aa in self.atom[:]:
                if aa.acceptor: self.acceptors.append(aa.index)
                if aa.donor: self.donors.append(aa.index)
            self.acceptors=array(self.acceptors,order='Fortran')
            self.donors=array(self.donors,order='Fortran')
            self.num_residues=len(self.residue)
            self.num_waters=len(self.water)


            ### Print info:
            self.info()

        ####################################



        elif self.select:

            '''Bad type of file'''



############################################################################
############################################################################
#>>>>>>>>>> FILE.PDB

    def read_pdb (self,name_file):

        for line in open(name_file,'r'):
            ss=line.split()

            if ss[0] in ['HELIX','SHEET','TURN']:
                self.ss_pdb.append(line)
            if ss[0] in ['ATOM','HETATM']:

                temp_atom=cl_unit()
                temp_atom.type_pdb=line[0:6].replace(' ', '')
                temp_atom.pdb_index=int(line[6:11])
                temp_atom.name=(line[12:16].split())[0]
                temp_atom.alt_loc=line[16]
                temp_atom.resid_name=(line[17:20])
                temp_atom.chain=line[21]
                temp_atom.resid_pdb_index=int(line[22:26])
                temp_atom.code_ins_res=line[26]
                temp_atom.occup=float(line[54:60])
                temp_atom.bfactor=float(line[60:66])
                temp_atom.seg_ident=line[72:76].replace(' ', '')
                temp_atom.elem_symb=line[76:78].replace(' ', '')
                temp_atom.charge=line[78:80].replace(' ', '')

                temp_atom.index=len(self.atom)
                self.atom.append(temp_atom)

#>>>>>>>>>> FILE.GRO

    def read_gro (self,name_file):
        '''Reading a gro file'''

        f=open(name_file,'r')

        line=f.readline()                                          # Header of the gro file

        line=f.readline()                                        
        self.num_atoms=int(line)

        for i in range(self.num_atoms):           
            
            temp_atom=cl_unit()

            line=f.readline().split()
            temp_atom.pdb_index=int(line[2])
            temp_atom.name=line[1]
            temp_atom.resid_name=line[0][-3:]
            temp_atom.resid_pdb_index=int(line[0][:-3]) 

            temp_atom.index=i           

            self.atom.append(temp_atom)



#>>>>>>>>>> MAKING....


############################################################################
############################################################################


    def load_coors (self,input_file,frame=None,begin=None,end=None):

        self.coors_file=input_file
        self.traj_mark=0

        if self.coors_file.endswith('pdb'):
            self.read_coors_pdb(self.file)

        elif self.file.endswith('gro'):
            self.read_coors_gro(self.file)

        #elif self.file.endswith('xtc'):
        #    self.read_coors_xtc(self.file)

        elif self.coors_file.endswith('bin'):
            self.read_coors_bin(self.coors_file,frame,begin,end)


    def read_coors_pdb (self,name_file):

        temp_frame=cl_coors(name_file)
        self.coors.append(temp_frame)

    def read_coors_gro (self,name_file):

        temp_frame=cl_coors(name_file)
        self.coors.append(temp_frame)


    def read_coors_bin(self,name_file,frame,begin,end):
        if begin==None and frame==None and end==None:
            temp_frame=cl_coors(name_file,self.recent_frame)
            self.coors.append(temp_frame)
            self.recent_frame+=1
        elif begin==None and end==None and frame!=None:
            temp_frame=cl_coors(name_file,frame)
            self.coors.append(temp_frame)
            self.recent_frame=frame
        elif begin!=None and end!=None:
            for ii in range(begin,end):
                temp_frame=cl_coors(name_file,ii)
                self.coors.append(temp_frame)
            self.recent_frame=end
       


    def delete_coors (self,begin=None,end=None,frame=None):
        
        if frame==begin==end :
            print '#','deleting coors'
            del self.coors[:]


    def write_pdb (self,filename=None):
        
        if filename==None:
            print 'Enter filename: '
            print '      foo.write_pdb("foo.pdb")'
        else:
            write_pdb(filename,self)
    


#FUNCTIONS....

    def distance(self,pbc=False):
        
        
        for frame in self.coors:

            dist_frame=f.aux_funcs.dist(pbc,frame.xyz,frame.box,frame.xyz,self.num_atoms,self.num_atoms)
            
        self.dist_matrix.append(dist_frame)

    def info(self):

        print '#','System created from the file ',self.file,':'
        print '#',self.num_atoms,' atoms'
        print '#',self.num_residues,' residues'
        print '#',len(self.chains),' chains'
        print '#',self.num_waters,' waters'


###########################################################################
############################################################################
# GNM:

    def contact_map(self,cutoff=10.0):

        contact_map=f_gnm.contact_map(cutoff,self.coors[0].xyz,self.num_atoms)
        return array(contact_map,dtype=bool)


############################################################################
############################################################################
# Para los hbonds

    def neighbs(self,system2=None,limit=0,dist=0.0,pbc=False):


        if system2==None:
            system2=self
            ident=True
        else:
            ident=False

        if limit != 0:
            neighbs=f.aux_funcs.neighbs_limit(pbc,ident,limit,self.coors[0].xyz,self.coors[0].box,system2.coors[0].xyz,self.num_atoms,system2.num_atoms)
        
        else:
            print type(self.coors[0].xyz[1][:]),self.coors[0].xyz[1][:],system2.num_atoms
            neighbs=f.aux_funcs.neighbs_dist(pbc,ident,dist,self.coors[0].xyz,self.coors[0].box,system2.coors[0].xyz,system2.num_atoms)


        return neighbs

    



#    def donors_neighbs_accept(self,limit=None):

#        a=f.aux_funcs.don_neighbs_accept(limit,self.coors[0].xyz,self.coors[0].box,self.donors,self.acceptors,self.num_atoms,len(self.donors),len(self.acceptors))
#        return a



    
    def hbonds(self,syst=None,opt=None,neighbs=None):

        if syst==None:
            syst=self
            ident=True
        else:
            ident=False

        if opt=='Skinner':
            hbonds_out=f.aux_funcs.hbonds_skinner(self.coors[0].xyz,self.coors[0].box,self.donors,self.acceptors,self.num_atoms,len(self.donors),len(self.acceptors))
            
        if opt=='Roh':
            hbonds_out=f.aux_funcs.hbonds_roh(ident,self.coors[0].xyz,self.coors[0].box,self.donors,self.acceptors,
                                          syst.coors[0].xyz,syst.donors,syst.acceptors,neighbs[0],neighbs[1],neighbs[2],
                                          self.num_atoms,len(self.donors),len(self.acceptors),syst.num_atoms,
                                          len(syst.donors),len(syst.acceptors),len(neighbs[0][0]))
        


        hbonds={}

        if ident==True :
            for ii in self.donors:
                for jj in syst.acceptors:
                    if hbonds_out[0][ii,jj]==-1 and jj not in self.atom[ii].covalent_bond[:]:
                        
                        try: 
                            hbonds[ii][jj]=hbonds_out[1][ii,jj]
                        except:
                            hbonds[ii]={}
                            hbonds[ii][jj]=hbonds_out[1][ii,jj]

                        try: 
                            hbonds[jj][ii]=hbonds_out[1][jj,ii]
                        except:
                            hbonds[jj]={}
                            hbonds[jj][ii]=hbonds_out[1][jj,ii]  



        num_hbonds=0
        for aa in hbonds.values():
            num_hbonds+=len(aa)
        print num_hbonds

        return hbonds


    def rdf(self):

        

        delta_x=0.2
        t=ceil(self.coors[0].box[0][0]/delta_x)
        YY=zeros(t)
        XX=zeros(t)
       # dist_matrix=self.distance(pbc=True)
        for kk in range(len(self.coors)):

            self.distance(pbc=True)
            t=ceil(self.dist_matrix[kk]/delta_x)
            for ii in t:
                YY[t]+=2


        """
            for ii in range(self.num_atoms):
                for jj in range(ii+1,self.num_atoms):
                    t=ceil(self.dist_matrix[kk][ii][jj]/delta_x)
                    YY[t]+=2
            dens=(self.coors[kk].box[0][0]*self.coors[kk].box[1][1]*self.coors[kk].box[2][2])/(self.num_atoms)
        """
       



        RDF=[]
        for ii in range(len(YY)):
            XX[ii]=delta_x*ii
            YY[ii]=YY[ii]+YY[ii]*dens
            
            YY[ii]=YY[ii]/((XX[ii]**2)*delta_x)
            YY[ii]=YY[ii]/(4*pi*self.num_atoms*len(self.coors))
            A=[XX[ii],YY[ii]]
            RDF.append(A)
        
        RDF=0
        #pylab.plot(XX,YY)
        #pylab.show()
        return RDF
        

#######################################################
#######################################################
### Analysis of the entire water trajectory:

    def water_analysis (self,traj_name=None,write_out=None,init_frame=0,last_frame=-1):

        if traj_name == None:
            print 'input the traj_name'
            return None
        if traj_name.endswith('.xtc'):
            if not path.exists(traj_name[:-3]+'bin'):
                xtc2bin(traj_name,traj_name[:-3]+'bin')
            traj_name=traj_name[:-3]+'bin'

        if write_out==None:
            write_out='No'
            

        command='./pyn_anw '+traj_name+' '+write_out+' '+self.water_model+' '+str(self.num_residues)+' '+str(init_frame)+' '+str(last_frame)
        system(command)


        self.file_hbonds='aux_hbs.bin'
        self.file_mss='aux_mss.bin'
        self.file_shell='aux_shell.bin'
        self.file_net='aux_net.oup'
        self.file_key_mss='aux_key_mss.oup'

    def get_hbonds(self,frame=None,molecule=None):

        L=[]
        HB=open(self.file_hbonds,'rb')
        HB.seek((2*4*4*(frame-1)*self.num_waters)+(4*(frame)))
        for ii in range(self.num_waters*2):    
            B=(stc.unpack('3i',HB.read(12)))
            S=(stc.unpack('1f',HB.read(4)))
            L.append([B,S])
        if molecule!=None:
            
            print L[(molecule-1)*2][0][0],'H1', L[(molecule-1)*2][0][2],'O',L[(molecule-1)*2][1][0]
            print L[(molecule-1)*2+1][0][0],'H2', L[(molecule-1)*2+1][0][2],'O',L[(molecule-1)*2+1][1][0]

            
            for jj in range(self.num_waters*2):
                if molecule==L[jj][0][2]:
                    if L[jj][0][1]==1:
                        print L[jj][0][0],'O', L[jj][0][2],'H1',L[jj][1][0]
                    else:
                        print L[jj][0][0],'O', L[jj][0][2],'H2',L[jj][1][0]
        else:
            print L
    

    def get_mss(self,frame=None,molecule=None):
        L=[]    
       
        BB=open(self.file_mss,'rb')

        BB.seek((17*4*(frame-1)*self.num_waters)+(4*(frame)))
        for ii in range(self.num_waters):    
            B=(stc.unpack('17i',BB.read(68)))
            L.append(B)
        if molecule!=None:
            print L[molecule-1]
        else:
            print L

    def get_shell(self,frame=None,molecule=None):
        L=[]    
       
        BB=open(self.file_shell,'rb')

        BB.seek((17*4*(frame-1)*self.num_waters)+(4*(frame)))
        for ii in range(self.num_waters):    
            B=(stc.unpack('17i',BB.read(68)))
            L.append(B)
        if molecule!=None:
            print L[molecule-1]
        else:
            print L

    def get_network(self):

        return cl_net(self.file_net,self.file_key_mss)






#######################################################
#######################################################
#######################################################

def plot_contact_map(contact_map):

    pylab.gray()
    pylab.imshow(contact_map==False,origin='lower',interpolation=None) # igual seria mejor no interpolar
    #pylab.matshow(contact_map==False)
    return pylab.show()



def norm_vect(vector):
    return sqrt(dot(vector,vector.conj()))


#######################################################
#######################################################
#######################################################
#### Handling files and trajs:

#### write pdb

def xtc2bin (xtc=None,bin=None):

    if xtc==None:
        print 'Enter the name of the xtc file:'
        print '      xtc2bin (xtc="traj.xtc",bin="traj.bin")'
    else:
        if bin==None:
            bin=xtc[:-3]+'bin'
            print 'writing the file ',bin
        
        f.aux_funcs.xtc2bin(xtc,bin)
        

def write_pdb(file_name,sel):
        file=open(file_name,'w')
        for ii in sel.ss_pdb:
            file.write(str(ii))
        for ii in range(sel.num_atoms):
            a='ATOM  '                                  # 1-6
            a+="%5d" % (ii+1)                           # 7-11
            #a+="%5d" % sel.atom[ii].pdb_index          # 7-11
            a+=' '                                      # 12
            a+="%-4s" % sel.atom[ii].name           # 13-16
            a+=' '                                      # 17
            a+="%3s" % sel.atom[ii].resid_name          # 18-20
            a+=' '                                      # 21
            a+="%1s" % sel.atom[ii].chain               # 22
            a+="%4d" % sel.atom[ii].resid_pdb_index     # 23-26
            a+=' '                                      # 27
            a+='   '                                    # 28-30
            a+="%8.3f" % float(sel.coors[0].xyz[ii][0]) # 31-38
            a+="%8.3f" % float(sel.coors[0].xyz[ii][1]) # 39-46
            a+="%8.3f" % float(sel.coors[0].xyz[ii][2]) # 47-54
            a+="%6.2f" % sel.atom[ii].occup             # 55-60
            a+="%6.2f" % sel.atom[ii].bfactor           # 61-66
            a+='          '                             # 67-76
            a+="%2s" % sel.atom[ii].elem_symb           # 77-78
            a+="%2s" % sel.atom[ii].charge              # 79-80
            a+='\n' 
            file.write(str(a))         
        file.close()
        return None

#######################################################
#######################################################









#######################################################
#### Selection algorithm:
#######################################################

#############Subfunction for make parenthesis in any way
def prolong_string(condition):
    b=' '
    for ii in range(len(condition)):
        if condition[ii]=='(' or condition[ii]==')':
           b+=' '
           b+=condition[ii]
           b+=' '
        else:
            b+=condition[ii]
    return b

################Beginning of reading the string with logical operators(or,and)##############
########Block for "and"########

def make_selection(system,condition):

    if type(condition)==str:

        condition=prolong_string(condition)
        if 'and' in condition:
                st=condition.split()  
                st2=st
                aux_2=[]
                for ii in st:
                    try:
                        aux_1=st.index('and')
                        aux_2.append(aux_1)
                        st[aux_1]=''
                           
                    except:
                          a=0
                a=aux_2[0]
                
                b=aux_2[len(aux_2)-1]
                begini= st[:a]
                endi= st[b:]
                BEG=' '
                for jj in range(len(begini)):
                    BEG+=' '
                    BEG+=begini[jj]
                ENDI=' '
                for jj in range(len(endi)):
                    ENDI+=' '
                    ENDI+=endi[jj]
                BEG_LIST=good_select(system,BEG)
                END_LIST=good_select(system,ENDI)
                pp=[]                            #additional loops
                # len(aux_2)-1         - test to count how many options are inside the condition
         
                lyst=[]
         
                for ii in range(len(aux_2)-1):  #########conditions which are not first neither last
                     for jj in range(len( st[aux_2[ii]:aux_2[ii+1]])):
                         pp=st[aux_2[ii]:aux_2[ii+1]]
                         PIPO=' '
                         for fff in range(len(pp)):               #making string from pieces
                             PIPO+=' '
                             PIPO+=pp[fff]
                         lex=good_select(system,PIPO)
                         for yy in range(len(lex)):
                             lyst.append(lex[yy])
                SHY=[]
                for ii in range(system.num_atoms):
                   # print  lyst.count(ii), len(aux_2)
                    if lyst.count(ii)>=(len(aux_2)-1):
                        SHY.append(ii)
                       
                LAST_LIST=[]
                for jj in range(system.num_atoms):
                    if jj in BEG_LIST and jj in END_LIST and jj in SHY:
                        LAST_LIST.append(jj)
                sux=appending_sel(system,LAST_LIST)
        else:
            LIST=good_select(system,condition)
            sux=appending_sel(system,LIST)

    elif type(condition) in [list,ndarray,tuple]:
        LIST=condition
        sux=appending_sel(system,LIST)
    else:
        print "ERROR sel01"

    #for kk in range(sux.num_atoms):
    #    sux.atom[kk].index=kk
   
   
    return sux

####################Function for make selection from piece of ye string

def good_select(system,condition):


    list_of_ind=[]
    temp=condition.split()
    or_index=[]
    and_index=[]
    par_index_op=[]
    par_index_cl=[]
    for ii in range(len(temp)):
    
        if temp[ii]=='(':
            par_index_op.append(ii)
        elif temp[ii]==')':
            par_index_cl.append(ii)
    cond_ar=[]
    for ii in range(len(par_index_op)):
   
        a=(par_index_op[ii],par_index_cl[ii])
        cond_ar.append(a)

   # print cond_ar   - pairs of parenthesis


    for ii in range(len(cond_ar)):
        aaa=[]                                   # for handling place(s) of 'or'
        possib= temp[cond_ar[ii][0]+1:cond_ar[ii][1]]    #possible value of parameter
        param= temp[cond_ar[ii][0]-1]                    # parameter

    
        Numb_or=possib.count(' or ')
        
        if Numb_or==0:          # if only one possible value of parameter - seems to be easiest case
           
            if param=='atom_name':   
                list_of_ind=g_atom_name(system,possib,list_of_ind) 
              
            elif param=='resid_name':
                list_of_ind=g_resid_name(system,possib,list_of_ind)
            elif param=='donors':
                possib='donor'
                list_of_ind=g_donors(system,possib,list_of_ind)
            elif param=='acceptors':
                possib='acceptor'
                list_of_ind=g_acceptors(system,possib,list_of_ind)
            elif param=='atom_index':
                list_of_ind=g_atom_index(system,possib,list_of_ind)
            else:
                print 'ERROR sel02: unknown parameter ',param
             
          #  print len(list_of_ind)
        for jjj in range(Numb_or):  # if there is a lot of possible options
                      
            j=possib.index(' or ')
            aaa.append(j)
            possib[j]=''
        
        aaa.append(len(possib))
        for jjj in range(Numb_or+1):
            possib2=possib[aaa[jjj-1]-1]# extractin' possible value of parameters from list like ' X or Y or Z ' - previous line for this piece of function  also
           
            if param=='atom_name':
                list_of_ind=g_atom_name(system,possib2,list_of_ind)
                
            elif param=='resid_name':   
                list_of_ind=g_resid_name(system,possib2,list_of_ind)
            elif param=='donors':
              
                list_of_ind=g_donors(system,possib2,list_of_ind)
            elif param=='acceptors':
                list_of_ind=g_acceptors(system,possib2,list_of_ind)
            elif param=='atom_index':
                list_of_ind=g_atom_index(system,possib2,list_of_ind)
            else:
                 print 'ERROR : unknown parameter ',param
    #######after all -  output is  list of index ##############
    return(list_of_ind)



### additional func for all parameters ###

def g_atom_name(system,possib,list_of_ind):    
        for ii in range(system.num_atoms):
           
            if system.atom[ii].name==possib:
               
                list_of_ind.append(ii)
                
        return list_of_ind

def g_resid_name(system,possib,list_of_ind):    
        for ii in range(system.num_atoms):

            if system.atom[ii].resid_name==possib:
               
                list_of_ind.append(ii)
             
        return list_of_ind

def g_donors(system,possib,list_of_ind):    
        for ii in range(system.num_atoms):
            if system.atom[ii].donor==True:
                list_of_ind.append(ii)
             
        return list_of_ind

def g_acceptors(system,possib,list_of_ind):    
        for ii in range(system.num_atoms):

            if system.atom[ii].acceptor==True:
               
                list_of_ind.append(ii)
             
        return list_of_ind

def g_atom_index(system,possib,list_of_ind):
    
    if ':' in possib:
        a=possib
        i = a[a.index('[')+1:a.index(':')]
        j = a[a.index(':')+1:a.index(']')]
       
        if int(j)>system.num_atoms:
            j=system.num_atoms
            print 'Total number of atoms in your system is ', system.num_atoms
        for kk in range(int(i),int(j)):
            list_of_ind.append(kk)
        
    return list_of_ind

##############appending from list of atom indexes###########################

def appending_sel(syst,list_of_ind):   
    temp_set=cl_set()

    ############Appending atoms to selection#################
    for ii in range(syst.num_atoms):
        if ii in list_of_ind:
           temp_set.atom.append(syst.atom[ii])
      

    # Building global attributes:
    temp_set.num_atoms=len(temp_set.atom)

 
    ##########Appending coordinates to selection#########################
    for jj in syst.coors:
        temp_frame=cl_coors()
        temp_frame.box=jj.box
  
        for ii in range(syst.num_atoms):         
                    if ii in list_of_ind:
                     
                        temp_frame.xyz.append(jj.xyz[ii])
        temp_set.coors.append(temp_frame)

   
         


    return temp_set

def write_set_to_file(system,name_of_file):
    file=open(name_of_file,'w')
    pic.dump(system,file)
    file.close()

def read_set_from_file(self,name_of_file):
        file=open(name_of_file,'r')
        A=pic.load(file)
        file.close()
        return A

#### Romans's new stuff
#######################################################
#######################################################
#######################################################

#######################################################
#######################################################
#######################################################
#### To build the oscillations with the anm modes


def build_fluct_anm(system,anm,mode='all',output='None',amplitude=8.0,steps=60):

    prov_list=[]
    for ii in system.atom:
        if ii.name in ['N','CA','C','O'] and ii.type_pdb =='ATOM':
            prov_list.append(ii.index)

    prov_system=make_selection(system,prov_list)

    for aa in prov_system.atom[:]:
        if aa.type_pdb in ['ATOM']:
            if aa.chain not in prov_system.chains:
                prov_system.chains.append(aa.chain)


    num_nodes=len(anm.eigenvects[:])
    list_modes=[]

    if mode=='all':
        num_modes=num_modes
        for ii in range(1,num_modes+1):
            list_modes.append(ii)
    elif type(mode)==int:
        num_modes=1
        list_modes.append(mode)
    elif type(mode) in [list,tuple]:
        num_modes=len(mode)
        list_modes=list(mode)
    
    osc=zeros(shape=(prov_system.num_atoms,3))

    prefix=system.name
    if prefix[-1]=='.':
        prefix=prefix[:-1]

    in_net=[]
    for ii in anm.system.atom:
        in_net.append(ii.index)
    in_syst=[]
    for ii in prov_system.atom:
        in_syst.append(ii.index)

    tt=zeros(shape=(prov_system.num_atoms))


    jj=-1

    for chch in prov_system.chains :
        interr=-1
        for ii in anm.system.atom:
            if ii.chain == chch:
                extreme=ii.index
                extreme=in_net.index(extreme)

        for ii in prov_system.atom:
            if ii.chain == chch:
                jj+=1
        
                if ii.index in in_net:
                    interr+=1
                    net_initial=in_net.index(ii.index)
                    net_end=net_initial+1
                    if net_end>extreme:
                        interr=-1
                    else:
                        initial=jj
                        end=in_syst.index(in_net[net_end])

                if interr==-1:
                    tt[jj]=1.0
                else:
                    tt[jj]=dot((prov_system.coors[0].xyz[jj]-prov_system.coors[0].xyz[initial]),(prov_system.coors[0].xyz[end]-prov_system.coors[0].xyz[initial]))
                    tt[jj]=tt[jj]/(dot((prov_system.coors[0].xyz[end]-prov_system.coors[0].xyz[initial]),(prov_system.coors[0].xyz[end]-prov_system.coors[0].xyz[initial])))


    for ind_mode in list_modes :
        kk=0
        jj=-1
        for chch in prov_system.chains :
            interr=-1
            for ii in anm.system.atom:
                if ii.chain == chch:
                    extreme=ii.index
                    extreme=in_net.index(extreme)
            ant=anm.eigenvects[ind_mode-1][kk]
            for ii in prov_system.atom:
                if ii.chain == chch:
                    jj+=1
        
                    if ii.index in in_net:
                        interr+=1
                        net_initial=in_net.index(ii.index)
                        net_end=net_initial+1
                        if net_end>extreme:
                            interr=-1
                            ant[:]=anm.eigenvects[ind_mode-1][kk]
                        else:
                            aaa=anm.eigenvects[ind_mode-1][net_initial]
                            bbb=anm.eigenvects[ind_mode-1][net_end]
                        kk+=1
                    if interr==-1:
                        osc[jj][:]=ant[:]
                    else:
                        osc[jj][:]=aaa[:]+(bbb[:]-aaa[:])*tt[jj]
            

        file_name=prefix+'_anm_'+str(ind_mode)+'.pdb'
        file=open(file_name,'w')
        a='HEADER    '+prefix+'     ANM: Mode '+str(ind_mode)+'\n'
        file.write(str(a))

        delta_f=2.0*pi/(steps*1.0)

        for frame in range(0,steps):

            a='MODEL '+str(frame)+'\n'
            file.write(str(a))

            for ii in system.ss_pdb:
                file.write(str(ii))

            for ii in range(prov_system.num_atoms):
                a='ATOM  '                                 # 1-6
                a+="%5d" % (ii+1)                          # 7-11
                #a+="%5d" % prov_system.atom[ii].pdb_index  # 7-11
                a+=' '                                     # 12
                a+=' '+"%-3s" % prov_system.atom[ii].name  # 13-16
                a+=' '                                     # 17
                a+="%3s" % prov_system.atom[ii].resid_name # 18-20
                a+=' '                                     # 21
                a+="%1s" % prov_system.atom[ii].chain      # 22
                a+="%4d" % prov_system.atom[ii].resid_pdb_index # 23-26
                a+=' '                                     # 27
                a+='   '                                   # 28-30
                a+="%8.3f" % float(prov_system.coors[0].xyz[ii][0]+amplitude*sin(delta_f*frame)*osc[ii][0]) # 31-38
                a+="%8.3f" % float(prov_system.coors[0].xyz[ii][1]+amplitude*sin(delta_f*frame)*osc[ii][1]) # 39-46
                a+="%8.3f" % float(prov_system.coors[0].xyz[ii][2]+amplitude*sin(delta_f*frame)*osc[ii][2]) # 47-54
                a+="%6.2f" % prov_system.atom[ii].occup    # 55-60
                a+="%6.2f" % prov_system.atom[ii].bfactor  # 61-66
                a+='          '                            # 67-76
                a+="%2s" % prov_system.atom[ii].elem_symb  # 77-78
                a+="%2s" % prov_system.atom[ii].charge     # 79-80
                a+='\n' 
                file.write(str(a))

            a='ENDMDL \n'
            file.write(str(a))

        file.close() 


    return


def xtc2bin(xtc_name,bin_name):
    command='./xtc2bin %s %s'%(xtc_name,bin_name)

    if path.exists(bin_name):
        command2='mv %s %s#'%(bin_name,bin_name)
        print 'file',bin_name,' was moved to ', bin_name+'#'
        system(command2)

    system(command)









#######################################################
#######################################################
#######################################################



"""
    def hb()

    def dist()

    def neighb()

    def com()







        #self.molid= #parent molecular object id

#--------------------------------------------------------------
#    def com( self ):
#        '''Returning the center of mass'''

#--------------------------------------------------------------
#    def neighbors( self, radii=None, system=None ):
#        '''Calculating neighbors.
#           It returns a list of residue objects'''

#--------------------------------------------------------------
#    def hb( self, criteria=None, system=None, residue=None ):
#        '''Calculating hydrogen bonds:
#               criteria example: ('sk',0.0085)
#           NB: it calculates all the hydrogen-bonds of a residue
#               when 'system' is provided, or the hydrogen-bonds
#               between residues when 'residue' is given'''


"""



#################################################################
#################################################################
#################################################################
#################################################################
#################################################################

class cl_water(cl_set):
    '''Water specific attributes and methods'''
#--------------------------------------------------------------
#    def __init__( self ):
#        '''The Water object inherits from Residue'''
#
#        pass

    O=-1
    H1=-1
    H2=-1
    
    uvect_norm=[]


#    def get_uvect_norm(self):
        


class cl_residue(cl_set):
#    '''Water specific attributes and methods'''
#--------------------------------------------------------------
#    def __init__( self ):
#        '''The Water object inherits from Residue'''
#
    pass





""" Roman's way to implement residues """
"""
if self.file:
            self.residues=res_sel(self)

####################SELECTION_of_RES#####################
def res_sel(system):
    list_of_ind=[]
    residues=[]
    for ii in range(system.num_atoms):
        a=system.atom[ii].resid_pdb_index
     
        if a==(system.atom[ii-1].resid_pdb_index):
            list_of_ind.append(ii)
            
        else:
            if list_of_ind!=[]:
                res=cl_set()
                res=appending_sel(system,list_of_ind)
               
                residues.append(res)

            list_of_ind=[]
            list_of_ind.append(ii)
  
    return residues


"""
