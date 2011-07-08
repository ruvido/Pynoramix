import numpy as npy
import struct as stc #To read binary files

class cl_coors:

    def __init__(self,file_input=None,frame=None):
        self.file=file_input
        self.time=0
        self.frame=0
        self.step=0
        self.precision=0


        self.xyz=[]
        self.box=npy.zeros(shape=(3,3),order='Fortran')
        if file_input!=None:

            if self.file.endswith('inp'):
                self.read_inp(self.file)
                
            if self.file.endswith('pdb'):
                self.read_pdb(self.file)
            
            if self.file.endswith('gro'):
                self.read_gro(self.file)

            if self.file.endswith('bin'):
                self.read_bin(self.file,frame)

            self.xyz=npy.array(self.xyz,order='Fortran')
        
            if self.file.endswith('bin') or self.file.endswith('gro') or self.file.endswith('xtc'):
                self.xyz=10.0*self.xyz
                self.box=10.0*self.box


#>>>>>>>>>>#>>>>>>>>>>
#>>>>>>>>>>#>>>>>>>>>>


    def read_pdb (self,name_file):


     for line in open(name_file,'r'):
            ii=line.split()

            if ii[0]=='CRYST1':                               # Reading the size of the box:
                self.box[0][0]=float(ii[1])                # Using the global variable (pyn_var_glob.py)
                self.box[1][1]=float(ii[2])                # for the size of the box vg.box
                self.box[2][2]=float(ii[3])               


            if ii[0] in ['ATOM','HETATM']:

                aux=(float(line[30:38]),float(line[38:46]),float(line[46:54]))
                self.xyz.append(aux)


#>>>>>>>>>> FILE.GRO

    def read_gro (self,name_file):


        f=open(name_file,'r')

        line=f.readline()                                          # Header of the gro file

        line=f.readline()                         
        num_atoms=int(line)                                   # Number of atoms

        for i in range(num_atoms): 
            line=f.readline().split()
            self.xyz.append(map(float,line[3:6]))

        line=f.readline().split()              # Reading the size of the box
        
        self.box[0][0]=float(line[0])       # Using the global variable (pyn_var_glob.py) 
        self.box[1][1]=float(line[1])       # for the size of the box vg.box      
        self.box[2][2]=float(line[2])       

        
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    READING THE TRAJECTORY   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    def update_coors ( self, name_file,frame=None):     
        '''Updating coordinates'''

        self.traj_file=name_file                     # Name of the trajectory file
                                                     # self.frame was defined at the beginning
                                                     # If frame==None, the next frame will be read

        if self.traj_file.endswith('bin'):            # Reading binary traj.
            self.read_bin(name_file,frame)
            
        #if self.name_aux.endswith('xtc'):            # Readint binary traj.
        #    self.read_traj_bin(name_file,frame)      <--- To be inplemented


#>>>>>>>>>> TRAJ.BIN (Binary)
    def read_bin (self, name_file,frame):
        self.traj_file=name_file  
        FF=file(self.traj_file)

      
       # if FF.closed==False:                     # Checking if the file was not opened before
        self.f_traj=open(name_file,'rb')      
           
           
            
        N_A=stc.unpack('i', self.f_traj.read(4))[0]
	self.box[0][:]=stc.unpack('3f', self.f_traj.read(3*4))[0:3]    # Using the global variable (pyn_var_glob.py)
        self.box[1][:]=stc.unpack('3f', self.f_traj.read(3*4))[0:3]    # for the size of the box vg.box        
        self.box[2][:]=stc.unpack('3f', self.f_traj.read(3*4))[0:3]         
        
        if frame!= None: 
            self.frame=frame            # Going to the choosen frame
            bytes_frame=(13+N_A*3)*4       # Bytes per frame
            self.f_traj.seek(bytes_frame*self.frame,1)     # Jumping to the choosen frame
            
        else:
            self.frame+=1
     
            #print self.frame
        #temp_coors.num_atoms=stc.unpack('i', self.f_traj.read(4))[0]            # Reading Attributes
        N_A=stc.unpack('i', self.f_traj.read(4))[0]                              # N_A  = number of atoms

        self.step=stc.unpack('i', self.f_traj.read(4))[0]
        self.time=stc.unpack('f', self.f_traj.read(4))[0]
        self.box[0][:]=stc.unpack('3f', self.f_traj.read(3*4))[0:3]    # Using the global variable (pyn_var_glob.py)
        self.box[1][:]=stc.unpack('3f', self.f_traj.read(3*4))[0:3]    # for the size of the box vg.box        
        self.box[2][:]=stc.unpack('3f', self.f_traj.read(3*4))[0:3]                                                 
        
        format=str(3*N_A)+'f'                                   # Format of system's coordinates
        temp=stc.unpack(format, self.f_traj.read(N_A*4*3))[0:N_A*3]
        
        for ii in range(0,3*N_A,3):
            aux=temp[ii:ii+3]
            
            self.xyz.append(aux)
        self.precision=stc.unpack('f',self.f_traj.read(4))[0]             # Precision of the trajectory
       
        self.f_traj.close()
       # return temp_coors



"""
    def read_bin(self,name_file,frame):                    

        if (self.traj_mark == 0):                     # Checking if the file was not opened before
            self.f_traj=open(name_file,'rb')      
            self.traj_mark=1                          # Marking the file as opened
            self.frame=0
        
        if (frame != None):                           # Going to the choosen frame
            bytes_frame=(13+self.num_atoms*3)*4       # Bytes per frame
            self.f_traj.seek(bytes_frame*self.frame,1)     # Jumping to the choosen frame
            self.frame=frame
        else:
            self.frame+=1
        
        self.num_atoms=stc.unpack('i', self.f_traj.read(4))[0]            # Reading Attributes
        self.step=stc.unpack('i', self.f_traj.read(4))[0]
        self.time=stc.unpack('f', self.f_traj.read(4))[0]
        self.box[0][:]=stc.unpack('3f', self.f_traj.read(3*4))[0:3]    # Using the global variable (pyn_var_glob.py)
        self.box[1][:]=stc.unpack('3f', self.f_traj.read(3*4))[0:3]    # for the size of the box vg.box        
        self.box[2][:]=stc.unpack('3f', self.f_traj.read(3*4))[0:3]                                                 
        self.box=vg.box                                         # Size of the box for system                 
        format=str(self.dim_system)+'f'                                   # Format of system's coordinates
        self.coors=stc.unpack(format, self.f_traj.read(self.num_atoms*4*3))[0:self.num_atoms*3]     
        self.precision=stc.unpack('f',self.f_traj.read(4))[0]             # Precision of the trajectory
"""
      
      

      
