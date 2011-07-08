class cl_unit():
#--------------------------------------------------------------
    def __init__(self):
        '''Initialize an atom object'''


        # > Topological properties
        self.name=''
        self.index=0
        self.pdb_index=0
        self.resid_name=''
        self.resid_pdb_index=0
        self.resid_index=0
        self.chain=''
        self.covalent_bond=[]
        self.alt_loc=0
        self.code_ins_res=0
        self.seg_ident=''
        self.elem_symb=''
        self.type_pdb=''

        # > Physical and Chemical properties
        self.bfactor=0.0
        self.acceptor=False          # True or false
        self.donor=False            # True or false
        self.polar_class='Nothing'       # A or D or nothing
        self.polarizability=False    # True of falsel
        self.mass=0.0
        self.charge=0.0
        self.vdw=0.0
        self.occup=0.0



#        self.coors=[0.0,0.0,0.0]

        # > Topological properties
        #self.index=
        #self.resid=
        #self.molid=




