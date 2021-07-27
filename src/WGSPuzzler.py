import re
from sys import prefix
from numpy.random import MT19937, RandomState, SeedSequence

def list_rc(seq):
    r=[]
    for c in seq[::-1]:
        if c=='A': r.append('T')
        elif c=='C': r.append('G')
        elif c=='G': r.append('C')
        else: r.append('A')
    return r

class WGSPuzzler:
    def __init__(self) -> None:
        self.ref_seqs=[]
        self.operations=[]
        self.planed_sequencing=[]
        self.saved_seed=None

    def load_fasta(self,filename,molecule_list=None) -> None:
        self.operations=[("load",filename,molecule_list)]   
    
    def make_mutated_copies(self,ratio=.01,molecule_list=None):
        self.operations.append(("mutated_copies",ratio,molecule_list))

    def insert_repeats(self,repeat_sequence,count,mut_ratio=0,molecule_list=None):
        pass

    def make_inversions(self,count,min_size=10000,max_size=1000000,molecule_list=None):
        pass

    def pe_sequencing(self,read_size=250,coverage=50,min_fragment=400,max_fragment=700,error_ratio=.02):
        self.operations.append(("pe_sequencing",read_size,coverage,min_fragment,max_fragment,error_ratio))
    
    def write_reference(self,filename):
        self.operations.append(("write_reference",filename))

    def se_sequencing(self,min_read_size=2000,max_read_size=40000,coverage=50,error_ratio=.07):
        self.operations.append(("se_sequencing",min_read_size,max_read_size,coverage,error_ratio))
    
    def select_molecules(self,molecule_list):
        self.operations.append(("select_molecules",molecule_list))

    def execute(self,seed=None,prefix="puzzle"):
        if seed is not None:
            self.saved_seed=seed
        elif self.saved_seed is None:
            self.saved_seed=RandomState(MT19937(None)).randint(100000)
        rangen=RandomState(MT19937(SeedSequence(self.saved_seed)))
        filenumber=1
        for op in self.operations:
            
            if op[0]=="load":
                optxt,filename,molecule_list=op
                self.ref_seqs=[]
                with open(filename) as infile:
                    seq=""
                    for line in infile.readlines():
                        line=line.strip().upper()
                        if line=="": continue
                        if line[0]==">":
                            if seq: self.ref_seqs.append(seq)
                            seq=""
                        else:
                            seq+=re.sub('[^ACGT]','A',line)
                    if seq: self.ref_seqs.append(seq)
                if molecule_list:
                    self.ref_seqs=[m for i,m in enumerate(self.ref_seqs) if i in molecule_list]

            elif op[0]=="mutated_copies":
                optxt,ratio,molecule_list=op
                if molecule_list is None: molecule_list=range(len(self.ref_seqs))
                for mi in molecule_list:
                    new_seq=[x for x in self.ref_seqs[mi]]
                    for si in rangen.randint(0,len(new_seq),int(len(new_seq)*ratio)):
                        new_seq[si]=rangen.choice(['A','C','G','T'])
                    self.ref_seqs.append(''.join(new_seq))

            elif op[0]=="select_molecules":
                self.ref_seqs=[x for i,x in enumerate(self.ref_seqs) if i in op[1]]
                
            elif op[0]=="pe_sequencing":
                optxt,read_size,coverage,min_fragment,max_fragment,error_ratio=op
                quals="#"*read_size
                r1f=open(f"{prefix}_{filenumber:02}_pe_r1.fastq","w")
                r2f=open(f"{prefix}_{filenumber:02}_pe_r2.fastq","w")
                mol_lens=[len(x) for x in self.ref_seqs]
                gen_size=sum(mol_lens)
                numpairs=int(coverage*gen_size/(2*read_size))
                print(f"simulating {numpairs} pairs")
                errors_per_pair=2*read_size*error_ratio*4/3 # multiply by 4/3 because random base can be the right base
                for rid in range(1,numpairs+1):
                    #choose where to start r1 and r2, and get their base sequence
                    r1start=rangen.randint(0,gen_size)
                    rmol=0
                    while r1start>mol_lens[rmol]:
                        r1start-=mol_lens[rmol]
                        rmol+=1
                    if mol_lens[rmol]-r1start<min_fragment:
                        r1start=0
                    if mol_lens[rmol]-r1start<max_fragment:
                        r2end=mol_lens[rmol]-read_size
                    else:
                        r2end=r1start+rangen.randint(min_fragment,max_fragment+1)-read_size
                    r1seq=[x for x in self.ref_seqs[rmol][r1start:r1start+read_size]]
                    r2seq=list_rc([x for x in self.ref_seqs[rmol][r2end:r2end+read_size]])
                    #add errors
                    for error_i in range(rangen.randint(0,errors_per_pair+1)):
                        r1seq[rangen.randint(read_size)]=rangen.choice(['A','C','G','T'])
                    for error_i in range(rangen.randint(0,errors_per_pair+1)):
                        r2seq[rangen.randint(read_size)]=rangen.choice(['A','C','G','T'])
                    
                    #half of the time, flip r1 and r2 and their directions
                    if rangen.random()>=.5:
                        r1seq,r2seq=list_rc(r2seq),list_rc(r1seq)
                    #write to disk
                    r1f.write(f"@read_{rid:08}/1\n{''.join(r1seq)}\n+\n{quals}\n")
                    r2f.write(f"@read_{rid:08}/1\n{''.join(r2seq)}\n+\n{quals}\n")
                r1f.close()
                r2f.close()
                filenumber+=1
            elif op[0]=="se_sequencing":
                optxt,min_fragment,max_fragment,coverage,error_ratio=op
                rf=open(f"{prefix}_{filenumber:02}_se.fastq","w")
                mol_lens=[len(x) for x in self.ref_seqs]
                gen_size=sum(mol_lens)
                numreads=int(2*coverage*gen_size/(min_fragment+max_fragment))
                print(f"simulating {numreads} single reads")
                for rid in range(1,numreads+1):
                    #choose where to start the read
                    rstart=rangen.randint(0,gen_size)
                    rmol=0
                    while rstart>mol_lens[rmol]:
                        rstart-=mol_lens[rmol]
                        rmol+=1
                    while mol_lens[rmol]-rstart<min_fragment:
                        rstart-=rangen.randint(min_fragment)
                    rend=rstart+rangen.randint(min_fragment,max_fragment+1)
                    
                    rseq=[x for x in self.ref_seqs[rmol][rstart:rend]]
                    #add errors (so far, this is completely random)
                    for error_i in range(rangen.randint(0,2*error_ratio*len(rseq)*4/3)):
                        rseq[rangen.randint(len(rseq))]=rangen.choice(['A','C','G','T'])
                    
                    #half of the time, flip read directions
                    if rangen.random()>=.5:
                        list_rc(rseq)
                    #write to disk
                    rf.write(f"@read_{rid:08}/1\n{''.join(rseq)}\n+\n{'#'*len(rseq)}\n")
                rf.close()
                filenumber+=1
            elif op[0]=="write_reference":
                optxt,filename=op
                print("writing reference to",filename)
                with open(filename,"w") as reff:
                    for rsi, rs in enumerate(self.ref_seqs,start=1):
                        reff.write(f">refseq_{rsi:04}\n{rs}\n")

            else:
                print(f"Operation {op[0]} not implemented!")
                raise NotImplementedError
    
    def write_recipe(self,filename):
        """just the operations, all random values will change once executed"""
        pass

    def load_recipe(self,filename):
        """just the operations, all random values will change once executed"""
        pass

    def write_blueprint(self,filename):
        """all random values will be kept, including position of sequences and their errors"""
        pass

    def load_blueprint(self,filename):
        """all random values will be kept, including position of sequences and their errors"""
        pass

