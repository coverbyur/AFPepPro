import numpy as np
import time
from tqdm import tqdm
from numba import jit, njit, prange
from numba_progress import ProgressBar
from Bio import SeqIO

class PeptideContainer:
    #This is a container object for everything we need later
    def __init__(self, input_file):
        self.n_pos, self.n_neg, self.n_neu, self.p_pos, self.p_neg, \
        self.p_neu, self.pos_aa, self.pos_aa_p, self.neg_aa, \
        self.neg_aa_p, self.neu_aa, self.neu_aa_p, self.mer_len, \
        self.n_chsq, self.n_subsq = read_param(input_file)
        self.charge_seq = []
        self.pos_aa_idx = [] 
        self.neg_aa_idx = [] 
        self.neu_aa_idx = []
        self.term_aa = []
        self.term_aa_idx = []

    def add_terminal_AA(self, term_aa):
        '''Add the terminal amino acid, usually resin-defined'''
        self.term_aa = term_aa

    def convert_map(self, aa_map):
        '''Convert the index map of amino acid lookups '''
        if len(self.pos_aa_idx) + len(self.neg_aa_idx) + len(self.neu_aa_idx) == 0:        
            for idx in range(len(self.pos_aa)):
                self.pos_aa_idx.append(aa_map[self.pos_aa[idx]])
            for idx in range(len(self.neg_aa)):
                self.neg_aa_idx.append(aa_map[self.neg_aa[idx]])
            for idx in range(len(self.neu_aa)):
                self.neu_aa_idx.append(aa_map[self.neu_aa[idx]])
            self.term_aa_idx = (aa_map[self.term_aa])
        else:
            raise ValueError('Internal maps already contain values!')

    def prep_chsq(self):
        '''Produce arg list for make_pep_from_chsq for acceleration purposes'''
        return (self.mer_len, self.n_pos, self.pos_aa_p, self.pos_aa_idx, \
        self.n_neg, self.neg_aa_p, self.neg_aa_idx,self. n_neu, self.neu_aa_p, \
        self.neu_aa_idx, self.term_aa_idx)

def read_param(input_file):
    '''Read the parameters for the coming simulation in the following format:
    Order: Number of Positively charged AAs, Number of Negatively charged AAs, Number of Neutral AAs,
    Probability of positively/negatively/neutrally charged AAs,
    Positively charged AAs | probabilities |Negatively charged AAs | probabilities | Neutral AAs | probabilities '''
    #Also read the interaction energy files
    n_pos = []
    n_neg = []
    n_neu = []
    p_pos = []
    p_neg = []
    p_neu = []
    pos_aa = []
    pos_aa_p = []
    neg_aa = []
    neg_aa_p = []
    neu_aa = []
    neu_aa_p = []
    mer_len = []
    n_chsq = []
    n_subsq = []

    index = 0 #Just to simplify the code a bit
    
    data = np.genfromtxt(input_file, dtype=str, delimiter=",",autostrip=True, comments='#')
    n_pos.append(int(data[0]))
    n_neg.append(int(data[1]))
    n_neu.append(int(data[2]))
    p_pos.append(float(data[3]))
    p_neg.append(float(data[4]))
    p_neu.append(float(data[5]))
    #This looks strange at the moment b/c it will eventually work on multiple rows
    #Positive parameters
    index = 6+n_pos[0]
    pos_aa.extend(data[6:index])
    index = index
    pos_aa_p.extend(map(float, data[index:index+n_pos[0]]))
    #Negative parameters
    index = index + n_pos[0]
    neg_aa.extend(data[index:index+n_neg[0]])
    index = index + n_neg[0]
    neg_aa_p.extend(map(float,data[index:index+n_neg[0]]))
    #Neutral parameters
    index = index + n_neg[0]
    neu_aa.extend(data[index:index+n_neu[0]])
    index = index + n_neu[0]
    neu_aa_p.extend(map(float,data[index:index+n_neu[0]]))
    index = index+n_neu[0]
    mer_len.append(int(data[index]))
    index = index + 1
    n_chsq.append(int(data[index]))
    index = index + 1
    n_subsq.append(int(data[index]))

    return n_pos, n_neg, n_neu, p_pos, p_neg, p_neu, pos_aa, \
    pos_aa_p, neg_aa, neg_aa_p, neu_aa, neu_aa_p, mer_len, n_chsq, n_subsq

def read_aa_tables_pasta(input_file_parallel, input_file_antiparallel, \
    input_file_AA_map):
    '''Read the PASTA-style AA interaction data tables and organization scheme from files and 
    generate 2d arrays along with a AA-to-index dictionary'''
    par_table = np.genfromtxt(input_file_parallel, dtype=float, delimiter=',', autostrip=True, comments='#')
    antipar_table = np.genfromtxt(input_file_antiparallel, dtype=float, delimiter=',', autostrip=True, comments='#')
    aa_map = {}
    with open(input_file_AA_map) as f:
        for line in f:
            (key, val) = line.split()
            aa_map[str(key)] = int(val)
    return par_table, antipar_table, aa_map

def convert_map_pro(pro_seq, aa_map):
    '''Convert the index map of amino acid lookups '''
    pro_seq_mapped = np.zeros([len(pro_seq),1], dtype=int)
    for idx, aa in enumerate(pro_seq):
        pro_seq_mapped[idx] = int(aa_map[pro_seq[idx]])
    return pro_seq_mapped

def read_protein_seq(pro_file): #Read the FASTA file containing protein sequences
    fasta_sequences = list(SeqIO.parse(open(pro_file), 'fasta'))
    pro_names = {}
    pro_seq = {}
    for i in range(len(fasta_sequences)):
            record = fasta_sequences[i]
            #print(record.id)
            pro_names[i], pro_seq[i] = record.id, str(record.seq)
    return pro_names, pro_seq

def gen_charge_seq(pcontainer):
    '''Generate a charge sequence in the container with the given input parameters'''
    n_ml = pcontainer.mer_len[0]
    r_list = np.random.rand(n_ml)
    pcontainer.charge_seq = [] #For list comprehension assignment later
    for i in range(n_ml):
        if (r_list[i] < pcontainer.p_pos[0]):
            pcontainer.charge_seq.append('+')
            continue
        if (r_list[i] < (pcontainer.p_neg[0] + pcontainer.p_pos[0])):
            pcontainer.charge_seq.append('-')
            continue
        else:
            pcontainer.charge_seq.append('0')
    return

def gen_charge_seq_list(pcontainer, n_chsq):
    '''Generate a n_chsq charge sequences from the container parameters'''
    n_ml = pcontainer.mer_len[0]
    r_list = []
    seq_list = []
    charge_tol = 1 #Hard coded for now
    j = 0
    while j < n_chsq:
        net_charge = 0
        temp_seq = []
        r_list = np.random.rand(n_ml)
        pcontainer.charge_seq = [] #For list comprehension assignment later
        for i in range(n_ml):
            if (r_list[i] < pcontainer.p_pos[0]):
                temp_seq.append('+')
                net_charge += 1
                continue
            if (r_list[i] < (pcontainer.p_neg[0] + pcontainer.p_pos[0])):
                temp_seq.append('-')
                net_charge += -1
                continue
            else:
                temp_seq.append('0')
        if np.abs(net_charge) < charge_tol:
            seq_list.append(temp_seq)
            j += 1
    return seq_list

def gen_charge_seq_list_fixed(pcontainer, n_chsq):
    '''Generate a n_chsq charge sequences from the container parameters'''
    n_ml = pcontainer.mer_len[0]
    r_list = []
    seq_list = []
    charge_tol = 1 #Hard coded for now
    j = 0
    while j < n_chsq:
        net_charge = 0
        temp_seq = []
        r_list = np.random.rand(n_ml)
        aan_pos = int(n_ml*pcontainer.p_pos[0])
        aan_neg = int(n_ml*pcontainer.p_neg[0])
        aan_neu = int(n_ml*pcontainer.p_neu[0])
        if n_ml > (aan_pos + aan_neg + aan_neu):
            aan_rand = n_ml - (aan_pos + aan_neg + aan_neu)
        else:
            aan_rand = 0
        pcontainer.charge_seq = [] #For list comprehension assignment later
        for i in range(aan_rand):
            if (r_list[i] < pcontainer.p_pos[0]):
                temp_seq.append('+')
                net_charge += 1
                continue
            if (r_list[i] < (pcontainer.p_neg[0] + pcontainer.p_pos[0])):
                temp_seq.append('-')
                net_charge += -1
                continue
            else:
                temp_seq.append('0')
        aa_list = []
        for i in range(aan_pos):
            aa_list.append('+')
        for i in range(aan_neg):
            aa_list.append('-')
        for i in range(aan_neu):
            aa_list.append('0')
        #print(aa_list)
        temp_seq = temp_seq + aa_list
        np.random.shuffle(temp_seq)
        #print(temp_seq)
        if np.abs(net_charge) < charge_tol:
            seq_list.append(temp_seq)
            j += 1
    return seq_list

def compare_sub_seq(pcontainer, e_table_p, e_table_ap):
    '''Generate discrete subsequences from charge sequences; all parameters from input file '''
    n_chsq = pcontainer.n_chsq[0]
    ch_sq_list = gen_charge_seq_list_fixed(pcontainer, n_chsq)
    n_subsq = pcontainer.n_subsq[0] #Is this cheaper than just asking the object?
    p_args = pcontainer.prep_chsq()
    return comp_seq(ch_sq_list, n_subsq, p_args, e_table_p, e_table_ap)

def compare_sub_seq_pro(pcontainer, e_table_p, e_table_ap, pro_mapped_list):
    '''Generate discrete subsequences from charge sequences; all parameters from input file '''
    n_chsq = pcontainer.n_chsq[0]
    ch_sq_list = gen_charge_seq_list_fixed(pcontainer, n_chsq)
    n_subsq = pcontainer.n_subsq[0] #Is this cheaper than just asking the object?
    p_args = pcontainer.prep_chsq()
    return comp_seq(ch_sq_list, n_subsq, p_args, e_table_p, e_table_ap)
    
@jit()
def comp_seq(ch_sq_list, n_subsq, p_args, e_table_p, e_table_ap):
    '''Broken out from compare_sub_seq to allow numba JIT acceleration'''
    e_list = np.ones(n_subsq)*np.inf #Initialize containers
    e_mean = np.ones(len(ch_sq_list))*np.inf
    e_std = np.ones(len(ch_sq_list))*np.inf
    #Parameters hardcoded for now
    delta_s = -.2
    min_len = 6
    for ch_sq in tqdm(range(len(ch_sq_list))):
        for ch_idx in range(n_subsq):
            #The majority of computational work is done in these 3 lines
            pep1 = make_pep_from_chsq(*p_args, ch_sq_list[ch_sq])
            pep2 = make_pep_from_chsq(*p_args, ch_sq_list[ch_sq])
            e_list[ch_idx] = pasta(pep1, pep2, e_table_p, e_table_ap, min_len, delta_s)
            #print('e output:', e_list[ch_idx])
        e_mean[ch_sq] = np.mean(e_list)
        e_std[ch_sq] = np.std(e_list)
    return ch_sq_list, e_mean, e_std

@jit()
def make_pep_from_chsq(mer_len, n_pos, pos_aa_p, pos_aa_idx, \
    n_neg, neg_aa_p, neg_aa_idx, n_neu, neu_aa_p, neu_aa_idx, term_aa_idx, chsq):
    '''Make a "real" peptide from a charge sequence and its associated data
    Except I'm using the map index instead of the actual AA sequence'''
    pep = []
    for i in range(mer_len[0]):
        r_int = np.random.rand()
        if chsq[i] == '+':
            for j in range(n_pos[0]):
                if r_int <= pos_aa_p[j]:
                    pep.append(pos_aa_idx[j])
                    break
                else:
                    r_int = r_int - pos_aa_p[j]
        elif chsq[i] == '-':
            for j in range(n_neg[0]):
                if r_int <= neg_aa_p[j]:
                    pep.append(neg_aa_idx[j])
                    break
                else:
                    r_int = r_int - neg_aa_p[j]
        elif chsq[i] == '0':
            for j in range(n_neu[0]):
                if r_int <= neu_aa_p[j]:
                    pep.append(neu_aa_idx[j])
                    break
                else:
                    r_int = r_int - neu_aa_p[j]
        else:
            print(pep)
            print(chsq[i])
            raise ValueError('Unrecognized charge sequence!')
    pep.append(term_aa_idx) #Add the resin-bound AA
    return pep

@jit(nopython=True, parallel=True)
def pasta(pep1, pep2, e_table_p, e_table_ap, min_len, delta_s):
    '''Perform the pasta algorithm on the provided peptides'''
    p_len = len(pep1)
    min_score = np.inf
    pep2r = pep2[::-1] #Make a reverse copy of pep2
    for l_cur in range(min_len,p_len):
        for k_idx in range(p_len - l_cur):
            for v_idx in range(p_len-l_cur):
                sum_score_p = 0
                sum_score_ap = 0
                for p_idx in range(l_cur):
                    sum_score_p += e_table_p[pep1[k_idx+p_idx], pep2[p_idx+v_idx]]
                    sum_score_ap += e_table_ap[pep1[k_idx+p_idx], pep2r[p_idx+v_idx]]
                    min_score_temp = np.minimum(sum_score_p - l_cur*delta_s, \
                    sum_score_ap- l_cur*delta_s)
                min_score = np.minimum(min_score_temp, min_score)
    return min_score

#@jit(parallel=True) Numba update broke this module
def pasta_pro(pep, pro, e_table_p, e_table_ap, min_len, delta_s):
    '''Perform the pasta algorithm on the provided peptide and protein, assuming protein is longer'''
    p_len = len(pep)
    pro_len = len(pro)
    pro_run_len = pro_len - p_len + min_len
    pepr = pep[::-1] #Make a reverse copy of peptide
    #For this next step: calculate the interaction of each AA in peptide with all protein AAs
    pepEF = np.full((p_len, pro_run_len), np.inf) #Forward interaction table
    pepER = np.full((p_len, pro_run_len), np.inf) #Reversed interaction table

    for j in range(pro_run_len):
        for i in range(p_len): #This is the "no memory" unoptomized version that does many redundant lookups
            idx1 = int(i)
            idx2 = int(j)
            if (j+i > (pro_len - 1)):
                continue
            ipep = int(pep[i][0]) #Numba wants this or it complains. The speed is worth it.
            ipro = int(pro[j+i][0])
            ipepr = int(pepr[i][0])
            pepEF[idx1, idx2] = e_table_p[ipep, ipro] - delta_s
            pepER[idx1, idx2] = e_table_ap[ipepr, ipro] - delta_s

    #Now to do my version of Kadane's algorithm
    min_sum = sum(pepEF[0:min_len,0]) #Set minimum to very first possible value
    fstart = 0
    fend = 0
    #Iterate for each position along protein
    for j in range(pro_run_len):
        min_ending_here = np.sum(pepEF[0:min_len,j])
        k = 0
        for i in range(p_len):
            cur_val = pepEF[i, j]
            if cur_val > min_ending_here:
                k = 0 #Reset tracking position
            min_ending_here = min_ending_here + cur_val
            k = k + 1
            if k >= min_len: #Only consider if above the minimum length
                min_sum = min(min_sum, min_ending_here)
            #Once this finishes, carry min_sum through all possible arrangements

    for j in range(pro_run_len):
        min_ending_here = np.sum(pepER[0:min_len,j])
        k = 0
        for i in range(p_len):
            cur_val = pepER[i, j]
            if cur_val > min_ending_here:
                k = 0 #Reset tracking position
            min_ending_here = min_ending_here + cur_val
            k = k + 1
            if k >= min_len: #Only consider if above the minimum length
                min_sum = min(min_sum, min_ending_here)

    return min_sum

@jit(parallel=True)
def run_with_proteins(pep_list, pro_list, e_table_p, e_table_ap, min_len, delta_s): #Supply peptides and FASTA style proteins, operate pasta_pro, return interaction energies

    num_pep = len(pep_list)
    num_pro = len(pro_list)

    energy_list = np.zeros((num_pep, num_pro))
    #Suddenly tqdm seems to be causing problems, temporarily removed
    for i in tqdm(range(num_pep)):
        for j in range(num_pro):
            energy_list[i,j] = pasta_pro(pep_list[i], pro_list[j], e_table_p, e_table_ap, min_len, delta_s)
    return energy_list

def main_handler_pro(pepfile, p_tablef, ap_tablef, aa_mapf, pro_file, outfile):
    '''Operate the other functions in AFPep after initial file location designations'''
    #Read input files and get relevent parameters
    min_len = 6
    delta_s = -0.2

    p_table, ap_table, aa_map = read_aa_tables_pasta(p_tablef, ap_tablef, aa_mapf)
    pep_names, pep_seqs = read_protein_seq(pepfile)
    pro_names, pro_seqs = read_protein_seq(pro_file)

    pep_list_mapped = {}
    for i in range(len(pep_seqs)):
        pep_list_mapped[i] = convert_map_pro(pep_seqs[i], aa_map)

    pro_list_mapped = {} #np.ones(len(pro_seqs), dtype='int')
    for i in range(len(pro_seqs)):
        pro_list_mapped[i] = convert_map_pro(pro_seqs[i], aa_map)

    p_table, ap_table, aa_map = read_aa_tables_pasta(p_tablef, ap_tablef, aa_mapf)

    energy_list = run_with_proteins(pep_list_mapped, pro_list_mapped, p_table, ap_table, min_len, delta_s)
    outdata = energy_list
    headtxt = pro_names
    #headtxt = ''.join(pro_names) #+ "\nCharge Sequence\t Mean Interaction Energy (Pasta Units)\t Std Dev"
    np.savetxt(outfile, outdata, delimiter='\t', fmt = '%3.3f')

def main_handler(pepfile, p_tablef, ap_tablef, aa_mapf, outfile):
    '''Operate the other functions in AFPep after initial file location designations'''
    #Read input files and get relevant parameters
    pcont = PeptideContainer(pepfile)
    pcont.add_terminal_AA('G')
    p_table, ap_table, aa_map = read_aa_tables_pasta(p_tablef, ap_tablef, aa_mapf)
    pcont.convert_map(aa_map)
    ch_sq_list, e_mean, e_std = compare_sub_seq(pcont, p_table, ap_table)
    outdata = np.array([ch_sq_list,e_mean,e_std]) #Bunch them up
    fmt = '%s\t %3.3f\t %3.3f'
    input_file = open(pepfile)
    input_data = input_file.readlines()
    headtxt = ''.join(input_data) + "\nCharge Sequence\t Mean Interaction Energy (Pasta Units)\t Std Dev"
    np.savetxt(outfile,np.transpose(outdata), fmt, header=headtxt)

def fast_pep_pro(pep_seq, pro_seq, p_tablef, ap_tablef, aa_mapf):
    '''Operate the other functions in AFPep after initial file location designations'''
    #Read input files and get relevent parameters
    p_table, ap_table, aa_map = read_aa_tables_pasta(p_tablef, ap_tablef, aa_mapf)
    pep = convert_map_pro(pep_seq, aa_map)
    pro = convert_map_pro(pro_seq, aa_map)
    energy = pasta_pro(pep, pro, p_table, ap_table, 6, -0.2)
    print(energy)

def start(): #pragma: no cover
    '''Give the user a prompt to enter in the input file location and desired output file name'''
    p_tablef = 'AFPeptide/data/ParAAT.txt'
    ap_tablef = 'AFPeptide/data/AParAAT.txt'
    aa_mapf = 'AFPeptide/data/AAMap.txt'
    print("This program will generate semi-randomized zwitterionic peptides based on an input")
    pep_q1 = input("Enter Z for ZIP mode and P for protein mode:")
    if pep_q1 == "P":
        protein_flag = 0
    else:
        protein_flag = 1
    if protein_flag == 1:
        pep_q2 = input("Do you want to enter a discrete sequence for analysis?")
        if pep_q2 == 'y':
            pep_seq = input("Enter peptide sequence now: ")
            fast_pep(pep_seq, p_tablef, ap_tablef, aa_mapf)
            exit()
        
        parafile = input("Enter the name of the file containing parameter information (format in readme): ")
        outfile = input("Enter the name of the output file to save output to: ")
        outfile = "output/"+outfile
        t0 = time.time()
        main_handler(parafile, p_tablef, ap_tablef, aa_mapf, outfile)
        t1 = time.time()
        print("Execution time was: ", t1-t0)
    else:
        pep_q = input("Do you want to enter a discrete sequence for analysis?")
        if pep_q == 'y':
            pep_seq = input("Enter peptide sequence now: ")
            pro_seq = input("Enter protein sequence now: ")
            fast_pep_pro(pep_seq, pro_seq, p_tablef, ap_tablef, aa_mapf)
        pepfile = input("Enter the name of the file containing peptide sequences in FASTA format: ")
        protfile = input("Enter name of file containing proteins in FASTA format: ")
        outfile = input("Enter the name of the output file to save output to: ")
        outfile = "output/"+outfile
        t0 = time.time()
        main_handler_pro(pepfile, p_tablef, ap_tablef, aa_mapf, protfile, outfile)
        t1 = time.time()
        print("Execution time was: ", t1-t0)
