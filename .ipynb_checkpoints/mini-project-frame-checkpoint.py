#!/bin/python

#importing the porgrams
import os
import sys

pdb_file = None
global file_name
lines = None


### THE MENU FUNCTIONS
def print_menu():
    print(
"""
********************************************************************************
* PDB FILE ANALYZER                                                            *
********************************************************************************
* Select an option from below:                                                 *
*                                                                              *
*      1) Open a PDB File                      (O)                             *
*      2) Information                          (I)                             *
*      3) Show histogram of amino acids        (H)                             *
*      4) Display Secondary Structure          (S)                             *
*      5) Export PDB File                      (X)                             *
*      6) Exit                                 (Q)                             *
*                                                                              *
*                                                       Current PDB: %s        *
********************************************************************************
"""%pdb_file)
    
    option = input(": ")
    #return option  

    # """Takes user's option following main menu display, assesses it and calls the respective funtions"""
    if option.lower() in ('o','i','h','s','x','q'):
        if option.lower() == 'o':
            open_file()
            print_menu()
            
        if option.lower() == 'i':
            print("Information")
            main_info_function(file_name)
            print_menu()
            
        if option.lower() == 'h':
            print("Histogram of Amino Acids")
            main_histogram(file_name)
            print_menu()
            
        if option.lower() == 's':
            print("Display Secondary Structure")
            main_sec_structure(file_name)
            print_menu()
            
        if option.lower() == 'x':
            export_func(file_name)
            print_menu()
            
        if option.lower() == 'q':
            print("Exit")
            main_quit()
           
    else:
        print("\n" + "unsupported option: %s, please enter valid option" %(option))
        print_menu()


#------------------------------------------------------------------------------------------------------------------------
#                                            OPEN SECTION                                                        #
#------------------------------------------------------------------------------------------------------------------------

def path_check(file):
    """ Checks if file exists. if it does, then call function that checks whether file is open """
    if os.path.exists(file): 
        check_open(file)
    else:
        print("file not found")
        print_menu()

def check_open(file):
    """checks if file is open or not. If file is already open, it calls function for replacement enquiry"""
    if pdb_file:
        file_name = file
    else:
        replacement_check(file)

def replacement_check(file):
    """checks if user wants to replace file or not, incase file was already loaded"""
    if True:
        x = input("Would you like to replace the current file, Y/N: ")
        if x.lower() == 'y':
            file_name = file
            print("The File 3AYU.pdb has been successfully loaded")
        else:
            print("proceed to select options to work with", pdb_file)
            print_menu()

def read_pdb(pdb_file):
    with open(file_name,"r") as myfile:
        lines = myfile.readlines()
    return lines    
#------------------------------------------------------------------------------------------------------------------------
#                                            INFORMATION SECTION                                                        #
#------------------------------------------------------------------------------------------------------------------------
def title_print(lines):
    """print title of the protein"""
    title_string = ""
    for line in lines:
        if line.startswith("TITLE"):
            title_string = title_string + line[9:].strip()
    print("Title : " , title_string)
       
def extract_chain_sequences(lines):
    """extract all the sequence residue lines from the file"""
    seq = []
    for line in lines:
        if line.startswith('SEQRES'):
            seq.append(line[0:])
    return seq

def collect_chain_ids(all_sequences):
    """identify chains in protein"""
    chains = []
    for line in all_sequences:
        if line[11] not in chains:
            chains.append(line[11])
    return chains

def print_chains(chains_in_prot):
    """print chains in protein"""
    x = '' .join(chains_in_prot)
    print("- Chains:", x[0], "and", x[1])
    
def pdb_info(all_sequences, chains_in_prot, lines):
    title_print(lines)
    for chain in chains_in_prot:
        print("- Chain ", chain)
        
        residues = []
        for line in all_sequences:
            if line[11] == chain:
                one_letter_code = {'GLY':'G', 'ALA':'A', 'VAL':'V', 'CYS':'C', 'PRO':'P', 'LEU':'L',\
                                   'ILE':'I', 'MET':'M', 'TRP':'W', 'PHE':'F', 'SER':'S', 'THR':'T',\
                                   'TYR':'Y', 'ASN':'N', 'GLN':'Q', 'LYS':'K', 'ARG':'R', 'HIS':'H',\
                                   'ASP':'D', 'GLU':'E'}
                residues.extend(line[18:].split()) #splits the string into a list of residues after appending to the list of residues
                chain_seq = '' .join([one_letter_code[i] for i in residues])#converts the 3 code residues to their corresponding 1 letter denotation
                
                helix = []
                for line in lines:
                    if line.startswith('HELIX') and line[19] == chain:
                        helix.append(line[0:])
                numb = len(helix)
                
                sheet = []
                for line in lines:
                    if line.startswith('SHEET') and line[21] == chain:
                        sheet.append(line[0:])
                num = len(sheet)
                
        print("Number of amino acids:", len(chain_seq))
        print("Number of helix: ", numb)
        print("Number of sheet: ", num)
        print("Sequence:", '\n'.join(''.join(chain_seq[i:i+50]) for i in range(0, len(chain_seq), 50)))

        

#------------------------------------------------------------------------------------------------------------------------
#                                            HISTOGRAM SECTION                                                          #
#------------------------------------------------------------------------------------------------------------------------
def chain_sequence(all_sequences):
    """general print sequences in all chains"""
    residues = []
    for line in all_sequences:
        one_letter_code = {'GLY':'G', 'ALA':'A', 'VAL':'V', 'CYS':'C', 'PRO':'P', 'LEU':'L', 'ILE':'I',\
                           'MET':'M', 'TRP':'W', 'PHE':'F', 'SER':'S', 'THR':'T', 'TYR':'Y', 'ASN':'N', \
                           'GLN':'Q', 'LYS':'K', 'ARG':'R', 'HIS':'H', 'ASP':'D', 'GLU':'E'}
        residues.extend(line[18:].split()) #splits the string into a list of residues after appending to the list of residues
        chain_seq = '' .join([one_letter_code[i] for i in residues]) #converts the 3 code residues to their corresponding 1 letter denotation
    return chain_seq

def ordering_option():
    """prints the ordering options and gives the user an input cell. returns the input"""
    print("""
     Choose an option to order by:
         number of amino acids - ascending (an)
         number of amino acids - descending (dn)
         alphabetically - ascending (aa)
         alphabetically - descending (da)
    """)
    choice = input("order by: ")
    return choice

def an_order(chain_seq):
    """ generates a list of the 20 amino acids ordered in ascending order of their abundance in the sequence"""
    # availing list of all the possible (20) amino acids
    aa_residues = ['G', 'A', 'V', 'C', 'P', 'L', 'I', 'M', 'W', 'F', 'S', 'T', 'Y', 'N', 'Q', 'K', 'R', 'H', 'D', 'E']
    aa_number = [] # initializing an empty list of count for every amino acid in the sequence
    
    # Create a dictionary of each amino acid paired with its count in the sequence
    for residue in aa_residues:
        aa_number.append(chain_seq.count(residue))
    aa_count_dict = dict((residue,aa) for residue,aa in zip(aa_residues, aa_number))
    
     # create a new list of the 20 amino acids ordered according to user option(an)
    residues_list= []
    for k,v in sorted(aa_count_dict.items(), key=lambda p:p[1]):
        residues_list.append(k)
    return residues_list

def dn_order(chain_seq):
    """generates a list of the 20 amino acids ordered in descending order of their abundance in the sequence"""
    # availing list of all the possible (20) amino acids
    aa_residues = ['G', 'A', 'V', 'C', 'P', 'L', 'I', 'M', 'W', 'F', 'S', 'T', 'Y', 'N', 'Q', 'K', 'R', 'H', 'D', 'E']
    aa_number = []  # initializing an empty list of count for every amino acid in the sequence
    
    # Create a dictionary of each amino acid paired with its count in the sequence
    for residue in aa_residues:
        aa_number.append(chain_seq.count(residue))
    aa_count_dict = dict((residue,aa) for residue,aa in zip(aa_residues, aa_number))
    
    # create a new list of the 20 amino acids ordered according to user option(dn)
    residues_list = []
    for k,v in sorted(aa_count_dict.items(), key=lambda p:p[1], reverse=True):
        residues_list.append(k)
        
    return residues_list

def aa_order(chain_seq):
    """generates a list of the 20  amino acids ordered in ascending alphabetical order"""
    residue_list = []
    one_letter_code = {'G':'Gly', 'A':'Ala', 'V':'Val', 'C':'Cys', 'P':'Pro', 'L':'Leu', 'I':'Ile', 'M':'Met', 'W':'Trp', 'F':'Phe', 'S':'Ser', 'T':'Thr', 'Y':'Tyr', 'N':'Asn', 'Q':'Gln', 'K':'Lys', 'R':'Arg', 'H':'His', 'D':'Asp', 'E':'Glu'}
    for k,v in sorted(one_letter_code.items(), key=lambda p:p[1]):
        residue_list.append(k)
    return residue_list


def da_order(chain_seq):
    """generates a list of the 20  amino acids ordered in ascending alphabetical order"""
    residue_list = []
    one_letter_code = {'G':'Gly', 'A':'Ala', 'V':'Val', 'C':'Cys', 'P':'Pro', 'L':'Leu', 'I':'Ile', 'M':'Met', 'W':'Trp', 'F':'Phe', 'S':'Ser', 'T':'Thr', 'Y':'Tyr', 'N':'Asn', 'Q':'Gln', 'K':'Lys', 'R':'Arg', 'H':'His', 'D':'Asp', 'E':'Glu'}
    for k,v in sorted(one_letter_code.items(), key=lambda p:p[1], reverse=True):
        residue_list.append(k)
    return residue_list

def draw_hist(chain_seq, residues):
    """generates a histogram of amino acids in the sequence"""
    one_letter_code = {'G':'Gly', 'A':'Ala', 'V':'Val', 'C':'Cys', 'P':'Pro', 'L':'Leu', 'I':'Ile', 'M':'Met', 'W':'Trp', 'F':'Phe', 'S':'Ser', 'T':'Thr', 'Y':'Tyr', 'N':'Asn', 'Q':'Gln', 'K':'Lys', 'R':'Arg', 'H':'His', 'D':'Asp', 'E':'Glu'}
    #residues = ['G', 'A', 'V', 'C', 'P', 'L', 'I', 'M', 'W', 'F', 'S', 'T', 'Y', 'N', 'Q', 'K', 'R', 'H', 'D', 'E']
    for residue in residues:
        freq = []
        for i in chain_seq:
            if residue == i:
                l = one_letter_code[i] 
                freq.append(l)
        amino_acid = "".join(set(freq))
        if residue in chain_seq:
            print(amino_acid, "(", len(freq),"):", "*" * len(freq))
        else:
            pass

def summary_caller(order, chain_seq):
    """takes theorder option given by the user and calls appropriate function to excecute task"""
    
    if order.lower() in ('an','dn','aa','da'):
        if order.lower() == 'an': 
            residues = an_order(chain_seq)
            draw_hist(chain_seq, residues)
            
        elif order.lower() == 'dn':
            residues = dn_order(chain_seq)
            draw_hist(chain_seq, residues)
            
        elif order.lower() == 'aa':
            residues = aa_order(chain_seq)
            draw_hist(chain_seq, residues)
            
        else:
            residues = da_order(chain_seq)
            draw_hist(chain_seq, residues)
    else:
        print("invalid input")
        ordering_option()       
        

#------------------------------------------------------------------------------------------------------------------------
#                                            SECONDARY STRUCTURE SECTION                                                #
#------------------------------------------------------------------------------------------------------------------------

def sec_structure_generate(chains_in_prot, all_sequences,lines):
    for chain in chains_in_prot: # Adress a chain at a time from my unique chains list
        print('\n' + "Chain",chain,":")
    
        residues = [] # Initiate empty list for my residues per chain
        count = 0
        structure = [] # Initiate an empty list to create the secondary structure symbols as we move over our residue list once generated
        tag =[]
        for line in all_sequences: # Generate primary sequence of the chain
             if line[11] == chain:
                count = count+1
                one_letter_code = {'GLY':'G', 'ALA':'A', 'VAL':'V', 'CYS':'C', 'PRO':'P', 'LEU':'L', 'ILE':'I',\
                    'MET':'M', 'TRP':'W', 'PHE':'F', 'SER':'S', 'THR':'T', 'TYR':'Y', 'ASN':'N', \
                    'GLN':'Q', 'LYS':'K', 'ARG':'R', 'HIS':'H', 'ASP':'D', 'GLU':'E'}
                x = line[18:].split()
                residues.extend(x) #splits the string into a list of residues after appending to the list of residues
        chain_seq = '' .join([one_letter_code[i] for i in residues]) #converts the 3 code residues to their corresponding 1 letter denotation
        for i in residues: #Fill the structure list with dashes ('-') as place holders per residue
            structure.append("-")
            tag.append(" ")
        for line in lines: # Identify where each structure starts and ends using the secondary structure info in the pdb file
        
            if line.startswith('SHEET'): #Process for sheet part of the chain
                new_line = line.split()
                if new_line[5] == chain:
                    start = int(new_line[6])
                    stop = int(new_line[9])
                    num = (stop - start) +1
                    update_structure = num * "|" 
            
                    update_tag = (new_line[1] + new_line[2])
                    tag[start-1:start+1] = update_tag
                    structure[start - 1 : stop] = update_structure  
            
            if line.startswith('HELIX'): #process for helix part of chain
                new_line = line.split()
                if new_line[4] == chain:
                    start = int(new_line[5])
                    stop = int(new_line[8])
                    num = (stop - start) +1
                    update_structure = num * "/" 
                
                    update_tag = (new_line[1])
                    tag[start-1:start+1] = update_tag
                    structure[start - 1 : stop] = update_structure
                   
        for i in range(0, len(chain_seq),80):
        
            print('\n' + ''.join(chain_seq[i:i+80]) + \
                  '\n' + ''.join(structure[i:i+80]) +\
                  '\n' + ''.join(tag[i:i+80]))
            print("(%d)" %(len(chain_seq)))
        
#------------------------------------------------------------------------------------------------------------------------
#                                            EXPORT SECTION                                                             #
#------------------------------------------------------------------------------------------------------------------------

def export_func(file_name):
    out_file_name = file + 'export'
    with open(out_file_name,"w") as myfile:
        lines_to_transfer = open(file_name,"r")
        for line in lines_to_transfer:
            myfile.write(line)
        file.close()

#------------------------------------------------------------------------------------------------------------------------
#                                            EXIT SECTION                                                              #
#------------------------------------------------------------------------------------------------------------------------

def quit():
    """Takes users input of whether to exit program or go back to the main menu"""
    option = input("Do you want to exit (E) or do you want go back to the menu (M)")
    return option

def quit_options(option): 
    """Executes user option"""
    if option == "E" or option == "e":
        sys.exit()
    elif option == "M" or option == "m":
        print_menu()
    else:
        quit()
        
        
#------------------------------------------------------------------------------------------------------------------------
#                                            THE MASTER FUNCTIONS FOR EACH OPTION                                       #
#------------------------------------------------------------------------------------------------------------------------
def open_file():
    """ Asks the user to put in file name and once received, it calls the function that checks whether file path exists"""
    file = input("Enter a Valid PATH for a PDB File: ")
    path_check(file)
    return file

def main_info_function(file_name):
    """The information summary boss function"""
    lines = read_pdb(file)
    all_sequences = extract_chain_sequences(lines)  #extracts the sequence residues lines
    sequence = chain_sequence(all_sequences)  # extracts the sequence of amino acids in the protein
    chains_in_prot = collect_chain_ids(all_sequences)
    pdb_info(all_sequences, chains_in_prot, lines)

def main_histogram(file_name):
    """Main function. Calls the rest of the functions within the histogram option"""
    lines = read_pdb(file)   
    all_sequences = extract_chain_sequences(lines)  #extracts the sequence residues lines
    sequence = chain_sequence(all_sequences)  # extracts the sequence of amino acids in the protein
    order = ordering_option()  # displays order options for the user and records the choice input
    summary_caller(order,sequence)  # sieves through the order options and calls the appropriate functions based on order selected by user

def main_sec_structure(file_name):
    """This is master secondary structure function"""
    lines = read_pdb(file)
    all_sequences = extract_chain_sequences(lines)
    chains_in_prot = collect_chain_ids(all_sequences)
    sec_structure_generate(chains_in_prot, all_sequences,lines)

def main_quit():
    """The master quit function"""
    option = quit()
    quit_options(option)
    
def software_funct():
    print_menu()
software_funct()

