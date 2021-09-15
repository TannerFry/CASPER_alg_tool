from Bio.Seq import Seq

def reverse_complement(seq):
    seq = seq[::-1]
    reverse_comp_seq = ""
    for char in seq:
        if char == "A":
            reverse_comp_seq += "T"
        elif char == "T":
            reverse_comp_seq += "A"
        elif char ==  "C":
            reverse_comp_seq += "G"
        elif char == "G":
            reverse_comp_seq += "C"    
    return reverse_comp_seq

def get_pams(pam):
    if pam == "NGG":
        return  ["AGG", "GGG", "TGG", "CGG"]
    else:
        return []

class file_operations():
    def __init__(self) -> None:
        pass

    def read_CRISPRSCAN(self, CASPER_info_file_path, on_target_matrix):
        #CRISPRSCAN variable
        CRISPRSCAN_data = {}

        #create file object
        with open(CASPER_info_file_path, "r") as f:
            for line in f:
                if line.startswith(on_target_matrix):
                    line = f.readline()
                    while "--" not in line:
                        tmp = line.split()
                        if int(tmp[1]) - 1 in CRISPRSCAN_data.keys():
                            CRISPRSCAN_data[int(tmp[1]) - 1][tmp[0]] = float(tmp[2])
                        else:
                            tmp_dict = dict()
                            tmp_dict[tmp[0]] = float(tmp[2])
                            CRISPRSCAN_data[int(tmp[1]) - 1] = tmp_dict
                        line = f.readline()

        return CRISPRSCAN_data

    def read_fna_file(self, fna_path):
        organism_data = {}
        current_chromosome = ""
        with open(fna_path, "r") as f:
            for line in f:
                if line.find(">") != -1:
                    line = line.replace(">","")
                    line = line.strip("\n")
                    current_chromosome = line
                    organism_data[current_chromosome] = ""
                else:
                    line = line.strip("\n")
                    organism_data[current_chromosome] += line

        return organism_data
    

# on targeting scoring algorithm class
class on_target():
    def __init__(self):
        #pam penalty data
        self.pam_scores = []
        self.pam_scores.append([1, 1, 1, 1, 1.66, 2.5, 5])
        self.pam_scores.append([1, 1, 1, 1, 1.66, 2.5, 5])
        self.pam_scores.append([1, 1, 1, 1.66, 3.33, 3.33, 0])
        self.pam_scores.append([2, 2, 2.5, 2.5, 3.33, 0, 0])
        self.pam_scores.append([2, 2, 2.5, 2.5, 0, 0, 0, 0])
        self.pam_scores.append([8, 8, 10, 0, 0, 0, 0])
        self.pam_scores.append([10, 10, 0, 0, 0, 0, 0])

    def get_sc(self, seq, CRISPRSCAN_data, three_prime):
        """ Calculates the CRISPRscan score for the on-target sequence. Seq should be 35nt in length.
            CRISPRscan assumes index 1 is the position farthest from the PAM. """
        seq = seq.upper()
        score = 0
        """ Tally the CRISPRscan score for the sequence """
        if three_prime:
            for i, nt in enumerate(seq):
                if i == len(seq) - 1:
                    dnt = ""
                dnt = seq[i:i + 2]
                if i in CRISPRSCAN_data.keys():
                    if (nt + "x") in CRISPRSCAN_data[i].keys():
                        score += CRISPRSCAN_data[i][nt + "x"]
                        #print(i,nt+"x",CRISPRSCAN_data[i][nt + "x"])
                    if dnt in CRISPRSCAN_data[i].keys():
                        score += CRISPRSCAN_data[i][dnt]
                        #print(i,dnt,CRISPRSCAN_data[i][dnt])
        else:
            seq = seq[::-1] ### Reverse sequence so PAMs are in same location relative to CRISPRscan data
            for i, nt in enumerate(seq):
                if i == len(seq) - 1:
                    dnt = ""
                dnt = seq[i:i + 2]
                if i in CRISPRSCAN_data.keys():
                    if (nt + "x") in CRISPRSCAN_data[i].keys():
                        score += CRISPRSCAN_data[i][nt + "x"]
                        #print(i,nt+"x",CRISPRSCAN_data[i][nt + "x"])
                    if dnt in CRISPRSCAN_data[i].keys():
                        score += CRISPRSCAN_data[i][dnt]
                        #print(i,dnt,CRISPRSCAN_data[i][dnt])

        # Return the CRISPRscan score
        score += 0.183930944
        return (score * 100)

    def get_sij(self, seq, three_prime, gRNA_len, pam):
        if three_prime:
            seq = seq[6:(6+gRNA_len)] ###This works
        else:
            seq = seq[6+len(pam):6+len(pam)+gRNA_len] ###This works
        seq = Seq(seq.upper())
        rev_seq = seq.reverse_complement()
        pam_list = ["AGG", "TGG", "CGG", "GGG"]
        count = 0
        rev_count = 0
        pam_penalty = 0
        for i in range(len(seq)):
            tmp = seq[i:i + 3]
            rev_tmp = rev_seq[i:i + 3]
            if tmp in pam_list:
                count += 1
            if rev_tmp in pam_list:
                rev_count += 1
        if rev_count <= 6 and count <= 6:
            pam_penalty = self.pam_scores[count][rev_count]

        return pam_penalty

    def get_sg(self, seq, three_prime, gRNA_len, pam):
        if three_prime:
            seq = seq[6:(6+gRNA_len)]
        else:
            seq = seq[6+len(pam):6+len(pam)+gRNA_len] ###This works
        score = 0
        for nt in seq:
            if nt == "G":
                score += 1.0
            elif nt == "C":
                score += 0.5
            elif nt == "A":
                score -= 0.1
        return score / gRNA_len

    def get_p(self, sij_score, sg_score):
        if sij_score > 1:
            p = sij_score * sg_score
        elif sij_score == 1:
            p = 1 - (sg_score / 5)
        else:
            p = 0
        return p

    def ggg_penalty(self, seq):
        ggg_count = 0
        for i in range(0, len(seq)-2):
            if seq[i:i+3] == "GGG":
                ggg_count += 1

        if ggg_count < 2:
            return 1
        elif ggg_count == 2:
            return 0.85
        elif ggg_count == 3:
            return 0.7
        else:
            return 0.5

    #wrapper to score a sequence using the scoring algorithms above
    def score_sequence(self, gRNA, full_seq, CRISPRSCAN_data, endo, three_prime, gRNA_len, pam):
        sc_score = self.get_sc(full_seq, CRISPRSCAN_data,three_prime)
        sij_score = self.get_sij(full_seq,three_prime,gRNA_len,pam)
        sg_score = self.get_sg(full_seq,three_prime,gRNA_len,pam)
        p_score = self.get_p(sij_score, sg_score)

        if endo != "spCas9":
            p_score += self.ggg_penalty(gRNA)

        if p_score == 0:
            score = sc_score
        else:
            score = sc_score / p_score
        return score


# seq finder object
class seq_finder():
    def __init__(self):
        pass

    def find_pams(self, organism_data, pam):
        #hardcoding spCas9 pams for now
        pams = get_pams(pam)

        #pam_locations will store the indexes of the first letter of the pams found for each chromosome section
        pam_locations = {}
        pam_locations_reverse_comp = {}
        
        #formard passthrough
        for chromosome in organism_data.keys():
            pam_locations[chromosome]= []
            for i in range(0, len(organism_data[chromosome])-2):
                if organism_data[chromosome][i:i+3] in pams:
                    pam_locations[chromosome].append(i)
        
        #reverse complement passthrough
        for chromosome in organism_data.keys():
            pam_locations_reverse_comp[chromosome]= []
            reverse_comp_data = reverse_complement(organism_data[chromosome])

            for i in range(0, len(reverse_comp_data)-2):
                if reverse_comp_data[i:i+3] in pams:
                    pam_locations_reverse_comp[chromosome].append(i)

        return pam_locations, pam_locations_reverse_comp
    
    def extract_sequences(self, organism_data, pam_locations, pam_locations_reverse_comp, sequence_length, pam_length):
        sequence_data = {}
        
        for chromosome in organism_data.keys():
            sequence_data[chromosome] = []
            #process forward passthrough sequences
            for pam_location in pam_locations[chromosome]:
                #calculate leftover padding needed
                leftover_padding = 35 - 6 - sequence_length - pam_length
                full_seq = organism_data[chromosome][pam_location - sequence_length - leftover_padding: pam_location + pam_length + 6]
                gRNA = organism_data[chromosome][pam_location - sequence_length: pam_location]
                
                #get rid of entries too close to edge
                if len(gRNA) == 20 and len(full_seq) == 35:
                    sequence_data[chromosome].append([pam_location - 1, gRNA, full_seq])
                
            #process reverse complement passthrough sequences
            reverse_comp_data = reverse_complement(organism_data[chromosome])
            for pam_location in pam_locations_reverse_comp[chromosome]:
                leftover_padding = 35 - 6 - sequence_length - pam_length
                full_seq = reverse_comp_data[pam_location - 6: pam_location + pam_length + sequence_length + leftover_padding]
                gRNA = reverse_comp_data[pam_location + pam_length: pam_location + pam_length + sequence_length]

                #get rid of entries too close to edge
                if len(gRNA) == 20 and len(full_seq) == 35:
                    sequence_data[chromosome].append([-1 * (pam_location + len(pam)), gRNA, full_seq])
                
        return sequence_data


if __name__ == "__main__":
    #input args: endo, pam, multi-threading, directionality, repeats, five length, seed length, three length, org code, output path, CASPERinfo path, fna path, org name, notes, on target matrix
    endo = "spCas9"
    pam = "NGG"
    mt = False
    directionality = True
    repeats = True
    five_length = 4
    seed_length = 16
    three_length = 0
    org_code = "bsu"
    output_path = ""
    CASPERinfo_path = "CASPERinfo"
    fna_path = "staphylococcus_aureus.fna"
    org_name = "staph"
    notes = ""
    on_target_matrix = "DATA:CRISPRSCAN"

    sequence_length = five_length + seed_length + three_length

    #class objects
    file_ops_obj = file_operations()
    on_target_obj = on_target()
    seq_finder_obj = seq_finder()

    #read on-target matrix
    CRISPRSCAN_data = file_ops_obj.read_CRISPRSCAN(CASPERinfo_path, on_target_matrix)
    #print(CRISPRSCAN_data)

    #read fna file
    organism_data = file_ops_obj.read_fna_file(fna_path)

    #analyze fna file
    pam_locations, pam_locations_reverse_comp = seq_finder_obj.find_pams(organism_data, pam)

    #extract sequences from pam locations in each chromosome
    sequence_data = seq_finder_obj.extract_sequences(organism_data, pam_locations, pam_locations_reverse_comp, sequence_length, pam_length=len(pam))

    #write out the sequencing data
    with open("seq_finder_output.csv", "w") as f:
        for chromosome in sequence_data.keys():
            f.write("Location, gRNA, full sequence, on-target score\n")
            for sequence in sorted(sequence_data[chromosome], key=lambda x: abs(x[0])):
                f.write(f"{sequence[0]}, {sequence[1]}, {sequence[2]}, {round(on_target_obj.score_sequence(sequence[1], sequence[2], CRISPRSCAN_data, endo, directionality, sequence_length, pam))}\n")