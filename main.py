from Bio.Seq import Seq
import matplotlib.pyplot as plt
import math


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
    elif pam == "TTTV":
        return ["TTTA", "TTTC","TTTG"] 
    else:
        return []


class file_operations():
    def __init__(self):
        pass

    def read_CRISPRSCAN(self, CASPER_info_file_path):
        #CRISPRSCAN variable
        CRISPRSCAN_data = {}

        #create file object
        with open(CASPER_info_file_path, "r") as f:
            for line in f:
                if line.startswith("DATA:CRISPRSCAN"):
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

    def read_HSU_matrix(self, CASPER_info_file_path):
        #create HSU matrix object
        HSU_matrix_data = []

        #create file object
        with open(CASPER_info_file_path) as fin:
            #loop line by line looking for HSU matrix
            for line in fin:
                if line.find("MATRIX:HSU MATRIX-spCas9-2013") != -1:
                    for line in fin:
                        #check if end of HSU matrix
                        if line == "\n":
                            break
                        else:
                            #parse and typecast the HSU data to floats and store
                            line = line.split()
                            line = [float(val) for val in line]
                            HSU_matrix_data.append(line)
                    break

                # if line.find("MATRIX:HSU MATRIX-asCas12-2016") != -1:
                #     for line in fin:
                #         #check if end of HSU matrix
                #         if line == "\n":
                #             break
                #         else:
                #             #parse and typecast the HSU data to floats and store
                #             line = line.split()
                #             line = [float(val) for val in line]
                #             HSU_matrix_data.append(line)
                #     break

        HSU_dict = {}
        HSU_dict["GT"] = HSU_matrix_data[0]
        HSU_dict["AC"] = HSU_matrix_data[1]
        HSU_dict["GG"] = HSU_matrix_data[2]
        HSU_dict["TG"] = HSU_matrix_data[3]
        HSU_dict["TT"] = HSU_matrix_data[4]
        HSU_dict["CA"] = HSU_matrix_data[5]
        HSU_dict["CT"] = HSU_matrix_data[6]
        HSU_dict["GA"] = HSU_matrix_data[7]
        HSU_dict["AA"] = HSU_matrix_data[8]
        HSU_dict["AG"] = HSU_matrix_data[9]
        HSU_dict["TC"] = HSU_matrix_data[10]
        HSU_dict["CC"] = HSU_matrix_data[11]

        return HSU_dict

    def read_off_target_sequences(self, input_sequences_file_path):
        #create input sequences object
        query_seqs = []
        ref_seqs = []

        #create file object
        with open(input_sequences_file_path) as fin:
            #loop through each line to parse sequences
            temp = []
            for line in fin:
                line = line.strip("\n")
                #check for end of section
                if line.find(">") != -1:
                    for line in fin:
                        line = line.strip("\n")
                        if line.find(">") != -1:
                            temp = []
                            break
                        else:
                            line = line.split(",")
                            ref_seqs.append(temp)
                            query_seqs.append(line)
                else:
                    line = line.split(",")
                    temp.append(line)

        return query_seqs, ref_seqs

    def read_on_target_sequences(self, input_on_target_sequences_file_path):
        on_target_sequences = []
        with open(input_on_target_sequences_file_path, "r") as fin:
            for line in fin:
                line = line.strip("\n")
                on_target_sequences.append(line)
        return on_target_sequences

    def write_OT_results(self, query_seqs, ref_seqs, OT_scores, sh_scores,casper_scores):
        with open("OT_results.txt", "w") as fout:
            for i in range(len(query_seqs)):
                query_seq = query_seqs[i][0]
                # fout.write(f"{OT_scores[i]}\n\n")
                for j in range(len(ref_seqs[i])):
                    ref_seq = ref_seqs[i][j][0]
                    fout.write(f"{sh_scores[i][j]},{casper_scores[i][j]}\n")
                fout.write("\n")


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
        if score <=0:
            return 0
        else:
            return score
        # return (((score + 0.65546) / 1.94947)) 

    def get_sij(self, seq, three_prime, gRNA_len, pam):
        seq = Seq(seq.upper())
        rev_seq = seq.reverse_complement()
        pam_list = get_pams(pam)
        count = 0
        rev_count = 0
        pam_penalty = 0
        for i in range(len(seq)):
            tmp = seq[i:i + len(pam)]
            rev_tmp = rev_seq[i:i + len(pam)]
            if tmp in pam_list:
                count += 1
            if rev_tmp in pam_list:
                rev_count += 1
        if rev_count <= 6 and count <= 6:
            pam_penalty = self.pam_scores[count][rev_count]

        return pam_penalty

    def get_sg(self, seq, three_prime, gRNA_len, pam):
        score = 0
        for nt in seq:
            if nt == "G":
                score += 1.0
            elif nt == "C":
                score += 0.5
            elif nt == "A":
                score -= 0.1
        return (score / gRNA_len)
        # return (score + 0.1*gRNA_len) / gRNA_len

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
        sij_score = self.get_sij(gRNA,three_prime,gRNA_len,pam)
        sg_score = self.get_sg(gRNA,three_prime,gRNA_len,pam)
        p_score = self.get_p(sij_score, sg_score)

        # if gRNA == "TTTAATTGATCTAAGTTTGC":
        #     print(f"sc score: {sc_score}")
        #     print(f"p score: {p_score}")
        #     print(f"sg score: {sg_score}")
        #     print(f"sij score: {sij_score}")

        if endo != "spCas9" and three_prime:
            p_score += self.ggg_penalty(gRNA)

        if p_score == 0:
            score = sc_score
        else:
            score = sc_score / p_score

        if three_prime:
            score /= 1 #Scalar modifier to get scores between 0 and 100
            if (score * 100) <= 100:
                return (score * 100)
            else:
                return 100
        else:
            score /= 1 #Scalar modifier to get scores between 0 and 100
            if (score * 100) <= 100:
                return (score * 100)
            else:
                return 100


class off_target():
    def __init__(self):
        pass

    def sh_score(self, mismatches, hsu_keys, hsu_matrix,gRNA_len):
        sh_score = 1.0
        for i in range(len(mismatches)):
            sh_score *= hsu_matrix[hsu_keys[i]][gRNA_len-mismatches[i]]
            # print(mismatches[i])
            # print(hsu_matrix[hsu_keys[i]][gRNA_len-mismatches[i]])
        return sh_score

    def ss_score(self, mismatches,gRNA_len):
        ss_score = 1.0
        # print(mismatches)
        for i in range(len(mismatches)):
            if gRNA_len == 24: 
                if mismatches[i] <= 8:
                    ss_score -= 0.1
                elif mismatches[i] <= 20:
                    ss_score -= 0.0125
                else:
                    ss_score -= 0
            else: 
                if mismatches[i] <= 6:
                    ss_score -= 0.1
                elif mismatches[i] <= 12:
                    ss_score -= 0.05
                else:
                    ss_score -= 0.0125
        return ss_score

    def st_score(self, mismatches):
        st_score = 3.547
        for mismatch in mismatches:
            st_score -= (1.0 / mismatch)
        return st_score / 3.5477

    def reverse_comp(self, c):
        if c == "A":
            return "T"
        elif c == "T":
            return "A"
        elif c == "G":
            return "C"
        elif c == "C":
            return "G"
        else:
            return "N"

    def get_mismatches(self, ref_seq, query_seq, three_prime, gRNA_len):
        """ Returns location and type of mismatch, relative to bp from PAM. Index next to PAM is set as 1. """
        mismatch_locations = []
        mismatch_keys = []
        # print(ref_seq,query_seq)

        for i in range(gRNA_len-1, -1, -1):
            # print(i)
            # print(ref_seq[i],query_seq[i])
            if ref_seq[i] != query_seq[i]:
                if three_prime:
                    mismatch_locations.append(gRNA_len-i)
                    # mismatch_keys.append(query_seq[i] + self.reverse_comp(ref_seq[i]))
                    mismatch_keys.append(ref_seq[i] + self.reverse_comp(query_seq[i]))
                else:
                    mismatch_locations.append(i+1) ### Account for indexing convention
                    # mismatch_keys.append(query_seq[i] + self.reverse_comp(ref_seq[i]))
                    mismatch_keys.append(ref_seq[i] + self.reverse_comp(query_seq[i]))

        return mismatch_locations, mismatch_keys

    def find_similars(self, query_seq, ref_seqs, seqs_on_target_score, HSU_matrix,endo,three_prime,gRNA_len):
        target_scores = []
        sh_scores = []
        for i in range(len(ref_seqs)):
            #print(f"ref on-score: {seqs_on_target_score[ref_seqs[i][0]]}")
            #print(f"query on-score: {seqs_on_target_score[query_seq[0]]}")
            # r_ratio = seqs_on_target_score[query_seq[0]]/ seqs_on_target_score[ref_seqs[i][0]]
            r_ratio = seqs_on_target_score[ref_seqs[i][0]]/seqs_on_target_score[query_seq[0]] # Off/On

            ref_seq = ref_seqs[i][0]
            mismatch_locations, mismatch_keys = self.get_mismatches(ref_seq, query_seq[0],three_prime,gRNA_len)
            if len(mismatch_locations) <= 10:
                value = (math.sqrt(self.sh_score(mismatch_locations, mismatch_keys, HSU_matrix,gRNA_len)) + self.st_score(mismatch_locations)) * (r_ratio ** 2) * (self.ss_score(mismatch_locations,gRNA_len)** 6)
                # value = (self.st_score(mismatch_locations)) * (r_ratio ** 2) * (self.ss_score(mismatch_locations,endo)** 6)
                # value = (math.sqrt(1.2984) + self.st_score(mismatch_locations)) * (r_ratio ** 2) * (self.ss_score(mismatch_locations)** 6)
                value /= 4

                target_scores.append(value)
                sh_scores.append(self.sh_score(mismatch_locations,mismatch_keys,HSU_matrix,gRNA_len)) # Append Hsu scores to list

                """ Print statements to check values """
                if ref_seq == "ATGCGAGTCCGGGCAGGAGAAAAA":
                    print(f"mismatch_locations: {mismatch_locations}")
                    print(f"mismatch_keys: {mismatch_keys}")
                    print(f"r_ratio: {r_ratio}")
                    print(f"st_score: {self.st_score(mismatch_locations)}")
                    print(f"ss_score: {self.ss_score(mismatch_locations,gRNA_len)}")
                    print(f"sh_score: {self.sh_score(mismatch_locations, mismatch_keys, HSU_matrix,gRNA_len)}")
                    print(f"score: {value}")

        avg_score = sum(target_scores)
        if len(target_scores) != 0:
            avg_score /= len(target_scores)

        return avg_score, sh_scores, target_scores

    def score_sequences(self, query_seqs, ref_seqs, sequence_scores, HSU_matrix,endo,three_prime,gRNA_len):
        all_OT_scores = []
        all_sh_scores = []
        all_casper_scores = []

        for i in range(len(query_seqs)):
            off_target_scores, sh_scores, casper_scores = self.find_similars(query_seqs[i], ref_seqs[i], sequence_scores, HSU_matrix,endo,three_prime,gRNA_len)
            all_OT_scores.append(off_target_scores)
            all_casper_scores.append(casper_scores)
            all_sh_scores.append(sh_scores)

        return all_OT_scores, all_sh_scores, all_casper_scores


if __name__ == "__main__":
    #input/output file variables
    CASPER_info_file_path = "CASPERinfo"
    input_on_target_sequences_file_path = "on_target_sequences.txt"
    input_off_target_sequences_file_path = "5_seqs.txt"

    output_file_path = "results.txt"

    #create class objects
    file_operations_object = file_operations()
    on_target_object = on_target()
    off_target_object = off_target()

    #parse all input data: CRISPR scan, HSU matrix, and input sequences
    CRISPRSCAN_data = file_operations_object.read_CRISPRSCAN(CASPER_info_file_path)
    HSU_matrix = file_operations_object.read_HSU_matrix(CASPER_info_file_path)
    input_on_target_sequences = file_operations_object.read_on_target_sequences(input_on_target_sequences_file_path)
    query_seqs, ref_seqs = file_operations_object.read_off_target_sequences(input_off_target_sequences_file_path)

    #get on-target scores for off-target data - score based on long sequence
    OT_on_scores = {}

    for seq in query_seqs:
        # OT_on_scores[seq[0]] = on_target_object.score_sequence(seq[0], seq[1], CRISPRSCAN_data, endo="spCas9",three_prime=True,gRNA_len=20,pam="NGG")
        OT_on_scores[seq[0]] = on_target_object.score_sequence(seq[0], seq[1], CRISPRSCAN_data, endo="asCas12",three_prime=False,gRNA_len=24,pam="TTTC")

    for section in ref_seqs:
        for seq in section:
            # OT_on_scores[seq[0]] = on_target_object.score_sequence(seq[0], seq[1], CRISPRSCAN_data, endo="spCas9", three_prime=True,gRNA_len=20,pam="NGG")
            OT_on_scores[seq[0]] = on_target_object.score_sequence(seq[0], seq[1], CRISPRSCAN_data, endo="asCas12", three_prime=False,gRNA_len=24,pam="TTTC")


    #calcualte off-target scores
    # OT_scores, sh_scores, casper_scores = off_target_object.score_sequences(query_seqs, ref_seqs, OT_on_scores, HSU_matrix,endo="spCas9",three_prime=True,gRNA_len=20)
    OT_scores, sh_scores, casper_scores = off_target_object.score_sequences(query_seqs, ref_seqs, OT_on_scores, HSU_matrix,endo="asCas12",three_prime=False,gRNA_len=24)

    # print(OT_scores)
    # print(sh_scores)

    #write out off-target results
    file_operations_object.write_OT_results(query_seqs, ref_seqs, OT_scores, sh_scores,casper_scores)

    # Graph stuff
    my_figure = plt.Figure()
    plt.hist(casper_scores,bins=10,range=(0,1.0))
    plt.show()
    