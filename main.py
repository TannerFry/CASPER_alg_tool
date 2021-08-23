from Bio.Seq import Seq
import math

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
        input_sequences = {}

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
                            input_sequences[line] = temp
                else:
                    temp.append(line)

        return input_sequences

    def read_on_target_sequences(self, input_on_target_sequences_file_path):
        on_target_sequences = []
        with open(input_on_target_sequences_file_path, "r") as fin:
            for line in fin:
                line = line.strip("\n")
                on_target_sequences.append(line)
        return on_target_sequences

    def write_OT_results(self):
        pass


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

    def get_sc(self, seq, CRISPRSCAN_data):
        seq = seq.upper()
        score = 0
        """ Tally the CRISPRscan score for the sequence """
        for i, nt in enumerate(seq):
            if i == len(seq) - 1:
                dnt = ""
            dnt = seq[i:i + 2]
            if i in CRISPRSCAN_data.keys():
                if (nt + "x") in CRISPRSCAN_data[i].keys():
                    score += CRISPRSCAN_data[i][nt + "x"]
                if dnt in CRISPRSCAN_data[i].keys():
                    score += CRISPRSCAN_data[i][dnt]

        # Return the CRISPRscan score
        score += 0.183930944
        return (score * 100)

    def get_sij(self, seq):
        seq = seq[6:len(seq) - 9]
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

    def get_sg(self, seq):
        seq = seq[6:len(seq) - 9]
        score = 0
        for nt in seq:
            if nt == "G":
                score += 1.0
            elif nt == "C":
                score += 0.5
            elif nt == "A":
                score -= 0.1
        return score / 20

    def get_p(self, sij_score, sg_score):
        if sij_score > 1:
            p = sij_score * sg_score
        elif sij_score == 1:
            p = 1 - (sg_score / 5)
        else:
            p = 0
        return p

    def score_sequence(self, seq, CRISPRSCAN_data):
        sc_score = self.get_sc(seq, CRISPRSCAN_data)
        sij_score = self.get_sij(seq)
        sg_score = self.get_sg(seq)
        p_score = self.get_p(sij_score, sg_score)
        if p_score == 0:
            score = sc_score
        else:
            score = sc_score / p_score
        return score


class off_target():
    def __init__(self):
        pass

    def sh_score(self, mismatches, hsu_keys, hsu_matrix):
        sh_score = 1.0
        for i in range(len(mismatches)):

            sh_score *= hsu_matrix[hsu_keys[i]][mismatches[i]]

        return sh_score

    def ss_score(self, mismatches):
        ss_score = 1.0
        for i in range(len(mismatches)):
            if mismatches[i] < 6:
                ss_score -= 0.1
            elif mismatches[i] < 12:
                ss_score -= 0.05
            else:
                ss_score -= 0.0125
        return ss_score

    def st_score(self, mismatches):
        st_score = 3.5477
        for i in range(len(mismatches)):
            st_score -= (1.0 / (mismatches[i] + 1))
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

    def get_mismatches(self, ref_seq, query_seq, seq_length=20):
        mismatch_locations = []
        mismatch_keys = []

        for i in range(seq_length-1, -1, -1):
            if ref_seq[i] != query_seq[i]:
                mismatch_locations.append(i)
                mismatch_keys.append(query_seq[i] + self.reverse_comp(ref_seq[i]))

        return mismatch_locations, mismatch_keys

    def find_similars(self, query_seq, ref_seqs, seqs_on_target_score, HSU_matrix):
        target_scores = []
        sh_scores = []

        for i in range(len(ref_seqs)):
            r_ratio = seqs_on_target_score[ref_seqs[i]] / seqs_on_target_score[query_seq]
            ref_seq = ref_seqs[i]
            mismatch_locations, mismatch_keys = self.get_mismatches(ref_seq, query_seq)
            if len(mismatch_locations) <= 5:
                value = (math.sqrt(self.sh_score(mismatch_locations, mismatch_keys, HSU_matrix)) + self.st_score(mismatch_locations)) * (r_ratio ** 2) * (self.ss_score(mismatch_locations)** 6)
                value /= 4
                sh_scores.append(self.sh_score(mismatch_locations, mismatch_keys, HSU_matrix))
                target_scores.append(value)

        avg_score = sum(target_scores)
        if len(target_scores) != 0:
            avg_score /= len(target_scores)

        return avg_score

    def score_sequences(self, sequences, sequence_scores, HSU_matrix):
        off_target_scores = []

        for query_seq in sequences.keys():
            off_target_scores.append(self.find_similars(query_seq, sequences[query_seq], sequence_scores, HSU_matrix))

        return off_target_scores

if __name__ == "__main__":
    #input/output file variables
    CASPER_info_file_path = "CASPERinfo"
    input_on_target_sequences_file_path = "on_target_sequences.txt"
    input_off_target_sequences_file_path = "off_target_sequences.txt"
    output_file_path = "results.txt"

    #create class objects
    file_operations_object = file_operations()
    on_target_object = on_target()
    off_target_object = off_target()

    #parse all input data: CRISPR scan, HSU matrix, and input sequences
    CRISPRSCAN_data = file_operations_object.read_CRISPRSCAN(CASPER_info_file_path)
    HSU_matrix = file_operations_object.read_HSU_matrix(CASPER_info_file_path)
    input_on_target_sequences = file_operations_object.read_on_target_sequences(input_on_target_sequences_file_path)
    input_off_target_sequences = file_operations_object.read_off_target_sequences(input_off_target_sequences_file_path)

    #get on-target scores for off-target data
    OT_on_scores = {}
    for query_seq in input_off_target_sequences.keys():
        OT_on_scores[query_seq] = on_target_object.score_sequence(query_seq, CRISPRSCAN_data)
        for ref_seq in input_off_target_sequences[query_seq]:
            OT_on_scores[ref_seq] = on_target_object.score_sequence(ref_seq, CRISPRSCAN_data)

    #calcualte off-target scores
    OT_scores = off_target_object.score_sequences(input_off_target_sequences, OT_on_scores, HSU_matrix)
    print(OT_scores)

    #write out off-target results
    #file_operations_object.write_OT_results()