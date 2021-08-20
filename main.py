from Bio.Seq import Seq

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

        return HSU_matrix_data

    def read_off_target_sequences(self, input_sequences_file_path):
        #create input sequences object
        input_sequences = []

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
                            input_sequences.append(temp)
                            temp = []
                            break
                        else:
                            temp.append(line)
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

    def write_results(self):
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

    def score_sequences(self, sequences, CRISPRSCAN_data):
        sc_scores = []
        scores = []
        for i, seq in enumerate(sequences):
            sc_score = self.get_sc(seq, CRISPRSCAN_data)
            sc_scores.append(sc_score)
            sij_score = self.get_sij(seq)
            sg_score = self.get_sg(seq)
            p_score = self.get_p(sij_score, sg_score)
            if p_score == 0:
                scores.append(sc_score)
            else:
                scores.append(sc_score / p_score)
        return scores


class off_target():
    def __init__(self):
        pass


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
    HSU_matrix_data = file_operations_object.read_HSU_matrix(CASPER_info_file_path)
    input_on_target_sequences = file_operations_object.read_on_target_sequences(input_on_target_sequences_file_path)

    #test on-scoring algorithm with on-score sequences
    on_target_scores = on_target_object.score_sequences(input_on_target_sequences, CRISPRSCAN_data)

    print(on_target_scores)