import numpy as np
import csv
import main
import pandas as pd
import random
import math

class random_sampler():
    def __init__(self) -> None:
        self.Off_target = main.off_target()
        self.On_target = main.on_target()
        self.file_operations_object = main.file_operations()

        #read CASPERinfo data for algorithms
        self.CRISPRSCAN_data = self.file_operations_object.read_CRISPRSCAN("CASPERinfo")
        self.HSU_matrix = self.file_operations_object.read_HSU_matrix("CASPERinfo")
    
    def read_input_data(self):
        #read input data
        self.OT_data = []
        with open("off_target_sampling_data.csv", "r") as f:
            f_csv = csv.reader(f)
            for line in f_csv:
                self.OT_data.append(line)
        
        #get On target scores for all sequences
        for i in range(len(self.OT_data)):
            on_target_score_1 = self.On_target.score_sequence(self.OT_data[i][0], self.OT_data[i][1], self.CRISPRSCAN_data, "spCa9", three_prime=True, gRNA_len=20, pam="NGG")
            on_target_score_2 = self.On_target.score_sequence(self.OT_data[i][3], self.OT_data[i][4], self.CRISPRSCAN_data, "spCa9", three_prime=True, gRNA_len=20, pam="NGG")
            self.OT_data[i] = [self.OT_data[i][0], self.OT_data[i][1], on_target_score_1, self.OT_data[i][3], self.OT_data[i][4], on_target_score_2, float(self.OT_data[i][2])]
        
        #print(self.OT_data[0])
        #self.OT_data sample structure: [OffTarget gRNA, OffTarget Full Seq, OffTarget On-Target score, OnTarget gRNA, OnTarget Full Seq, OnTarget On-Target Score, indel score]

    def run(self):
        self.r_squared_sh_scores = []
        self.r_squared_target_scores = []

        for i in range(1000):
            #get 90 random samples from the OT data
            random_samples = random.sample(self.OT_data, 90)

            #get s_h and target scores for each sequence sample
            sh_scores = []
            target_scores = []
            relative_indels = []
            for sample in random_samples:
                off_target_grna = sample[0]
                on_target_grna = sample[3]
                
                r_ratio = sample[2] / sample[5]
                
                mismatch_locations, mismatch_keys = self.Off_target.get_mismatches(off_target_grna, on_target_grna)
                
                value = (math.sqrt(self.Off_target.sh_score(mismatch_locations, mismatch_keys, self.HSU_matrix)) + self.Off_target.st_score(mismatch_locations)) * (r_ratio ** 2) * (self.Off_target.ss_score(mismatch_locations)** 6)
                value /= 4
                
                sh_scores.append(self.Off_target.sh_score(mismatch_locations, mismatch_keys, self.HSU_matrix))
                target_scores.append(value)

                #store relative indel
                relative_indels.append(sample[-1])


            # #calc sh r squared value
            sh_r_squared = self.calc_r_squared(relative_indels, sh_scores)
            self.r_squared_sh_scores.append(sh_r_squared)

            # #calc target r squared value
            target_r_squared = self.calc_r_squared(relative_indels, target_scores)
            self.r_squared_target_scores.append(target_r_squared)


    def calc_r_squared(self, x, y):
        correlation_matrix = np.corrcoef(x, y)
        correlation_xy = correlation_matrix[0,1]
        r_squared = correlation_xy**2
        return r_squared

    def write_results(self):
        with open("ot_random_sampling_results.csv","w") as f:
            f.write("HSU Score R^2,Target Score R^2\n")
            for i in range(len(self.r_squared_sh_scores)):
                f.write(f"{self.r_squared_sh_scores[i]},{self.r_squared_target_scores[i]}\n")

if __name__ == "__main__":
    #create sampling object
    random_sampler_object = random_sampler()
    random_sampler_object.read_input_data()
    random_sampler_object.run()
    random_sampler_object.write_results()