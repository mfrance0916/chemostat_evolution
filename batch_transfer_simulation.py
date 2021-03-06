#!/usr/bin/python
import random
import numpy as np
import os

def dist_fit_eff():

    selective_effect = random.uniform(-1.025, 0.05)

    if -1.25 < selective_effect < 0:
        
        selective_effect = selective_effect * selective_effect * -1.0

    if selective_effect < -1.0:

        selective_effect = -1.0

    if selective_effect > 0.0:

        selective_effect = np.random.exponential(scale=0.025) - 0.07

    return selective_effect

class Population(object):
    
    ## defines the following parameters that are needed for the simulation
    ## pop_dens = density of bacterial population in cells/ml
    ## volume = volume of the chemostat in mL
    ## generations per day specicfies an integer number of generations to run each day of batch culture 
    ## total_generations = numerical number of generations to run the simulation for
    ## batch length = amount of real time each transfer takes in hours, usually 24
    ## generation_length = amount of real time each generation of growth takes, usually ~1 hr
    ## gd_mut_rate = growth dependent mutation rate
    ## gi_mut_rate = growth independent mutation rate
    
    def __init__(self, pop_dens, volume, generations_per_day, total_generations, batch_length,generation_length, gd_mut_rate, gi_mut_rate, gd_ab_mut_rate, gi_ab_mut_rate,target):
        self.pop_dens = pop_dens
        self.volume = volume
        self.generations_per_day = generations_per_day
        self.total_generations = total_generations
        self.batch_length = batch_length
        self.generation_length = generation_length
        self.gd_mut_rate = gd_mut_rate
        self.gi_mut_rate = gi_mut_rate
        self.gd_ab_mut_rate = gd_ab_mut_rate
        self.gi_ab_mut_rate = gi_ab_mut_rate
        self.target = target
        
    def initialize(self):
        
        #calculate population size
        self.pop_size = self.pop_dens*self.volume
        #calculate dilution based on desired number of generations per day 
        self.dilution = 2 ** self.generations_per_day

        #create list of individuals with a relative fitness of 1, turns out it's too memory intensive, because competition is global we don't really need to track individual cells anyway, instead we will track the subpopulations
        #self.fitness = list([1.0]*self.pop_size)    
        #self.ab_status = list([1]*self.pop_size)

        #creating a csv file that will keep track of all of the subpopulations

        subpopulations = open('subpopulations.csv','w')
        subpopulations.write("1,1.0,%s,1,0,0,0" %(self.pop_size))
        subpopulations.close()
        
        #setting the current generation to 0
        self.current_generation = 0
        
        #creating a counting variable for the time in hours
        self.time = 0

        #determining the number to kill at each serial transfer 
        self.survive_prob = 1.0/self.dilution
        #print self.survive_prob

        self.subpop_id_num = 1

        #opening file to store output data
        pop_stat_out = open("population_statistics.csv",'w')
        pop_stat_out.write('gen,pop_size,lineage_count,ave_fit,maxfit,abund_lin,freq_res,fit_res\n')
        pop_stat_out.write('0,0,' + str(self.pop_size) + ',1,1,1.0,1,0.0,0\n') 
        pop_stat_out.close()

        top_lineages = open("top_lineages.csv","w")
        top_lineages.write('gen,id,rel_abund,total_abund,fitness,rel_fitness,gen_form,cycleform,barcode,abstatus\n')
        top_lineages.write('0,0,1.0,' + str(self.pop_size) + ',1.0,1.0,0,0,0,0\n')

        top_lineages.close()




    def update_time(self):

        #setting the starting generation to 0
        stationary_phase_time = int(self.batch_length) - int(self.generations_per_day)*int(self.generation_length)
        current_generation = 0

        #looping from 0 generations to desired number of generations set by an input parameter self.total_generations
        while current_generation < self.total_generations:

            #implementing batch transfer at the start of the day

            #creating the post transfer file to write subpopulations to
            post_transfer = open('post_transfer.csv','w')

            #looping through the current subpopulations to select how many will survive
            for subpop in open('subpopulations.csv'):

                #stripping off new line character
                subpop = subpop.rstrip('\n') 
                
                #pulling out the different subpopulation variables from the csv text
                subpop_id = subpop.split(',')[0]
                subpop_fitness = float(subpop.split(',')[1])
                subpop_count = int(subpop.split(',')[2])
                subpop_cycleform = int(subpop.split(",")[3])
                subpop_genform = int(subpop.split(",")[4])
                subpop_barcode = str(subpop.split(",")[5])
                subpop_ab_stat = (subpop.split(",")[6])


                #implementing random death process only if there is an extant member of the lineage, otherwise drop 
                if subpop_count > 0:
                    #random death process simulated on the subpopulation as a draw from a binomial distribution with p = probability of being transfered and n = subpopulation count
                    num_survivors = np.random.poisson(subpop_count*self.survive_prob)
                    #print num_survivors

                    #writing out the lineage to the post transfer file as long as the 
                    if num_survivors > 0:    
                        post_transfer.write(subpop_id + ',' + str(subpop_fitness) + ',' + str(num_survivors) + ',' + str(subpop_cycleform) + ',' + str(subpop_genform) + ',' + str(subpop_barcode) + ',' + str(subpop_ab_stat) + '\n')

            #closing the post transfer file so that it can be used in a loop later
            post_transfer.close()

            cycle_generation = 0



            #looping through the number of generations in each day, cultures will be serial transfer at end of day
            for generation_step in range(1,(int(self.generations_per_day) + 1)):

                #CALCULATING THE AVERAGE FITNESS OF THE POPULATION
                
                #setting up a bunch of reporter variables
                total_pop_size = 0
                total_fitness = 0.0
                currentfit = 0.0
                total_ab_size = 0
                total_ab_fitness = 0

                #looping through lineages in the post transfer file
                for lineage in open('post_transfer.csv'):
                
                    #stripping off new line from lineage data
                    lineage = lineage.rstrip('\n')

                    #adding up counts to determine the total population size
                    total_pop_size = total_pop_size + int(lineage.split(',')[2])

                    #adding up fitness as the multiple of individual fitness and the number of individuals
                    total_fitness = total_fitness + (self.target - (self.target- float(lineage.split(',')[1])))*int(lineage.split(',')[2])


                average_fitness = total_fitness/total_pop_size

                if lineage.split(',')[6] == '1':
                    total_ab_size = int(total_ab_size) + int(lineage.split('')[2])
                    total_ab_fitness = total_ab_fitness + (self.target - (self.target- float(lineage.split(',')[1])))*int(lineage.split(',')[2])

                if int(total_ab_size) > 0:
                    average_ab_fitness = total_ab_fitness/total_ab_size 
                else:
                    average_ab_fitness = 0.0

                num_progeny = total_pop_size * 2

                #print average_fitness
                
                current_generation += 1
                
                



                #looping through to implement growth independent mutations 
                post_gi_mutation = open('post_gi_mutation.csv','w')

                for subpop_post in open('post_transfer.csv'):

                    subpop_post = subpop_post.rstrip('\n') 

                    #print subpop_post

                    #storing variables about the ancestral lineages for future use
                    subpop_id = subpop_post.split(',')[0]
                    subpop_fitness = float(subpop_post.split(',')[1])
                    subpop_count = int(subpop_post.split(',')[2])
                    subpop_cycleform = int(subpop_post.split(',')[3])
                    subpop_genform = int(subpop_post.split(',')[4])
                    subpop_barcode = str(subpop_post.split(',')[5])
                    subpop_ab_stat = (subpop_post.split(",")[6])


                    #determining the number of new growth independent mutants to generate
                    new_gi_mutants = np.random.poisson(subpop_count*self.gi_mut_rate)

                    #adjusting the ancestral subpopulation size
                    subpop_count = subpop_count - new_gi_mutants
                
                    #looping through the new growth independent mutants and storing their information in a file
                    for mutant in range(0, new_gi_mutants):

                        #incrementing the new subpop id number by 1
                        self.subpop_id_num += 1

                        #store the new subpop id in a variable for output later
                        new_mutant_id = self.subpop_id_num

                        #calculating the fitness effect of the new mutatation, this is a simple model for the fitness of effect of new mutations and will need to be improved by a real model
                        selective_eff = dist_fit_eff()

                        #determining the new mutant fitness as the addition of the selective effect to the ancestors fitness
                        new_mutant_fitness = subpop_fitness + selective_eff

                        #taking care of negative fitness cases
                        if new_mutant_fitness < 0.0:
                            #print('gi  ' + str(new_mutant_fitness))
                            new_mutant_fitness = 0.0

                        #setting the count of the new mutant equal to 1
                        new_mutant_count = 1

                        new_mutant_cycleform = cycle_generation
                        new_mutant_genform = current_generation
                        new_mutant_barcode = str(subpop_barcode) + ":" + str(new_mutant_id)


                        #writing out the information for the new mutants to a temporary file
                        post_gi_mutation.write(str(new_mutant_id) + ',' + str(new_mutant_fitness) + ',' + str(new_mutant_count) + "," + str(new_mutant_cycleform) + "," + str(new_mutant_genform) + "," + str(new_mutant_barcode) + ',' + str(subpop_ab_stat) + '\n')
                    
                    #writing the original lineages data to the file
                    if subpop_count <= 0:
                        continue

                    if subpop_ab_stat == '0':

                        new_gi_ab_mutants = np.random.poisson(subpop_count*self.gi_ab_mut_rate)
                        subpop_count = subpop_count - new_gi_ab_mutants

                        if new_gi_ab_mutants > 0:

                            self.subpop_id_num += 1
                            new_mutant_id = self.subpop_id_num
                            new_mutant_genform = current_generation
                            new_mutant_cycleform = cycle_generation
                            new_mutant_barcode = str(subpop_barcode) + ':' + str(new_mutant_id)
                            post_gi_mutation.write(str(new_mutant_id) + ',' + str(subpop_fitness) + ',' + str(new_gi_ab_mutants) + "," + str(new_mutant_cycleform) + "," + str(new_mutant_genform) + "," + str(new_mutant_barcode) + ',1,\n')

                    if subpop_count <= 0:
                        continue

                    post_gi_mutation.write(str(subpop_id) + ',' + str(subpop_fitness) + ',' + str(subpop_count) + "," + str(subpop_cycleform) + "," + str(subpop_genform) + "," + str(subpop_barcode) + ',' + str(subpop_ab_stat) + '\n')


                post_gi_mutation.close()

                cycle_generation += 1

                #ENACTING GROWTH AND MUTATION PROCESS THAT IS DEPENDENT ON GROWTH
                os.remove('post_transfer.csv')

                post_growth = open('post_transfer.csv','w')

                for lineage in open('post_gi_mutation.csv'):

                    lineage = lineage.rstrip('\n')
                
                    #storing variables for use
                    subpop_id = lineage.split(',')[0]
                    subpop_fitness = float(lineage.split(',')[1])
                    subpop_count = int(lineage.split(',')[2])
                    subpop_cycleform = int(lineage.split(',')[3])
                    subpop_genform = int(lineage.split(',')[4])
                    subpop_barcode = str(lineage.split(',')[5])
                    subpop_ab_stat = (lineage.split(',')[6])

                
                    #SELECTION THINGS WITH HIGHER FITNESS DO BETTER

                    #calculating the relative fitness as the division of the lineage fitness by the average fitness
                    #no target version
                    #rel_fitness = subpop_fitness/average_fitness
                    #target version
                    rel_fitness = (self.target-abs(self.target - subpop_fitness))/(average_fitness)


                    #determining the number of cellular divisions generated by the lineage 
                    divisions = np.random.poisson(subpop_count*rel_fitness)

                    subpop_count = subpop_count + divisions

                    if divisions > 0:

                        daughters = divisions * 2

                        #GROWTH DEPENDENT MUTATION
                        new_gd_mutants = np.random.poisson(daughters*self.gd_mut_rate)

                        subpop_count = subpop_count - new_gd_mutants

                        for mutant in range(0, new_gd_mutants):

                            #incrementing the new subpop id number by 1
                            self.subpop_id_num += 1

                            #store the new subpop id in a variable for output later
                            new_mutant_id = self.subpop_id_num

                            #calculating the fitness effect of the new mutatation, this is a simple model for the fitness of effect of new mutations and will need to be improved by a real model
                            selective_eff = dist_fit_eff()

                            #determining the new mutant fitness as the addition of the selective effect to the ancestors fitness
                            new_mutant_fitness = subpop_fitness + selective_eff

                            #taking care of negative fitness cases
                            if new_mutant_fitness < 0.0:
                                #  print('gd  ' + str(new_mutant_fitness))
                                new_mutant_fitness = 0.0

                            #setting the count of the new mutant equal to 1
                            new_mutant_count = 1

                            new_mutant_cycleform = cycle_generation
                            new_mutant_genform = current_generation
                            new_mutant_barcode = str(subpop_barcode) + ":" + str(new_mutant_id)
    
                            #writing out the data to the post growth file
                            post_growth.write(str(new_mutant_id) + ',' + str(new_mutant_fitness) + ',' + str(new_mutant_count) + "," + str(new_mutant_cycleform) + "," + str(new_mutant_genform) + "," + str(new_mutant_barcode) + ',' + str(subpop_ab_stat) + '\n')

                    
                    #writing the original lineages data to the file
                    
                    if subpop_count <= 0:
                        continue

                    if subpop_ab_stat == '0':

                        new_gd_ab_mutants = np.random.poisson((daughters - new_gd_mutants)*self.gd_ab_mut_rate)
                        subpop_count = subpop_count - new_gd_ab_mutants

                        if new_gi_ab_mutants > 0:

                            self.subpop_id_num += 1
                            new_mutant_id = self.subpop_id_num
                            new_mutant_genform = current_generation
                            new_mutant_cycleform = cycle_generation
                            new_mutant_barcode = str(subpop_barcode) + ':' + str(new_mutant_id)
                            post_growth.write(str(new_mutant_id) + ',' + str(subpop_fitness) + ',' + str(new_gi_ab_mutants) + "," + str(new_mutant_cycleform) + "," + str(new_mutant_genform) + "," + str(new_mutant_barcode) + ',1,\n')

                    if subpop_count <= 0:
                        continue

                    post_growth.write(str(subpop_id) + ',' + str(subpop_fitness) + ',' + str(subpop_count) + "," + str(subpop_cycleform) + "," + str(subpop_genform) + "," + str(subpop_barcode) + ',' + str(subpop_ab_stat) + '\n')

                post_growth.close()

                


            #looping through the final post transfer file to initialize a growth independent step which will be then be used in the transfer step at the start of the loops
            new_subpopulations = open('subpopulations.csv','w')

            for subpop_post in open("post_transfer.csv"):

                subpop_post = subpop_post.rstrip('\n') 

                #print subpop_post

                #storing variables about the ancestral lineages for future use
                subpop_id = subpop_post.split(',')[0]
                subpop_fitness = float(subpop_post.split(',')[1])
                subpop_count = int(subpop_post.split(',')[2])
                subpop_cycleform = int(subpop_post.split(',')[3])
                subpop_genform = int(subpop_post.split(',')[4])
                subpop_barcode = str(subpop_post.split(',')[5])
                subpop_ab_stat = (subpop_post.split(',')[6])

                #calculating amount of time in stationary phase
                
                #print stationary_phase_time
                #determining the number of new growth independent mutants to generate
                #print subpop_count
                new_gi_mutants = np.random.poisson(subpop_count*self.gi_mut_rate*stationary_phase_time)
                #print new_gi_mutants

                #adjusting the ancestral subpopulation size
                subpop_count = subpop_count - new_gi_mutants
                
                #looping through the new growth independent mutants and storing their information in a file
                for mutant in range(0, new_gi_mutants):

                    #incrementing the new subpop id number by 1
                    self.subpop_id_num += 1

                        #store the new subpop id in a variable for output later
                    new_mutant_id = self.subpop_id_num

                    #calculating the fitness effect of the new mutatation, this is a simple model for the fitness of effect of new mutations and will need to be improved by a real model
                    selective_eff = dist_fit_eff()

                    #determining the new mutant fitness as the addition of the selective effect to the ancestors fitness
                    new_mutant_fitness = subpop_fitness + selective_eff

                    #taking care of negative fitness cases
                    if new_mutant_fitness < 0.0:
                        #print('gi  ' + str(new_mutant_fitness))
                        new_mutant_fitness = 0.0

                    #setting the count of the new mutant equal to 1
                    new_mutant_count = 1

                    new_mutant_cycleform = 11
                    new_mutant_genform = current_generation
                    new_mutant_barcode = str(subpop_barcode) + ":" + str(new_mutant_id)


                    #writing out the information for the new mutants to a temporary file
                    new_subpopulations.write(str(new_mutant_id) + ',' + str(new_mutant_fitness) + ',' + str(new_mutant_count) + "," + str(new_mutant_cycleform) + "," + str(new_mutant_genform) + "," + str(new_mutant_barcode) + ',' + str(subpop_ab_stat) + '\n')
                    
                #writing the original lineages data to the file
                if subpop_count <= 0:
                        continue
                
                if subpop_ab_stat == '0':

                    new_gi_ab_mutants = np.random.poisson(subpop_count*self.gi_ab_mut_rate*stationary_phase_time)
                    subpop_count = subpop_count - new_gi_ab_mutants

                    if new_gi_ab_mutants > 0:

                        self.subpop_id_num += 1
                        new_mutant_id = self.subpop_id_num
                        new_mutant_genform = current_generation
                        new_mutant_cycleform = 11
                        new_mutant_barcode = str(subpop_barcode) + ':' + str(new_mutant_id)
                        new_subpopulations.write(str(new_mutant_id) + ',' + str(subpop_fitness) + ',' + str(new_gi_ab_mutants) + "," + str(new_mutant_cycleform) + "," + str(new_mutant_genform) + "," + str(new_mutant_barcode) + ',1,\n')

                    if subpop_count <= 0:
                        continue


                new_subpopulations.write(str(subpop_id) + ',' + str(subpop_fitness) + ',' + str(subpop_count) + "," + str(subpop_cycleform) + "," + str(subpop_genform) + "," + str(subpop_barcode) + ',' + str(subpop_ab_stat) + '\n')

            new_subpopulations.close()

                      
            pop_stat_out = open("population_statistics.csv",'a')

            top_lineages = open("top_lineages.csv","a")


            total_pop_size = 0
            total_fitness = 0.0
            currentfit = 0.0
            maxfit = 0.0
            lineage_count = 0
            total_ab_size = 0
            total_ab_fitness = 0.0

            print current_generation

            for subpop in open('subpopulations.csv'):

                subpop = subpop.rstrip('\n')
                    
                #pulling out the different subpopulation variables
                subpop_id = subpop.split(',')[0]
                subpop_fitness = float(subpop.split(',')[1])
                subpop_count = int(subpop.split(',')[2])

                total_pop_size = total_pop_size + int(subpop.split(',')[2])

                #adding up fitness as the multiple of individual fitness and the number of individuals
                total_fitness = total_fitness + (self.target-abs(self.target - subpop_fitness))*subpop_count

                #determining the best organism at the time
                if (self.target-abs(self.target - subpop_fitness)) > maxfit:
                    maxfit = currentfit

                if lineage.split(',')[6] == '1':
                    total_ab_size = int(total_ab_size) + int(lineage.split('')[2])
                    total_ab_fitness = total_ab_fitness + (self.target - (self.target- float(lineage.split(',')[1])))*int(lineage.split(',')[2])

            average_ab_fitness = total_ab_fitness/total_ab_size 

            average_fitness = float(total_fitness)/float(total_pop_size)

            num_abund_lineages = 0

            for subpop in open('subpopulations.csv'):

                subpop = subpop.rstrip('\n')

                lineage_count += 1
                
                subpop_id = subpop.split(',')[0]
                subpop_fitness = float(subpop.split(',')[1])
                subpop_count = int(subpop.split(',')[2])
                subpop_cycleform = int(subpop.split(',')[3])
                subpop_genform = int(subpop.split(',')[4])                
                subpop_barcode = str(subpop.split(',')[5])
                subpop_ab_stat = str(subpop.split(',')[6])

                if int(subpop_count) > int(total_pop_size/1000):

                    num_abund_lineages += 1

                    subpop_rel_fitness = (self.target-abs(self.target - subpop_fitness))/average_fitness
                    subpop_rel_abund = float(subpop_count)/float(total_pop_size)

                    top_lineages.write(str(current_generation) + "," + str(subpop_id) + "," + str(subpop_rel_abund) + "," + str(subpop_count) + "," + str(subpop_fitness) + "," + str(subpop_rel_fitness) + "," + str(subpop_cycleform) + "," + str(subpop_genform) + "," + str(subpop_barcode) + ',' + str(subpop_ab_stat) + "\n")


            #print average_fitness
            #print total_pop_size

            freq_ab = total_ab_size / float(total_pop_size)

            pop_stat_out.write(str(current_generation) + "," + str(total_pop_size)  + ',' + str(lineage_count) + "," + str(average_fitness) + "," + str(maxfit) + "," + str(num_abund_lineages) + ',' + str(freq_ab) + ',' + str(average_ab_fitness) + "\n")

            pop_stat_out.close()

            top_lineages.close()

#def __init__(self, pop_dens, volume, generations_per_day, total_generations, batch_length,generation_length, gd_mut_rate, gi_mut_rate):   

exp_evol = Population(3000000000,5,10,1000,24,1,0.0006,0.00015,0.00000012,0.0000000055,2.0)

exp_evol.initialize()

exp_evol.update_time()
