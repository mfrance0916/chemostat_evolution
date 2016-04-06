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

        selective_effect = np.random.exponential(scale=0.03) - 0.06 

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
    
    def __init__(self, pop_dens, volume, generations_per_day, total_generations, batch_length,generation_length, gd_mut_rate, gi_mut_rate):
        self.pop_dens = pop_dens
        self.volume = volume
        self.generations_per_day = generations_per_day
        self.total_generations = total_generations
        self.batch_length = batch_length
        self.generation_length = generation_length
        self.gd_mut_rate = gd_mut_rate
        self.gi_mut_rate = gi_mut_rate
        
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
        subpopulations.write("1,1.0,%s" %(self.pop_size))
        subpopulations.close()
        
        #setting the current generation to 0
        self.current_generation = 0
        
        #creating a counting variable for the time in hours
        self.time = 0

        #determining the number to kill at each serial transfer 
        self.survive_prob = 1.0/self.dilution
        #print self.survive_prob

        self.subpop_id_num = 1.0/self.dilution

        #opening file to store output data
        pop_stat_out = open("population_statistics.csv",'w')
        pop_stat_out.write('time,gen,pop_size,num_lineages,ave_fit,abund_lin\n')
        pop_stat_out.write('0,0,' + str(self.pop_size) + ',1,1.0,1\n') 
        pop_stat_out.close()

    def update_time(self):

        #setting the starting generation to 0
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

                #implementing random death process only if there is an extant member of the lineage, otherwise drop 
                if subpop_count > 0:
                    #random death process simulated on the subpopulation as a draw from a binomial distribution with p = probability of being transfered and n = subpopulation count
                    num_survivors = np.random.binomial(subpop_count,self.survive_prob)
                    print num_survivors

                    #writing out the lineage to the post transfer file as long as the 
                    if num_survivors > 0:    
                        post_transfer.write(subpop_id + ',' + str(subpop_fitness) + ',' + str(num_survivors) + '\n')

            #closing the post transfer file so that it can be used in a loop later
            post_transfer.close()



            #looping through the number of generations in each day, cultures will be serial transfer at end of day
            for generation_step in range(1,(int(self.generations_per_day) + 1)):

                #CALCULATING THE AVERAGE FITNESS OF THE POPULATION
                
                #setting up a bunch of reporter variables
                total_pop_size = 0
                total_fitness = 0.0
                currentfit = 0.0
                maxfit = 0.0

                #looping through lineages in the post transfer file
                for lineage in open('post_transfer.csv'):
                
                    #stripping off new line from lineage data
                    lineage = lineage.rstrip('\n')

                    #adding up counts to determine the total population size
                    total_pop_size = total_pop_size + int(lineage.split(',')[2])

                    #adding up fitness as the multiple of individual fitness and the number of individuals
                    total_fitness = total_fitness + float(lineage.split(',')[1])*int(lineage.split(',')[2])

                    #setting current fitness which will be used to identify the best individual lineage in hte population
                    currentfit = float(lineage.split(',')[1])


                    #determining the best organism at the time
                    if currentfit > maxfit:
                        maxfit = currentfit


                print total_pop_size
                print total_fitness

                average_fitness = total_fitness/total_pop_size

                print average_fitness
                
                current_generation += 1
                print current_generation

               


        
        #pop_stat_out = open("population_statistics.csv",'a')

        #subpop = subpop.rstrip('\n')
                
                #pulling out the different subpopulation variables
        #        subpop_id = subpop.split(',')[0]
        #        subpop_fitness = float(subpop.split(',')[1])
        #        subpop_count = int(subpop.split(',')[2])

                #generation = open('post_kill.csv','w')

                #count up total fitness, total population size, and then average fitness


                #print subpop_count
                #if subpop_count > 0:

exp_evol = Population(500000,10,10,100,24,1,0.0006,0.00015)

exp_evol.initialize()

exp_evol.update_time()
