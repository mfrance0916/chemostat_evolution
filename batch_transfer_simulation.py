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
    ## batch length = amount of real time each transfer takes
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
        self.dilution = self.generations_per_day/

        #create list of individuals with a relative fitness of 1, turns out it's too memory intensive, because competition is global we don't really need to track individual cells anyway, instead we will track the subpopulations
        #self.fitness = list([1.0]*self.pop_size)    
        #self.ab_status = list([1]*self.pop_size)

        #creating a csv file that will keep track of all of the subpopulations

        subpopulations = open('subpopulations.csv','w')
        subpopulations.write("1,1.0,%s,0" %(self.pop_size))
        subpopulations.close()
        
        #calculating the generation time as the natural log of 2 over the dilution rate (flow over volume)
        self.current_generation = 0
        
        #creating a counting variable for the time in hours
        self.time = 0

        #determining the number to kill at each serial transfer 
        self.survive_prob = 
        #print self.survive_prob

        self.subpop_id_num = 1.0/dilution

        #opening file to store output data
        pop_stat_out = open("population_statistics.csv",'w')
        pop_stat_out.write('time,gen,pop_size,num_lineages,ave_fit,abund_lin\n')
        pop_stat_out.write('0,0,' + str(self.pop_size) + ',1,1.0,1\n') 
        pop_stat_out.close()

    def update_time(self):

        for generation_step in range(1,)
        
        pop_stat_out = open("population_statistics.csv",'a')

        subpop = subpop.rstrip('\n')
                
                #pulling out the different subpopulation variables
                subpop_id = subpop.split(',')[0]
                subpop_fitness = float(subpop.split(',')[1])
                subpop_count = int(subpop.split(',')[2])

                generation = open('post_kill.csv','w')

                #count up total fitness, total population size, and then average fitness


                #print subpop_count
                if subpop_count > 0:

