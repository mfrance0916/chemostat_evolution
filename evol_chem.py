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

        selective_effect = np.random.exponential(scale=0.02) - 0.075 

    return selective_effect

class Population(object):
    
    ## defines the following parameters that are needed for the simulation
    ## pop_dens = density of bacterial population in cells/ml
    ## n_generations = numberical number of generations to run the simulation for
    ## gd_mut_rate = growth dependent mutation rate
    ## gi_mut_rate = growth independent mutation rate
    ## volume = volume of the chemostat in mL
    ## flow_rate = rate of flow of new media into the chemostat
    ## gd_ab_mut_rate = rate of mutation towards antibiotic resistance dependent on growth
    ## gi_ab_mut_rate = rate of mutation towards antibiotic resistance independent of growth

    def __init__(self, pop_dens, n_generations, gd_mut_rate, gi_mut_rate, volume, flow_rate, gd_ab_mut_rate, gi_ab_mut_rate):
        self.pop_dens = pop_dens
        self. n_generations = n_generations
        self.gd_mut_rate = gd_mut_rate
        self.gi_mut_rate = gi_mut_rate
        self.volume = volume
        self.flow_rate = flow_rate
        self.gd_ab_mut_rate = gd_ab_mut_rate
        self.gi_ab_mut_rate = gi_ab_mut_rate

    def initialize(self):
        
        #calculate population size
        self.pop_size = self.pop_dens*self.volume

        #create list of individuals with a relative fitness of 1, turns out it's too memory intensive, because competition is global we don't really need to track individual cells anyway, instead we will track the subpopulations
        #self.fitness = list([1.0]*self.pop_size)    
        #self.ab_status = list([1]*self.pop_size)

        #creating a csv file that will keep track of all of the subpopulations

        subpopulations = open('subpopulations.csv','w')
        subpopulations.write("1,1.0,%s,0" %(self.pop_size))
        subpopulations.close()
        
        #calculating the generation time as the natural log of 2 over the dilution rate (flow over volume)
        self.generation_time = np.log(2.0)/(self.flow_rate/self.volume)
        
        #creating a counting variable for the time in hours
        self.time = 0
        self.stop_time = self.n_generations * self.generation_time

        #print self.generation_time
        #print self.time

        #determining the number to kill each timestep
        self.survive_prob = 1.0 - ((self.pop_size*(self.flow_rate/self.volume))/self.pop_size)
        #print self.survive_prob

        self.subpop_id_num = 1

        #opening file to store output data
        pop_stat_out = open("population_statistics.csv",'w')
        pop_stat_out.write('time,gen,pop_size,num_lineages,ave_fit,num_abres,abund_abres,ave_fit_abres,shannon_diversity\n')
        pop_stat_out.write('0,0,' + str(self.pop_size) + ',1,1.0,0,0.0,0,0\n') 
        pop_stat_out.close()



    def update_time(self):
        pop_stat_out = open("population_statistics.csv",'a')

        #looping through until time is greater than or equal to the stop time as determiend by the generation time adn the number of generations requested

        while self.time < self.stop_time:
            
            self.time += 1
            #print self.time
            
            #ENACTING THE KILL THROUGH RANDOM DEATH PROCESS using a draw from a binomial distribution and tallying up the remaining fitness
            #looping through subpopulations file

            post_kill = open('post_kill.csv','w')

            num_kills = 0

            for subpop in open('subpopulations.csv'):
                
                
                subpop = subpop.rstrip('\n') 
                
                #pulling out the different subpopulation variables
                subpop_id = subpop.split(',')[0]
                subpop_fitness = float(subpop.split(',')[1])
                subpop_count = int(subpop.split(',')[2])
                subpop_abstatus = subpop.split(',')[3]

                #print subpop_count

                
                if subpop_count > 0:
                    #random death process simulated on the subpopulation as a draw from a binomial distribution with p = probability of surviving and n = subpopulation count
                    num_survivors = np.random.binomial(subpop_count,self.survive_prob)
                    
                    #incrimenting the total number of kills
                    num_kills = num_kills + (subpop_count - num_survivors)

                    #writing out to the kill file
                    if num_survivors > 0:    
                        post_kill.write(subpop_id + ',' + str(subpop_fitness) + ',' + str(num_survivors) + ',' + subpop_abstatus + '\n')

            
            #ENACTING MUTATIONAL PROCESS THAT IS INDEPENDENT OF GROWTH

            #creating file that will contain the output of this step
            post_kill.close()
            post_gi_mutation = open('post_gi_mutation.csv','w')

            for subpop_post in open('post_kill.csv'):

                subpop_post = subpop_post.rstrip('\n') 

                #print subpop_post

                #storing variables about the ancestral lineages for future use
                subpop_id = subpop_post.split(',')[0]
                subpop_fitness = float(subpop_post.split(',')[1])
                subpop_count = int(subpop_post.split(',')[2])
                subpop_abstatus = subpop_post.split(',')[3]


                #determining the number of new growth independent mutants to generate
                new_gi_mutants = np.random.binomial(subpop_count,self.gi_mut_rate)

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

                    #new mutant inherits the ancestors antibiotic status
                    new_mutant_abstatus = subpop_abstatus

                    #writing out the information for the new mutants to a temporary file
                    post_gi_mutation.write(str(new_mutant_id) + ',' + str(new_mutant_fitness) + ',' + str(new_mutant_count) + ',' + new_mutant_abstatus + '\n')


                #determining the number of new ANTIBIOTIC resistant mutants generated by GI mutations
                #only accounting for forward mutations i.e. resistant things can't generate new resistant things via mutation
                if subpop_abstatus == '0':

                    #new mutants drawn from a binomial with mutation rate as prob and number of trials equal to the subpop size
                    new_gi_ab_mutants = np.random.binomial(subpop_count,self.gi_ab_mut_rate)
                    #print new_gi_ab_mutants

                    #adjusting the count for the ancestral population
                    subpop_count = subpop_count - new_gi_ab_mutants
    
                    #only adding a line to the file if the binomial draw came up with new mutants
                    if new_gi_ab_mutants > 0:

                        #incrementing the lineage tracking number by 1
                        self.subpop_id_num += 1

                        #assigning the new number to the ab resistant population
                        new_mutant_id = self.subpop_id_num

                        #writing out the new lineages data into the file, because all mutants are identical in this model they can 
                        post_gi_mutation.write(str(new_mutant_id) + ',' + str(subpop_fitness) + ',' + str(new_gi_ab_mutants) + ',' + '1' + '\n')


                #writing the original lineages data to the file
                post_gi_mutation.write(str(subpop_id) + ',' + str(subpop_fitness) + ',' + str(subpop_count) + ',' + str(subpop_abstatus) + '\n')


            post_gi_mutation.close()

            

            #CALCULATING TOTAL POPULATION SIZE, TOTAL FITNESS AND AVERAGE FITNESS and the MAXIMUM FITNESS
            total_fitness = 0.0
            total_pop_size = 0

            maxfit = 0.0


            for lineage in open('post_gi_mutation.csv'):
                
                lineage = lineage.rstrip('\n') 

                total_pop_size = total_pop_size + int(lineage.split(',')[2])

                total_fitness = total_fitness + float(lineage.split(',')[1])*int(lineage.split(',')[2])

                average_fitness = total_fitness/total_pop_size

                #print average_fitness

                currentfit = float(lineage.split(',')[1])

                #determining the best organism at the time
                if currentfit > maxfit:
                    maxfit = currentfit


            #ENACTING GROWTH AND MUTATION PROCESS THAT IS DEPENDENT ON GROWTH
            os.remove('subpopulations.csv')

            post_growth = open('subpopulations.csv','w')

            for lineage in open('post_gi_mutation.csv'):

                lineage = lineage.rstrip('\n')
                
                #storing variables for use
                subpop_id = lineage.split(',')[0]
                subpop_fitness = float(lineage.split(',')[1])
                subpop_count = int(lineage.split(',')[2])
                subpop_abstatus = lineage.split(',')[3]

                
                #SELECTION THINGS WITH HIGHER FITNESS DO BETTER

                #calculating the relative fitness as the division of the lineage fitness by the average fitness
                rel_fitness = subpop_fitness/average_fitness
                rel_abundance = float(subpop_count)/float(total_pop_size)
                

                #calculating the probability to divide as the multiplication of the number of kills divided by the total population size times the relative fitness
                prob_divide = rel_fitness*float(num_kills)/total_pop_size

                if prob_divide < 0.0:
                    print prob_divide
                    prob_divide = 0.0

                #determining the number of cellular divisions generated by the lineage
                divisions = np.random.binomial(subpop_count,prob_divide)

                subpop_count = subpop_count + divisions

                if divisions > 0:
                    daughters = divisions * 2

                    #GROWTH DEPENDENT MUTATION
                    new_gd_mutants = np.random.binomial(daughters,self.gd_mut_rate)

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
                            print('gd  ' + str(new_mutant_fitness))
                            new_mutant_fitness = 0.0

                        #setting the count of the new mutant equal to 1
                        new_mutant_count = 1
    
                        #new mutant inherits the ancestors antibiotic status
                        new_mutant_abstatus = subpop_abstatus

                        #writing out the data to the post growth file
                        post_growth.write(str(new_mutant_id) + ',' + str(new_mutant_fitness) + ',' + str(new_mutant_count) + ',' + new_mutant_abstatus + '\n')

                    if subpop_abstatus == '0':

                        #new mutants drawn from a binomial with mutation rate as prob and number of trials equal to the subpop size
                        new_gd_ab_mutants = np.random.binomial(subpop_count,self.gd_ab_mut_rate)
            
                        #adjusting the count for the ancestral population
                        subpop_count = subpop_count - new_gd_ab_mutants
    
                        #only adding a line to the file if the binomial draw came up with new mutants
                        if new_gd_ab_mutants > 0:

                            #incrementing the lineage tracking number by 1
                            self.subpop_id_num += 1

                            #assigning the new number to the ab resistant population
                            new_mutant_id = self.subpop_id_num

                            #writing out the new lineages data into the file, because all mutants are identical in this model they can 
                            post_growth.write(str(new_mutant_id) + ',' + str(subpop_fitness) + ',' + str(new_gd_ab_mutants) + ',' + '1' + '\n')

                #writing the original lineages data to the file
                post_growth.write(str(subpop_id) + ',' + str(subpop_fitness) + ',' + str(subpop_count) + ',' + str(subpop_abstatus) + '\n')

            post_growth.close()

            #storing output information for plotting later
            
            population_size = 0
            total_pop_fitness = 0.0
            abres_pop_size = 0
            abres_fitness = 0.0
            lineage_count = 0

            for lineage in open('subpopulations.csv'):

                lineage = lineage.rstrip('\n')

                lineage_count += 1

                population_size = population_size + int(lineage.split(',')[2])

                if lineage.split(',')[3] == '1':
                    abres_pop_size = abres_pop_size + int(lineage.split(',')[2])
                    abres_fitness = abres_fitness + float(lineage.split(',')[1])*int(lineage.split(',')[2])

                total_pop_fitness = total_pop_fitness + float(lineage.split(',')[1])*int(lineage.split(',')[2])

            average_pop_fitness = float(total_pop_fitness)/population_size
            abres_rel_abund = float(abres_pop_size)/population_size
            
            if abres_pop_size > 0:
                average_ab_fitness = float(abres_fitness)/abres_pop_size
            else:
                average_ab_fitness = 0.0

            shannon_diversity = 0.0
            
            for lineage in open('subpopulations.csv'):

                fin_rel_abund = float(lineage.split(',')[2])/population_size

                shannon_diversity = shannon_diversity + (fin_rel_abund*np.log(fin_rel_abund))

            shannon_diversity = shannon_diversity * -1

            pop_stat_out.write(str(self.time) + ',' + str(self.time/self.generation_time) + ',' + str(population_size) + ',' + str(lineage_count) + ',' + str(average_pop_fitness) + ',' + str(abres_pop_size) + ',' + str(abres_rel_abund) + ',' + str(average_ab_fitness) + ',' + str(shannon_diversity) + '\n')




exp_evol = Population(50000000,1000,0.0000005,0.0000005,10,1.5,0.00000006,0.000000006)

exp_evol.initialize()

exp_evol.update_time()

