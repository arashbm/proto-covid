import sys
import os
import csv
import logging

import pandas as pd
import numpy as np

from covid import Config, Compartment, AgeGroup, Patch

# returns a dictionary where keys are region names and values are patches
def create_region_patches(population_filename):
    patches = {}

    # Read region population data:
    with open(population_filename, newline='') as csvfile:
        reader = csv.reader(csvfile)

        next(reader) # first row is the header

        for row in reader:
            name = row[0]

            young = sum(float(s) for s in row[2:9])
            adults = sum(float(s) for s in row[9:17])
            elderly = sum(float(s) for s in row[17:22])

            area = row[22] # surfce area, not used for now

            # for simplicity we assume everyone has equal chance of being
            # infected (and symptomatic) and twice that of being asymptomatic


            presymptomatic_p = 4000/5.51e6
            infection_p = 2000/5.51e6
            asymptomatic_p = 4000/5.51e6
            hospitalized_p = 200/5.51e6
            dead_p = 34/5.51e6
            susceptible_p = 1.0 - presymptomatic_p - infection_p \
                    - asymptomatic_p - hospitalized_p - dead_p

            # population can be anything (dict, dataframe etc.) as long as
            # something like p[AgeGroup.YOUNG][Compartment.SUSCEPTIBLE] is
            # valid for all age/compartment combinations
            pop = pd.DataFrame(data=0.0, columns=AgeGroup, index=Compartment)

            for g, group_pop in [(AgeGroup.YOUNG, young),
                                    (AgeGroup.ADULTS, adults),
                                    (AgeGroup.ELDERLY, elderly)]:
                pop[g][Compartment.SUSCEPTIBLE] = group_pop*susceptible_p
                pop[g][Compartment.PRESYMPTOMATIC] = group_pop*presymptomatic_p
                pop[g][Compartment.INFECTED] = group_pop*infection_p
                pop[g][Compartment.ASYMPTOMATIC] = group_pop*asymptomatic_p
                pop[g][Compartment.HOSPITALIZED] = group_pop*hospitalized_p
                pop[g][Compartment.DEAD] = group_pop*dead_p

            # Each row is a region. We create a patch for each
            patches[name] = Patch(population=pop)

    logging.info(f"Created {len(patches)} patches.")
    return patches



if __name__ == '__main__':
    contacts = pd.DataFrame([[0.6, 0.2, 0.2],
                             [0.2, 0.6, 0.2],
                             [0.2, 0.2, 0.6]],
                            columns=AgeGroup, index=AgeGroup)

    conf = Config(
            pi=0.5, eta=1/2.34, theta=0.05, nu=1/2.86, alpha=1/3.00,
            rho=1/2.0, chi=1/10.0, delta=1/7.0, beta_infected=0.06,
            beta_asymptomatic=0.06, beta_presymptomatic=0.06,
            kappa=0.4, contact=contacts)

    if os.environ.get('VERBOSE', 'FALSE').upper() == 'TRUE':
        logging.basicConfig(level='DEBUG')
    else:
        logging.basicConfig(level='INFO')

    patches = create_region_patches('data/regions.csv')

    regions = list(patches)

    random_mobility = pd.DataFrame(
            0.0, columns=regions, index=regions)

    for i in regions:
        pop_i = patches[i].population().values.sum()
        for j in regions:
            pop_j = patches[j].population().values.sum()
            random_mobility[i][j] =  pop_i*pop_j if i!= j else 0.0

    trips_per_day = 5e5 # for sake of argument
    random_mobility = trips_per_day*random_mobility/random_mobility.values.sum()

    t = 0.0  # time in days
    t_max = 100
    dt = 0.2
    print("t total_hospitalized total_infected total_dead")
    while t < t_max:
        t += dt

        total_hospitalized = 0.0
        total_infected = 0.0
        total_dead = 0.0

        for i in patches:
            logging.debug(f"delta for i={i}")

            delta = patches[i].delta_population(conf)

            pop_i = patches[i].population()
            pop_i_norm = pop_i/pop_i.values.sum()

            for j in patches:
                pop_j = patches[j].population()
                pop_j_norm = pop_j/pop_j.values.sum()
                delta += random_mobility[j][i]*pop_j_norm
                delta -= random_mobility[j][i]*pop_i_norm

            patches[i].apply_delta(delta*dt)

            for g in AgeGroup:
                total_hospitalized += pop_i[g][Compartment.HOSPITALIZED]
                total_infected += pop_i[g][Compartment.INFECTED]
                total_dead += pop_i[g][Compartment.DEAD]

        print(t, total_hospitalized, total_infected, total_dead)
