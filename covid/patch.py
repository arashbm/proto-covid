import enum
from collections import namedtuple

import numpy as np
import pandas as pd

class Compartment(enum.Enum):
    SUSCEPTIBLE = enum.auto()
    EXPOSED = enum.auto()
    ASYMPTOMATIC = enum.auto()
    PRESYMPTOMATIC = enum.auto()
    INFECTED = enum.auto()
    HOSPITALIZED = enum.auto()
    DEAD = enum.auto()
    RECOVERED = enum.auto()

class AgeGroup(enum.Enum):
    YOUNG = enum.auto()
    ADULTS = enum.auto()
    ELDERLY = enum.auto()

Config = namedtuple('Config',
        ['pi', 'eta', 'alpha', 'theta', 'nu', 'rho', 'chi', 'delta',
            'beta_infected', 'beta_presymptomatic', 'beta_asymptomatic',
            'kappa', 'contact'],
        defaults=[pd.DataFrame(0.0, columns=AgeGroup, index=AgeGroup)])

class Patch:
    def __init__(self, population):
        self.__population = pd.DataFrame(0.0,
                columns=AgeGroup, index=Compartment)

        for g in AgeGroup:
            for c in Compartment:
                self.__population[g][c] = population[g][c]

    def delta_population(self, c: Config) -> pd.DataFrame:
        delta = pd.DataFrame(0.0, columns=AgeGroup, index=Compartment)

        for g in AgeGroup:
            lambdaa = self.__force_of_infection(g, c)

            delta[g][Compartment.SUSCEPTIBLE] = \
                -lambdaa*self.__population[g][Compartment.SUSCEPTIBLE]

            delta[g][Compartment.EXPOSED] = \
                +lambdaa*self.__population[g][Compartment.SUSCEPTIBLE] \
                -c.eta*self.__population[g][Compartment.EXPOSED]

            delta[g][Compartment.ASYMPTOMATIC] = \
                +(1-c.pi)*c.eta*self.__population[g][Compartment.EXPOSED] \
                -c.rho*self.__population[g][Compartment.ASYMPTOMATIC]

            delta[g][Compartment.PRESYMPTOMATIC] = \
                +c.pi*c.eta*self.__population[g][Compartment.EXPOSED] \
                -c.alpha*self.__population[g][Compartment.PRESYMPTOMATIC]

            delta[g][Compartment.INFECTED] = \
                +c.alpha*self.__population[g][Compartment.PRESYMPTOMATIC] \
                -(c.theta+c.nu)*self.__population[g][Compartment.INFECTED]

            delta[g][Compartment.HOSPITALIZED] = \
                +c.theta*self.__population[g][Compartment.INFECTED] \
                -(c.chi+c.delta)*self.__population[g][Compartment.HOSPITALIZED]

            delta[g][Compartment.RECOVERED] = \
                +c.rho*self.__population[g][Compartment.ASYMPTOMATIC] \
                +c.nu*self.__population[g][Compartment.INFECTED] \
                +c.chi*self.__population[g][Compartment.HOSPITALIZED]

            delta[g][Compartment.DEAD] = \
                +c.delta*self.__population[g][Compartment.HOSPITALIZED]

        return delta

    def population(self):
        return self.__population

    def apply_delta(self, delta):
        self.__population += delta

    def __contact_making_population(self, g: AgeGroup, c: Config):
        return sum(self.__population[g][[
            Compartment.SUSCEPTIBLE,
            Compartment.EXPOSED,
            Compartment.PRESYMPTOMATIC,
            Compartment.ASYMPTOMATIC,
            Compartment.RECOVERED,
            Compartment.INFECTED
            ]]) + c.kappa*self.__population[g][Compartment.INFECTED]

    def __force_of_infection(self, g: AgeGroup, c: Config) -> float:
        lambdaa = 0.0
        contact = self.__contact_function(c)

        for h in AgeGroup:
            x_h = self.__contact_making_population(h, c)
            a_h = self.__population[g][Compartment.ASYMPTOMATIC]
            p_h = self.__population[g][Compartment.PRESYMPTOMATIC]
            i_h = self.__population[g][Compartment.INFECTED]
            lambdaa += contact[g][h]*(
                    c.beta_asymptomatic*a_h/x_h +
                    c.beta_presymptomatic*p_h/x_h +
                    c.beta_infected*c.kappa*i_h/x_h)

        return lambdaa

    def __contact_function(self, c: Config) -> pd.DataFrame:
        contact = pd.DataFrame(0.0, columns=AgeGroup, index=AgeGroup)

        for g in AgeGroup:
            x_g = self.__contact_making_population(g, c)
            for h in AgeGroup:
                x_h = self.__contact_making_population(h, c)
                contact[g][h] = min(
                        c.contact[g][h]*x_g,
                        c.contact[h][g]*x_h)/x_g
        return contact
