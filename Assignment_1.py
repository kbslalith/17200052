import math
import random
import numpy as np
import scipy.stats
import random


def f(x):
    """The objective is defined as the cost + a per-demographic penalty
    for each demographic not reached."""

    n = len(x)
    assert n == n_venues

    # initially we have reached no categories
    reached = np.zeros(n_demographics, dtype=int)
    # and incurred no cost
    cost = 0.0

    # we look at each bit in x
    for xi, ri, ci in zip(x, r, c):
        if xi:
            # xi = 1, so we reach some categories and incur a cost
            reached = reached | ri #
            cost += ci # incur the cost

    for ri, pi in zip(reached, p):
        # for each demographic
        if ri == 0:
            # failed to reach a demographic? pay a penalty specific to that demographic
            cost += pi
    return cost

def nbr(x):
    """Generate a neighbour in a bitstring space"""
    x = x[:] # make a copy, so that we don't overwrite the original
    i = random.choice(range(len(x))) # choose an index in the bitstring
    x[i] = 1 - x[i] # "flip" the bit at that index, from 0 to 1 or 1 to 0.
    return x

def init():
    """Generate a random point in a bitstring space"""
    return [random.choice([0, 1]) for i in range(len(c))]


########################################################################################
########################## Definitions for various Algorithms ##########################

# Hill CLimbing method for an advertising problem

def hill_climb(f, nbr, init, maxits, direction="max"):
    """ Hill Climbing Algorithm"""
    assert direction in ("min", "max")
    x = init() # generate an initial candidate solution
    fx = f(x) #generate the cost for initial candidate solution

    #Check extreme points but comment above x and fx
    #x = [1 for i in range(30)]   # check for all [1,1,1,....,1]
    #x = [0 for i in range(30)]   # check for all [0,0,0,....,0]
    #print(x)
    #fx = f(x)

    for i in range(maxits):
        xdash = nbr(x) # generate a neighbour of x
        fxdash = f(xdash) # generate cost for the neighbour
        if ((direction == "max" and fxdash > fx)
            or (direction == "min" and fxdash < fx)):
            # if the neighbour is better, "go there", else stay at x
            x = xdash
            fx = fxdash
        #print(fx)
    return x, fx

# Simulated Annealing method for an advertising problem

def anneal(f,nbr,init,maxits):
    """ Simulated Annealing Algorithm"""
    x = init()                          # generate an initial candidate solution
    fx = f(x)                           # generate cost
    T = 1.0                             # set the initial value for the control parameter
    T_min = 0.0001                      # set minimum value to stop the loop
    alpha = 0.99                        # rate at which control parameter value decrease

    for i in range(maxits):
        if T > T_min:
            xdash = nbr(x)              # generate neighbour
            fxdash = f(xdash)
            if (fxdash < fx) or (math.exp((fx-fxdash)/T) > random.random()):
                # if new value is less or acceptance condition is true
                #move to a new value
                x = xdash
                fx = fxdash
        T = T*alpha
    return x, fx

# Iterated Hill Climbing method for an advertising problem

def iterated_hill_climb(f, nbr, init, maxits, direction="max"):
    """ Iterated Hill CLimbing Algorithm"""
    assert direction in ("min", "max")

    n_restarts = 50                    # Total number of restarts

    best = init()                       # initial best value randomly selected
    fbest = f(best)

    for j in range(n_restarts):

        x = init()                      # generate an initial random solution
        fx = f(x)

        for i in range(maxits // n_restarts):
            xdash = nbr(x)              # generate a neighbour of x
            fxdash = f(xdash)
            if ((direction == "max" and fxdash > fx)
                or (direction == "min" and fxdash < fx)):
            # if the neighbour is better, "go there", else stay at x
                x = xdash
                fx = fxdash
            if ((direction == "max" and fx > fbest)
                or (direction == "min" and fx < fbest)):
            # if the neighbour is the best-ever, save it
                best = x
                fbest = fx
    return best, fbest


# Late Acceptance Hill Climbing for an advertising problem

def LAHC(f,nbr,init,maxits,L):
    """ Late Acceptance Hill Climbing Algorithm"""
    x = init()                          # generate an initial candidate solution
    fx = f(x)                           # find the cost of the solution
    best = x                            # keep best as the initial candidate
    fbest = f(best)                     # find the cost of the best candidate
    ls = [fx]*L                         # Create a list of length L

    for i in range(maxits):
        xdash = nbr(x)
        fxdash = f(xdash)
        if fxdash < fx:
            best = xdash
            fbest = fxdash
        v = i % L                       # find the index of the list
        if fxdash < ls[v] or fxdash <= fx:
            x = xdash
            fx = fxdash
        ls[v] = fx
    return best, fbest



r = np.genfromtxt("reach.dat").astype(int)# Read data from Reach.dat
c = np.genfromtxt("cost.dat")             # Read data from Cost.dat
p = np.genfromtxt("penalty.dat")          # Read data from Penalty.dat

n_venues = c.shape[0]                     # Number of venues
n_demographics = p.shape[0]               # Number of demographics to reach
assert n_venues == r.shape[0]
assert n_demographics == r.shape[1]

maxits = 5000                             # maximum number of iterations

#x_, fx_ = hill_climb(f, nbr, init, maxits, direction="min")            # Call Hill climbing algorithm
#x_, fx_ = anneal(f, nbr, init, maxits)                                 # Call Simulated Annealing algorithm
#x_, fx_ = iterated_hill_climb(f, nbr, init, maxits, direction="min")   # Call Iterated Hill Climbing
#x_, fx_ = LAHC(f, nbr, init, maxits,L=20)                              # Call LAHC algorithm

print(x_)
print(fx_)
