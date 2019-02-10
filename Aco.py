
import numpy as np
from scipy.stats import norm
import Functions
class ACO:

    def __init__(self):

        self.verbosity = True

        # Initial algorithm parameters
        self.maximum_iteration = 100
        self.population_size = 5
        self.s = 50
        self.t = 0.1
        self.u = 0.85

        # Initial (NULL) problem definition
        self.num_variable = 2
        self.variable_ranges = [[0, 1],
                           [0, 1]]
        self.cost_function = None

        # Optimization results
        self.SA = None
        self.best_solution = None
        # end def


    def variables(self, nvariable, ranges):

        if len(ranges) != nvariable:
            print "Error, Ranges of variables and numbers does not match"
        else:
            self.num_variable = nvariable
            self.variable_ranges = ranges
            self.SA = np.zeros((self.s, self.num_variable))
    # end def


    def cost(self, costf):

        self.cost_function = costf
    # end def


    def parameters(self, maximum_iteration, population_size, s, t, u):

        self.maximum_iteration = maximum_iteration
        self.population_size = population_size
        self.s = s
        self.t = t
        self.u = u
    # end def


    def verbosity(self, status):

        if type(status) is bool:
            self.verbosity = status
        else:
            print "Error, received verbosity parameter is not boolean"
    # end def


    def _biased_selection(self, probabilities):

        r = np.random.uniform(0, sum(probabilities))
        for p, f in enumerate(probabilities):
            r -= f
            if r <= 0:
                return p
    # end def


    def optimize(self):

        if self.num_var == 0:
            print "Error, first set the number of variables and their boundaries"
        elif self.cost_function == None:
            print "Error, first define the function cost  to be used"
        else:

            if self.verbosity:   print "Initialization to the achrive"
            # Initialize the archive by random sampling, respecting each variable's constraints
            pop = np.zeros((self.population_size, self.num_var))
            w = np.zeros(self.s)

            for i in xrange(self.s):
                for j in xrange(self.num_variable):
                    self.SA[i, j] = np.random.uniform(self.var_ranges[j][0], self.var_ranges[j][1])
                self.SA[i, -1] = self.cost_function(self.SA[i, 0:self.num_var])
            self.SA = self.SA[self.SA[:, -1].argsort()]

            x = np.linspace(1 ,self.s ,self.k)
            w = norm.pdf(x ,1 ,self. t *self.k)
            p = w/ sum(w)

            if self.verbosity:   print
            "Algo for main loop"

            # Algorithm runs until it reaches maximum number of iterations
            for iteration in xrange(self.maximum_iter):
                if self.verbosity:
                    print
                    "[%d]" % iteration
                    print
                    self.SA[0, :]

                Mi = self.SA[:, 0:self.num_variable]
                for ant in xrange(self.population_size):
                    l = self._biased_selection(p)

                    for var in xrange(self.num_variable):
                        sigma_sum = 0
                        for i in xrange(self.k):
                            sigma_sum += abs(self.SA[i, var] - self.SA[l, var])
                        sigma = self.xi * (sigma_sum / (self.k - 1))

                        pop[ant, var] = np.random.normal(Mi[l, var], sigma)

                        # Deals with search space violation using the random position strategy
                        if pop[ant, var] < self.var_ranges[var][0] or pop[ant, var] > self.var_ranges[var][1]:
                            pop[ant, var] = np.random.uniform(self.var_ranges[var][0], self.var_ranges[var][1])

                    pop[ant, -1] = self.cost_function(pop[ant, 0:self.num_var])

                self.SA = np.append(self.SA, pop, axis=0)
                self.SA = self.SA[self.SA[:, -1].argsort()]
                self.SA = self.SA[0:self.s, :]

            self.best_solution = self.SA[0, :]
            # sol=Functions.Sphere(self.best_solution)
            # print "cost"
            # print sol
            return self.best_solution
            # end def

# end class 
