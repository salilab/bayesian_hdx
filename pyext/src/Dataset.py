


class Timepoint(object):
    '''
    Timepoint objects are bound to a specific Fragment, where Timepoint represents
    a regular observation time of the HDX experiment.
    The class contains both the experimental data (as a list of Replicate objects)
    as well as a list of calculated values from the forward model (self.models)
    '''
    def __init__(self, time, sigma0=5.0):
        '''
        @param time - Time in seconds
        @param sigma0 - Initial estimate of timepoint error sigma in pctD units. 
        '''
        self.time = time
        self.models = []
        self.replicates = []
        self.sigma = sigma0

    def number_of_replicates(self):
        return len(self.replicates)

    def add_replicate(self, deut, sat=1.0, recovery=1.0, temp=293.15):
        self.replicates.append(Replicate(deut, sat, recovery, temp))

    def calculate_tp_score(self, freq_grid, exp_grid, sig, num_amides, force_calc=False):
        """ given the exp_grid of exponential rates and number of observations
        in each grid, calculates the Bayesian score """
        score = 0
        for r in self.replicates:
            if force_calc:
                self.calc_model_deut(freq_grid, exp_grid, num_amides)
                rep_score = r.calculate_replicate_score(self.model_deut, sig)
            else:               
                try:
                    rep_score = r.calculate_replicate_score(self.model_deut, sig)
                except:
                    self.calc_model_deut(freq_grid, exp_grid, num_amides)
                    rep_score = r.calculate_replicate_score(self.model_deut, sig)

            # This if statement prevents log overflows when score ~0
            if r.calculate_replicate_score(self.model_deut, sig) == 0.0:
                score += 10000 
            else:
                score += -1.0*numpy.log(rep_score)
        return score / len(self.replicates)


    def calc_model_deut(self, freq_grid, exp_grid, num_amides):
        # Given an exp_frequency grid and the exp_grid, calculate the deuteration of the model
        # at this timepoint. 
        deut=0
        for n in range(len(exp_grid)):
            #print(n, freq_grid)
            exchange_rate = 10**exp_grid[n]
            num_amides_at_this_rate = freq_grid[n]
            deut += num_amides_at_this_rate*(1-numpy.exp(-exchange_rate*self.time))
        self.model_deut = deut / num_amides * 100
        #print(self.time, grid, exp_grid, self.model_deut)
        return self.model_deut

    def get_avg_sd(self):
        if len(self.replicates) < 1:
            #print "No replicates in timepoint, ", self.time
            self.avg=None
            self.sd=None
        elif len(self.replicates) <= 2:
            #print "Not enough replicates, ", self.time
            self.avg=sum(r.deut for r in self.replicates)/len(self.replicates)
            self.sd=None
        else:
            sum_deut=0
            sumsq_deut=0
            for r in self.replicates:
                sum_deut=sum_deut+float(r.deut)
                sumsq_deut=sumsq_deut+r.deut**2
            self.avg=sum_deut/self.num_replicates
            self.sd=numpy.sqrt((self.num_replicates*sumsq_deut-sum_deut**2)/self.num_replicates**2)
        return self.avg, self.sd

    def get_model_avg(self):
        if len(self.models) < 1:
            print("No models imported into timepoint ", self.time)
            self.model_avg=None
            return None
        else:
            #print "Not enough replicates, ", self.time
            self.model_avg=numpy.average(self.models)
        return self.model_avg

    def get_model_sd(self):
        if len(self.models) <= 2:
            print("Not enough models to calculate SD ", self.time)
            self.model_sd=None
        else:
            #print "Not enough replicates, ", self.time
            self.model_sd=numpy.std(self.models)
        return self.model_sd

    def add_model_deuteration_value(self, deut):
        self.models.append(deut)

    def clear_model_deuteration_values(self):
        self.models=[]

    def remove_outliers(self, sigma, numsig=3):
        mean,sd=self.get_avg_sd()
        outliers=0
        for d in self.deut:
            if abs(d-mean)/sigma>numsig:
                self.deut.remove(d)
                print("outlier removed: ", d, mean, sigma)
                outliers=outliers+1
                self.num_replicates=self.num_replicates-1
        return outliers

    def set_sigma(self, sigma):
        '''
        @param sigma - float to set timepoint sigma; or "exp_sd" to use SD from 
             experimental data.
        '''
        if sigma=="exp_sd":
            if self.num_replicates > 2:
                avg, self.sigma=self.get_avg_sd()
            else:
                raise Exception("Not enough models to calculate SD")
        else:
            self.sigma=sigma


class Replicate(object):
    """ Replicate objects store the individual data observations
        @param deut - deuteration level in percentD
        @param saturation - the deuterium concentration
        @param recovery - the 2D recovery estimate
    """
    def __init__(self, deut, deuterium_concentration=1.0, recovery=1.0):
        self.deut = deut
        self.sat = sat
        self.nearest_gridpoint = -1
        self.recovery = recovery

    def calculate_replicate_score(self, model_deut, sig):
        return self.gaussian_model(self.deut, model_deut, sig)

    def lognormal_model(self,exp,model,sig):
        return numpy.exp(-( (numpy.log(model)-numpy.log(exp) )**2)/(2*sig**2))/(sig*numpy.sqrt(numpy.pi)*model)

    def gaussian_model(self,exp,model,sig):
        return numpy.exp(-((model-exp)**2)/(2*sig**2))/(sig*numpy.sqrt(2*numpy.pi))

    def calc_nearest_gridpoint(self, grid):
        self.nearest_gridpoint = (numpy.abs(grid-self.deut)).argmin()

    def get_nearest_gridpoint(self, grid):
        if self.nearest_gridpoint==-1:
            self.calc_nearest_gridpoint(grid)

        return self.nearest_gridpoint