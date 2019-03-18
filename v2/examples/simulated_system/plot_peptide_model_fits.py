import sys
sys.path.append( "../../pyext/src" )
import plots
import numpy
import math

sigma=0.1

model_deuteration = {
    10 : [0.2, 0.34, 0.3, 0.24, 0.24],
    30 : [0.33, 0.33, 0.5, 0.44, 0.36],
    90 : [0.5, 0.4, 0.51, 0.52, 0.55, 0.61],
    300 : [0.65, 0.42, 0.7, 0.71, 0.5, 0.66],
    900 : [0.68, 0.54, 0.72, 0.66, 0.63, 0.67],
    3600 : [0.9, 0.85, 0.83, 0.96, 0.85, 0.87]
}

exp_deuteration = {
    10 : [23, 25, 27],
    30 : [35, 37, 39],
    90 : [48, 55, 58],
    300 : [62, 67, 70],
    900 : [68, 69, 70],
    3600 : [80, 81, 85]
}

'''
exp_deuteration = {
    10 : [0.02],
    30 : [0.05],
    90 : [0.5],
    300 : [0.3],
    900 : [0.4],
    3600 : [0.5]
}

model_deuteration = {
    10 : [0.02, 0.02002, 0.2001, 0.2001, 0.1999, 0.19999],
    30 : [0.05, 0.05],
    90 : [0.5, 0.51],
    300 : [0.3, 0.31],
    900 : [0.4, 0.41],
    3600 : [0.5, 0.51]
}
'''

model_deuteration = {
    10 : [99, 96, 90, 92, 65],
    30 : [99, 66, 55, 88],
    90 : [99],
    300 : [30, 31],
    900 : [40, 41],
    3600 : [50, 51]
}
fig, ax = plots.plot_incorporation_curve(model_deuteration, exp_deuteration, plot=True)

fig.savefig("output.png")

sumstat = 0
for i in model_deuteration.keys():
    models = model_deuteration[i]
    exp = exp_deuteration[i]

    welchs_t = (numpy.average(models)- numpy.average(exp)) / math.sqrt(numpy.var(models)/len(models)+sigma**2/len(exp))

    print i, welchs_t

    sumstat+=welchs_t

print sumstat / len(model_deuteration.keys())


