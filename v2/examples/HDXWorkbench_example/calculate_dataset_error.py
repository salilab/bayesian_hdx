import hxio
from matplotlib import pyplot as plt
# Script for plotting timepoint avg vs. std from a HDXWorkbench file

datafiles = ["Data_Export_HDX_Workbench.csv"]


error_dict = {}

i = 0
pctd_error = []
d2_error = []
for df in datafiles:
    datasets = hxio.import_HDXWorkbench(df)
    for d in datasets:
        j = 0
        d_error = {}
        for pep in d.get_peptides():
            #print(pep.sequence)
            nam = pep.get_number_of_observable_amides()
            for tp in pep.get_timepoints():
                (avg, sd) = tp.get_avg_sd()
                #print(tp.time, avg, sd)
                d_error[j] = {"sd" : sd, "avg" : avg, "amides" : nam, "time" : tp.time}
                if sd is not None and avg is not None:
                    pctd_error.append((avg, sd))
                    d2_error.append(((avg*nam)/100, sd*(avg)/100))

        error_dict[i] = d_error


# Error in percent D
#plt.scatter([x[0] for x in pctd_error], [y[1] for y in pctd_error])
plt.scatter([x[0] for x in d2_error], [y[1] for y in d2_error])

plt.show()


