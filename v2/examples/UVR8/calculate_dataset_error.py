import hxio
from matplotlib import pyplot as plt


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
            nam = pep.get_number_of_observable_amides()
            for tp in pep.get_timepoints():
                (avg, sd) = tp.get_avg_sd()
                d_error[j] = {"sd" : sd, "avg" : avg, "amides" : nam, "time" : tp.time}
                pctd_error.append((avg, sd))
                d2_error.append((avg*nam, sd))

        error[i] = d_errors



plt.scatter([x[0] for x in pctd_error], [y[1] for y in pctd_error])

plt.show()


