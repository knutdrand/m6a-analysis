import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from bdgtools.plotter import plot
from bdgtools.aggregateplot import MetaGenePlot
ip, control = [pd.read_pickle(f) for f in snakemake.input]
assert np.all(ip["x"]==control["x"])
ip["y"] /= control["y"]
plot(ip, MetaGenePlot)
ip.to_pickle(snakemake.output[0])
