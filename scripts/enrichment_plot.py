import numpy as np
import matplotlib.pyplot as plt

ip, control = [np.load(f) for f in snakemake.input]
enrichment = (ip/np.sum(ip))/(control/np.sum(control))
plt.plot(enrichment)
plt.vlines([124, 1124], ymin=0.5, ymax=1.5)
plt.hlines([1], xmin=0, xmax=enrichment.size)
plt.savefig(snakemake.output[0])
