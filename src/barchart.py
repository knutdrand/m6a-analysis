import matplotlib
import matplotlib.pyplot as plt
import numpy as np
def autolabel(rects, ax):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

def barchart(value_array, x_labels, y_labels):
    x = np.arange(len(x_labels))  # the label locations
    N = len(value_array)
    width = 0.9/N
    fig, ax = plt.subplots(figsize=(16.0, 10.0))
    rects = [ax.bar(x-width/2+i*width, values, width, label=label)
             for i, (values, label) in enumerate(zip(value_array, y_labels))]
    
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('N reads')
    ax.set_title('N reads by sample and type')
    ax.set_xticks(x)

    ax.set_xticklabels(x_labels)
    ax.legend()
    for r in rects:
        autolabel(r, ax)
    fig.tight_layout()

if __name__ == "__main__":
    import sys
    ylabels = sys.argv[1].split(",")
    xlabels = open(sys.argv[2]).read().strip().split()
    files = sys.argv[3:-1]
    values_array = [[int(v) for v in open(name).read().strip().split()]
                    for name in files]
    barchart(values_array, xlabels, ylabels)
    plt.savefig(sys.argv[-1])
