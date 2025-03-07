import os
import json
import argparse
from matplotlib import pyplot as plt

def plot(figure_data, save):
    plt.clf()
    for plot_data in figure_data["plots"]:
        plt.plot(plot_data["xs"], plot_data["ys"], label=plot_data["label"])
    plt.title(figure_data["title"])
    plt.grid(True)
    plt.legend()
    if save:
        plt.savefig(f"""asset/{figure_data["title"]}.png""")
    else:
        plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser("plot")
    parser.add_argument("-s", "--save", help="Save plots as images instead of displaying in GUI", action="store_true")
    args = parser.parse_args()
    data_path = "data"
    for filename in os.listdir(data_path):
        with open(os.path.join(data_path, filename)) as f:
            data = json.load(f)
        plot(data, args.save)
