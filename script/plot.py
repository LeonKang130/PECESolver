import json
from matplotlib import pyplot as plt

def plot(data):
    plt.clf()
    ts = data["ts"]
    for i in range(len(data["ys"])):
        plt.plot(ts, data["ys"][i], label=f"y{i}")
    plt.legend()
    plt.grid(True)
    plt.title(data["title"])
    plt.show()

def plot_curve(data):
    plt.clf()
    y1, y2 = data["ys"][:2]
    plt.plot(y1, y2)
    plt.grid(True)
    plt.title(data["title"])
    plt.show()

if __name__ == '__main__':
    with open('part-1.json') as f:
        data = json.load(f)
    plot(data)
    with open("part-2-1-1.json") as f:
        data = json.load(f)
    plot(data)
    plot_curve(data)
    with open("part-2-1-2.json") as f:
        data = json.load(f)
    plot(data)
    plot_curve(data)
    with open("part-2-2-1.json") as f:
        data = json.load(f)
    plot(data)
    with open("part-2-2-2.json") as f:
        data = json.load(f)
    plot(data)