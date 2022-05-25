import h5py, pandas as pd, numpy as np
from matplotlib import pyplot as plt


def load_stats(filename):
    f = h5py.File(filename, "r")

    linear = {}
    model = {}

    if "linear/POF" in f:
        linear["POF"] = f["linear/POF"]
    if "linear/POU" in f:
        linear["POU"] = f["linear/POU"]
    if "linear/ideal_distance" in f:
        linear["ideal/ideal_distance"] = f["linear/ideal_distance"]

    closef(f)


def aggregate_data(filename):
    schemes = ["IF", "Rawls", "Aristotle", "Nash"]
    stats = ["time", "support_size", "ideal_distance", "nadir_distance", "POU", "POF"]
    tuples = list(zip([scheme for scheme in schemes for i in range(len(stats))], stats * 4))
    tuples.append(("Graph", "|P|"))
    tuples.append(("Graph", "|N"))
    index = pd.MultiIndex.from_tuples(tuples, names = ["scheme", "stat"])
    nswp_data = pd.DataFrame(columns = index) 
    linear_data = pd.DataFrame(columns = index)

    f = h5py.File(filename, "r")

    for scheme in schemes:
        for stat in stats:

            if f"model/{stat}/{scheme}" in f:
                nswp_data[scheme, stat] = f[f"model/{stat}/{scheme}"]

            if f"linear/{stat}/{scheme}" in f:
                linear_data[scheme, stat] = f[f"linear/{stat}/{scheme}"]

    # add the sizes of the sets P and N
    nswp_data["Graph", "|P|"] = f["model/|P|"]
    nswp_data["Graph", "|N|"] = f["model/|N|"]
    linear_data["Graph", "|P|"] = f["linear/|P|"]
    linear_data["Graph", "|N|"] = f["linear/|N|"]

    f.close()
    return nswp_data, linear_data


def plot_POF(df):
    if "IF" in df:
        mask = df["IF", "POF"] >= 0
        x = df["IF", "POF"][mask]
        print("POF (IF): ", np.mean(x), " \u00B1 ", np.std(x)) 
    if "Rawls" in df:
        mask = df["Rawls", "POF"] >= 0
        x = df["Rawls", "POF"][mask]
        print("POF (Rawls): ", np.mean(x), " \u00B1 ", np.std(x)) 
    if "Aristotle" in df:
        mask = df["Aristotle", "POF"] >= 0
        x = df["Aristotle", "POF"][mask]
        print("POF (Aristotle): ", np.mean(x), " \u00B1 ", np.std(x)) 
    if "Nash" in df:
        mask = df["Nash", "POF"] >= 0
        x = df["Nash", "POF"][mask]
        print("POF (Nash): ", np.mean(x), " \u00B1 ", np.std(x))
   

def plot_POU(df):
    if "IF" in df:
        mask = df["IF", "POU"] >= 0
        x = df["IF", "POU"][mask]
        print("POU (IF): ", np.mean(x), " \u00B1 ", np.std(x)) 
    if "Rawls" in df:
        mask = df["Rawls", "POU"] >= 0
        x = df["Rawls", "POU"][mask]
        print("POU (Rawls): ", np.mean(x), " \u00B1 ", np.std(x)) 
    if "Aristotle" in df:
        mask = df["Aristotle", "POU"] >= 0
        x = df["Aristotle", "POU"][mask]
        print("POU (Aristotle): ", np.mean(x), " \u00B1 ", np.std(x)) 
    if "Nash" in df:
        mask = df["Nash", "POU"] >= 0
        x = df["Nash", "POU"][mask]
        print("POU (Nash): ", np.mean(x), " \u00B1 ", np.std(x))


def plot_time(filename, nswp_df, linear_df):
    fig, axs = plt.subplots(2, 4, sharex = False, sharey = False)
    fig.text(0.5, 0.00, "Size of |P|", ha = "center")
    fig.text(0.02, 0.5, "Seconds", va = "center", rotation = "vertical")
    nswp_sizes = sorted(set(nswp_df["Graph", "|P|"]))
    linear_sizes = sorted(set(linear_df["Graph", "|P|"]))

    # plot the elapsed times for the NSWP
    if "IF" in nswp_df:
        data = []
        for size in nswp_sizes:
            mask = np.logical_and(nswp_df["IF", "time"] >= 0, nswp_df["Graph", "|P|"] == size) 
            x = nswp_df["IF", "time"][mask]
            data.append(x)
        axs[0, 0].xaxis.set_tick_params(rotation = 30, labelsize = 8)
        axs[0, 0].boxplot(data, labels = nswp_sizes)

    if "Rawls" in nswp_df:
        data = []
        for size in nswp_sizes:
            mask = np.logical_and(nswp_df["Rawls", "time"] >= 0, nswp_df["Graph", "|P|"] == size) 
            x = nswp_df["Rawls", "time"][mask]
            data.append(x)
        axs[0, 1].xaxis.set_tick_params(rotation = 30, labelsize = 8)
        axs[0, 1].boxplot(data, labels = nswp_sizes)

    if "Aristotle" in nswp_df:
        data = []
        for size in nswp_sizes:
            mask = np.logical_and(nswp_df["Aristotle", "time"] >= 0, nswp_df["Graph", "|P|"] == size) 
            x = nswp_df["Aristotle", "time"][mask]
            data.append(x)
        axs[0, 2].xaxis.set_tick_params(rotation = 30, labelsize = 8)
        axs[0, 2].boxplot(data, labels = nswp_sizes)

    if "Nash" in nswp_df:
        data = []
        for size in nswp_sizes:
            mask = np.logical_and(nswp_df["Nash", "time"] >= 0, nswp_df["Graph", "|P|"] == size) 
            x = nswp_df["Nash", "time"][mask]
            data.append(x)
        axs[0, 3].xaxis.set_tick_params(rotation = 30, labelsize = 8)
        axs[0, 3].boxplot(data, labels = nswp_sizes)

    
    # plot the elapsed times for the linear combination of objectives 
    if "IF" in linear_df:
        data = []
        for size in nswp_sizes:
            mask = np.logical_and(linear_df["IF", "time"] >= 0, linear_df["Graph", "|P|"] == size) 
            x = linear_df["IF", "time"][mask]
            data.append(x)
        axs[1, 0].xaxis.set_tick_params(rotation = 30, labelsize = 8)
        axs[1, 0].boxplot(data, labels = nswp_sizes)

    if "Rawls" in linear_df:
        data = []
        for size in nswp_sizes:
            mask = np.logical_and(linear_df["Rawls", "time"] >= 0, linear_df["Graph", "|P|"] == size) 
            x = linear_df["Rawls", "time"][mask]
            data.append(x)
        axs[1, 1].xaxis.set_tick_params(rotation = 30, labelsize = 8)
        axs[1, 1].boxplot(data, labels = nswp_sizes)

    if "Aristotle" in linear_df:
        data = []
        for size in nswp_sizes:
            mask = np.logical_and(linear_df["Aristotle", "time"] >= 0, linear_df["Graph", "|P|"] == size) 
            x = linear_df["Aristotle", "time"][mask]
            data.append(x)
        axs[1, 2].xaxis.set_tick_params(rotation = 30, labelsize = 8)
        axs[1, 2].boxplot(data, labels = nswp_sizes)

    if "Nash" in linear_df:
        data = []
        for size in nswp_sizes:
            mask = np.logical_and(linear_df["Nash", "time"] >= 0, linear_df["Graph", "|P|"] == size) 
            x = linear_df["Nash", "time"][mask]
            data.append(x)
        axs[1, 3].xaxis.set_tick_params(rotation = 30, labelsize = 8)
        axs[1, 3].boxplot(data, labels = nswp_sizes)

    fig.tight_layout()
    fig.savefig(filename, dpi = 250, bbox_inches = "tight")


def plot_support(df):
    if "IF" in df:
        mask = df["IF", "support_size"] > 0
        x = df["IF", "support_size"][mask]
        print("Support size (IF): ", np.mean(x), " \u00B1 ", np.std(x)) 
    if "Rawls" in df:
        mask = df["Rawls", "support_size"] > 0
        x = df["Rawls", "support_size"][mask]
        print("Support size (Rawls): ", np.mean(x), " \u00B1 ", np.std(x)) 
    if "Aristotle" in df:
        mask = df["Aristotle", "support_size"] > 0
        x = df["Aristotle", "support_size"][mask]
        print("Support size (Aristotle): ", np.mean(x), " \u00B1 ", np.std(x)) 
    if "Nash" in df:
        mask = df["Nash", "support_size"] > 0
        x = df["Nash", "support_size"][mask]
        print("Support size (Nash): ", np.mean(x), " \u00B1 ", np.std(x))


def plot_solved(filename, nswp_df, linear_df):
    fig, ax = plt.subplots()
    ax.set_xlabel("Time elapsed")
    ax.set_ylabel("Percentage of instances solved")

    # plot the percentage of instances solved
    if "IF" in nswp_df:
        mask = nswp_df["IF", "time"] >= 0
        time = nswp_df["IF", "time"][mask]
        x = np.sort(time)
        y = np.array([np.mean(time <= x[i]) for i in range(len(x))])
        ax.plot(x, y, color = "black", linestyle = "-", label = "NSWP+IF")
    
    if "Rawls" in nswp_df:
        mask = nswp_df["Rawls", "time"] >= 0
        time = nswp_df["Rawls", "time"][mask]
        x = np.sort(time)
        y = np.array([np.mean(time <= x[i]) for i in range(len(x))])
        ax.plot(x, y, color = "black", linestyle = "--", label = "NSWP+Rawls")

    if "Aristotle" in nswp_df:
        mask = nswp_df["Aristotle", "time"] >= 0
        time = nswp_df["Aristotle", "time"][mask]
        x = np.sort(time)
        y = np.array([np.mean(time <= x[i]) for i in range(len(x))])
        ax.plot(x, y, color = "black", linestyle = "-.", label = "NSWP+Aristotle")

    if "Nash" in nswp_df:
        mask = nswp_df["Nash", "time"] >= 0
        time = nswp_df["Nash", "time"][mask]
        x = np.sort(time)
        y = np.array([np.mean(time <= x[i]) for i in range(len(x))])
        ax.plot(x, y, color = "black", linestyle = ":", label = "NSWP+Nash")

    # plot the the percentage of instances solved for the linear combination of objectives
    if "IF" in linear_df:
        mask = linear_df["IF", "time"] >= 0
        time = linear_df["IF", "time"][mask]
        x = np.sort(time)
        y = np.array([np.mean(time <= x[i]) for i in range(len(x))])
        ax.plot(x, y, color = "orange", linestyle = "-", label = "IF")
    
    if "Rawls" in linear_df:
        mask = linear_df["Rawls", "time"] >= 0
        time = linear_df["Rawls", "time"][mask]
        x = np.sort(time)
        y = np.array([np.mean(time <= x[i]) for i in range(len(x))])
        ax.plot(x, y, color = "orange", linestyle = "--", label = "Rawls")

    if "Aristotle" in linear_df:
        mask = linear_df["Aristotle", "time"] >= 0
        time = linear_df["Aristotle", "time"][mask]
        x = np.sort(time)
        y = np.array([np.mean(time <= x[i]) for i in range(len(x))])
        ax.plot(x, y, color = "orange", linestyle = "-.", label = "Aristotle")

    if "Nash" in linear_df:
        mask = linear_df["Nash", "time"] >= 0
        time = linear_df["Nash", "time"][mask]
        x = np.sort(time)
        y = np.array([np.mean(time <= x[i]) for i in range(len(x))])
        ax.plot(x, y, color = "orange", linestyle = ":", label = "Nash")

    ax.legend()
    ax.set_xscale("log")
    fig.savefig(filename, dpi = 250)


def plot_ideal(df):
    if "IF" in df:
        mask = df["IF", "ideal_distance"] >= 0
        x = df["IF", "ideal_distance"][mask]
        print("Distance to ideal (IF): ", np.mean(x), " \u00B1 ", np.std(x)) 
    if "Rawls" in df:
        mask = df["Rawls", "ideal_distance"] >= 0
        x = df["Rawls", "ideal_distance"][mask]
        print("Distance to ideal (Rawls): ", np.mean(x), " \u00B1 ", np.std(x)) 
    if "Aristotle" in df:
        mask = df["Aristotle", "ideal_distance"] >= 0
        x = df["Aristotle", "ideal_distance"][mask]
        print("Distance to ideal (Aristotle): ", np.mean(x), " \u00B1 ", np.std(x)) 
    if "Nash" in df:
        mask = df["Nash", "ideal_distance"] >= 0
        x = df["Nash", "ideal_distance"][mask]
        print("Distance to ideal (Nash): ", np.mean(x), " \u00B1 ", np.std(x))


def plot_nadir(df):
    if "IF" in df:
        mask = df["IF", "nadir_distance"] >= 0
        x = df["IF", "nadir_distance"][mask]
        print("Distance to nadir (IF): ", np.mean(x), " \u00B1 ", np.std(x)) 
    if "Rawls" in df:
        mask = df["Rawls", "nadir_distance"] >= 0
        x = df["Rawls", "nadir_distance"][mask]
        print("Distance to nadir (Rawls): ", np.mean(x), " \u00B1 ", np.std(x)) 
    if "Aristotle" in df:
        mask = df["Aristotle", "nadir_distance"] >= 0
        x = df["Aristotle", "nadir_distance"][mask]
        print("Distance to nadir (Aristotle): ", np.mean(x), " \u00B1 ", np.std(x)) 
    if "Nash" in df:
        mask = df["Nash", "nadir_distance"] >= 0
        x = df["Nash", "nadir_distance"][mask]
        print("Distance to nadir (Nash): ", np.mean(x), " \u00B1 ", np.std(x))


def main():
    nswp_data, linear_data = aggregate_data("stats.hdf")
    
    # Plot the POF
    print("### NSWP ###")
    plot_POF(nswp_data)
    print("\n\n")
    print("### Linear ###")
    plot_POF(linear_data)
    print("\n\n")

    # Plot the POU
    print("### NSWP ###")
    plot_POU(nswp_data)
    print("\n\n")
    print("### Linear ###")
    plot_POU(linear_data)
    print("\n\n")

    # Plot the distance to ideal
    print("### NSWP ###")
    plot_ideal(nswp_data)
    print("\n\n")
    print("### Linear ###")
    plot_ideal(linear_data)
    print("\n\n")

    # Plot the distance to nadir 
    print("### NSWP ###")
    plot_nadir(nswp_data)
    print("\n\n")
    print("### Linear ###")
    plot_nadir(linear_data)
    print("\n\n")

    # Plot the size of the support
    print("### NSWP ###")
    plot_support(nswp_data)
    print("\n\n")
    print("### Linear ###")
    plot_support(linear_data)
    print("\n\n")

    # Plot the times
    plot_time("time.png", nswp_data, linear_data)

    # Plot the instances solved
    plot_solved("solved.png", nswp_data, linear_data)

main()
