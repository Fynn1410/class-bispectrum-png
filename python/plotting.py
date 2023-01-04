from email.header import Header
from pickle import FALSE, TRUE
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cmx
import csv
import numpy as np
from random import randint

def shape(a):
    if not (type(a) == list or type(a) == np.ndarray):
        return []
    return [len(a)] + shape(a[0])

def dim(a):
    return len(shape(a))

def read_bib(path, skip_rows = 0, columns = None, idx0 = 0, idx1 = -1):
    head = []
    bib = {}

    with open(path, "r", newline='\n') as csvfile:
        reader = csv.reader(csvfile, delimiter=' ')

        # reading header
        head = next(reader)
            
        # init arrays in bib
        for h in head:
            bib[h] = []

        # skipping lines
        for i in range(skip_rows): next(reader)

        # reading data
        for row in reader:
            for i, h in enumerate(head):
                bib[h].append(float(row[i]))

    csvfile.close()

    for h in head:
        bib[h] = np.array(bib[h][idx0:idx1])

    return bib


def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

def negative_selector(x):
    idx = [] #transition of positive-negative and vice-versa
    
    idx.append(0)
    for i in range(1,len(x)-2):
        if (x[i] > 0 and x[i+1] <= 0):
            idx.append(i+1)
        elif (x[i] < 0 and x[i+1] >= 0):
            idx.append(i+1)
    idx.append(len(x)-1)
    return idx

def pn(x):
    idx = negative_selector(x)
    pos = np.zeros(len(x))
    neg = np.zeros(len(x))
    if(x[0]>=0): 
        for i in range(0,len(idx)-1):
            if(i%2 == 0):
                pos[idx[i]:idx[i+1]] = x[idx[i]:idx[i+1]]
            else:
                neg[idx[i]:idx[i+1]] = x[idx[i]:idx[i+1]]

    else: 
        for i in range(0,len(idx)-1):
            if(i%2 == 1):
                pos[idx[i]:idx[i+1]] = x[idx[i]:idx[i+1]]
            else:
                neg[idx[i]:idx[i+1]] = x[idx[i]:idx[i+1]]

    return pos, neg, idx

def negate_plot(plot, x, y, color, label):
    p_y, n_y, id_y = pn(y)
    if(n_y[0] == 0):
        for i in range(len(id_y)-1):
            if(i==0): plot.plot(x[id_y[i]:id_y[i+1]], p_y[id_y[i]:id_y[i+1]], color = color, linestyle = "-", label = label)

            if(i%2 == 0):
                plot.plot(x[id_y[i]:id_y[i+1]], p_y[id_y[i]:id_y[i+1]], color = color, linestyle = "-")
            else:
                plot.plot(x[id_y[i]:id_y[i+1]], abs(n_y[id_y[i]:id_y[i+1]]), color = color, linestyle = "--")
            
    else:
        for i in range(len(id_y)-1):
            if(i==0): plot.plot(x[id_y[i]:id_y[i+1]], abs(n_y[id_y[i]:id_y[i+1]]), color = color, linestyle = "--", label = label)
            if(i%2 == 1):
                plot.plot(x[id_y[i]:id_y[i+1]], p_y[id_y[i]:id_y[i+1]], color = color, linestyle = "-")
            else:
                plot.plot(x[id_y[i]:id_y[i+1]], abs(n_y[id_y[i]:id_y[i+1]]), color = color, linestyle = "--")


def negate_plot_bib(plot, x_key, data, y_meta = None, y_keys = None):
    if(y_keys == None):
        for key in data.keys():
            if (key != x_key):
                if(y_meta == None):
                    label = key
                    color = list(mcolors.CSS4_COLORS.keys())[randint(0, len(mcolors.CSS4_COLORS.keys()))]
                else:
                    if(key in y_meta):
                        if("label" in y_meta[key]): label = y_meta[key]["label"]
                        else: label = key

                        if("color" in y_meta[key]): color = y_meta[key]["color"]
                        else: color = list(mcolors.CSS4_COLORS.keys())[randint(0, len(mcolors.CSS4_COLORS.keys()))]

                    else:
                        label = key
                        color = list(mcolors.CSS4_COLORS.keys())[randint(0, len(mcolors.CSS4_COLORS.keys()))]

                negate_plot(plot, data[x_key], data[key], color, label)
    else: 
        for y_key in y_keys:
            if (y_key in data):
                if(y_meta == None):
                    label = y_key
                    color = list(mcolors.CSS4_COLORS.keys())[randint(0, len(mcolors.CSS4_COLORS.keys()))]
                else:
                    if(y_key in y_meta):
                        
                        if("label" in y_meta[y_key]): 
                            label = y_meta[y_key]["label"]
                        else: label = y_key

                        if("color" in y_meta[y_key]): color = y_meta[y_key]["color"]
                        else: color = list(mcolors.CSS4_COLORS.keys())[randint(0, len(mcolors.CSS4_COLORS.keys()))]
                        
                    else:
                        label = y_key
                        color = list(mcolors.CSS4_COLORS.keys())[randint(0, len(mcolors.CSS4_COLORS.keys()))]

                negate_plot(plot, data[x_key], data[y_key], color, label)

def plot_bib(x_key, data, title = '', y_meta = None, y_key = None, save = False, path = ""):
    fig = plt.figure(figsize=(10, 8))
    sub = fig.add_subplot(1,1,1)

    negate_plot_bib(sub, x_key, data, y_meta)

    sub.set_xscale("log")
    sub.set_yscale("log")
    sub.set_title(title, fontsize = 14)
    sub.set_ylabel("P [$h^{-3}$$Mpc^3$]", fontsize = 12)
    sub.set_xlabel("k [h/Mpc]", fontsize = 12)
    sub.legend()
    sub.grid()
    if(save): 
        plt.savefig(path+".pdf")
    plt.show()

def plot_data_bib(path, x_key, columns = None, y_meta = None, title = "", idx0 = 0, idx1 = -1, skip_rows = 0):
    bib = {}
    bib = read_bib(path, skip_rows = skip_rows, columns = columns, idx0 = idx0, idx1 = idx1)

    fig = plt.figure(figsize=(10, 8))
    sub = fig.add_subplot(1,1,1)

    negate_plot_bib(sub, x_key, bib, y_meta)

    sub.set_xscale("log")
    sub.set_yscale("log")
    sub.set_title(title, fontsize = 14)
    sub.set_ylabel("P [$h^{-3}$$Mpc^3$]", fontsize = 12)
    sub.set_xlabel("k [h/Mpc]", fontsize = 12)
    sub.legend()
    sub.grid()
    plt.show()


def plot_comparison(x, y1, y2, label1, label2, title = None, plot_meta = {"x-label": "k [h/Mpc]", "y-label": "P [$h^{-3}$$Mpc^3$]", "x_scale": "log", "y_scale": "log", "x_lim": [], "y_lim": [], "y_lim_res": []}, save = False, path = ""):
    fig = plt.figure(figsize=(10, 12))
    sub = fig.add_axes([0.15,0.35,0.75,0.55])
    sub_res = fig.add_axes([0.15,0.15,0.75,0.2])
    if(dim(y1) == 1):
        negate_plot(sub, x, y1, label = label1, color = "b")
        negate_plot(sub, x, y2, label = label2, color = "r")
        sub_res.plot(x, np.subtract(np.divide(y2,y1)*100., 100.), label = "({y2}/{y1})-1".format(y2=label2, y1=label1), color = "b", linestyle = "-")  
    elif(dim(y1) == 2):
        cmp = plt.get_cmap("jet")
        cNorm = mcolors.Normalize(vmin = 0, vmax = 2*shape(y1)[0])
        scalarMap = cmx.ScalarMappable(norm = cNorm, cmap = cmp)
        for i in range(shape(y1)[0]):
            negate_plot(sub, x, y1[i], label = label1[i], color = scalarMap.to_rgba(2*i))
            negate_plot(sub, x, y2[i], label = label2[i], color = scalarMap.to_rgba(2*i+1))
            sub_res.plot(x, np.subtract(np.divide(y2[i],y1[i])*100., 100.), label = "({y2}/{y1})-1".format(y2=label2[i], y1=label1[i]), color = scalarMap.to_rgba(2*i+1), linestyle = "-")  
    if "x_scale" in plot_meta:
        sub.set_xscale(plot_meta["x_scale"])
        sub_res.set_xscale(plot_meta["x_scale"])
    else:
        sub.set_xscale("log")
        sub_res.set_xscale("log")
    if "y_scale" in plot_meta:
        sub.set_yscale(plot_meta["y_scale"])
    else:
        sub.set_yscale("log")
    if "x_lim" in plot_meta:
        if plot_meta["x_lim"]: 
            sub.set_xlim(plot_meta["x_lim"])
            sub_res.set_xlim(plot_meta["x_lim"])
    if "y_lim" in plot_meta:
        if plot_meta["y_lim"]: 
            sub.set_ylim(plot_meta["y_lim"])
    if "y_lim_res" in plot_meta:
        if plot_meta["y_lim_res"]: 
            sub_res.set_ylim(plot_meta["y_lim_res"])      
    sub.set_xticks([])
    if(title != None):
        sub.set_title(title, fontsize = 16)
    if "x-label" in plot_meta:
        sub.set_xlabel(plot_meta["x-label"], fontsize = 14)
        sub_res.set_xlabel(plot_meta["x-label"], fontsize = 14)
    else:
        sub.set_xlabel("k [h/Mpc]", fontsize = 14) 
        sub_res.set_xlabel("k [h/Mpc]", fontsize = 14) 
    if "y-label" in plot_meta:
        sub.set_ylabel(plot_meta["y-label"], fontsize = 14)
    else:
        sub.set_ylabel("P [$h^{-3}$$Mpc^3$]", fontsize = 14) 
    sub.legend()
    sub.grid()
    sub_res.set_ylabel("rel. Deviation [%]", fontsize = 14)
    sub_res.legend(fontsize = 14)
    sub_res.grid()
    if(save): 
        plt.savefig(path, facecolor = "white")
    plt.show()
    plt.close()

def plot_comparison_bib(x_key, y1_data, y2_data, y_meta, y1_label, y2_label, save = False, path = ""):
    keys = intersection(list(y1_data.keys()),list(y2_data.keys()))
    y_keys = intersection(keys,list(y_meta.keys()))
    for y_key in y_keys:
        if(y_key != x_key):
            plot_comparison(y1_data[x_key], y1_data[y_key], y2_data[y_key], y1_label, y2_label, y_meta[y_key]["label"]+" comparison", save = save , path = path+y_key+"_comp")

def plot_comparison_dif(x, y1, y2, ref, label1, label2, label_ref, title, save = False, path = ""):
    fig = plt.figure(figsize=(10, 12))
    sub = fig.add_subplot(1,1,1)
    negate_plot(sub, x, y1, label = label1, color = "b")
    negate_plot(sub, x, y2, label = label2, color = "r")
    negate_plot(sub, x, y2-y1, label = label2+"-"+label1, color = "g")
    if(len(ref) != 0):
        negate_plot(sub, x, ref, label = label_ref, color = "k")

    sub.set_xscale("log")
    sub.set_yscale("log")
    sub.set_xticks([])
    sub.set_title(title, fontsize = 16)
    sub.set_ylabel("P [$h^{-3}$$Mpc^3$]", fontsize = 14)
    sub.set_xlabel("k [h/Mpc]", fontsize = 14)
    sub.legend()
    sub.grid()
    if(save): 
        plt.savefig(path+".pdf")
    plt.show()
    plt.close()

def plot_comparison_bib_dif(x_key, y1_data, y2_data, y_meta, y1_label, y2_label, ref = [], ref_label = None, save = False, path = ""):
    keys = intersection(list(y1_data.keys()),list(y2_data.keys()))
    y_keys = intersection(keys,list(y_meta.keys()))
    for y_key in y_keys:
        if(y_key != x_key):
            plot_comparison_dif(y1_data[x_key], y1_data[y_key], y2_data[y_key], ref, y1_label, y2_label, ref_label, y_meta[y_key]["label"]+" comparison", save = save , path = path+y_key+"_comp")


def plot_comparison_bib_fast(x_key, y1_data, y2_data, y1_label, y2_label, save = False, path = ""):
    keys = intersection(list(y1_data.keys()),list(y2_data.keys()))
    for y_key in keys:
        if(y_key != x_key):
            plot_comparison(y1_data[x_key], y1_data[y_key], y2_data[y_key], y1_label, y2_label, y_key+" comparison", save = save , path = path+"_"+y_key)


def moment_comparison(x_key, y1_data, y2_data, y_meta, y1_label, y2_label, title, z, plot_meta, mu = None, save = False, path = ""):
    keys = intersection(list(y1_data.keys()),list(y2_data.keys()))
    y_keys = intersection(list(keys),list(y_meta.keys()))
    fig = plt.figure(figsize=(10, 12))
    sub = fig.add_axes([0.15,0.3,0.75,0.5])
    negate_plot_bib(sub, x_key, y1_data, y_meta, y_keys=y_keys)
    sub.set_xscale("log")
    sub.set_yscale("log")
    sub.set_xticks([])
    # sub.set_title(title, fontsize = 14)
    if(len(keys)>6): height = 0.14
    else: height = 0.1
    if(mu==None):
        ylabel = "P(k, z={z:.0f}) [$(Mpc/h)^3$]".format(z=z)
    else:
        ylabel = "P(k, $\mu$={mu}, z={z:.0f}) [$(Mpc/h)^3$]".format(mu=mu,z=z)
    sub.set_ylabel(ylabel, fontsize = 12)
    sub.set_xlabel("k [h/Mpc]", fontsize = 12)
    legend = sub.legend(bbox_to_anchor=(0.,1.02,1.,height), loc = "upper center", ncol = 6, mode = "expand", title = title)
    plt.setp(legend.get_title(), fontsize="x-large")
    sub.grid()
    sub.set_xlim(plot_meta["xlim"])
    sub.set_ylim(plot_meta["ylim1"])
    plt.axvline(0.25, 0, 1, linestyle = "--", color="k",linewidth=1.2)

    sub_res= fig.add_axes([0.15,0.15,0.75,0.15])
    for key in y_keys:
            if (key != x_key):
                sub_res.plot(y1_data[x_key], np.subtract(np.divide(y1_data[key],y2_data[key])*100., 100.), color = y_meta[key]["color"], linestyle = "-")
    sub_res.set_xscale("log")
    sub_res.set_ylabel("({label1}/{label2})-1 [%]".format(label1 = y1_label, label2 = y2_label), fontsize = 10)
    sub_res.set_xlabel("k [h/Mpc]", fontsize = 12)
    sub_res.grid()

    sub_res.set_xlim(plot_meta["xlim"])
    sub_res.set_ylim(plot_meta["ylim2"])

    plt.axvline(0.25, 0, 1, linestyle = "--", color="k",linewidth=1.2)

    if(save): 
        plt.savefig(path+".pdf")
    plt.show()
    plt.close()

if __name__ == "__main__":
    out_oneloop = "/home/dennis/Software/class/python/W_NW_OL.txt"
    wnw_oneloop = read_bib(out_oneloop)

    out_pt = "/home/dennis/Software/CLASS-PT/python/CLASS-PT.txt"
    wnw_pt = read_bib(out_pt)

    plot_meta = {"x-label": "1/Mpc", "x_lim": [1e-3, 1e-1], "y_lim": [1e2, 1e5]}
    plot_comparison(wnw_oneloop["k"], [wnw_oneloop["Plin"], wnw_oneloop["Pnw"]], [wnw_pt["Plin"],wnw_pt["Pnw"]], ["Oneloop Plin", "Oneloop Pnw"], ["PT Plin", "PT Pnw"], plot_meta = plot_meta, save = False)

