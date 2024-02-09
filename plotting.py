import pandas as pd
from helper_mouse_pars import *
import matplotlib.pyplot as plt


def plot_process(ls_ode_time, d):
    Plasma = d["MBV"] / (W_Blood * 0.8)
    Liver = (d["M_Liver_b"] + d["M_Liver_t"] + d["M_Liver_pc"]) / W_Liver
    Kidney = (d["M_Kidney_b"] + d["M_Kidney_t"] + d["M_Kidney_pc"]) / W_Kidney
    Lung = (d["M_Lung_b"] + d["M_Lung_t"] + d["M_Lung_pc"]) / W_Lung
    Spleen = (d["M_Spleen_b"] + d["M_Spleen_t"] + d["M_Spleen_pc"]) / W_Spleen
    GI = (d["M_GI_b"] + d["M_GI_t"] + d["M_GI_pc"]) / W_GI
    Brain = (d["M_Brain_b"] + d["M_Brain_t"]) / W_Brain
    Rest = (d["M_Rest_b"] + d["M_Rest_t"]) / W_Brain

    obs = pd.read_csv("C:/switchdriver/dataset/tk/mouse/R_input_mouse_study1_100nm.csv")
    label = "100nm"

    plt.figure()
    plt.plot(
        obs.Time,
        obs.CL,
        "o",
        alpha=0.7,
        color="blue",
        label="measured Liver-Au ",
    )  # ng/g, gold
    plt.plot(
        ls_ode_time, Liver * 1000, label="calculated Liver", color="black"
    )  # Liver mg/kg transfer to ng/g
    plt.xlabel("Time (h) ")
    plt.ylabel("concentration of PEG-Au (ng / g per gram body weight)")
    plt.legend()
    plt.title(label + "_" + "Liver model and experiment data comparison")
    plt.savefig("plots/" + label + "_Liver.png", dpi=600)

    plt.figure()
    plt.plot(
        obs.Time,
        obs.CK,
        "o",
        alpha=0.7,
        color="blue",
        label="measured Kidney",
    )  # ng/g
    plt.plot(ls_ode_time, Kidney * 1000, label="Kidney", color="black")
    plt.xlabel("Time (h) ")
    plt.ylabel("concentration of PEG-Au (ng / g per gram body weight)")
    plt.legend()
    plt.title(label + "_" + "Kidney model and experiment data comparison")
    plt.savefig("plots/" + label + "_Kidney.png", dpi=600)
    # plt.show()

    plt.figure()
    plt.plot(
        obs.Time,
        obs.Clung,
        "o",
        alpha=0.7,
        color="blue",
        label="measured Lung",
    )  # ng/g
    plt.plot(ls_ode_time, Lung * 1000, label="Lung", color="black")
    plt.xlabel("Time (h) ")
    plt.ylabel("concentration of PEG-Au (ng / g per gram body weight)")
    plt.legend()
    plt.title(label + "_" + "Lung model and experiment data comparison")
    plt.savefig("plots/" + label + "_Lung.png", dpi=600)
    # plt.show()

    plt.figure()
    plt.plot(
        obs.Time,
        obs.CS,
        "o",
        alpha=0.7,
        color="blue",
        label="measured Spleen",
    )  # ng/g
    plt.plot(ls_ode_time, Spleen * 1000, label="Spleen", color="black")
    plt.xlabel("Time (h) ")
    plt.ylabel("concentration of PEG-Au (ng / g per gram body weight)")
    plt.legend()
    plt.title(label + "_" + "Spleen model and experiment data comparison")
    plt.savefig("plots/" + label + "_Spleen.png", dpi=600)
    # plt.show()

    # plt.figure()
    # plt.plot(
    #     obs.Time,
    #     obs.CB,
    #     "o",
    #     alpha=0.7,
    #     color="blue",
    #     label="measured Brain",
    # )  # ng/g
    # plt.plot(ls_ode_time, Brain * 1000, label="Brain", color="black")
    # plt.xlabel("Time (h) ")
    # plt.ylabel("concentration of PEG-Au (ng / g per gram body weight)")
    # plt.legend()
    # plt.title(label + "_" + "Brain model and experiment data comparison")
    # plt.savefig("plots/" + label + "_Brain.png", dpi=600)
    # plt.show()

    plt.figure()
    plt.plot(ls_ode_time, GI * 1000, label="GI", color="black")
    plt.xlabel("Time (h) ")
    plt.ylabel("concentration of PEG-Au (ng / g per gram body weight)")
    plt.legend()
    plt.title(label + "_" + "GI model results")
    plt.savefig("plots/" + label + "_GI.png", dpi=600)

    plt.figure()
    plt.plot(ls_ode_time, Rest * 1000, label="Rest", color="black")
    plt.xlabel("Time (h) ")
    plt.ylabel("concentration of PEG-Au (ng / g per gram body weight)")
    plt.legend()
    plt.title(label + "_" + "Rest model results")
    plt.savefig("plots/" + label + "_Rest.png", dpi=600)

    plt.figure()
    plt.plot(ls_ode_time, d["M_bile"], label="bile")
    plt.plot(ls_ode_time, d["M_Urine"], label="urine")
    plt.plot(ls_ode_time, d["M_fecal"], label="fecal")
    plt.xlabel("Time (h) ")
    plt.ylabel("amount of PEG-Au (mg)")
    plt.legend()
    plt.title(label + "_" + "Excrection model results")
    plt.savefig("plots/" + label + "_Excrection.png", dpi=600)
