import numpy as np
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from helper_mouse_pars import *
from helper_mouse_model import *
from tqdm import tqdm


# ------------------------------------------set up the exposure time & initial condition--------------------------------------------------------------

# every 24 hours for 1 total doses
N_doses = 1  # number of doses,single dose
# DOSE = 4.2  # Single oral dose from Luccisano et al. (2012) ug/kg/day
Dose_test = 0.85  # once injected 0.85 mg/kg from https://www.sciencedirect.com/science/article/pii/S0041008X10000724?via%3Dihub#fig3
# seting the time range when solving the ode equation
ode_t_step = 0.5  # h
ode_t_end = 24 * 180  # h, in total 6 months

# Dose_test = 3  # 3 mg/kg intravenously https://www.sciencedirect.com/science/article/pii/S0927776518302091?via%3Dihub
# ode_t_step = 1  # h
# ode_t_end = 2160  # h, in total 90 days
# Dose_test = (
#     250  # mg/kg https://link.springer.com/article/10.1007/s11095-010-0068-z#Sec2
# )
# ode_t_step = 1  # h
# ode_t_end = 72  # min, in total 72 hours
# Dose_test = 1  # mg/kg https://www.nature.com/articles/s41598-019-48748-3#Sec12
# ode_t_step = 1  # h
# ode_t_end = 336  # min, in total 72 hours

Dose = Dose_test * Body_weight  # Amount of oral dose, mg


ls_ode_time = np.arange(0, ode_t_end + ode_t_step, ode_t_step)


# sequence of concentrations needs to be calculated
ls_result = [
    # "MIV",
    # "M_Oral",
    "MBA",
    "MBV",
    "M_Lung_b",
    "M_Lung_t",
    "M_Lung_pc",
    "M_Brain_b",
    "M_Brain_t",
    "M_Rest_b",
    "M_Rest_t",
    "M_Kidney_b",
    "M_Kidney_t",
    "M_Kidney_pc",
    "M_Spleen_b",
    "M_Spleen_t",
    "M_Spleen_pc",
    "M_Liver_b",
    "M_Liver_t",
    "M_Liver_pc",
    "M_Urine",
    "M_bile",
    "M_GI_b",
    "M_GI_t",
    "M_GI_pc",
    "M_fecal",
]

## Get a prediction
idx_dose_loc = 19  # if oral dose
idx_dose_loc = 1  # injection

# initial condition for each part
y0 = [0] * idx_dose_loc + [Dose] + [0] * (len(ls_result) - idx_dose_loc - 1)


# ------------------------------------------------------solve the differential equation--------------------------------------------------
# y0 = [0] * (len(ls_result))
result = [y0]


sum_has_run = False
for i in tqdm(range(1, len(ls_ode_time))):
    sol = odeint(
        mouse_pbpk,
        y0,
        ls_ode_time[i - 1 : i + 1],
        rtol=1e-8,
        atol=1e-08,
        tfirst=True,
    )

    y0 = sol[1, :]
    result.append(y0)

result = np.array(result)
# ------------------------------------------------------solve the differential equation--------------------------------------------------


# save the results into a dictionary
d = dict()
for i in range(len(ls_result)):
    d[ls_result[i]] = result[:, i]

for j in range(1000, 2000):
    sum = 0
    for i in d.keys():
        sum = sum + d[i][j]
    print(sum)


# ------------------------------------------------------plotting--------------------------------------------------

from plotting import *

plot_process(ls_ode_time, d)

# ------------------------------------------------------plotting--------------------------------------------------
# from lmfit import minimize, Parameters, Parameter, report_fit
from scipy.integrate import odeint

plt.plot(ls_ode_time, d["M_Liver_b"], label="liver_b")
plt.plot(ls_ode_time, d["M_Liver_pc"], label="liver_pc")
plt.plot(ls_ode_time, d["M_Liver_t"], label="liver_t")
plt.legend()


Lung = (d["M_Lung_b"] + d["M_Lung_t"] + d["M_Lung_pc"]) / W_Lung
plt.plot(ls_ode_time, Lung, label="Lung_tot")
plt.plot(ls_ode_time, d["M_Lung_b"], label="Lung_b")
plt.plot(ls_ode_time, d["M_Lung_pc"], label="Lung_pc")
plt.plot(ls_ode_time, d["M_Lung_t"], label="Lung_t")

plt.legend()

Brain = (d["M_Brain_b"] + d["M_Brain_t"]) / W_Brain
plt.plot(ls_ode_time, d["M_Brain_b"], label="brain_b")
plt.plot(ls_ode_time, d["M_Brain_t"], label="brain_t")
plt.legend()

Rest = (d["M_Rest_b"] + d["M_Rest_t"]) / W_Rest
plt.plot(ls_ode_time, d["M_Rest_b"], label="rest_b")
plt.plot(ls_ode_time, d["M_Rest_t"], label="rest_t")
plt.legend()

Kidney = (d["M_Kidney_b"] + d["M_Kidney_t"] + d["M_Kidney_pc"]) / W_Kidney
plt.plot(ls_ode_time, Kidney, label="Kidney_tot")
plt.plot(ls_ode_time, d["M_Kidney_b"], label="Kidney_b")
plt.plot(ls_ode_time, d["M_Kidney_pc"], label="Kidney_pc")
plt.plot(ls_ode_time, d["M_Kidney_t"], label="Kidney_t")
plt.legend()

Spleen = (d["M_Spleen_b"] + d["M_Spleen_t"] + d["M_Spleen_pc"]) / W_Spleen
plt.plot(ls_ode_time, Spleen, label="Spleen_tot")
plt.plot(ls_ode_time, d["M_Spleen_b"], label="Spleen_b")
plt.plot(ls_ode_time, d["M_Spleen_pc"], label="Spleen_pc")
plt.plot(ls_ode_time, d["M_Spleen_t"], label="Spleen_t")
plt.legend()


plt.plot(ls_ode_time, d["M_bile"], label="bile")
plt.plot(ls_ode_time, d["M_Urine"], label="urine")
plt.plot(ls_ode_time, d["M_fecal"], label="fecal")
plt.legend()


plt.plot(ls_ode_time, d["MBV"], label="MBV")
plt.plot(ls_ode_time, d["MBA"], label="MBA")
plt.legend()

GI = (d["M_GI_b"] + d["M_GI_t"] + d["M_GI_pc"]) / W_GI
plt.plot(ls_ode_time, GI, label="GI_tot")
plt.plot(ls_ode_time, d["M_GI_b"], label="GI_b")
plt.plot(ls_ode_time, d["M_GI_pc"], label="GI_pc")
plt.plot(ls_ode_time, d["M_GI_t"], label="GI_t")
plt.legend()
# AUC_CA = d["AUCCA_free"]
# AUC_CL = d["AUCCL"]
# AUC_CK = d["AUCKb"]
# plt.plot(
#     [0.5, 4, 24, 168, 720, 2160, 4320],
#     [
#         592.8853755,
#         1581.027625,
#         16942.625,
#         19819.875,
#         5691.699605,
#         5138.339921,
#         4308.300395,
#     ],
#     "o",
#     label="measured Liver-gold 4nm",
# )  # ng/g, gold 4nm
# plt.plot(
#     [0.5, 4, 24, 168, 720, 2160, 4320],
#     [
#         471.6981132,
#         717.8014766,
#         1579.163249,
#         7824.036095,
#         7454.88105,
#         6408.941756,
#         3240.360952,
#     ],
#     "o",
#     label="measured Liver-gold 13nm",
# )  # ng/g, gold 4nm


# plt.plot(
#     [4, 24, 48, 144, 240, 480, 720, 2160],
#     [
#         3617.647059,
#         4676.470588,
#         5823.529412,
#         2647.058824,
#         2294.117647,
#         2294.117647,
#         1852.941176,
#         1941.176471,
#     ],
#     "o",
#     label="measured Liver-PEG coated Au NPs 6.2nm",
# )  # ng/g, gold 4nm
# plt.plot(
#     [4, 24, 48, 144, 240, 480, 720, 2160],
#     [
#         7323.529412,
#         11029.41176,
#         14735.29412,
#         3882.352941,
#         6088.235294,
#         2823.529412,
#         2294.117647,
#         2470.588235,
#     ],
#     "o",
#     label="measured Liver-PEG coated Au NPs 24.3nm",
# )  # ng/g, gold 4nm
# plt.plot(
#     [4, 24, 48, 144, 240, 480, 720, 2160],
#     [
#         9441.176471,
#         32470.58824,
#         22500,
#         9705.882353,
#         6705.882353,
#         5470.588235,
#         6176.470588,
#         3088.235294,
#     ],
#     "o",
#     label="measured Liver-PEG coated Au NPs 42.5nm",
# )  # ng/g, gold 4nm
# plt.plot(
#     [4, 24, 48, 144, 240, 480, 720, 2160],
#     [
#         30617.64706,
#         40411.76471,
#         29205.88235,
#         14911.76471,
#         13411.76471,
#         13500,
#         8294.117647,
#         6352.941176,
#     ],
#     "o",
#     label="measured Liver-PEG coated Au NPs 61.2nm",
# )  # ng/g, gold 4nm

# plt.plot(
#     [4, 24, 48, 72],
#     [500000, 756250, 700000, 662500],
#     "o",
#     label="measured Liver-111In-NT-BCMS 57nm",
# )  # ng/g, gold 4nm
# plt.plot(
#     [4, 24, 48, 72],
#     [593750, 575000, 600000, 556250],
#     "o",
#     label="measured Liver-111In-T-BCMS 61nm",
# )  # ng/g, gold 4nm
# --------------------------------------------------------------------------
y0 = [0] * idx_dose_loc + [Dose] + [0] * (len(ls_result) - idx_dose_loc - 1)


def func(
    x,
    # K_max_Liver,
    # K_release_Liver,
    # K_50_Liver,
    # n_Liver,
    P_Liver,
    # K_max_Spleen,
    # K_release_Spleen,
    # K_50_Spleen,
    # n_Spleen,
    P_Spleen,
    # K_max_Kidney,
    # K_release_Kidney,
    # K_50_Kidney,
    # n_Kidney,
    P_Kidney,
):
    ls_ode_time = np.arange(0, max(x) + ode_t_step, ode_t_step)
    # sequence of concentrations needs to be calculated
    ls_result = [
        # "MIV",
        # "M_Oral",
        "MBA",
        "MBV",
        "M_Lung_b",
        "M_Lung_t",
        "M_Lung_pc",
        "M_Brain_b",
        "M_Brain_t",
        "M_Rest_b",
        "M_Rest_t",
        "M_Kidney_b",
        "M_Kidney_t",
        "M_Kidney_pc",
        "M_Spleen_b",
        "M_Spleen_t",
        "M_Spleen_pc",
        "M_Liver_b",
        "M_Liver_t",
        "M_Liver_pc",
        "M_Urine",
        "M_bile",
        # "M_GI_b",
        # "M_GI_t",
        # "M_GI_pc",
        "M_fecal",
    ]

    ## Get a prediction
    idx_dose_loc = 19  # if oral dose
    idx_dose_loc = 1  # injection

    # initial condition for each part
    y0 = [0] * idx_dose_loc + [Dose] + [0] * (len(ls_result) - idx_dose_loc - 1)

    result = [y0]

    for i in tqdm(range(1, len(x))):
        sol = odeint(
            mouse_pbpk,
            y0,
            x[i - 1 : i + 1],
            args=(
                # K_max_Liver,
                # K_release_Liver,
                # K_50_Liver,
                # n_Liver,
                P_Liver,
                # K_max_Spleen,
                # K_release_Spleen,
                # K_50_Spleen,
                # n_Spleen,
                P_Spleen,
                # K_max_Kidney,
                # K_release_Kidney,
                # K_50_Kidney,
                # n_Kidney,
                P_Kidney,
            ),
            rtol=1e-8,
            atol=1e-08,
            tfirst=True,
        )

        y0 = sol[1, :]
        result.append(y0)

    result = np.array(result)
    d = dict()
    for i in range(len(ls_result)):
        d[ls_result[i]] = result[:, i]

    Plasma = d["MBV"] / (W_Blood * 0.8)
    Liver = (d["M_Liver_b"] + d["M_Liver_t"] + d["M_Liver_pc"]) / W_Liver
    # idx_t = [np.where(ls_ode_time == i)[0][0] for i in x]
    # return Liver[idx_t] * 1000
    return Liver * 100


from scipy.optimize import curve_fit

# popt, pcov = curve_fit(
#     func,
#     [0.5, 4, 24, 168],
#     [
#         471.6981132,
#         717.8014766,
#         1579.163249,
#         7824.036095,
#     ],
# )
# plt.plot(
#     [0.5, 4, 24, 168],
#     func([0.5, 4, 24, 168], *popt),
#     "r-",
#     # label="fit: a=%5.3f, b=%5.3f, c=%5.3f" % tuple(popt),
# )
# plt.plot(
#     [0.5, 4, 24, 168],
#     #  , 720, 2160, 4320],
#     [
#         471.6981132,
#         717.8014766,
#         1579.163249,
#         7824.036095,
#         # 7454.88105,
#         # 6408.941756,
#         # 3240.360952,
#     ],
#     "o",
#     label="measured Liver-Gold 13nm",
# )  # ng/g, gold 13nm
# plt.show()
# --------------------------------------------------------------------------
