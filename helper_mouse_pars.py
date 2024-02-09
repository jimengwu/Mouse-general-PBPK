# -------------------1. Physiological parameters for weight of organs, constant---------
Body_weight = 0.02  # kg, Bodyweight (Davies 1993)

F_W_Liver = 0.055  # Unitless, Fractional liver tissue (Brown 1997) table 4 https://journals.sagepub.com/doi/epdf/10.1177/074823379701300401
F_W_Brain = 0.0165  # Unitless, Fractional brain tissue (Brown 1997)
F_W_Lung = 0.0073  # Unitless, Fractional brain tissue (Brown 1997)
F_W_Kidney = 0.017  # Unitless, Fractional kidney tissue (Brown 1997)
F_W_Spleen = 0.0035  # Unitless, Fractional spleen tissue (Brown 1997)
F_W_GI = 0.0422  # Unitless, Fractional GI tract tissue (Brown 1997, Davies 1993 https://link.springer.com/article/10.1023/A:1018943613122)
F_W_Blood = 0.049  # Unitless, Fractional kidney tissue (Brown 1997)
F_W_Plasma = 0.029  # Unitless, Fractional plasma (Davies 1993)
# -------------------1. Physiological parameters for weight of organs, constant---------


# -------------------2. Physiological parameters for blood flow of organs, constant---------
QC = 13.98  # mL/min, Cardiac output table 22
QC = 13.98 * 60 / 1000  # L/h
# QC = 0.8775
F_Q_Brain = (
    0.033  # Unitless, Fraction blood flow to brain of the cardiac output (Brown 1997)
)
F_Q_Liver = (
    0.162  # Unitless, Fraction blood flow to liver of the cardiac output (Brown 1997)
)
F_Q_Kidney = 0.091  # Unitless, Fraction blood flow to kidney (Brown 1997)
F_Q_Spleen = 0.013  # Unitless, Fraction blood flow to spleen
F_Q_GI = 0.215  # Fraction blood flow to GI https://analyticalsciencejournals.onlinelibrary.wiley.com/doi/epdf/10.1002/jat.2550030607
# -------------------2. Physiological parameters for blood flow of organs, constant---------

# -------------------3. Physiological parameters for blood volume fraction in tissues---------------
F_BV_Liver = (
    0.31  # Unitless, Fraction blood volume in liver, percentage of liver table 30
)
F_BV_GI = 0.04  # Unitless, needs to be changed
F_BV_Brain = 0.03  # Unitless, Fraction blood volume in brain, percentage of brain
F_BV_Kidney = 0.24  # Unitless, Fraction blood volume in kidney, percentage of kidney
F_BV_Spleen = 0.17  # Unitless, Fraction blood volume in spleen, percentage of spleen
F_BV_Lung = 0.5  # Unitless, Fraction blood volume in Lung, percentage of Lung
F_BV_Rest = 0.04  # Unitless, Fraction blood volume in rest of body
# -------------------3. Physiological parameters for blood volume fraction in tissues---------------


# ------------------4. Partition coefficient between tissue and blood---------------
P_Liver = 0.08  # https://www.tandfonline.com/doi/full/10.3109/17435390.2013.863406
# P_Liver = 1  # for 100 nm
P_GI = 0.147
P_Brain = 0.147
# P_Instetine = 0.147
P_Kidney = 0.147  # https://www.tandfonline.com/doi/full/10.3109/17435390.2013.863406
P_Spleen = 0.147  # https://www.tandfonline.com/doi/full/10.3109/17435390.2013.863406
P_Lung = 0.147  # https://www.tandfonline.com/doi/full/10.3109/17435390.2013.863406
P_Rest = 0.147  # https://www.tandfonline.com/doi/full/10.3109/17435390.2013.863406
# ------------------4. Partition coefficient between tissue and blood---------------


# -------------------5. Diffusion limitation coefficient constants (Permeability coefficient)------------------------
DLC_Liver = 0.1  # unitless
DLC_Brain = 0.000001  # due to blood brain barrier
DLC_Kidney = 0.001
DLC_Spleen = 0.03
DLC_Lung = 0.001
DLC_GI = 0.001
DLC_Rest = 0.000001
# -------------------5. Diffusion limitation coefficient constants (Permeability coefficient)------------------------


# -------------------6. Endocytosis-related parameters------------------------
K_release_Liver = 0.01  # h-1
K_max_Liver = 20  # h-1
K_50_Liver = 48  # h
n_Liver = 1  # Unitless


K_release_GI = 30
K_max_GI = 0.075
K_50_GI = 24
n_GI = 5


K_release_Spleen = 0.001
K_max_Spleen = 40
K_50_Spleen = 48
n_Spleen = 5


K_release_Kidney = 0.0004
K_max_Kidney = 0.075
K_50_Kidney = 24
n_Kidney = 5

K_release_Lung = 0.003
K_max_Lung = 0.075
K_50_Lung = 24
n_Lung = 5


# -------------------6. Endocytosis-related parameters------------------------


# -------------------7. Excretion----------------------------
Kbile = 0.000001  # Biliary clearance (L/hr)
Kurine = 0.000003  # Urine clearance (L/hr)
Kfecal = 0.000003  # fitting h-1


# -------------------7. Excretion----------------------------


# -------------------Blood flow------------

Q_Brain = F_Q_Brain * QC  # L/h, Blood flow to brain
Q_Kidney = F_Q_Kidney * QC  # L/h, Blood flow to kidney
Q_Spleen = F_Q_Spleen * QC  # L/h, Blood flow to Spleen
Q_GI = F_Q_GI * QC  # L/h, Blood flow to GI
Q_Liver = (
    F_Q_Liver * QC
)  # L/h, Blood flow to liver, including the blood to GI and spleen

Q_Rest = QC - Q_Liver - Q_Brain - Q_Kidney  # L/h, Blood flow to the rest of body

# ------------------Blood flow--------------


# --------------Tissue Volume-------------
W_Liver = Body_weight * F_W_Liver  # Liver, kg
W_Brain = Body_weight * F_W_Brain  # Brain
W_Kidney = Body_weight * F_W_Kidney  # Kidney
W_Spleen = Body_weight * F_W_Spleen  # Spleen
W_Lung = Body_weight * F_W_Lung  # Lungs
W_GI = Body_weight * F_W_GI  # GI
W_Blood = Body_weight * F_W_Blood  # blood
W_Plasma = Body_weight * F_W_Plasma  # plasma


W_Rest = (
    Body_weight
    - W_Liver
    - W_Brain
    - W_Kidney
    - W_Spleen
    - W_Lung
    - W_Plasma
    # - W_GI
)

W_balance = (
    Body_weight
    - W_Liver
    - W_Brain
    - W_Kidney
    - W_Spleen
    - W_Lung
    - W_Plasma
    - W_Rest
    - W_GI
)

W_Liver_b = (
    W_Liver * F_BV_Liver
)  # Weight/volume of capillary blood in liver compartment
W_Liver_t = W_Liver - W_Liver_b  # Weight/volume of tissue in liver compartment

W_GI_b = W_GI * F_BV_GI  # Weight/volume of capillary blood in GI compartment
W_GI_t = W_GI - W_GI_b  # Weight/volume of tissue in GI compartment


W_Brain_b = (
    W_Brain * F_BV_Brain
)  # Weight/volume of capillary blood in Brain compartment
W_Brain_t = W_Brain - W_Brain_b  # Weight/volume of tissue in Brain compartment


W_Kidney_b = (
    W_Kidney * F_BV_Kidney
)  # Weight/volume of capillary blood in Kidney compartment
W_Kidney_t = W_Kidney - W_Kidney_b  # Weight/volume of tissue in Kidney compartment

W_Spleen_b = (
    W_Spleen * F_BV_Spleen
)  # Weight/volume of capillary blood in Spleen compartment
W_Spleen_t = W_Spleen - W_Spleen_b  # Weight/volume of tissue in Spleen compartment

W_Lung_b = W_Lung * F_BV_Lung  # Weight/volume of capillary blood in Lung compartment
W_Lung_t = W_Lung - W_Lung_b  # Weight/volume of tissue in Lung compartment

W_Rest_b = W_Rest * F_BV_Rest  # Weight/volume of capillary blood in Rest compartment
W_Rest_t = W_Rest - W_Rest_b  # Weight/volume of tissue in Rest compartment

# --------------Tissue Volume-------------


# fITTING WITH ONLY LIVER DATA
import pandas as pd

params = pd.read_excel(
    "C:/switchdriver/dataset/tk/mouse/mouse_tk_pars.xlsx", sheet_name=1
)
# K_release_Liver = 0.01
# K_50_Liver = 1387
# K_release_GI = 0.08
# K_max_GI = 0.236
# K_50_GI = 179.38
# K_release_Spleen = 0.0002
# K_max_Spleen = 0.449
# K_50_Spleen = 13900
# n_Spleen = 0.0007
# K_release_Kidney = 6.8e-5
# K_50_Kidney = 10.64
# K_release_Lung = 0.0018
# K_50_Lung = 0.2404
# Kbile = 1.8938e-5  # Biliary clearance (L/hr)
# Kurine = 7.9e-6  # Urine clearance (L/hr)
# Kfecal = 0.00186  # fitting h-1

col_name = "100nm_liver"
K_release_Liver = params[col_name][0]
K_max_Liver = params[col_name][1]
K_50_Liver = params[col_name][2]
K_release_GI = params[col_name][3]
K_50_GI = params[col_name][4]
K_release_Spleen = params[col_name][5]
K_max_Spleen = params[col_name][6]
K_50_Spleen = params[col_name][7]
K_release_Kidney = params[col_name][8]
K_50_Kidney = params[col_name][9]
K_release_Lung = params[col_name][10]
K_50_Lung = params[col_name][11]
DLC_Brain = params[col_name][12]
DLC_Kidney = params[col_name][13]
DLC_Spleen = params[col_name][14]
DLC_Lung = params[col_name][15]
DLC_GI = params[col_name][16]
DLC_Rest = params[col_name][17]
Kbile = params[col_name][18]  # Biliary clearance (L/hr)
Kurine = params[col_name][19]  # Urine clearance (L/hr)
Kfecal = params[col_name][20]  # fitting h-1
