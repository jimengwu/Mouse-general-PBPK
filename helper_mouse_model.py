from helper_mouse_pars import *

sum_has_run = False


# TODO: HOW to merge the inject into; stomach not included?; lymph node included
def mouse_pbpk(
    t,
    y,
):
    (
        MBA,
        MBV,
        M_Lung_b,
        M_Lung_t,
        M_Lung_pc,
        M_Brain_b,
        M_Brain_t,
        M_Rest_b,
        M_Rest_t,
        M_Kidney_b,
        M_Kidney_t,
        M_Kidney_pc,
        M_Spleen_b,
        M_Spleen_t,
        M_Spleen_pc,
        M_Liver_b,
        M_Liver_t,
        M_Liver_pc,
        M_Urine,
        M_bile,
        M_GI_b,
        M_GI_t,
        M_GI_pc,
        M_fecal,
    ) = y

    K_up_Liver = (K_max_Liver * (t ** (n_Liver))) / (
        K_50_Liver ** (n_Liver) + t ** (n_Liver)
    )
    K_up_Spleen = (K_max_Spleen * (t ** (n_Spleen))) / (
        K_50_Spleen ** (n_Spleen) + t ** (n_Spleen)
    )
    K_up_Kidney = (K_max_Kidney * (t ** (n_Kidney))) / (
        K_50_Kidney ** (n_Kidney) + t ** (n_Kidney)
    )
    K_up_Lung = (K_max_Lung * (t ** (n_Lung))) / (K_50_Lung ** (n_Lung) + t ** (n_Lung))

    K_up_GI = (K_max_GI * (t ** (n_GI))) / (K_50_GI ** (n_GI) + t ** (n_GI))

    # no brain compartment included in the endocytosis part

    # # ---------------------permeability area cross product between the capillary blood and the tissue of the organ
    # PA_Liver = DLC_Liver * Q_Liver  # unitless
    # PA_Brain = DLC_Brain * Q_Brain  # unitless
    # PA_Kidney = DLC_Kidney * Q_Kidney  # unitless
    # PA_Spleen = DLC_Spleen * Q_Spleen  # unitless
    # PA_Lung = DLC_Lung * QC  # unitless
    # PA_GI = DLC_GI * Q_GI  # unitless
    # PA_Rest = DLC_Rest * Q_Rest  # unitless

    # -----------------------1. Blood: Concentrations in the blood compartment---------------

    # --------------------------- 1.1 Arterial Blood part ---------------------------------------
    CA = MBA / (W_Plasma * 0.2)

    RBA = QC * (M_Lung_b / W_Lung_b) - QC * (
        MBA / (W_Plasma * 0.2)
    )  # rate of change in Arterial blood
    dxdt_MBA = RBA

    # --------------------------- 1.1 Arterial Blood part ---------------------------------------

    # --------------------------- 1.2 Venous Blood part ---------------------------------------
    CV = MBV / (W_Plasma * 0.8)

    # dxdt_MIV = MIV / t

    def dose(t):
        global sum_has_run
        if sum_has_run:
            return 0
        sum_has_run = True
        if t < 0.005:
            return 0.017
        else:
            return 0

    RBV = (
        Q_Liver * (M_Liver_b / W_Liver_b)
        + Q_Brain * (M_Brain_b / W_Brain_b)
        + Q_Kidney * (M_Kidney_b / W_Kidney_b)
        + Q_Rest * (M_Rest_b / W_Rest_b)
        # + dose(t) / 0.005
        - QC * (MBV / (W_Plasma * 0.8))
    )
    dxdt_MBV = RBV

    # --------------------------- 1.2 Venous Blood part ---------------------------------------

    # --------------------------------- 2. Lung compartment------------------------------------

    # ----------------------------------2.1 lung blood compartment ----------------------------
    R_Lung_b = (
        QC * ((MBV / (W_Plasma * 0.8)) - (M_Lung_b / W_Lung_b))
        - (DLC_Lung * QC) * (M_Lung_b / W_Lung_b)
        + ((DLC_Lung * QC) * (M_Lung_t / W_Lung_t)) / P_Lung
    )
    dxdt_M_Lung_b = R_Lung_b
    # -----------------------------------2.2 in the lung tissue compartment---------------------
    R_Lung_t = (
        (DLC_Lung * QC) * (M_Lung_b / W_Lung_b)
        - ((DLC_Lung * QC) * (M_Lung_t / W_Lung_t)) / P_Lung
        - K_up_Lung * M_Lung_t
        + K_release_Lung * M_Lung_pc
    )
    dxdt_M_Lung_t = R_Lung_t

    # -------------------------------------2.3 in the lung pc compartment----------------------------
    R_lung_pc = K_up_Lung * M_Lung_t - K_release_Lung * M_Lung_pc
    dxdt_M_Lung_pc = R_lung_pc

    # -------------------------------------3. brain compartment--------------------------------------

    # ------------------------------------3.1 brain blood compartment--------------------------------------
    # PA_Brain is a factor of 10 decreased compared to other organs
    R_Brain_b = (
        Q_Brain * (MBA / (W_Plasma * 0.2) - (M_Brain_b / W_Brain_b))
        - (DLC_Brain * Q_Brain) * (M_Brain_b / W_Brain_b)
        + ((DLC_Brain * Q_Brain) * (M_Brain_t / W_Brain_t)) / P_Brain
    )
    dxdt_M_Brain_b = R_Brain_b

    # ------------------------------------3.2 brain tissue compartment--------------------------------------

    R_Brain_t = (DLC_Brain * Q_Brain) * (M_Brain_b / W_Brain_b) - (
        (DLC_Brain * Q_Brain) * (M_Brain_t / W_Brain_t)
    ) / P_Brain
    dxdt_M_Brain_t = R_Brain_t

    C_Brain = (M_Brain_b + M_Brain_t) / W_Brain

    # --------------------------------------------4. Rest of the body -----------------
    # --------------------------------------------4.1 Rest of body Blood------------
    R_Rest_b = (
        Q_Rest * ((MBA / (W_Plasma * 0.2)) - (M_Rest_b / W_Rest_b))
        - (DLC_Rest * Q_Rest) * (M_Rest_b / W_Rest_b)
        + ((DLC_Rest * Q_Rest) * (M_Rest_t / W_Rest_t)) / P_Rest
    )
    dxdt_M_Rest_b = R_Rest_b

    # --------------------------------------------4.2 Rest of body Tissue------------
    R_Rest_t = (DLC_Rest * Q_Rest) * (M_Rest_b / W_Rest_b) - (
        (DLC_Rest * Q_Rest) * (M_Rest_t / W_Rest_t)
    ) / P_Rest
    dxdt_M_Rest_t = R_Rest_t

    C_Rest = (M_Rest_b + M_Rest_t) / W_Rest

    # ------------------------------------5. Kidney ---------------
    # ------------------------------------5.1 Kidney Blood------------

    R_urine = Kurine * (M_Kidney_b / W_Kidney_b)
    R_Kidney_b = (
        Q_Kidney * ((MBA / (W_Plasma * 0.2)) - (M_Kidney_b / W_Kidney_b))
        - (DLC_Kidney * Q_Kidney) * (M_Kidney_b / W_Kidney_b)
        + ((DLC_Kidney * Q_Kidney) * (M_Kidney_t / W_Kidney_t)) / P_Kidney
        - R_urine
    )
    dxdt_M_Kidney_b = R_Kidney_b
    dxdt_M_urine = R_urine

    # ------------------------------------5.2 Kidney Tissue------------
    R_Kidney_t = (
        (DLC_Kidney * Q_Kidney) * (M_Kidney_b / W_Kidney_b)
        - ((DLC_Kidney * Q_Kidney) * (M_Kidney_t / W_Kidney_t)) / P_Kidney
        - K_up_Kidney * M_Kidney_t
        + K_release_Kidney * M_Kidney_pc
    )
    dxdt_M_Kidney_t = R_Kidney_t

    # ------------------------------------5.2 Kidney pc------------
    R_Kidney_pc = K_up_Kidney * M_Kidney_t - K_release_Kidney * M_Kidney_pc
    dxdt_M_Kidney_pc = R_Kidney_pc
    C_Kidney = (M_Kidney_b + M_Kidney_t + M_Kidney_pc) / W_Kidney

    # -------------------------------------6. Spleen------
    # -------------------------------------6.1 Spleen blood ------

    R_Spleen_b = (
        Q_Spleen * ((MBA / (W_Plasma * 0.2)) - (M_Spleen_b / W_Spleen_b))
        - (DLC_Spleen * Q_Spleen) * (M_Spleen_b / W_Spleen_b)
        + ((DLC_Spleen * Q_Spleen) * ((M_Spleen_t) / W_Spleen_t)) / P_Spleen
    )
    dxdt_M_Spleen_b = R_Spleen_b

    # -------------------------------------6.2 Spleen Tissue ------
    R_Spleen_t = (
        (DLC_Spleen * Q_Spleen) * (M_Spleen_b / W_Spleen_b)
        - ((DLC_Spleen * Q_Spleen) * ((M_Spleen_t) / W_Spleen_t)) / P_Spleen
        - K_up_Spleen * M_Spleen_t
        + K_release_Spleen * M_Spleen_pc
    )
    dxdt_M_Spleen_t = R_Spleen_t

    # -------------------------------------6.2 Spleen pc ------
    R_Spleen_pc = K_up_Spleen * M_Spleen_t - K_release_Spleen * M_Spleen_pc
    dxdt_M_Spleen_pc = R_Spleen_pc
    C_Spleen = (M_Spleen_b + M_Spleen_t + M_Spleen_pc) / W_Spleen

    # -----------------------------------7. Liver--------------------
    # -----------------------------------7.1 Liver  blood --------------------

    R_Liver_b = (
        (Q_Liver - Q_GI - Q_Spleen) * (MBA / (W_Plasma * 0.2))
        - (Q_Liver) * (M_Liver_b / W_Liver_b)
        + Q_GI * (M_GI_b / W_GI_b)
        + Q_Spleen * (M_Spleen_b / W_Spleen_b)
        - (DLC_Liver * (Q_Liver)) * (M_Liver_b / W_Liver_b)
        + ((DLC_Liver * (Q_Liver)) * (M_Liver_t / W_Liver_t)) / P_Liver
    )
    dxdt_M_Liver_b = R_Liver_b

    # ---------------------------------------7.2 Liver Tissue------------
    R_Liver_t = (
        (DLC_Liver * (Q_Liver)) * (M_Liver_b / W_Liver_b)
        - ((DLC_Liver * (Q_Liver)) * (M_Liver_t / W_Liver_t)) / P_Liver
        - K_up_Liver * M_Liver_t
        + K_release_Liver * M_Liver_pc
        - Kbile * (M_Liver_t / W_Liver_t)
    )
    R_bile = Kbile * (M_Liver_t / W_Liver_t)
    dxdt_M_bile = R_bile
    dxdt_M_Liver_t = R_Liver_t

    # -------------------------------------7.3 Liver pc compartment----------------
    R_Liver_pc = K_up_Liver * M_Liver_t - K_release_Liver * M_Liver_pc
    dxdt_M_Liver_pc = R_Liver_pc
    C_Liver = (M_Liver_b + M_Liver_t + M_Liver_pc) / W_Liver

    # # --------------------------8. GI-----------------------
    # # --------------------------8. 1 GI blood -----------------------
    R_GI_b = (
        Q_GI * ((MBA / (W_Plasma * 0.2)) - (M_GI_b / W_GI_b))
        - (DLC_GI * Q_GI) * (M_GI_b / W_GI_b)
        + ((DLC_GI * Q_GI) * (M_GI_t / W_GI_t)) / P_GI
        - Kfecal * M_GI_b
    )
    dxdt_M_fecal = Kfecal * M_GI_b
    dxdt_M_GI_b = R_GI_b

    # --------------------------8. 1 GI tissue -----------------------
    # dxdt_M_Oral = M_Oral / t
    R_GI_t = (
        # dxdt_M_Oral +
        (DLC_GI * Q_GI) * (M_GI_b / W_GI_b)
        - ((DLC_GI * Q_GI) * (M_GI_t / W_GI_t)) / P_GI
        - K_up_GI * M_GI_t
        + K_release_GI * M_GI_pc
    )
    dxdt_M_GI_t = R_GI_t

    # --------------------------8.2 GI pc -----------------------
    R_GI_pc = K_up_GI * M_GI_t - K_release_GI * M_GI_pc
    dxdt_M_GI_pc = R_GI_pc
    C_GI = (M_GI_b + M_GI_t + M_GI_pc) / W_GI

    dydt = [
        # dxdt_MIV,
        # dxdt_M_Oral,
        dxdt_MBA,
        dxdt_MBV,
        dxdt_M_Lung_b,
        dxdt_M_Lung_t,
        dxdt_M_Lung_pc,
        dxdt_M_Brain_b,
        dxdt_M_Brain_t,
        dxdt_M_Rest_b,
        dxdt_M_Rest_t,
        dxdt_M_Kidney_b,
        dxdt_M_Kidney_t,
        dxdt_M_Kidney_pc,
        dxdt_M_Spleen_b,
        dxdt_M_Spleen_t,
        dxdt_M_Spleen_pc,
        dxdt_M_Liver_b,
        dxdt_M_Liver_t,
        dxdt_M_Liver_pc,
        dxdt_M_urine,
        dxdt_M_bile,
        dxdt_M_GI_b,
        dxdt_M_GI_t,
        dxdt_M_GI_pc,
        dxdt_M_fecal,
    ]
    return dydt
