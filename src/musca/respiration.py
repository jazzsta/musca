def m_respiration(dryweight, T2, organ_type):
    """
    Compute mantainance respiration of a component

    :Parameters:
        - 'dryweight' - dry weight of the component (g)
        - 'organ_type' - organ type of the component
        - T2 - (C)

        - MRR - mantainance respiration rate (gC gDM-1 s-1)
        - Q10 - parameter for the impact of temperature on respiration (-)
        - CC - carbon contene (gC gDM-1)
        - T0 - base temperature (C)
        - T1 - (C)

    :Returns:
        value of mantainance respiration for the component.

    :Remark:
        parameters from:
        Supplementary Table S1 in Pallas et al., 2016. Simulation of carbon allocation and organ growth variability in apple tree
        by connecting architectural and source-sink models. Annals of Botany. pp1-14.

    """
    MRR = {}
    MRR["fruit"] = MRRfruit = 0.53         # gC gDM-1 s-1
    MRR["leaf"] = MRRleaf = 0.32
    MRR["nl"] = MRRnl = 0.13
    MRR["current_year_shoot"] = MRRstem = 0.16
    MRR["old_wood"] = MRRoldwood = 0.39
    MRR["new_root"] = MRRnewroot = 0.99
    MRR["root"] =MRRoldroot = 0.39

    Q10 = {}

    Q10["fruit"] = Q10fruit = 1.731
    Q10["leaf"] = Q10leaf = 2.46
    Q10["nl"] = Q10nl = 2.34
    Q10["current_year_shoot"] = Q10stem = 2.34
    Q10["old_wood"] = Q10oldwood = 2.34
    Q10["new_root"] = Q10newroot = 3.1
    Q10["root"] = Q10oldroot = 2.34

    CC = {}
    CC["fruit"] = CCf = 0.40   # gC gDM-1
    CC["leaf"] = CCl = 0.48
    CC["nl"] = CCnl = 0.47
    CC["current_year_shoot"] = CCs = 0.47
    CC["old_wood"] = CCow = 0.47
    CC["new_root"] = CCnr = 0.48
    CC["root"] = CCor = 0.47

    T1 = 20. # C
    T0 = 10. # C

    secs_in_day = 86400. # s day-1
    dryweight = dryweight * 1000. # dry weights from Kg to g
    #  gDM day-1           =   gDM     * gC gDM-1 s-1     / gC gDM-1
    mantainanceRespiration = dryweight * (MRR[organ_type] / CC[organ_type]) * Q10[organ_type]**((T2-T1)/T0)
    mantainanceRespiration *= secs_in_day

    if mantainanceRespiration and organ_type == "fruit":
        print "m_respiration of ", organ_type, "of dry mass ", dryweight, " is ", mantainanceRespiration

    mantainanceRespiration = mantainanceRespiration / 1000. # from g to Kg
    return mantainanceRespiration

def g_respiration(demand, organ_type):
    """
    Compute growth respiration of a component

    :Parameters:
        - demand - dry matter demand, given by the sum of: a function of thermal time and measured growth, and mantainance respiration (gDM)
        - organ type - organ type of the component
        - GRC - the growth respiration coefficient (gC gDM-1)
        - 'CC' - carbon content in an organ type (gC gDM-1)

    :Returns:
        growthRespiration: growth respiration value of the component

    :Remark:
        parameters from:
        Supplementary Table S1 in Pallas et al., 2016. Simulation of carbon allocation and organ growth variability in apple tree
        by connecting architectural and source-sink models. Annals of Botany. pp1-14.

    """
    GRC["fruit"] = 0.06         # gC gDM-1
    GRC["leaf"] = 0.101
    GRC["nl"] = 0.082
    GRC["current_year_shoot"] = 0.089
    GRC["old_wood"] = 0.081
    GRC["new_root"] = 0.082
    GRC["root"] = 0.081

    CC = {}
    CC["fruit"] = CCf = 0.40   # gC gDM-1
    CC["leaf"] = CCl = 0.48
    CC["nl"] = CCnl = 0.47
    CC["current_year_shoot"] = CCs = 0.47
    CC["old_wood"] = CCow = 0.47
    CC["new_root"] = CCnr = 0.48
    CC["root"] = CCor = 0.47

    #    gDM             gDM   *  gC gDM-1        / gC gDM-1
    growthRespiration = demand * (GRC[organ_type] / CC[organ_type])

    if growthRespiration and organ_type == "fruit":
        print("g_respiration of ", organ_type, " is ", growthRespiration)

    return growthRespiration
