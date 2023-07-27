def ev_to_jkmol(data_ev):
    """
    This function performs the change of units from eV to kJ/mol
    """
    data_kjmol = (2625.5/27.2107)*data_ev

    return data_kjmol

