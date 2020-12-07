# -*- coding: utf-8 -*-

""" describe real gas behaviour with pressures > 0 bar using the compressibility factor Z """


import pandas as pd
import os as os


# points between who linear relationship is assumed, [[pressure values], [Z factor]]
compr_factors_data = {"H2": [[0, 1000], [1, 1.717]],
                      "C2H4": [[0, 60, 100, 1000], [1, 0.163, 0.304, 2.283]]}

def calc_pressure(m_kg, v_m3, t_k=293.15, fluid="H2", precision = 0.1):
    """
    Calculate the pressure of a gas tank looking up the compressibilty factor Z.

    INPUT:

        **m_kg** (float) -  Gas mass in the tank (kg)

        **v_m3** (float) - Tank's volume = volume of the gas (m^3)

        **t_k** (float) - The gas' temperature which is assumed to be constant (Kelvin)

        **fluid** (str, "H2") - gas's name

    OPTIONAL:

        **precision** (float, 0.01) - maximal error caused by break of calculation loop relative to final pressure

    OUTPUT:

        Calculated pressure (bar)
    """
    r_j_per_mol_k = 8.314472                                                    # universal gas constant
    m_kg_per_mol = calc_molar_mass(fluid) / 1000                                # calculate molar mass, convert to kg
    p_prel = (m_kg * r_j_per_mol_k * t_k) / (m_kg_per_mol * v_m3) / 100000      # calculate preliminary pressure
    while 1:                                                                    # repeat until precise enough
        z = get_z(p_prel, fluid)
        p = (m_kg * r_j_per_mol_k * t_k) / (m_kg_per_mol * v_m3) * z / 100000   # calculate p with Z, convert Pa -> bar
        if abs(1 - p / p_prel) > precision:
            p_prel = p
        else:                                                                   # stop calculation when p precise enough
            return p


def calc_m_kg(p_bar, v_m3, t_k=293.15, fluid="H2", precision=0.1):
    """
       Calculate the gas mass contained in a gas tank with a given pressure looking up the compressibility factor Z.

       INPUT:

           **p_bar** (float) -  Pressure inside the gas tank (bar)

           **v_m3** (float) - Tank's volume = volume of the gas (m^3)

           **t_k** (float) - The gas' temperature which is assumed to be constant (Kelvin)

           **fluid** (str, "H2") - gas's name

       OPTIONAL:

           **precision** (float, 0.01) - maximal error caused by break of calculation loop relative to final pressure

       OUTPUT:

           Calculated State of Charge (%), Calculated gas mass (kg)
       """
    r_j_per_mol_k = 8.314472                                                    # universal gas constant
    m_kg_per_mol = calc_molar_mass(fluid) / 1000                                # calculate molar mass, convert to kg
    z = get_z(p_bar, fluid)                                                     # calculate compressibility factor
    return (p_bar * 100000 * m_kg_per_mol * v_m3) / (r_j_per_mol_k * t_k * z)   # calculate and return gas mass


def get_z(p_bar, fluid):
    """
    Get compressibilty factor Z from tabled data.

    INPUT:

        **p_bar** (float) - Pressure at which Z is to be determined (bar)

        **fluid** (str, "H2") - Fluid for which Z is to be determined

    OUTPUT:


    """
    # if pressure higher than highest predefined value, calculate Z using the last available value
    index = -1
    for i, p in enumerate(compr_factors_data[fluid][0]):
        if p > p_bar:
            if i < 1:
                # <-> p < 0 bar -> outside of range -> raise error
                raise IndexError("value out of range: fcn get_z: pressure can't be under 0 bar")
            index = i
            break
        elif p == p_bar:
            return compr_factors_data[fluid][1][i]
    # linear interpolation: calculate Z
    # 1.:       add to next lowest value
    # 2.,  3.:  calculated share from next lowest pressure to p_bar / to next highest pressure
    # 4.        multiply with difference between next lowest and next highest compressibility factor
    return compr_factors_data[fluid][1][index - 1] + \
           (p_bar - compr_factors_data[fluid][0][index - 1]) / \
           (compr_factors_data[fluid][0][index] - compr_factors_data[fluid][0][index - 1]) * \
           (compr_factors_data[fluid][1][index] - compr_factors_data[fluid][1][index - 1])


def calc_molar_mass(formula="H2"):
    """
    Calculate a basic substance's molar mass from its molecular formular.

    INPUT:

        **formula** (str, "H2") - Substance's molecular formula

    OUTPUT:

        Substance's molar mass (g/mol)
    """

    j = ""                                          # initialize string to store current element symbol
    m_mol = 0                                       # initialize molar mass
    pte = PTE()
    for index, i in enumerate(formula):
        if not (i.isdigit() or i.isupper()):        # at number or uppercase letter, next symbol begins
            j += i
        elif i.isdigit():
            i = int(i)
            if index == 0:
                raise KeyError("Error in calculation of molar mass: molecular formula can't start with a number")
            m_mol += round(pte.lookup(j, "AtomicMass")) * i
            j = ""                                  # reset string to store new symbol
        else:
            if len(j) > 0:
                m_mol += round(pte.lookup(j, "AtomicMass"))
            j = i                                   # reset string: current character i yet to be processed
    if len(j) > 0:                                  # process last element if not yet processed
         m_mol += round(pte.lookup(j, "AtomicMass"))
    return m_mol


class PTE:
    """
    Periodic Table of Elements.
    """

    def __init__(self):
        self.df = pd.read_csv(os.path.join("component_models", "auxiliaries", "PubChemElements_all.csv"),
                              index_col="Symbol")

    def lookup(self, symbol, value):
        """
        Look up values in PTE

        INPUT:

            **symbol** (str) - Symbol of the element to look up values for

            ** value** (str or str[]) - Value / values to be looked up

        OUTPUT:

            Value / values looked up in the PTE
        """
        if isinstance(value, list):
            try:
                return [self.df.loc[symbol, v] for v in value]
            except:
                raise KeyError("Error in calculation of molar mass: element not found in pse")
        try:
            return self.df.loc[symbol, value]
        except:
            msg = "Error in calculation of molar mass: element '{}' not found in pse".format(symbol)
            raise KeyError(msg)


def test(p_bar, t_c, rho_kg_per_m3):
    print("conditions: t = {} °C, rho = {} kg/m3".format(t_c, rho_kg_per_m3))
    p_calc = round(calc_pressure(rho_kg_per_m3, 1, t_k=273.15+t_c), 5)
    print("p = {} bar; expected p = {} bar (diff: {} bar / {} %)"
          .format(p_calc, p_bar, round(p_calc-p_bar, 5), round((1-p_calc/p_bar)*100)))


comp_data = pd.DataFrame([[1, 0, 0.0887],
                          [1, 100, 0.0649],
                          [100, 0, 8.3447],
                          [100, 100, 6.1840],
                          [300, 0, 22.151],
                          [300, 25, 20.537],
                          [300, 100, 16.883],
                          [1000, 0, 52.115],
                          [1000, 25, 49.424],
                          [1000, 100, 42.819]],
                         columns=["p_bar", "t_c", "rho_kg_per_m3"])

[test(x, y, z) for x, y, z in zip(comp_data["p_bar"], comp_data["t_c"], comp_data["rho_kg_per_m3"])]
