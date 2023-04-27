import os
from pymatgen.util.plotting import periodic_table_heatmap, pretty_plot
import matplotlib.pyplot as plt
import chemparse

# reference:
# https://pymatgen.org/pymatgen.util.plotting.html

item_path = 'data/infer_calculation'

def get_db_elecount(path):
    items = os.listdir(item_path)
    items = [i.split('-')[0] for i in items if os.path.isdir(os.path.join(item_path,i))]
    # manually add the heterobilayer 
    items.extend(['C2','BN','Bi2Se3','Bi2Te3','MoS2','MoSe2','MoTe2','WSe2','WS2','WSe2'])
    periodic_table = { 0: 'n', 1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl', 18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti', 23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co', 28: 'Ni', 29: 'Cu', 30: 'Zn', 31: 'Ga', 32: 'Ge', 33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb', 38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo', 43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag', 48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te', 53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La', 58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm', 63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho', 68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu', 72: 'Hf', 73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir', 78: 'Pt', 79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb', 83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr', 88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa', 92: 'U' }
    elemental_data={}
    for i in range(1,93):
        ele_name = periodic_table[i]
        elemental_data[ele_name] = 0
    for i in items:
        parsed_formula = chemparse.parse_formula(i)
        for j in parsed_formula.keys():
            elemental_data[j] += 1
    return elemental_data


def main():
    elemental_data = get_db_elecount(path = item_path)    
    # elemental_data={'Fe': 4.2, 'O': 5.0}
    plt.figure(figsize=(5,12))
    periodic_table_heatmap(elemental_data=elemental_data, cmap='tab20', cbar_label= 'Count', max_row=6, value_format='%d',show_plot=False)
    # to change ticks, go to the periodic_table_heatmap source code
    # change line 279 to cbar = fig.colorbar(heatmap,ticks=[0,2,4,6,8,10,12,14,16,18,20])
    pretty_plot(plt=plt,dpi=800)
    
    plt.savefig('figure/periodic_table.png',dpi = 1000)
    plt.savefig('figure/periodic_table.svg',dpi = 1000)
    
    print('finished')

if __name__=='__main__':
    main()
    