repeats = {'LTR': {'Gypsy': 7105, 'Bel-Pao': 6293, 'Copia': 1199, 'NA': 802}, 'Simple_repeat': {'NA': 139838}, 'DNA': {'NA': 9148, 'm3bp': 9075, 'PiggyBac': 1302, 'm8bp': 4352, 'hAT': 3343, 'CACTA': 923, 'Tc1-Mariner': 8745, 'PIF-Harbinger': 5052, 'mTA': 18084, 'Gambol': 2344, 'P': 3789, 'Transib': 229, 'm4bp': 1054, 'Helitron': 387}, 'SINE': {'tSINE': 11145}, 'LINE': {'CR1': 6284, 'R1': 2045, 'Jockey': 4954, 'I': 1022, 'RTE': 5048, 'Outcast': 560, 'L2': 784, 'Loner': 731, 'L1': 383, 'NA': 28}, 'Low_complexity': {'NA': 13141}, 'Undetermined': {'Maque': 191}}

class_colors = {}
for i, cls in enumerate(repeats.keys()):
    # if cls == 'NA':
    #     class_colors[cls] = (0.6, 0.6, 0.6, 1.0)  # RGBA for grey
    # else:
    #     
    class_colors[cls] = base_cmap(i)

print(class_colors)