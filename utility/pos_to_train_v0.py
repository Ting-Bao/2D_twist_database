import numpy as np
from pymatgen.core.structure import Structure
from string import digits

basis_path = '/home/xurz/POT/opmx_basis.txt'

def find_basis(x, b):
    # x = element name (str), b = basis accuracy (int, 1/2/3)
    # txt file format: Element_name  ZVAL  POT_name  NELECT  Quick  Standard  Precise
    pot_name = x + '_PBE19'
    # with open('opmx_basis.txt', 'r') as opmx_basis:
    with open(basis_path, 'r') as opmx_basis:
        lines = opmx_basis.readlines()
        for line in lines:
            if pot_name in line:
                ZVAL = int(line.split()[1])     # atomic number       (6)
                POT = str(line.split()[2])      # potential name      (C_PBE19)
                NELECT = float(line.split()[3]) # number of electrons (4.0)
                BASIS = str(line.split()[3+b])  # basis setting       (C6.0-s2p2d1)
    return ZVAL, POT, NELECT, BASIS

def pos_to_train_input(file, prec):  # generate key input for train.ini orbital setting
    # more option will be added to account for complex situation
    # prec = 1/2/3, type = 1, 2
    str_raw = Structure.from_file(file)  # read in structure
    formula = str_raw.formula.split()
    element = list(map(lambda formula: formula.translate(str.maketrans('','',digits)), formula))
    n_element = len(element)  # !!!

    atomic_number, n_orbital, n_orbital_each, orbitals_array = [], [], [], []
    basis_array = np.zeros([n_element, 4], dtype=int)
    orbital_info = np.array([1, 3, 5, 7])  # number for each orbitals, example: 3 for p, 5 for d
    
    for i in range(n_element):
        z_val, pot, n_elect, basis = find_basis(element[i], prec)
        atomic_number.append(z_val)  # !!!

        basis = basis.split('-')[1]
        for j in range(0, len(basis), 2):  # example: s3p2d1f1, output = [ns,np,nd,nf]
            basis_array[i, int(j/2)] = int(basis[j+1])
        
        n_orbital.append(np.dot(basis_array[i], orbital_info))  # tot num of orbitals for each element (list)
        n_orbital_each.append(basis_array[i] * orbital_info)  # num of each orbitals for each element (list)
        length = np.dot(basis_array[i], orbital_info)
        orbitals_array.append(list(range(length)))  # orbitals array for each element (list) !

    N_orbital = max(n_orbital)  # !!!
    n_orbital_each = np.array(n_orbital_each)

    basis_col_max = basis_array.max(axis=0)
    n_orbital_each_max = np.array(n_orbital_each).max(axis=0)
    
    # insert '-1' in orbitals_array
    for m in range(n_element):
        for n in range(4):
            if n_orbital_each[m, n] < n_orbital_each_max[n] and n == 0:
                for l in range(n_orbital_each_max[n] - n_orbital_each[m, n]):
                    orbitals_array[m].insert(n_orbital_each[m, n], -1)

            elif n_orbital_each[m, n] < n_orbital_each_max[n] and n != 0:
                for l in range(n_orbital_each_max[n] - n_orbital_each[m, n]):
                    orbitals_array[m].insert(n_orbital_each[m, n] + sum(n_orbital_each_max[0:n]), -1)

    return n_element, atomic_number, N_orbital, orbitals_array  # type: int, list, int, 2D list

def gen_orbital_set(n_element, atomic_number, N_orbital, orbitals, type):  # type = 1/2 (all/multi)
    # generate orbital setting for train.ini
    orbital_str = '['
    first_flag2 = True
    num = 0
    orbital_set, atomic_number_set = [], []

    if type == 1:  # training orbital hoppings in all-in-one style (all)
        for model_i, model_j in ((model_i, model_j) for model_i in range(N_orbital) for model_j in range(N_orbital)):
            first_flag = True
            if first_flag2: 
                orbital_str += ''
                first_flag2 = False
            else:
                orbital_str += ', '
            
            for element_i in range(n_element):
                for element_j in range(n_element):
                    orbital_i = orbitals[element_i][model_i]
                    orbital_j = orbitals[element_j][model_j]
                    # print(element_i,model_i,'|',element_j,model_j)
                    if orbital_i != -1 and orbital_j != -1:
                        if first_flag:
                            orbital_str += f'{{"{atomic_number[element_i]} {atomic_number[element_j]}": [{orbital_i}, {orbital_j}]}}'
                            first_flag = False
                        else:
                            orbital_str += f', {{"{atomic_number[element_i]} {atomic_number[element_j]}": [{orbital_i}, {orbital_j}]}}'
            
            num += 1
        orbital_str += ']'
        orbital_str = orbital_str.replace('{}, ', '')
        orbital_str = orbital_str.replace(', {}', '')
        orbital_set.append(orbital_str)
        
    if type == 2:  # training orbital hoppings in separate style (separate)
        for idx_left in range(n_element):
            for idx_right in range(n_element):
                atomic_number_set.append([atomic_number[idx_left], atomic_number[idx_right]])  # only meaningful for type 2
                
                for model_i, model_j in ((model_i, model_j) for model_i in range(N_orbital) for model_j in range(N_orbital)):
                    if first_flag2:
                        orbital_str += '{'
                        first_flag2 = False
                    else:
                        orbital_str += ', {'
                    first_flag = True

                    orbital_i = orbitals[idx_left][model_i]
                    orbital_j = orbitals[idx_right][model_j]
                    if orbital_i != -1 and orbital_j != -1:
                        if first_flag:
                            orbital_str += f'"{atomic_number[idx_left]} {atomic_number[idx_right]}": [{orbital_i}, {orbital_j}]'
                            first_flag = False
                        else:
                            orbital_str += f', "{atomic_number[idx_left]} {atomic_number[idx_right]}": [{orbital_i}, {orbital_j}]'
                    orbital_str += '}'
                    
                    num += 1
                orbital_str += ']'
                orbital_str = orbital_str.replace('{}, ', '')
                orbital_str = orbital_str.replace(', {}', '')
                orbital_set.append(orbital_str)

    return orbital_set, atomic_number_set  # type: both 2D list

##############################################################################################
def gen_train_ini(orbital_set, atomic_number_set, graph_dir, save_dir, raw_dir, dateset_name, cuda, eporchs, lr):
    # generate & write train.ini file
    if len(orbital_set) == 1:  # all-in-one case
        orbital_set_1 = orbital_set[0]
        in_file = fr"""[basic]
graph_dir = {graph_dir}
save_dir = {save_dir}
raw_dir = {raw_dir}
dataset_name = {dateset_name}
pyg_version =
;choices = ['h5', 'npz']
interface = h5
target = hamiltonian
disable_cuda = False
device = cuda:{cuda}
; -1 for cpu_count(logical=False) // torch.cuda.device_count()
num_threads = -1
save_to_time_folder = True
save_csv = True
tb_writer = True
seed = 42
multiprocessing = False
orbital = {orbital_set_1}
O_component = H
energy_component = summation
max_element = -1
statistics = False
normalizer = False
boxcox = False

[graph]
radius = -1.0
max_num_nbr = 0
create_from_DFT = True
if_lcmp_graph = True
separate_onsite = False
new_sp = False

[train]
epochs = {eporchs}
pretrained =
resume =
train_ratio = 0.6
val_ratio = 0.2
test_ratio = 0.2
early_stopping_loss = {0.0000005}
early_stopping_loss_epoch = [0.000000, 500]
revert_then_decay = True
revert_threshold = 100
revert_decay_epoch = [500, 2000, 6000]
revert_decay_gamma = [0.4, 0.5, 0.5]
clip_grad = True
clip_grad_value = 4.2
switch_sgd = False
switch_sgd_lr = 1e-4
switch_sgd_epoch = -1

[hyperparameter]
batch_size = 1
dtype = float32
;choices = ['sgd', 'sgdm', 'adam', 'lbfgs']
optimizer = adam
;initial learning rate
learning_rate = {lr}
;choices = ['', 'MultiStepLR', 'ReduceLROnPlateau', 'CyclicLR']
lr_scheduler =
lr_milestones = []
momentum = 0.9
weight_decay = 0
criterion = MaskMSELoss
retain_edge_fea = True
lambda_Eij = 0.0
lambda_Ei = 0.1
lambda_Etot = 0.0

[network]
atom_fea_len = 64
edge_fea_len = 128
gauss_stop = 9.0
;The number of angular quantum numbers that spherical harmonic functions have
num_l = 5
aggr = add
distance_expansion = GaussianBasis
if_exp = True
if_MultipleLinear = False
if_edge_update = True
if_lcmp = True
normalization = LayerNorm
;choices = ['CGConv', 'GAT', 'PAINN']
atom_update_net = CGConv
"""
        with open('train.ini', 'w') as f:
            f.write(in_file)

    if len(orbital_set) > 1:  # separate/multi case
        for i in range(len(orbital_set)):
            orbital_set_1 = orbital_set[i]
            file_name = 'train_' + str(atomic_number_set[i][0]) + '-' + str(atomic_number_set[i][1]) + '.ini'
            in_file = fr"""[basic]
graph_dir = {graph_dir}
save_dir = {save_dir}
raw_dir = {raw_dir}
dataset_name = {dateset_name}
pyg_version =
;choices = ['h5', 'npz']
interface = h5
target = hamiltonian
disable_cuda = False
device = cuda:{cuda}
; -1 for cpu_count(logical=False) // torch.cuda.device_count()
num_threads = -1
save_to_time_folder = True
save_csv = True
tb_writer = True
seed = 42
multiprocessing = False
orbital = {orbital_set_1}
O_component = H
energy_component = summation
max_element = -1
statistics = False
normalizer = False
boxcox = False

[graph]
radius = -1.0
max_num_nbr = 0
create_from_DFT = True
if_lcmp_graph = True
separate_onsite = False
new_sp = False

[train]
epochs = {eporchs}
pretrained =
resume =
train_ratio = 0.6
val_ratio = 0.2
test_ratio = 0.2
early_stopping_loss = {0.0000005}
early_stopping_loss_epoch = [0.000000, 500]
revert_then_decay = True
revert_threshold = 100
revert_decay_epoch = [500, 2000, 6000]
revert_decay_gamma = [0.4, 0.5, 0.5]
clip_grad = True
clip_grad_value = 4.2
switch_sgd = False
switch_sgd_lr = 1e-4
switch_sgd_epoch = -1

[hyperparameter]
batch_size = 1
dtype = float32
;choices = ['sgd', 'sgdm', 'adam', 'lbfgs']
optimizer = adam
;initial learning rate
learning_rate = {lr}
;choices = ['', 'MultiStepLR', 'ReduceLROnPlateau', 'CyclicLR']
lr_scheduler =
lr_milestones = []
momentum = 0.9
weight_decay = 0
criterion = MaskMSELoss
retain_edge_fea = True
lambda_Eij = 0.0
lambda_Ei = 0.1
lambda_Etot = 0.0

[network]
atom_fea_len = 64
edge_fea_len = 128
gauss_stop = 9.0
;The number of angular quantum numbers that spherical harmonic functions have
num_l = 5
aggr = add
distance_expansion = GaussianBasis
if_exp = True
if_MultipleLinear = False
if_edge_update = True
if_lcmp = True
normalization = LayerNorm
;choices = ['CGConv', 'GAT', 'PAINN']
atom_update_net = CGConv
"""
            with open(file_name, 'w') as f:
                f.write(in_file)