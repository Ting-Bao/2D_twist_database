from pymatgen.core.structure import Structure
import torch, random, torch.nn.functional as F, numpy as np
from torch_geometric.data import Data
from torch_geometric.nn import CGConv

def get_data(filename, data_expand):
    st = Structure.from_file(filename)
    st_lattice = st.lattice.matrix
    st_Z = st.atomic_numbers
    st_fc = st.frac_coords
    st_cc = st.cart_coords
    st_nb = st.get_all_neighbors(r = 7.4)

    # x = torch.tensor(np.vstack((st_Z,st_fc.T)).T) # graph nodes
    x = torch.tensor([st_Z]).T  # graph nodes (atomic numbers)

    # graph edge index & edge feature (distance)
    edge_1, edge_2, dist_12 = [], [], []
    for node_idx in range(len(st_nb)): # node atoms
        for nb_idx in range(len(st_nb[node_idx])): # neighbor atoms for each node
            
            nb_atom_idx = st_nb[node_idx][nb_idx].index # index of each neighbor
            edge_1.append(node_idx)
            edge_2.append(nb_atom_idx)

            nb_dist = st_nb[node_idx][nb_idx][1] # distance between neighbor and node
            dist_12.append(nb_dist)
    edge_12 = torch.tensor([edge_1, edge_2])     # edge index
    dist_12 = torch.tensor([dist_12]).T.float()  # edge feature (atom-pair distance)

    # expand edge feature by Gaussian basis
    u_k = np.linspace(0.0, 7.0, data_expand[1])
    dist_expand = torch.zeros([len(dist_12), data_expand[1]])
    for i in range(len(dist_12)):
        for j in range(data_expand[1]):
            e_k = np.exp(-0.1 * (float(dist_12[i]) - u_k[j])**2)
            dist_expand[i,j] = e_k

    # embedding node feature (torch.nn.embedding?)
    emb = torch.nn.Embedding(int(max(x)+8), data_expand[0])
    x_emb = emb(x)
    # print(x_emb)
    node_expand = torch.zeros([len(x), data_expand[0]])
    for m in range(len(x)):
        for n in range(data_expand[0]):
            node_expand[m,n] = x_emb[m,0,n]
    # x_final = torch.vstack([x_emb[0,0], x_emb[1,0], x_emb[2,0]])

    data = Data(x=node_expand, edge_index=edge_12, edge_attr=dist_expand)
    # print(data)
    return data

class Net(torch.nn.Module):
    def __init__(self):
        super(Net, self).__init__()
        self.conv1 = CGConv(channels = data_expand[0], # node feature length?
                            dim = data_expand[1], # edge feature length
                            aggr = 'add')
    def forward(self,data):
        atom_fea, edge_idx, edge_fea = data.x.float(), data.edge_index, data.edge_attr
        # print('1')
        atom_fea = self.conv1(x=atom_fea, edge_index=edge_idx, edge_attr=edge_fea)
        atom_fea = F.relu(atom_fea)
        # row, col = edge_idx
        # edge_fea = self.e_lin(torch.cat([atom_fea[row], atom_fea[col], edge_fea], dim=-1))
        return atom_fea#, edge_idx, edge_fea

np.random.seed(42)
torch.manual_seed(42)
torch.cuda.manual_seed_all(42)
random.seed(42)

pred_filename = 'POSCAR_12-6'
data_expand = [64, 128]
pred_data = get_data(pred_filename, data_expand)
model = Net().to('cpu')
pred_data = pred_data.to('cpu')
model.eval()
# print(model)
pred_out = model(pred_data)
print(np.shape(pred_out))
print(sum(torch.cosine_similarity(pred_out, pred_out))/np.shape(pred_out)[0])

similarity_list = []
for i in range(1):
    for j in range(3):
        # set_filename = './config/' + str(i) + '_' + str(j) + '/POSCAR_crystal'
        set_filename = './config/' + str(i) + '_' + str(j) + '/crystal.cif'
        set_data = get_data(set_filename, data_expand)
        set_data = set_data.to('cpu')
        model.eval()
        set_out = model(set_data)
        # cos_similarity = torch.cosine_similarity(set_out, pred_out)
        cos_similarity = sum(torch.cosine_similarity(pred_out, set_out))/np.shape(pred_out)[0]
        print(cos_similarity)
        similarity_list.append(cos_similarity)
print(sum(similarity_list))