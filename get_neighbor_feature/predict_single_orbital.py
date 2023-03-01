import csv
import json
import os
import random
import argparse
import sys
import time
import tqdm
from configparser import ConfigParser

import numpy as np
import torch
from pymatgen import Structure

from deeph import get_graph, DeepHKernal, collate_fn, MaskMSELoss

parser = argparse.ArgumentParser(description='Predict Hamiltonian')
parser.add_argument('--trained_model_dir', type=str,
                    # default='/home/lihe/DeepH-pack/result/graphene_1500_5_5/plain/2021-06-28_13-03-06',
                    default='/home/lihe/DeepH-pack/result/graphene_1500_5_5/LayerNorm/2021-06-28_17-15-58',
                    help='path of trained model')
parser.add_argument('--input_dir', type=str,
                    # default='/home/lihe/DeepH/graphene/graphene_0/4136',
                    default='/home/lihe/DeepH/graphene/graphene_0/0',
                    help='')
parser.add_argument('--output_dir', type=str,
                    # default='/home/lihe/DeepH/graphene/graphene_0/predict_output_4136',
                    default='/home/lihe/DeepH/graphene/graphene_0/predict_output',
                    help='')
parser.add_argument('--disable_cuda', action='store_true',help='Disable CUDA')
parser.add_argument(
    '--interface',
    type=str,
    default='h5',
    choices=['h5', 'npz'])
parser.add_argument('--huge_structure', type=bool, default=False, help='')
args = parser.parse_args()


if __name__ == '__main__':
    assert os.path.exists(os.path.join(args.trained_model_dir, 'best_model.pkl'))
    assert os.path.exists(os.path.join(args.trained_model_dir, 'config.ini'))
    assert os.path.exists(os.path.join(args.trained_model_dir, 'src'))

    os.makedirs(args.output_dir, exist_ok=True)

    config = ConfigParser()
    config.read(os.path.join(args.trained_model_dir, 'config.ini'))
    config.set('basic', 'save_dir', os.path.join(args.output_dir))
    config.set('basic', 'disable_cuda', str(args.disable_cuda))
    config.set('basic', 'save_to_time_folder', 'True')
    kernal = DeepHKernal(config, args.trained_model_dir)


    # orbital = json.loads(config.get('basic', 'orbital'))
    # num_orbital = len(orbital)
    # assert num_orbital == 1
    #
    # args.dtype = config.get('hyperparameter', 'dtype')
    # args.num_l = config.getint('network', 'num_l')
    # args.radius = config.getfloat('graph', 'radius')
    # args.max_num_nbr = config.getint('graph', 'max_num_nbr')
    # args.seed = config.getint('basic', 'seed')
    #
    # if not args.disable_cuda:
    #     args.device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')
    # else:
    #     args.device = torch.device('cpu')
    # if args.dtype == 'float32':
    #     default_dtype = np.float32
    #     default_dtype_torch = torch.float32
    # elif args.dtype == 'float64':
    #     default_dtype = np.float64
    #     default_dtype_torch = torch.float64
    # else:
    #     raise ValueError('Unknown dtype: {}'.format(args.dtype))
    # np.seterr(all='raise')
    # np.seterr(under='warn')
    # np.set_printoptions(precision=8, linewidth=160)
    # torch.set_default_dtype(default_dtype_torch)
    # torch.set_printoptions(precision=8, linewidth=160, threshold=np.inf)
    # if not args.seed:
    #     args.seed = np.random.randint(1, 10 ** 8)
    # np.random.seed(args.seed)
    # torch.manual_seed(args.seed)
    # torch.cuda.manual_seed_all(args.seed)
    # random.seed(args.seed)
    # torch.backends.cudnn.benchmark = False
    # torch.backends.cudnn.deterministic = True
    # print_args(args)
    #
    # model = HGNN(
    #     num_elements=config.getint('basic', 'max_element') + 1,
    #     in_atom_fea_len=config.getint('network', 'atom_fea_len'),
    #     in_edge_fea_len=config.getint('network', 'edge_fea_len'),
    #     num_orbital=num_orbital,
    #     distance_expansion=config.get('network', 'distance_expansion'),
    #     gauss_stop=config.getfloat('network', 'gauss_stop'),
    #     if_exp=config.getboolean('network', 'if_exp'),
    #     if_MultipleLinear=config.getboolean('network', 'if_MultipleLinear'),
    #     if_edge_update=config.getboolean('network', 'if_edge_update'),
    #     normalization=config.get('network', 'normalization')
    # )
    #
    # model_parameters = filter(lambda p: p.requires_grad, model.parameters())
    # params = sum([np.prod(p.size()) for p in model_parameters])
    # print("The model you built has: %d parameters" % params)
    # model.to(args.device)

    with torch.no_grad():
        input_dir = args.input_dir
        structure = Structure(np.loadtxt(os.path.join(args.input_dir, 'lat.dat')).T,
                              np.loadtxt(os.path.join(args.input_dir, 'element.dat')),
                              np.loadtxt(os.path.join(args.input_dir, 'site_positions.dat')).T,
                              coords_are_cartesian=True,
                              to_unit_cell=False)
        cart_coords = torch.tensor(structure.cart_coords)
        frac_coords = torch.tensor(structure.frac_coords)
        numbers = torch.tensor(structure.atomic_numbers)
        structure.lattice.matrix.setflags(write=True)
        lattice = torch.tensor(structure.lattice.matrix)
        inv_lattice = torch.inverse(lattice).type(torch.get_default_dtype())

        if os.path.exists(os.path.join(input_dir, 'graph.pkl')):
            data = torch.load(os.path.join(input_dir, 'graph.pkl'))
            print(f"Load processed graph from {os.path.join(input_dir, 'graph.pkl')}")
        else:
            begin = time.time()
            data = get_graph(cart_coords, frac_coords, numbers, 0,
                             r=kernal.config.getfloat('graph', 'radius'),
                             max_num_nbr=kernal.config.getint('graph', 'max_num_nbr'),
                             numerical_tol=1e-8, lattice=lattice, default_dtype_torch=torch.get_default_dtype(),
                             tb_folder=args.input_dir, interface=args.interface,
                             num_l=kernal.config.getint('network', 'num_l'),
                             huge_structure=args.huge_structure)
            torch.save(data, os.path.join(input_dir, 'graph.pkl'))
            print(f"Save processed graph to {os.path.join(input_dir, 'graph.pkl')}, cost {time.time() - begin} seconds")

        dataset_mask = kernal.make_mask([data])
        batch, subgraph = collate_fn(dataset_mask)
        sub_atom_idx, sub_edge_idx, sub_edge_ang, sub_index = subgraph

        checkpoint = torch.load(
            os.path.join(args.trained_model_dir, 'best_model.pkl'),
            map_location=torch.device('cpu')
        )
        print("=> load best checkpoint (epoch {})".format(checkpoint['epoch']))
        kernal.model.load_state_dict(checkpoint['state_dict'])
        kernal.model.eval()

        output = kernal.model(batch.x.to(kernal.device), batch.edge_index.to(kernal.device),
                              batch.edge_attr.to(kernal.device),
                              batch.batch.to(kernal.device),
                              sub_atom_idx.to(kernal.device), sub_edge_idx.to(kernal.device),
                              sub_edge_ang.to(kernal.device), sub_index.to(kernal.device),
                              huge_structure=args.huge_structure)

        label = batch.label
        mask = batch.mask
        output = output.cpu().reshape(label.shape)
        loss = MaskMSELoss()(output, label, mask)
        print(f'loss: {loss}')

        error = np.full((batch.edge_attr.shape[0]), np.nan)

        edge_stru_index = torch.squeeze(batch.batch[batch.edge_index[0]]).numpy()
        edge_slices = torch.tensor(batch.__slices__['x'])[edge_stru_index].view(-1, 1)
        atom_ids = torch.squeeze(batch.edge_index.T - edge_slices).tolist()
        atomic_numbers = torch.squeeze(batch.x[batch.edge_index.T]).tolist()
        edge_infos = torch.squeeze(batch.edge_attr[:, :7].detach().cpu()).tolist()

        with open(os.path.join(kernal.config.get('basic', 'save_dir'), 'error_distance.csv'), 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['index', 'atom_id', 'atomic_number', 'dist', 'atom1_x', 'atom1_y', 'atom1_z',
                                 'atom2_x', 'atom2_y', 'atom2_z']
                            + ['target'] * kernal.num_orbital + ['pred'] * kernal.num_orbital)
            for index_edge in range(batch.edge_attr.shape[0]):
                R = torch.round(batch.edge_attr[index_edge, 4:7] @ inv_lattice - batch.edge_attr[index_edge,
                                                                                 7:10] @ inv_lattice).int().tolist()
                i, j = batch.edge_index[:, index_edge]
                error[index_edge] = abs(output[index_edge].item() - label[index_edge].item())
                writer.writerow([
                    index_edge,
                    atom_ids[index_edge],
                    atomic_numbers[index_edge],
                    *(edge_infos[index_edge]),
                    label[index_edge].item(),
                    output[index_edge].item()
                ])




