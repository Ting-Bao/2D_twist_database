import csv
import os
import argparse
import time
from configparser import ConfigParser

import numpy as np
import torch
from pymatgen.core.structure import Structure

from src.deeph import get_graph, DeepHKernal, collate_fn


def main():
    parser = argparse.ArgumentParser(description='Predict Hamiltonian')
    parser.add_argument('--trained_model_dir', type=str,
                        # default='/home/lihe/DeepH-pack/result/graphene_1500_5_5/plain/2021-06-28_13-03-06',
                        # default='/home/lihe/DeepH-pack/result/graphene_1500_5_5/LayerNorm/2021-06-28_17-15-58',
                        default='/home/xurz/work/moire-twist/sn_vdw12_no_relax/4-5_7.34/H_predict/model_orbital_diag/',
                        help='path of trained model')
    parser.add_argument('--input_dir', type=str,
                        # default='/home/lihe/DeepH/graphene/graphene_0/4136',
                        # default='/home/lihe/DeepH/graphene/graphene_0/0',
                        default='/home/xurz/work/moire-twist/sn_vdw12_no_relax/4-5_7.34/processed/',
                        help='')
    parser.add_argument('--output_dir', type=str,
                        # default='/home/lihe/DeepH/graphene/graphene_0/predict_output_4136',
                        default='/home/xurz/work/moire-twist/sn_vdw12_no_relax/4-5_7.34/H_predict/model_orbital_diag/',
                        help='')
    parser.add_argument('--disable_cuda', action='store_true',help='Disable CUDA')
    parser.add_argument('--save_csv', action='store_true',help='Save the result for each edge in csv format')
    parser.add_argument(
        '--interface',
        type=str,
        default='h5',
        choices=['h5', 'npz'])
    parser.add_argument('--huge_structure', type=bool, default=False, help='')
    args = parser.parse_args()

    assert os.path.exists(os.path.join(args.trained_model_dir, 'best_model.pkl'))
    assert os.path.exists(os.path.join(args.trained_model_dir, 'config.ini'))

    os.makedirs(args.output_dir, exist_ok=True)

    config = ConfigParser()
    config.read(os.path.join(args.trained_model_dir, 'config.ini'))
    config.set('basic', 'save_dir', os.path.join(args.output_dir))
    config.set('basic', 'disable_cuda', str(args.disable_cuda))
    config.set('basic', 'save_to_time_folder', 'False')
    kernal = DeepHKernal(config)
    checkpoint = kernal.build_model(args.trained_model_dir)

    with torch.no_grad():
        input_dir = args.input_dir
        structure = Structure(np.loadtxt(os.path.join(args.input_dir, 'lat.dat')).T,
                              np.loadtxt(os.path.join(args.input_dir, 'element.dat')),
                              np.loadtxt(os.path.join(args.input_dir, 'site_positions.dat')).T,
                              coords_are_cartesian=True,
                              to_unit_cell=False)
        cart_coords = torch.tensor(structure.cart_coords, dtype=torch.get_default_dtype())
        frac_coords = torch.tensor(structure.frac_coords, dtype=torch.get_default_dtype())
        numbers = kernal.Z_to_index[torch.tensor(structure.atomic_numbers)]
        structure.lattice.matrix.setflags(write=True)
        lattice = torch.tensor(structure.lattice.matrix, dtype=torch.get_default_dtype())
        inv_lattice = torch.inverse(lattice)

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
                             create_from_DFT=kernal.config.getboolean('graph', 'create_from_DFT'),
                             if_lcmp_graph=kernal.config.getboolean('graph', 'if_lcmp_graph'),
                             target=kernal.config.get('basic', 'target'), huge_structure=args.huge_structure)
            torch.save(data, os.path.join(input_dir, 'graph.pkl'))
            print(f"Save processed graph to {os.path.join(input_dir, 'graph.pkl')}, cost {time.time() - begin} seconds")

        dataset_mask = kernal.make_mask([data])
        batch, subgraph = collate_fn(dataset_mask)
        sub_atom_idx, sub_edge_idx, sub_edge_ang, sub_index = subgraph

        output = kernal.model(batch.x.to(kernal.device), batch.edge_index.to(kernal.device),
                              batch.edge_attr.to(kernal.device),
                              batch.batch.to(kernal.device),
                              sub_atom_idx.to(kernal.device), sub_edge_idx.to(kernal.device),
                              sub_edge_ang.to(kernal.device), sub_index.to(kernal.device),
                              huge_structure=args.huge_structure)

        label = batch.label
        mask = batch.mask
        output = output.cpu().reshape(label.shape)

        assert label.shape == output.shape == mask.shape
        mse = torch.pow(label - output, 2)
        mae = torch.abs(label - output)

        for index_orb, orbital_single in enumerate(kernal.orbital):
            if index_orb != 0:
                print('================================================================')
            print('orbital:', orbital_single)
            if kernal.spinful == False:
                print(f'mse: {torch.masked_select(mse[:, index_orb], mask[:, index_orb]).mean().item()}, '
                      f'mae: {torch.masked_select(mae[:, index_orb], mask[:, index_orb]).mean().item()}')
            else:
                for index_soc, str_soc in enumerate([
                    'left_up_real', 'left_up_imag', 'right_down_real', 'right_down_imag',
                    'right_up_real', 'right_up_imag', 'left_down_real', 'left_down_imag',
                ]):
                    if index_soc != 0:
                        print('----------------------------------------------------------------')
                    print(str_soc, ':')
                    index_out = index_orb * 8 + index_soc
                    print(f'mse: {torch.masked_select(mse[:, index_out], mask[:, index_out]).mean().item()}, '
                          f'mae: {torch.masked_select(mae[:, index_out], mask[:, index_out]).mean().item()}')

        if args.save_csv:
            edge_stru_index = torch.squeeze(batch.batch[batch.edge_index[0]]).numpy()
            edge_slices = torch.tensor(batch.__slices__['x'])[edge_stru_index].view(-1, 1)
            atom_ids = torch.squeeze(batch.edge_index.T - edge_slices).tolist()
            atomic_numbers = torch.squeeze(batch.x[batch.edge_index.T]).tolist()
            edge_infos = torch.squeeze(batch.edge_attr[:, :7].detach().cpu()).tolist()

            with open(os.path.join(kernal.config.get('basic', 'save_dir'), 'error_distance.csv'), 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['index', 'atom_id', 'atomic_number', 'dist', 'atom1_x', 'atom1_y', 'atom1_z',
                                     'atom2_x', 'atom2_y', 'atom2_z']
                                + ['target'] * kernal.num_orbital + ['pred'] * kernal.num_orbital + ['mask'] * kernal.num_orbital)
                for index_edge in range(batch.edge_attr.shape[0]):
                    # R = torch.round(batch.edge_attr[index_edge, 4:7] @ inv_lattice - batch.edge_attr[index_edge,
                    #                                                                  7:10] @ inv_lattice).int().tolist()
                    # i, j = batch.edge_index[:, index_edge]
                    writer.writerow([
                        index_edge,
                        atom_ids[index_edge],
                        atomic_numbers[index_edge],
                        *(edge_infos[index_edge]),
                        *(label[index_edge].tolist()),
                        *(output[index_edge].tolist()),
                        *(mask[index_edge].tolist()),
                    ])

if __name__ == '__main__':
    main()
