import json
import pandas as pd
import os

file = 'find_spcaegroup/itemlist_formatted_1col.xlsx'
to_file = 'find_spcaegroup/itemlist_spcaegroup.xlsx'

data = pd.read_excel(file)

data['layer_group'] = [None] * len(data)
for i in range(len(data)):
    mat = data['material'][i]
    mat = mat.replace('$', '').replace('_', '')
    c2db_id = data['c2db_id'][i]
    name = '-'.join([mat, c2db_id])
    # print(name)
    info_path = 'data/c2dbdata/jsondata/{}.all_data.json'.format(name)
    if not os.path.exists(info_path):
        print(name, "not found")
        # data['space_group'][i] = None
    else:
        with open(info_path, 'r', encoding='utf-8') as fp:
            tmp = json.load(fp)

        tmp = tmp['results-asr.structureinfo.json']['kwargs']['data']
        spg = tmp['spacegroup']
        spg_num = tmp['spgnum']
        # use new greped data in c2db
        try:
            layergroup = tmp['layergroup']
            lgnum = tmp['lgnum']
        except:
            print('layer group not found')
        data['space_group'][i] = "{} ({})".format(spg, spg_num)
        data['layer_group'][i] = "{} ({})".format(layergroup, lgnum)
        print(name, "{} ({})".format(spg, spg_num))

data.to_excel(to_file)
