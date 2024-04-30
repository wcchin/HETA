"""
after
pip install -e .
"""

import os
import networkx as nx
import matplotlib.pyplot as plt

from heta import HETA as heta

test_data_dir = 'data/net/'
fs = sorted(os.listdir(test_data_dir))
fs = [ f for f in fs if f[-4:]=='.net' ]
#print(fs)
#f = fs[0]
for f in fs:
    print(f)
    fp = os.path.join(test_data_dir, f)
    g = nx.Graph(nx.read_pajek(fp))
    print(g.number_of_edges())
    g, ext_dic, int_dic = heta.bridge_or_bond(g)
    print(ext_dic, int_dic)
    print(heta.fingerprint(g))
    for u, v, d in g.edges(data=True):
        print(u, v, d['type'])
    heta.draw_result(g, layout='spring')
    plt.show()
    break
print('----------done----------')
