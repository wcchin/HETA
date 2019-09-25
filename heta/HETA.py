# -*- coding: utf-8 -*-

import random
import math
import os
import json
import time

import networkx as nx
import scipy
import numpy as np
#import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from tqdm import tqdm
import pathos
#from pathos.multiprocessing import ProcessingPool


EGO_NETWORK = 'ego'
GRAPH_KEY_COMMON_NODES_LIST = 'list_'
GRAPH_KEY_AVG_COMMON_NODES = 'avg'
GRAPH_KEY_STD_COMMON_NODES = 'std'
#TODO: check SUCCESSORS = 'suc'
#TODO: check PREDECESSORS = 'pre'
TIME = 0


def edge_type_identification(g, kmax, ext_dic, return_all_list=False):
    global silent
    if not(silent): epbar = tqdm(total=g.number_of_edges(), desc='identifying')
    for u,v in g.edges():
        g[u][v]['type'] = None

    silks = []
    bonds = []
    Lbridges = []
    Gbridges = []
    int_threshold = {}

    ## phase 1: identify silk links
    edges = list(g.edges(data=True))
    nextphase = []
    for e in edges:
        u,v,w = e
        if (g.degree(u) == 1) | (g.degree(v) == 1):
            g[u][v]['type'] = 'Silk'
            silks.append(e)
            if not(silent): epbar.update(1)
        else:
            nextphase.append(e)
    #print len(silks)

    ## phase 2: identify bond and local bridges
    for i in range(kmax):
        l = str(i+1)
        lindex = 'w'+l+'a'
        Boname = 'Bond'+l
        Bdname = 'Local_Bridge'+l
        T_outk = ext_dic[l]

        edges = nextphase
        nextphase = []

        nextstep = []
        Rnextstep = []
        for e in edges:
            u,v,w = e
            Re = w[lindex]
            if Re>=T_outk:
                g[u][v]['type'] = Boname
                bonds.append((Boname, e))
                if not(silent): epbar.update(1)
            else:
                nextstep.append(e)
                Rnextstep.append(Re)

        if len(Rnextstep)==0:
            T_ink = 0
        else:
            T_ink = scipy.mean(Rnextstep) - scipy.std(Rnextstep)
            if T_ink<0:
                T_ink = 0.0
            for e in nextstep:
                u,v,w = e
                Re = w[lindex]
                if Re>T_ink:
                    g[u][v]['type'] = Bdname
                    Lbridges.append((Bdname, e))
                    if not(silent): epbar.update(1)
                else:
                    nextphase.append(e)
        int_threshold[l] = T_ink
        ## for kmax loop end here

    ## phase 3: identify global bridge
    edges = nextphase
    #nextphase = []
    for e in edges:
        u,v,w = e
        g[u][v]['type'] = 'Global_Bridge'
        Gbridges.append(e)
        if not(silent): epbar.update(1)

    if not(silent): print('done identify edge types')
    if return_all_list:
        return g, bonds, Lbridges, Gbridges, silks, int_threshold
    else:
        return g, int_threshold


def get_ego_graph(g, s, t, l):
    index     = EGO_NETWORK + str(l - 1)
    node_list = set()
    for ng in nx.neighbors(g, s):
        if ng != t: node_list = node_list | g.node[ng][index]
    return node_list - set([s])


def processing_link_property(iter_item):
    g, c, sp, s, t = iter_item
    #debugmsg('analyze the edge (' + s + ', ' + t + ')...')
    ## from s to t

    base_st_nodes = set([s, t])
    c.node[s][0]  = set() ## for removing previously accessed neighbor nodes (0~(l-1) layer neighbors)
    c.node[t][0]  = set() ## same as above, for the other end
    s0 = set()
    t0 = set()

    for i in range(sp):
        l = i + 1

        #c.node[s][l] = get_outgoing_ego_graph(c, s, t, l) - c.node[s][0] - base_st_nodes
        #c.node[t][l] = get_incoming_ego_graph(c, t, s, l) - c.node[t][0] - base_st_nodes
        c.node[s][l] = get_ego_graph(c, s, t, l) - s0 - base_st_nodes
        c.node[t][l] = get_ego_graph(c, t, s, l) - t0 - base_st_nodes

        common_nodes = (c.node[s][l] & c.node[t][l]) | (c.node[s][l] & c.node[t][l-1]) | (c.node[s][l-1] & c.node[t][l])

        index1 = 'w'+str(l)+'a' ## same as article, from inferior view
        #index2 = 'w'+str(l)+'b' ## from superior view
        g[s][t][index1] = None
        #g[s][t][index2] = None
        if len(common_nodes)==0:
            g[s][t][index1] = 0
            #g[s][t][index2] = 0
        else:
            part1_a = min(len(c.node[s][l]  ), len(c.node[t][l])  )
            part2_a = min(len(c.node[s][l]  ), len(c.node[t][l-1]))
            part3_a = min(len(c.node[s][l-1]), len(c.node[t][l])  )
            denominator_a = float(part1_a + part2_a + part3_a)
            #part1_b = max(len(c.node[s][l]  ), len(c.node[t][l])  )
            #part2_b = max(len(c.node[s][l]  ), len(c.node[t][l-1]))
            #part3_b = max(len(c.node[s][l-1]), len(c.node[t][l])  )
            #denominator_b = float(part1_b + part2_b + part3_b)
            g[s][t][index1] = float(len(common_nodes)) / denominator_a
            #g[s][t][index2] = float(len(common_nodes)) / denominator_b

        c.graph[GRAPH_KEY_COMMON_NODES_LIST + str(l)].append(g[s][t][index1])

        #c.node[s][0] |= c.node[s][l]
        #c.node[t][0] |= c.node[t][l]
        s0 |= c.node[s][l]
        t0 |= c.node[t][l]
    #pair = { g[s][t][ind]:v for ind,v in g[s][t].items() }
    #gnote = { GRAPH_KEY_COMMON_NODES_LIST + str(l): c.graph[GRAPH_KEY_COMMON_NODES_LIST + str(l)] for l in range(1, sp+1) }
    #compute_link_prop_bar.update(1)
    # FIXME: review this; has to return something or multiprocessing will not work, it copied everythings
    # TODO: add showing current completed, replace pbar
    # REVIEW: not using multiprocessing here already
    return #pair, gnote


def generate_ego_graph(g, sp):
    for r in range(sp):
        for n in g.nodes(data = False):
            if r == 0:
                g.node[n][EGO_NETWORK + str(r)] = set([n])
            else:
                g.node[n][EGO_NETWORK + str(r)] = set(g.node[n][EGO_NETWORK + str(r - 1)])
                for ng in nx.neighbors(g, n):
                    g.node[n][EGO_NETWORK + str(r)] = g.node[n][EGO_NETWORK + str(r)] | g.node[ng][EGO_NETWORK + str(r - 1)]


def compute_link_property(g, sp):
    global silent, threads_no
    ## g = Graph
    ## sp = k_max layer
    #print 'computing link property R'
    """
    核心演算法：計算目標網絡的每一條連結的兩端節點在不同半徑下除了該連結之外的交集程度，以供稍後判斷 BOND/sink/local bridge/global bridge
    時間複雜度：O(m x l)
    m = 目標網絡的連結數目，在連通圖情況下，m 通常大於節點數目 n，卻遠小於節點數目的平方（n x n）
    l = 目標網絡的最短路徑，通常 l 遠小於 log(n)，可當作常數項 C 看待

    參數：g 目標網絡，必須是連通圖，若不是，函數 compute_link_property 將擷取目標網絡 g 最大的 component 來運算分析
    　　　sp 整數，通常是目標網絡的平均最短路徑，表示分析的階層數
    """
    c = g.copy()

    for i in range(sp): c.graph[GRAPH_KEY_COMMON_NODES_LIST + str(i + 1)] = []

    generate_ego_graph(c, sp)
    #if not silent: print('calculating link common rate')
    """
    if False:
        #global compute_link_prop_bar
        #compute_link_prop_bar = tqdm(total=g.number_of_edges(), desc='computing arc prop.')
        pool = pathos.multiprocessing.ProcessingPool(nodes=threads_no)
        iterlist = [ (g, c, sp, s, t) for s, t in g.edges(data=False) ]
        link_res = pool.imap(processing_link_property, iterlist)
        link_res = list(link_res)
        # FIXME: review this; something need to capture and write back to the graph here
        for pair, gnote in link_res:
            for ind1, v in pair.items():
                g[s][t][ind1] = v
            for k, v in gnote.items():
                c.graph[k] = v
        #compute_link_prop_bar.close()
    else:
        for s, t in g.edges(data = False):
            processing_link_property((g, c, sp, s, t))
    """
    for s, t in g.edges(data = False):
        processing_link_property((g, c, sp, s, t))

    for i in range(sp):
        l = str(i + 1)
        g.graph[GRAPH_KEY_AVG_COMMON_NODES + l] = scipy.mean(c.graph[GRAPH_KEY_COMMON_NODES_LIST + l])
        g.graph[GRAPH_KEY_STD_COMMON_NODES + l] = scipy.std( c.graph[GRAPH_KEY_COMMON_NODES_LIST + l])
    return g


def random_once(g, kmax, Q=100):
    #rg = nx.DiGraph(nx.directed_configuration_model(list(d for n, d in g.in_degree()), list(d for n, d in g.out_degree()), create_using = nx.DiGraph()))
    rg = g.copy()
    if g.number_of_edges() > 2:
        nx.connected_double_edge_swap(rg, (Q * g.number_of_edges()))
    rg = compute_link_property(rg, kmax)
    #rgs.append(rg)
    meas = {str(i+1): rg.graph[GRAPH_KEY_AVG_COMMON_NODES + str(i+1)] for i in range(kmax)}
    stds = {str(i+1): rg.graph[GRAPH_KEY_STD_COMMON_NODES + str(i+1)] for i in range(kmax)}
    return meas, stds


def randomizing(iter_item):
    c, g, kmax, random_pre = iter_item
    global output_random, random_dir#, random_pre
    not_before = True
    if output_random:
        # check if this c is processed, if so skip random and load result
        output_path = os.path.join(random_dir, random_pre+str(c)+'.json')
        if os.path.isfile(output_path):
            not_before = False
            with open(output_path, 'r') as fread:
                tmp = json.load(fread)
            meas = tmp['mean']
            stds = tmp['std']
        else:
            meas, stds = random_once(g, kmax)
    else:
        meas, stds = random_once(g, kmax)
    if output_random and not_before:
        output_path = os.path.join(random_dir, random_pre+str(c)+'.json')
        tmp = {'mean': meas, 'std': stds}
        with open(output_path, 'w') as fp_hand:
            json.dump(tmp, fp_hand, indent=2, sort_keys=True)
    if c%10==0:
        global TIME
        tempTIME = time.time()
        print('completed randomizing',c, 'approx. used time: ', (tempTIME-TIME)/60., 'mins for about 10 iter.')
        TIME = tempTIME
    return meas, stds


def get_external_threshold(g, kmax, times, random_pre):
    #global kmax
    global silent, threads_no, TIME#, output_random, random_dir, random_pre

    if g.number_of_edges()>2:
        #rgs     = []
        rgmeans = { str(k+1):[] for k in range(kmax) }
        rgstds  = { str(k+1):[] for k in range(kmax) }
        #Q       = 10 #int(math.log10(g.order()) * math.log10(g.size()))
        #print Q
        # 產生供比較對應用的 times 個隨機網絡
        random_results = []
        #global pbar_pool
        #pbar_pool = tqdm(total=times)#, desc='randomizing with no. thread: '+str(threads))
        tt0 = time.time()
        TIME = time.time()
        if g.number_of_edges()>100:
            #print 'start randomizing with no. thread: '+str(threads)
            pool = pathos.multiprocessing.ProcessingPool(nodes=threads_no)
            iterlist = [(c, g, kmax, random_pre) for c in range(times)]
            random_results = pool.imap(randomizing, iterlist)
            random_results = list(random_results)
        else:
            for c in range(times):
                random_results.append((randomizing((c, g, kmax, random_pre))))
        #pbar_pool.close()
        tt1 = time.time()
        print('randomizing used time: {} minutes'.format((tt1-tt0)/60.))

        for i in range(kmax):
            rgmeans[str(i + 1)] = [ meas[str(i + 1)] for meas, stds in random_results ]
            rgstds[str(i + 1)]  = [ stds[str(i + 1)] for meas, stds in random_results ]
        ext_threshold = {}
        for i in range(kmax):
            l = str(i+1)
            ext = scipy.mean(rgmeans[l]) + scipy.mean(rgstds[l])
            if ext>1:
                ext_threshold[l] = 1.0
            else:
                ext_threshold[l] = ext
        if not(silent): print('done randomized and calculate external threshold')
        return ext_threshold
    else:
        if not(silent): print('graph has less than 2 edges')
        return None


def average_shortest_path_length(g):
    nlist = list(node for node in g.nodes())
    total = 0
    count = 0
    for index in tqdm(range(1000), 'Calculating average shortest path length'):
        s, t = random.sample(nlist, k = 2)
        if nx.has_path(g, source = s, target = t):
           total += nx.shortest_path_length(g, source = s, target = t)
           count += 1
        elif nx.has_path(g, source = t, target = s):
           total += nx.shortest_path_length(g, source = t, target = s)
           count += 1
    if count<=0:
        aspl = 1
    else:
        aspl = (total / float(count))
    return aspl


def bridge_or_bond(ginput, times=99, external=None, threads=4, kmax=None,
                run_silent=False, output_random_res=False,
                random_dir_path='temp', random_prefix='rand_'):

    global silent, threads_no, output_random, random_dir#, random_pre
    silent = run_silent
    threads_no = threads
    output_random = output_random_res
    random_dir = random_dir_path
    #random_pre = random_prefix

    # use only the largest weakly connected component
    g = sorted(nx.connected_component_subgraphs(ginput, copy=False), key=len, reverse=True)[0]

    #kmax = max(1, int(nx.average_shortest_path_length(g) / 2.0)) # 決定每個節點要外看幾層，決定強弱連結

    if not(silent): print('no_nodes:', g.number_of_nodes(), ', no_edges:', g.number_of_edges())

    if kmax is None:
        if not(silent): print('calculating kmax')
        avg_sp  = average_shortest_path_length(g)
        kmax  = max(1, int(math.floor(avg_sp / 2.0)))
    if not(silent): print('max layer is '+str(kmax))

    if not(silent): print('computing link property R')
    g = compute_link_property(g, kmax)

    if not(silent): print('calculating external threshold')
    if output_random:
        if not os.path.exists(random_dir):
            os.makedirs(random_dir)
    if external is None:
        random_pre = random_prefix
        ext_dic = get_external_threshold(g, kmax, times, random_pre)
    else:
        if not(silent): print('external threshold is provided, skipped randomization')
        ext_dic = external

    #g, bonds, Lbridges, Gbridges, silks, int_dic = edge_type_identification(g, kmax, ext_dic, return_all_list=True)
    if not(ext_dic is None):
        g, int_dic = edge_type_identification(g, kmax, ext_dic, return_all_list=False)
    else:
        int_dic = None
    #print len(bonds), len(Lbridges), len(Gbridges), len(silks)

    return g, ext_dic, int_dic




"""
drawing
"""

def community_sorting(g):
    etypes = { (u,v):d['type'] for u,v,d in g.edges(data=True) }
    etypes_rev = {'Silk': [], 'Global_Bridge': []}
    for e,t in etypes.items():
        if t not in etypes_rev: etypes_rev[t] = []
        etypes_rev[t].append(e)
    LBs = [ k for k in etypes_rev.keys() if k[:5]=='Local' ]
    Bs = [ k for k in etypes_rev.keys() if k[:4]=='Bond' ]
    kmax = 0
    for lb in LBs:
        k = int(lb.replace('Local_Bridge',''))
        if k>kmax: kmax=k
    for b in Bs:
        k = int(b.replace('Bond',''))
        if k>kmax: kmax=k

    gcopy = g.copy()
    nodes_sorted = []

    #seps= []
    isolates_ori = list(nx.isolates(gcopy))[::-1] # reverse it
    if len(isolates_ori)>0:
        nodes_sorted.extend(isolates_ori)
        gcopy.remove_nodes_from(isolates_ori)
        #seps.append(len(nodes_sorted))

    isolates_bysilk = []
    if len(etypes_rev['Silk'])>0:
        gcopy.remove_edges_from(etypes_rev['Silk'])
        isolates_bysilk = list(nx.isolates(gcopy))[::-1]
        nodes_sorted.extend(isolates_bysilk)
        gcopy.remove_nodes_from(isolates_bysilk)
        #seps.append(len(nodes_sorted))

    if len(etypes_rev['Global_Bridge'])>0:
        res, coms, hcoms = part_this(gcopy, set(gcopy.nodes()), etypes_rev, kmax+1, isglobal=True)
    else:
        res, coms, hcoms = part_this(gcopy, set(gcopy.nodes()), etypes_rev, kmax)
    nodes_sorted.extend(res)

    flat_iso_global = isolates_ori + isolates_bysilk
    #for n in isolates_ori: flat_iso_global[n] = [n]
    #for n in isolates_bysilk: flat_iso_global[n] = [n]

    flat_cdic = {}
    #for n in isolates_ori: flat_cdic[n] = [n]
    #for n in isolates_bysilk: flat_cdic[n] = [n]
    for n,c in zip(res, coms):
        if not(c in flat_cdic): flat_cdic[c] = []
        flat_cdic[c].append(n)

    for c,ns in flat_cdic.items():
        flat_cdic[c] = ns[::-1]
    #print(hcoms)
    h_tree_dic = get_hcom_list(hcoms, flat_cdic)
    h_tree_dic.update({ n:n for n in flat_iso_global })

    return nodes_sorted[::-1], h_tree_dic, flat_cdic, kmax


def get_hcom_list(hcom, flat_cdic):

    level_name = {}
    if not isinstance(hcom, list):
        comlist = flat_cdic[hcom]
        if len(comlist)>1:
            level_name[hcom] = comlist
        else: #elif len(comlist)<=1:
            level_name[hcom] = comlist[0]
        return level_name # return a dict with one key, one value

    else:
        lead = get_lead(hcom)
        #diclist = []
        hdic = {} # a number of key, value pair, pair number same as length of hcom
        for c in hcom:
            cdic = get_hcom_list(c, flat_cdic) # get a dict with one key, one value
            #diclist.append( cdic )
            hdic.update( cdic )
        level_name[lead] = hdic
        #level_name[diclist[0].keys()[0]] = diclist
        return level_name # return a dict with one key, one value


def get_lead(alist):
    if not(isinstance(alist, list)):
        return alist
    else:
        return get_lead(alist[0])


def part_this(h, sublist, etypes_rev, k, isglobal=False):
    res_list = []
    com_list = []
    set_com = []
    if k<1:
        res_list = sublist # return sublist as-is
        com_list = [list(sublist)[0]]*len(sublist)
        set_com = com_list[0]
        return res_list, com_list, set_com

    else:
        #print('k:', k, isglobal)
        if isglobal:
            this_key = 'Global_Bridge'
        else:
            this_key = 'Local_Bridge'+str(k)
        if not(this_key in etypes_rev):
            res_list = sublist
            com_list = [list(sublist)[0]]*len(sublist)
            set_com = com_list[0]
            return res_list, com_list, set_com

        this_LB = etypes_rev[this_key]
        hsub = h.subgraph(sublist).copy()#nx.DiGraph(h.subgraph(sublist))
        hsub.remove_edges_from(this_LB)
        isolates_byLB = list(nx.isolates(hsub))[::-1]

        res_list.extend(isolates_byLB)
        com_list.extend(isolates_byLB)
        set_com.extend(isolates_byLB)
        #print(set_com)

        hsub.remove_nodes_from(isolates_byLB)
        communities = sorted(nx.connected_components(hsub), key=len)
        if len(communities)==0: # not network, all isolated
            return res_list, com_list, set_com

        elif len(communities)==1: # left one connected component, extend it with isolated and return
            c = communities[0]
            res_list.extend(c)
            com_list.extend([list(c)[0]]*len(c))
            set_com.extend([list(c)[0]])
            #print(com_list, set_com)
            return res_list, com_list, set_com

        else: # more than one connected component, process each CC
            for c in communities:
                res, coms, com_ids = part_this(hsub, c, etypes_rev, k-1 )
                res_list.extend(res)
                com_list.extend(coms)
                if isglobal:
                    if isinstance(com_ids, list):
                        set_com.extend(com_ids)
                    else:
                        set_com.append(com_ids)
                else:
                    set_com.append(com_ids)
            return res_list, com_list, set_com


def tree_to_level_count(h_tree):
    layer = 0
    levels = {}
    iso = {}
    ismax = False
    this_tree = h_tree
    while not ismax:
        next_trees = {}
        #this_level = 0
        levels[layer] = 0
        iso[layer] = 0
        for k,v in this_tree.items():
            if isinstance(v, dict):
                #this_level+=len(v)
                levels[layer]+=1
                next_trees.update(v)
            elif isinstance(v, list):
                if len(v)==1:
                    #print'something wrong'
                    iso[layer]+=1
                elif len(v)>1:
                    levels[layer]+=1
            else:
                iso[layer]+=1
        #levels[layer] = this_level
        if len(next_trees)<=0:
            ismax = True
        else:
            layer+=1
            this_tree = next_trees
    level_count = { 'level_{}'.format(str(lvl)):count for lvl,count in levels.items() }
    iso_count = { 'level_{}'.format(str(lvl)):count for lvl,count in iso.items() }
    return level_count, iso_count


def get_color(d, color_dic=None):
    if color_dic is None:
        color_dic = {1:'blue', 2:'red', 3:'green', 4:'yellow', -1:'black'}
    if 'Bond' in d['type']:
        return color_dic[1], 1
    elif 'Local' in d['type']:
        return color_dic[2], 2
    elif 'Global' in d['type']:
        return color_dic[3], 3
    elif 'Silk' in d['type']:
        return color_dic[4], 4
    else:
        return color_dic[-1], -1


def draw_mat(dg, ax=None, cmap=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(6,6))

    orders, h_tree_dic, flat_dic, kmax = community_sorting(dg)
    mat = np.zeros((len(orders), len(orders)))
    for u,v,d in dg.edges(data=True):
        #print u,v,d
        u2 = orders.index(u)
        v2 = orders.index(v)
        c,ind = get_color(d)
        mat[u2][v2] = ind
    if cmap is None: cmap = ListedColormap(['k', 'blue', 'red', 'green', 'yellow'])
    ax.matshow(mat, cmap=cmap, vmin=0, vmax=4)
    #return mat

def draw_net(dg, pos, ax=None, color_dic=None):
    if ax is None:
        fig, ax = plt.subplots(figsize=(6,6))
    f = nx.draw_networkx_nodes(dg, pos=pos, ax=ax, size=3, node_color='lightgrey')
    #f = nx.draw_networkx_labels(dg, pos=pos, ax=ax, size=3)

    for u,v,d in dg.edges(data=True):
        uxy = pos[u]
        vxy = pos[v]
        col = get_color(d, color_dic=color_dic)[0]
        ax.annotate('', xy=vxy, xytext=uxy,
                     arrowprops=dict(arrowstyle='-', color=col, connectionstyle='arc3,rad=-0.15')
                    )
    ax.set_aspect('equal')


def draw_result(dg, pos=None, layout='circular', color_dic=None, cmap=None):
    if pos is None:
        if layout=='circular':
            pos = nx.circular_layout(dg)
        elif layout=='spring':
            pos = nx.spring_layout(dg)
        else:
            pos = nx.spring_layout(dg)

    fig, axs = plt.subplots(1, 2, figsize=(14, 7))
    ax1, ax2 = axs
    for ax in axs:
        ax.axis('off')
    draw_mat(dg, ax=ax1, cmap=cmap)
    draw_net(dg, pos, ax=ax2, color_dic=color_dic)
    plt.tight_layout()
    return fig, axs

def get_type(d):
    if 'Bond' in d['type']:
        return 'bond', 1
    elif 'Local' in d['type']:
        return 'local', 2
    elif 'Global' in d['type']:
        return 'global', 3
    elif 'Silk' in d['type']:
        return 'silk', 4
    else:
        return 'unknown', -1


def fingerprint(dg, ebunch=None):
    counts = { 'bond':0., 'local':0., 'global':0., 'silk':0., 'unknown':0. }
    total = 0.

    if ebunch is None: ebunch = list(dg.edges())
    dic = { (u,v): d for u,v,d in dg.edges(data=True) }
    for u,v in ebunch:
        d = dic[(u,v)]
        typ = get_type(d)[0]
        counts[typ]+=1.
        total+=1.
    if total>0:
        proportions = { k:v/total for k,v in counts.items() }
        ords = ['bond', 'local', 'global', 'silk']
        proportions2 = [ proportions[k] for k in ords ]
        return proportions2
    else:
        return ['-','-','-','-']


def main():
    # some test
    test_data_dir = '../data/net/'
    fs = sorted(os.listdir(test_data_dir))
    fs = [ f for f in fs if f[-4:]=='.net' ]
    #print(fs)
    #f = fs[0]
    for f in fs:
        print(f)
        fp = os.path.join(test_data_dir, f)
        g = nx.Graph(nx.read_pajek(fp))
        print(g.number_of_edges())
        g, ext_dic, int_dic = bridge_or_bond(g)
        print(ext_dic, int_dic)
        print(fingerprint(g))
    print('----------done----------')

if __name__ == '__main__':
    main()
