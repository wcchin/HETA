# HETA

# Intro

Hierarchical Edge Type Analysis. An algorithm for detecting edge types using common neighbor concept.

This is a rewritten of: https://github.com/canslab1/Identification-Algorithm/ , into convenient function that take NetworkX Graph objects, and to be compatible with Python3.  

## Install


	git clone https://github.com/wcchin/HETA.git
	cd HETA
	pip install .
	# or 
	# pip install -e .


## Usage

Please check the test.py script.

In simple:

	
	### import necessary packages and heta package
	import networkx as nx
	import heta

	### read file as nx.Graph() object
	fp = 'data/net/14p.net'
	g = nx.Graph(nx.read_pajek(fp))

	### run algorithm using default parameter
	g, ext_dic, int_dic = heta.bridge_or_bond(g)

	### get all edge types
	for u, v, d in g.edges(data=True):
		print(u, v, d['type'])

	### draw the results
	import matplotlib.pyplot as plt
	heta.draw_result(g, layout='spring')
	plt.show()

	### external and internal threshold
	print(ext_dic, int_dic)

	### fingerprint analysis
	### proportion of ['bond', 'local', 'global', 'silk']
	print(heta.fingerprint(g))




# Article reference

## Beyond Bond Links in Complex Networks: Local Bridges, Global Bridges and Silk Links

## Highlights

- An identification algorithm for determining hierarchical link types is described.
- Multi-hierarchy level link types are identified using the common neighbor concept.
- Two applications are demonstrated: fingerprint analysis and network partitioning.
- Fingerprint analysis is used to compare topological structures among networks.
- Hierarchical link structures are used to partition networks into communities.

## Abstract

Many network researchers use intuitive or basic definitions when discussing the importance of strong and [weak links](https://www.sciencedirect.com/topics/physics-and-astronomy/weak-link) and their roles. Others use an approach best described as “if not strong, then weak” to determine the strengths and weaknesses of individual links, thus deemphasizing hierarchical network structures that allow links to express different strength levels. Here we describe our proposal for a hierarchical edge type analysis (HETA) algorithm for determining link types at multiple network hierarchy levels based on the common neighbor concept plus statistical factors such as bond links, *k*th-layer local bridges, global bridges, and silk links—all generated during long-term network development and evolution processes. Two sets of networks were used to validate our proposed algorithm, one consisting of 16 networks employed in multiple past studies, and one consisting of two types of one-dimensional small-world networks expressing different random rewiring or shortcut addition probabilities. Two applications with potential for developmental contributions are demonstrated: a network fingerprint analysis framework, and a hierarchical network community [partition method](https://www.sciencedirect.com/topics/mathematics/partition-method).

Keywords: Network topology, Hierarchy of links, Common neighbor concept, Fingerprint analysis, Hierarchical community partition, Edge type analysis


article url: https://www.sciencedirect.com/science/article/pii/S0378437119306375

cite this as:  
Huang C. Y., Chin, W. C. B., Fu, Y. H., & Tsai, Y. S. (in press) Beyond bond links in complex networks: Local bridges, global bridges and silk links. Physica A: Statistical Mechanics and its Applications. DOI: https://doi.org/10.1016/j.physa.2019.04.263 .
