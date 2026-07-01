from graphviz import Digraph
def de_bruijn_ize(st,k):
    edges=[]
    nodes = set()
    for i in range(len(st)-k+1):
        edges.append((st[i:i+k-1],st[i+1:i+k]))
        nodes.add(st[i:i+k-1])
        nodes.add(st[i+1:i+k])
    return nodes,edges



def visualize(st, k):
    nodes, edges = de_bruijn_ize(st, k)

    dot = Digraph()

    for node in nodes:
        dot.node(node, label=node)

    for src, dst in edges:
        dot.edge(src, dst)

    return dot


graph = visualize('ACTGAAGGGGTCGA', 3)

graph.render('debruijn_graph', view=True, format='png')