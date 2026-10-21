import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

G = nx.DiGraph()

tfs = {
    "SP1":   {"tier": 4, "label": "SP1\n(4-tier)"},
    "HIF1A": {"tier": 3, "label": "HIF1A\n(3-tier)"},
    "RELA":  {"tier": 3, "label": "RELA\n(3-tier)"},
}
for tf, info in tfs.items():
    G.add_node(tf, layer=0, kind="TF", label=info["label"], size=900 + info["tier"]*250)
    G.add_edge(tf, "MIF", kind="activation")

G.add_node("MIF", layer=1, kind="target", label="MIF", size=2200)

mirnas = {
    "hsa-miR-654-5p":   {"score": 4},
    "hsa-miR-451a":     {"score": 3},
    "hsa-miR-1228-5p":  {"score": 2},
    "hsa-miR-23a-5p":   {"score": 2},
    "hsa-miR-363-5p":   {"score": 2},
    "hsa-miR-23b-5p":   {"score": 2},
    "hsa-miR-541-3p":   {"score": 2},
    "hsa-miR-92a-1-5p": {"score": 2},
    "hsa-miR-608":      {"score": 2},
    "hsa-miR-1537-3p":  {"score": 2},
}
for mir, info in mirnas.items():
    G.add_node(mir, layer=2, kind="miRNA", label=mir, size=500 + info["score"]*180)
    G.add_edge(mir, "MIF", kind="repression")

lncrnas = ["AC020922.3", "AL355987.4", "AL356488.2", "LMNTD2-AS1"]
for lnc in lncrnas:
    G.add_node(lnc, layer=3, kind="lncRNA", label=lnc, size=900)

ceRNA_stats = {"AC020922.3": "p=1.9e-3", "AL355987.4": "p=3.1e-4",
               "AL356488.2": "p=1.1e-4", "LMNTD2-AS1": "p=1.9e-3"}
for lnc, stat in ceRNA_stats.items():
    G.add_edge(lnc, "MIF", kind="ceRNA_stat", label=stat)

clip_confirmed_edges = [
    ("AL356488.2", "hsa-miR-654-5p"), ("AL355987.4", "hsa-miR-654-5p"),
    ("LMNTD2-AS1", "hsa-miR-654-5p"), ("AC020922.3", "hsa-miR-654-5p"),
    ("AL356488.2", "hsa-miR-541-3p"), ("AL355987.4", "hsa-miR-541-3p"),
    ("LMNTD2-AS1", "hsa-miR-541-3p"), ("AC020922.3", "hsa-miR-541-3p"),
]
seed_match_edges = [
    ("AL355987.4", "hsa-miR-23a-5p"), ("AL355987.4", "hsa-miR-363-5p"),
    ("AL355987.4", "hsa-miR-23b-5p"), ("AL355987.4", "hsa-miR-1537-3p"),
    ("AC020922.3", "hsa-miR-92a-1-5p"), ("AL355987.4", "hsa-miR-92a-1-5p"),
]
for lnc, mir in clip_confirmed_edges:
    G.add_edge(lnc, mir, kind="sponge_clip")
for lnc, mir in seed_match_edges:
    G.add_edge(lnc, mir, kind="sponge_seed")

layer_x = {0: 0.0, 1: 1.0, 2: 2.0, 3: 3.0}
pos = {}
for layer_id in [3, 2, 1, 0]:
    nodes_in_layer = [n for n, d in G.nodes(data=True) if d["layer"] == layer_id]
    n = len(nodes_in_layer)
    for i, node in enumerate(sorted(nodes_in_layer)):
        y = (i - (n - 1) / 2) * 1.3
        pos[node] = (layer_x[layer_id], y)

color_map = {"TF": "#4472C4", "target": "#C00000", "miRNA": "#70AD47", "lncRNA": "#FFC000"}
node_colors = [color_map[G.nodes[n]["kind"]] for n in G.nodes]
node_sizes = [G.nodes[n]["size"] for n in G.nodes]
labels = {n: G.nodes[n]["label"] for n in G.nodes}

plt.figure(figsize=(16, 11))
nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, edgecolors="black", linewidths=0.8)
nx.draw_networkx_labels(G, pos, labels=labels, font_size=8)

edge_styles = {
    "activation":   dict(edge_color="#4472C4", width=2.0, style="solid"),
    "repression":   dict(edge_color="#C00000", width=1.2, style="dashed"),
    "ceRNA_stat":   dict(edge_color="#B8860B", width=1.5, style="solid"),
    "sponge_clip":  dict(edge_color="#2E7D32", width=2.0, style="solid"),
    "sponge_seed":  dict(edge_color="#9E9E9E", width=0.8, style="dotted"),
}
edge_curvature = {
    "activation": "arc3,rad=0.1", "repression": "arc3,rad=0.05",
    "ceRNA_stat": "arc3,rad=0.45", "sponge_clip": "arc3,rad=0.1", "sponge_seed": "arc3,rad=0.1",
}
for kind, style in edge_styles.items():
    edges = [(u, v) for u, v, d in G.edges(data=True) if d["kind"] == kind]
    nx.draw_networkx_edges(G, pos, edgelist=edges, edge_color=style["edge_color"],
        width=style["width"], style=style["style"], arrows=True, arrowsize=12,
        connectionstyle=edge_curvature[kind])

legend_elements = [
    Line2D([0], [0], color="#4472C4", lw=2, label="TF -> MIF (activation, Phase 1)"),
    Line2D([0], [0], color="#C00000", lw=1.2, linestyle="--", label="miRNA -| MIF (repression)"),
    Line2D([0], [0], color="#B8860B", lw=1.5, label="lncRNA <-> MIF (ceRNA statistic)"),
    Line2D([0], [0], color="#2E7D32", lw=2, label="lncRNA -> miRNA (CLIP-seq confirmed)"),
    Line2D([0], [0], color="#9E9E9E", lw=0.8, linestyle=":", label="lncRNA -> miRNA (seed-match only)"),
]
plt.legend(handles=legend_elements, loc="lower center", bbox_to_anchor=(0.5, -0.08), ncol=2, fontsize=8, frameon=False)
plt.title("MIF Regulatory Network: TFs, miRNAs, and lncRNA Sponges", fontsize=14)
plt.axis("off")
plt.tight_layout()
plt.savefig("MIF_regulatory_network_v2.png", dpi=300, bbox_inches="tight")
print("Saved MIF_regulatory_network_v2.png")
