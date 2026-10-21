import graphviz

# Create Digraph with document canvas styling
dot = graphviz.Digraph('MIF_Pipeline', comment='MIF Regulatory Discovery Pipeline')

# Global Graph Attributes
dot.attr(rankdir='TB', nodesep='0.4', ranksep='0.5', fontname='Helvetica')

# Global Node & Edge Styling
dot.attr('node', fontname='Helvetica', fontsize='10', shape='rect', 
         style='filled,rounded', fillcolor='#FFFFFF', color='#4A5568', penwidth='1.5')
dot.attr('edge', fontname='Helvetica', fontsize='8', color='#718096', 
         arrowhead='normal', penwidth='1.2')

# Top Root Node
dot.node('root', 'MIF Gene (Human)', fillcolor='#EDF2F7', style='filled,bold', shape='oval', fontsize='11')

# --- Stage I Cluster ---
with dot.subgraph(name='cluster_stage1') as c1:
    c1.attr(label='Stage I: Conservation-Guided Regulatory Discovery', 
            fontname='Helvetica-Bold', fontsize='11', style='dashed', 
            color='#A0AEC0', bgcolor='#F7FAFC', margin='16')
    
    c1.node('s1_promoter', 'Human + Chimpanzee + Mouse\nPromoter Conservation Analysis')
    c1.node('s1_jaspar', 'JASPAR TF Binding Prediction')
    c1.node('s1_trrust', 'TRRUST Regulatory Evidence')
    c1.node('s1_encode', 'ENCODE ChIP-seq Validation')
    c1.node('s1_cand', 'High-Confidence Candidates\n(SP1, RELA, HIF1A)', style='filled', fillcolor='#E2E8F0')
    c1.node('s1_clinical', 'Clinical Validation (6 TCGA Cohorts)\n• Expression  • Correlation  • Survival\n• Immune Context  • CNA / Mutation')
    c1.node('s1_prio', 'RELA Prioritized (KIRC)', style='filled,bold', fillcolor='#FEB2B2', color='#C53030')

    c1.edge('s1_promoter', 's1_jaspar')
    c1.edge('s1_jaspar', 's1_trrust')
    c1.edge('s1_trrust', 's1_encode')
    c1.edge('s1_encode', 's1_cand')
    c1.edge('s1_cand', 's1_clinical')
    c1.edge('s1_clinical', 's1_prio')

# --- Stage II Cluster ---
with dot.subgraph(name='cluster_stage2') as c2:
    c2.attr(label='Stage II: Expanded Human-Promoter Screen', 
            fontname='Helvetica-Bold', fontsize='11', style='dashed', 
            color='#A0AEC0', bgcolor='#F7FAFC', margin='16')
    
    c2.node('s2_promoter', 'Human MIF Promoter')
    c2.node('s2_jaspar', 'JASPAR Ranking (Top TF Candidates)')
    c2.node('s2_add', '11 Additional TFs')
    c2.node('s2_tcga', 'TCGA-TARGET-GTEx Validation\n• Tumor vs Normal\n• Survival\n• lncRNA Validation')
    c2.node('s2_prio', 'TFAP2A Prioritized', style='filled,bold', fillcolor='#FEEBC8', color='#C05621')

    c2.edge('s2_promoter', 's2_jaspar')
    c2.edge('s2_jaspar', 's2_add')
    c2.edge('s2_add', 's2_tcga')
    c2.edge('s2_tcga', 's2_prio')

# --- Integrated Modeling Section ---
dot.node('model', 'Integrated Modeling\n• LASSO-Cox\n• Stability Selection\n• Random Survival Forest')
dot.node('final', 'Integrated MIF Regulatory Map', style='filled,bold', 
         fillcolor='#C6F6D5', color='#22543D', shape='doubleoctagon', fontsize='11')

# --- Connecting Pipeline Sequence ---
dot.edge('root', 's1_promoter')
dot.edge('s1_prio', 's2_promoter', label=' Next Phase', fontcolor='#4A5568')
dot.edge('s2_prio', 'model')
dot.edge('model', 'final')

# Render to PDF and SVG (Vector outputs for manuscript submission)
dot.render('MIF_Pipeline', format='pdf', cleanup=True)
dot.render('MIF_Pipeline', format='svg', cleanup=True)

print("Diagram generated successfully as 'MIF_Pipeline.pdf' and 'MIF_Pipeline.svg'!")

