import csv,re,json,collections
D='data/'
hpa_look={r['hpa_location']:r for r in csv.DictReader(open('hpa_em_lookup.csv'))}
UNI=[('Nuclear pore complex','nuclear pore'),('nuclear pore','nuclear pore'),('Nucleus envelope','nuclear envelope'),('Nucleus membrane','nuclear envelope'),('Nucleus inner membrane','nuclear envelope'),('Nucleus lamina','nuclear lamina'),('nucleolus','nucleolus'),
 ('Mitochondrion outer membrane','mitochondrial outer membrane'),('Mitochondrion inner membrane','mitochondrial inner membrane'),('Mitochondrion matrix','mitochondrial matrix'),('Mitochondrion membrane','mitochondrial membrane'),('Mitochondrion intermembrane','mitochondrion'),('Mitochondrion','mitochondrion'),
 ('Endoplasmic reticulum membrane','ER membrane'),('Endoplasmic reticulum lumen','ER lumen'),('Endoplasmic reticulum','ER'),('Microsome membrane','ER membrane'),
 ('trans-Golgi network','Golgi (TGN)'),('Golgi apparatus membrane','Golgi membrane'),('Golgi apparatus','Golgi'),
 ('Lysosome membrane','lysosome membrane'),('Lysosome lumen','lysosome lumen'),('Lysosome','lysosome'),('Peroxisome membrane','peroxisome membrane'),('Peroxisome','peroxisome'),
 ('Early endosome','endosome'),('Late endosome','endosome'),('Recycling endosome','endosome'),('Endosome','endosome'),('Lipid droplet','lipid droplet'),('Autophagosome','autophagosome'),
 ('centriole','centriole'),('centrosome','centrosome'),('cilium basal body','basal body'),('kinetochore','kinetochore'),('centromere','centromere'),
 ('Apical cell membrane','plasma membrane'),('Basolateral cell membrane','plasma membrane'),('Cell membrane','plasma membrane'),
 ('tight junction','tight junction'),('adherens junction','adherens junction'),('focal adhesion','focal adhesion'),('Cell junction','cell junction'),('desmosome','desmosome'),('gap junction','gap junction'),
 ('Midbody','midbody'),('spindle','spindle'),('microtubule','microtubule'),('stress fiber','actin bundle'),('Nucleus speckle','nuclear speckle'),('Cajal body','Cajal body'),('PML body','PML body')]
UNI_PARTIAL=[('Chromosome','chromosome (broad)'),('cytoskeleton','cytoskeleton (broad)'),('Cytoplasmic vesicle','vesicle (unspecified)'),('Cell projection','cell projection')]
def uni_hits(s,table):
    out=[]
    for term,struct in table:
        if term.lower() in s.lower() and struct not in out: out.append(struct)
    return out
# --- load
uni=list(csv.DictReader(open(D+'uniprot_human_reviewed.tsv'),delimiter='\t'))
sub={r['Gene']:r for r in csv.DictReader(open(D+'subcellular_location.tsv'),delimiter='\t')}
full_by_ens={r['Ensembl']:r for r in csv.DictReader(open(D+'proteinatlas.tsv'),delimiter='\t')}
full_by_acc={}
for r in full_by_ens.values():
    for a in r['Uniprot'].split(','):
        a=a.strip()
        if a: full_by_acc.setdefault(a,r)
rna={r['Gene']:float(r['nTPM']) for r in csv.DictReader(open(D+'hpa_rna_U2OS.tsv'),delimiter='\t')}
def pc(s):
    s=s.strip().lstrip('><')
    try: return float(s)
    except: return None
beck_raw={r['uniprot']:r['copies_per_cell'] for r in csv.DictReader(open(D+'beck2011_U2OS_copies_per_cell.tsv'),delimiter='\t')}
mit={r['uniprot']:float(r['copy_number_control_over_synchronized']) for r in csv.DictReader(open(D+'beck2011_U2OS_mitotic_copy_ratio.tsv'),delimiter='\t')}
dep={r['gene']:r for r in csv.DictReader(open(D+'depmap_26Q3_gene_effect_U2OS_and_summary.tsv'),delimiter='\t')}
oc=collections.defaultdict(list)
for l in json.load(open(D+'opencell_lines.json')): oc[l['metadata']['ensg_id']].append(l)
allen=collections.defaultdict(list)
for r in csv.DictReader(open(D+'allen_cell_collection_lines.tsv'),delimiter='\t'): allen[r['gene']].append(r)
# --- universe: uniprot reviewed, plus HPA-only genes
seen_ens=set(); rows=[]
full_by_name={r['Gene']:r for r in csv.DictReader(open(D+'proteinatlas.tsv'),delimiter='\t')}
def ens_of(u):
    # UniProt 'Ensembl' column holds transcript IDs (ENST), so map accession -> HPA gene via proteinatlas.tsv, then gene name
    out=[]
    f=full_by_acc.get(u['Entry']) or full_by_name.get(u['Gene Names (primary)'])
    if f: out.append(f['Ensembl'])
    for e in u['Ensembl'].split(';'):
        e=e.strip().split(' ')[0].split('.')[0]
        if e.startswith('ENSG') and e not in out: out.append(e)
    return out
for u in uni:
    rows.append(('uniprot',u,ens_of(u)))
    seen_ens.update(ens_of(u))
for g,r in sub.items():
    if g not in seen_ens: rows.append(('hpa_only',None,[g]))
FAM=[('nuclear pore','nuclear envelope'),('nuclear envelope','nuclear envelope'),('nuclear lamina','nuclear envelope'),('nuclear membrane','nuclear envelope'),
 ('nucleol','nucleolus'),('speckle','nuclear bodies'),('cajal','nuclear bodies'),('pml','nuclear bodies'),('interchromatin','nuclear bodies'),('nuclear bod','nuclear bodies'),
 ('kinetochore','centromere/kinetochore'),('centromere','centromere/kinetochore'),('chromosome','chromosome'),
 ('mitochond','mitochondria'),('er ','ER'),('er','ER'),('endoplasmic','ER'),('golgi','Golgi'),('lysosom','endolysosomal'),('endosom','endolysosomal'),('vesicle','endolysosomal'),('autophag','endolysosomal'),
 ('peroxisom','peroxisome'),('lipid droplet','lipid droplet'),('centriol','centrosome'),('centrosom','centrosome'),('basal body','centrosome'),('cilium','centrosome'),
 ('plasma membrane','plasma membrane'),('junction','junctions/adhesions'),('focal adhesion','junctions/adhesions'),('desmosome','junctions/adhesions'),
 ('microtubule','cytoskeleton'),('actin','cytoskeleton'),('intermediate filament','cytoskeleton'),('spindle','cytoskeleton'),('midbody','cytoskeleton'),('cytoskeleton','cytoskeleton'),('cell projection','cytoskeleton'),('cytokinetic','cytoskeleton')]
def fam(s):
    t=s.lower().strip()
    for k,f in FAM:
        if k=='er': 
            if t=='er' or t.startswith('er ') : return f
        elif k in t: return f
    return 'other'
def key(s): return fam(s)
out=[]
for src,u,enss in rows:
    h=None; g=''
    for e in enss:
        if e in sub: h=sub[e]; g=e; break
    if not g and enss: g=enss[0]
    f=full_by_ens.get(g) or (full_by_acc.get(u['Entry']) if u else None) or (full_by_name.get(u['Gene Names (primary)']) if u else None) or {}
    name=(u['Gene Names (primary)'] if u else h['Gene name']) or (h['Gene name'] if h else '')
    acc=u['Entry'] if u else f.get('Uniprot','').split(',')[0].strip()
    main=[x for x in h['Main location'].split(';') if x] if h else []
    add=[x for x in h['Additional location'].split(';') if x] if h else []
    hpa_em=[hpa_look.get(l,{}).get('em_correlate','?') for l in main]
    hpa_id=2 if 'yes' in hpa_em else (1 if 'partial' in hpa_em else (0 if main else None))
    hpa_struct=[hpa_look[l]['structure'] for l in main if l in hpa_look and hpa_look[l]['em_correlate'] in('yes','partial')]
    uloc=u['Subcellular location [CC]'] if u else ''
    uclear=uni_hits(uloc,UNI); upart=uni_hits(uloc,UNI_PARTIAL)
    uni_id=2 if uclear else (1 if upart else (0 if uloc else None))
    identity=max([x for x in (hpa_id,uni_id) if x is not None], default=0)
    structures=sorted(set(hpa_struct+uclear+upart))
    families=sorted({fam(x) for x in structures if fam(x)!='other'})
    # agreement between HPA and UniProt
    if hpa_struct and uclear:
        agree=2 if {key(x) for x in hpa_struct}&{key(x) for x in uclear} else 0
    elif hpa_struct or uclear: agree=1
    else: agree=1
    rel=h['Reliability'] if h else ''
    ab=2 if rel in('Enhanced','Supported') else (1 if rel=='Approved' else 0)
    ntpm=rna.get(g); expr=2 if (ntpm is not None and ntpm>=10) else (1 if (ntpm is not None and ntpm>=1) else 0)
    cp=beck_raw.get(acc,''); quant=2 if cp else 1
    tag=[]
    for e in enss:
        for l in oc.get(e,[]): tag.append('OpenCell '+(l['metadata'].get('target_terminus') or '?')+'-term')
    for a in allen.get(name,[]): tag.append('Allen '+a['tag_location'])
    tag=sorted(set(tag)); taggable=2 if tag else 1
    d=dep.get(name,{}); ge=d.get('U2OS_gene_effect','')
    ess=('essential' if float(ge)<-1 else 'fitness' if float(ge)<-0.5 else 'neutral') if ge else ''
    cc=h['Cell cycle dependency'] if h else ''; mr=mit.get(acc)
    total=identity+agree+ab+expr+quant+taggable
    out.append(dict(gene_name=name,uniprot=acc,ensembl=g,source=('UniProt+HPA' if (u and h) else ('UniProt, HPA has antibody but no IF location' if (u and f.get('Antibody RRID')) else ('UniProt, in HPA but no antibody' if (u and f) else ('UniProt only, not in HPA' if u else 'HPA only')))),hpa_antibodies=len([a for a in f.get('Antibody RRID','').split(',') if a.strip()]) if f else 0,
        compartment_family=';'.join(families),em_structures=';'.join(structures),score_identity=identity,score_sources_agree=agree,score_antibody=ab,score_expressed_u2os=expr,score_copy_number=quant,score_tag_evidence=taggable,
        score_total=total,hpa_main=';'.join(main),hpa_additional=';'.join(add),hpa_reliability=rel,uniprot_location=uloc[:300],transmembrane='yes' if (u and u['Transmembrane']) else '',
        tag_lines=';'.join(tag),copies_per_cell_U2OS=cp,mitotic_ratio=round(mr,2) if mr else '',u2os_ntpm=ntpm if ntpm is not None else '',
        hpa_cell_cycle=cc,depmap_u2os_effect=ge,essentiality=ess,antibody_rrid=f.get('Antibody RRID',''),protein_class=f.get('Protein class','')))
out.sort(key=lambda o:(-o['score_total'],o['gene_name']))
with open('anchor_scores.tsv','w',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=list(out[0].keys()),delimiter='\t'); w.writeheader(); w.writerows(out)
C=collections.Counter(o['source'] for o in out); print(len(out),'proteins', dict(C))
print('score_total distribution:'); 
for s,n in sorted(collections.Counter(o['score_total'] for o in out).items(),reverse=True): print(f'  {s:2d}: {n}')
print('identity 2:',sum(o['score_identity']==2 for o in out),' 1:',sum(o['score_identity']==1 for o in out))
rows={o['gene_name']:o for o in out}
for n in ['NUP107','NUP98','LAMP1','CENPA','TOMM20','LMNB1','SEC61B','PEX14','ZYX','GOLGA2','LMNA','NUP96']:
    o=rows.get(n); print(n,'missing' if not o else f"total={o['score_total']} id={o['score_identity']} agree={o['score_sources_agree']} ab={o['score_antibody']} expr={o['score_expressed_u2os']} cp={o['score_copy_number']} tag={o['score_tag_evidence']} src={o['source']} struct={o['em_structures'][:45]}")
