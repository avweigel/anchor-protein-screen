"""Download the loss-of-function screens and write small per-gene extracts to data/screens/.

Run once, needs network plus pandas, openpyxl, pyarrow and h5py:
    python3 fetch_screens.py
Large raw files (Perturb-seq, BioGRID ORCS) go to data/screens/raw/ and are not committed.
score.py reads only the extracts (data/screens/extract_*.tsv).

Sources (downloaded 11 Nov 2026):
- JUMP Cell Painting, CRISPR knockouts in U2OS (Chandrasekaran et al., Nat Methods 2025).
  Per-gene phenotypic activity (corrected p of mAP) from Zenodo 14642166, crispr.parquet.
- PERISCOPE, genome-wide optical pooled Cell Painting in HeLa (DMEM, HPLM) and A549
  (Ramezani et al., Nat Methods 2025). Hit lists from github.com/broadinstitute/2022_PERISCOPE.
- Funk et al. 2022, Cell. Optical pooled screen of 5,072 fitness genes in HeLa (DNA, gH2AX, tubulin, actin).
  Table S2A, copied in the PERISCOPE repository (Extended_Data_9/inputs).
- MitoCheck, Neumann et al. 2010, Nature. Genome-wide siRNA live imaging of H2B-GFP HeLa.
  IDR idr0013 screen A annotation, copied in github.com/WayScience/mitocheck_data.
- Replogle et al. 2022, Cell. Perturb-seq, CRISPRi, genome-wide in K562 and essential genes in RPE1.
  Per-perturbation summary (obs) of the bulk h5ad files on plus.figshare.com (article 20029387).
- BioGRID ORCS 2.0.18, all curated human CRISPR screens, plus screen metadata from each screen page.
- HGNC complete set, for current gene symbols.
"""
import os, re, json, glob, gzip, tarfile, subprocess, collections
import pandas as pd

S = 'data/screens/'
RAW = S + 'raw/'
os.makedirs(RAW, exist_ok=True)

def get(url, out):
    if not os.path.exists(out) or os.path.getsize(out) == 0:
        subprocess.run(['curl', '-sL', '--http1.1', '-m', '1800', url, '-o', out], check=True)
    return out

# ---------- HGNC (compact)
h = pd.read_csv(get('https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt', RAW + 'hgnc.tsv'),
                sep='\t', dtype=str, low_memory=False)
h[['symbol', 'entrez_id', 'ensembl_gene_id', 'uniprot_ids', 'prev_symbol', 'alias_symbol']].fillna('').to_csv(
    S + 'extract_hgnc.tsv', sep='\t', index=False)

# ---------- JUMP CRISPR, U2OS
p = get('https://zenodo.org/records/14642166/files/crispr.parquet?download=1', RAW + 'jump_crispr.parquet')
d = pd.read_parquet(p, columns=['Gene/Compound', 'JCP2022 ID', 'Phenotypic activity'])
d = d.drop_duplicates('Gene/Compound')
meta = pd.read_csv(S + 'jump_crispr_metadata.csv.gz', dtype=str)
d = d.merge(meta, left_on='JCP2022 ID', right_on='Metadata_JCP2022', how='left')
d = d[~d['Gene/Compound'].isin(['no-guide', 'non-targeting'])]
pd.DataFrame({'gene': d['Gene/Compound'], 'entrez': d['Metadata_NCBI_Gene_ID'].fillna(''),
              'corrected_p': d['Phenotypic activity']}).to_csv(S + 'extract_jump_crispr_u2os.tsv', sep='\t', index=False)

# ---------- PERISCOPE
R = 'https://raw.githubusercontent.com/broadinstitute/2022_PERISCOPE/main/'
pre = S + 'periscope__'
for folder, scr in [('2_HeLa_Screens_Summary', 'HeLa_DMEM'), ('2_HeLa_Screens_Summary', 'HeLa_HPLM'), ('4_A549_Screen_Summary', 'A549')]:
    for kind in ['whole_cell', 'compartment_specific']:
        f = f'{folder}/outputs/{scr}_plate_level_median_per_feat_sig_genes_1_FDR_{kind}_hits.csv'
        get(R + f, pre + f.replace('/', '_'))
get(R + 'Extended_Data_9/inputs/Funk_et_al_Supplement-Table_S2.xlsx', S + 'funk2022_Funk_et_al_Supplement-Table_S2.xlsx')
get('https://raw.githubusercontent.com/jump-cellpainting/datasets/main/metadata/crispr.csv.gz', S + 'jump_crispr_metadata.csv.gz')
# genes tested: first column of the genome-wide HeLa DMEM gene-level profiles (same Brunello-based library in A549)
tl = S + 'periscope_HeLa_DMEM_genes_tested.txt'
if not os.path.exists(tl):
    gz = get('https://raw.githubusercontent.com/broadinstitute/2022_PERISCOPE/main/2_HeLa_Screens_Summary/outputs/'
             'HeLa_DMEM_PCA_gene_level_profiles.csv.gz', RAW + 'periscope_HeLa_DMEM_gene_profiles.csv.gz')
    pd.read_csv(gz, usecols=[0]).iloc[:, 0].to_csv(tl, index=False)
tested = pd.read_csv(tl).iloc[:, 0]
rows = {g: {'gene': g} for g in tested if not str(g).lower().startswith('nontargeting')}
for scr, folder in [('HeLa_DMEM', '2_HeLa_Screens_Summary'), ('HeLa_HPLM', '2_HeLa_Screens_Summary'), ('A549', '4_A549_Screen_Summary')]:
    for kind in ['whole_cell', 'compartment_specific']:
        f = f'{pre}{folder}_outputs_{scr}_plate_level_median_per_feat_sig_genes_1_FDR_{kind}_hits.csv'
        t = pd.read_csv(f)
        chans = [c for c in t.columns if c not in ('Gene', 'Sum')]
        for _, r in t.iterrows():
            o = rows.setdefault(r['Gene'], {'gene': r['Gene']})
            o[scr] = 1
            top = [c for c in chans if r[c] > 0]
            o.setdefault('channels', set()).update(top)
out = pd.DataFrame(rows.values())
for c in ['HeLa_DMEM', 'HeLa_HPLM', 'A549']: out[c] = out.get(c, 0).fillna(0).astype(int)
out['channels'] = out['channels'].apply(lambda s: ','.join(sorted(s)) if isinstance(s, set) else '')
out.to_csv(S + 'extract_periscope.tsv', sep='\t', index=False)

# ---------- Funk 2022
d = pd.read_excel(S + 'funk2022_Funk_et_al_Supplement-Table_S2.xlsx', sheet_name=0, header=1)
fdr = [c for c in d.columns if str(c).endswith('_FDR')]
d['n_sig_features'] = (d[fdr] < 0.05).sum(axis=1)
nt = d['Gene symbol'].astype(str).str.startswith('nontargeting')
m, s = d.loc[nt, 'Mitotic index (%)'].mean(), d.loc[nt, 'Mitotic index (%)'].std()
d['mitotic_index_high'] = (d['Mitotic index (%)'] > m + 3 * s).astype(int)
# hit: 3 or more of 32 image features at FDR < 0.05 (no non-targeting control reaches 3), or mitotic index above controls + 3 SD
d['hit'] = ((d['n_sig_features'] >= 3) | (d['mitotic_index_high'] == 1)).astype(int)
d = d[~nt]
pd.DataFrame({'gene': d['Gene symbol'], 'entrez': d['NCBI gene ID'], 'n_sig_features': d['n_sig_features'],
              'mitotic_index': d['Mitotic index (%)'], 'mitotic_index_high': d['mitotic_index_high'], 'hit': d['hit']}
             ).to_csv(S + 'extract_funk2022_hela.tsv', sep='\t', index=False)

# ---------- MitoCheck
d = pd.read_csv(get('https://raw.githubusercontent.com/WayScience/mitocheck_data/main/mitocheck_metadata/idr0013-screenA-annotation.csv.gz',
                    RAW + 'mitocheck_idr0013-screenA-annotation.csv.gz'), low_memory=False)
d = d[d['Control Type'].isna() & d['Gene Symbol'].notna()]
names = [c for c in d.columns if c.startswith('Phenotype') and 'Term Name' in c]
rows = []
for g, x in d.groupby('Gene Symbol'):
    ph = x[x['Has Phenotype'] == 'yes']
    terms = collections.Counter(t for t in ph[names].stack() if isinstance(t, str))
    rows.append({'gene': g, 'ensembl': x['Gene Identifier'].iloc[0],
                 'reproducible_phenotype': int(len(ph) > 0),
                 'validated_mitotic_hit': int((x['Validated Mitotic Hit At Gene Level'] == 'yes').any()),
                 'phenotypes': '; '.join(t.replace(' phenotype', '') for t, _ in terms.most_common(3))})
pd.DataFrame(rows).to_csv(S + 'extract_mitocheck_hela.tsv', sep='\t', index=False)

# ---------- Replogle 2022
import h5py
def obs(fn):
    o = h5py.File(fn)['obs']; out = {}
    for k in o.keys():
        v = o[k]
        if isinstance(v, h5py.Group): out[k] = v['categories'][:].astype(str)[v['codes'][:]]
        else:
            a = v[:]; out[k] = a.astype(str) if a.dtype.kind in 'SO' else a
    return pd.DataFrame(out)
for name, fid in [('k562_gwps', '35773217'), ('rpe1', '35775512')]:
    d = obs(get(f'https://ndownloader.figshare.com/files/{fid}', RAW + f'replogle_{name}_bulk.h5ad'))
    parts = d['gene_transcript'].str.split('_')
    d['gene'] = parts.str[1]; d['ensembl'] = parts.str[-1]
    d = d[~d['gene'].str.contains('non-targeting')]
    # strong transcriptional phenotype: energy-distance test p < 0.001 and 50 or more DE genes (Anderson-Darling)
    d['hit'] = ((d['energy_test_p_value'] < 0.001) & (d['anderson_darling_counts'] >= 50)).astype(int)
    d = d.sort_values('hit', ascending=False).drop_duplicates('gene')
    d[['gene', 'ensembl', 'energy_test_p_value', 'anderson_darling_counts', 'fold_expr', 'hit']].to_csv(
        S + f'extract_replogle_{name}.tsv', sep='\t', index=False)

# ---------- BioGRID ORCS
tgz = get('https://downloads.thebiogrid.org/Download/BioGRID-ORCS/Release-Archive/BIOGRID-ORCS-2.0.18/'
          'BIOGRID-ORCS-ALL-homo_sapiens-2.0.18.screens.tar.gz', RAW + 'orcs_hs.tar.gz')
os.makedirs(RAW + 'orcs_meta', exist_ok=True)
tf = tarfile.open(tgz)
members = [m for m in tf.getmembers() if m.name.endswith('.screen.tab.txt')]
def kw(kws, key):
    return '; '.join(k.split(' > ', 1)[1] for k in kws if k.startswith(key + ' > '))
screens = []
for m in members:
    sid = re.search(r'SCREEN_(\d+)-', m.name).group(1)
    mf = RAW + f'orcs_meta/{sid}.json'
    if not os.path.exists(mf) or os.path.getsize(mf) == 0:
        html = subprocess.run(['curl', '-s', '-m', '30', f'https://orcs.thebiogrid.org/Screen/{sid}'], capture_output=True, text=True).stdout
        j = re.search(r'(\{"@context".*?"citation":"[^"]*"\})', html, re.S)
        json.dump({'id': sid, 'kw': json.loads(j.group(1)).get('keywords') if j else None,
                   'name': json.loads(j.group(1)).get('name') if j else None}, open(mf, 'w'))
    meta = json.load(open(mf)); kws = meta.get('kw') or []
    screens.append({'screen_id': sid, 'name': meta.get('name') or '', 'cell_line': kw(kws, 'CELL LINE'),
                    'phenotype': kw(kws, 'PHENOTYPE'), 'library_type': kw(kws, 'LIBRARY TYPE'),
                    'methodology': kw(kws, 'METHODOLOGY'), 'screen_type': kw(kws, 'SCREEN TYPE'), 'member': m.name})
sc = pd.DataFrame(screens)
sc.drop(columns='member').to_csv(S + 'extract_orcs_screens.tsv', sep='\t', index=False)
lof = sc['methodology'].str.contains('KNOCKOUT|INHIBITION', case=False)   # loss of function only, not activation
groups = {'u2os': lof & sc['cell_line'].str.fullmatch('U-2 OS|U2OS|U-2OS', case=False),
          'hela': lof & sc['cell_line'].str.contains(r'^HELA', case=False, regex=True),
          'all': lof}
cnt = {c: collections.Counter() for c in ['u2os_hit', 'u2os_tested', 'hela_hit', 'hela_tested', 'all_hit', 'all_tested']}
sym = {}
for i, r in sc.iterrows():
    gs = [g for g, mask in groups.items() if mask[i]]
    if not gs: continue
    t = pd.read_csv(tf.extractfile(r['member']), sep='\t', dtype=str, usecols=['IDENTIFIER_ID', 'IDENTIFIER_TYPE', 'OFFICIAL_SYMBOL', 'HIT'])
    t = t[t['IDENTIFIER_TYPE'] == 'ENTREZ_GENE'].drop_duplicates('IDENTIFIER_ID')
    sym.update(zip(t['IDENTIFIER_ID'], t['OFFICIAL_SYMBOL']))
    ids = t['IDENTIFIER_ID'].tolist(); hits = t.loc[t['HIT'] == 'YES', 'IDENTIFIER_ID'].tolist()
    for g in gs:
        cnt[g + '_tested'].update(ids); cnt[g + '_hit'].update(hits)
rows = [{'entrez': e, 'gene': sym[e], **{c: cnt[c].get(e, 0) for c in cnt}} for e in sym]
pd.DataFrame(rows).to_csv(S + 'extract_orcs_genes.tsv', sep='\t', index=False)
print('ORCS loss-of-function screens:', {g: int(m.sum()) for g, m in groups.items()})
print('done')
