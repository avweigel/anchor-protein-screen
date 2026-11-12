"""Imaging literature per protein, from Europe PMC (full text of open-access papers plus abstracts of the rest).

A paper counts when the protein is its subject (gene symbol or HGNC full name in the title or abstract)
and the paper uses imaging (method terms anywhere in the text). Two searches per protein:
  any = any imaging method, light or electron microscopy  -> count + 5 most-cited papers
  em  = electron, correlative or cryo methods only         -> count
Short or word-like symbols (3 letters or fewer, or 4 letters with no digit) are ambiguous in text,
so for those only the full name is used.
Writes data/lit/europepmc_imaging.jsonl (one line per symbol); resumable. Run: python3 fetch_literature.py
"""
import csv, json, os, re, sys, time, threading, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor
OUT = 'data/lit/europepmc_imaging.jsonl'
IMG = '("immunofluorescence" OR "fluorescence microscopy" OR "confocal microscopy" OR "live-cell imaging" OR "live cell imaging" OR "super-resolution" OR "expansion microscopy" OR "electron microscopy" OR "cryo-electron tomography" OR "correlative light and electron" OR "FIB-SEM" OR "immunogold")'
EM = '("electron microscopy" OR "cryo-electron tomography" OR "cryo-ET" OR "correlative light and electron" OR "CLEM" OR "FIB-SEM" OR "volume electron microscopy" OR "immunogold" OR "subtomogram")'
API = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search'

def ambiguous(s): return len(s) <= 3 or (len(s) == 4 and s.isalpha())

def search(q, size):
    url = API + '?' + urllib.parse.urlencode({'query': q, 'format': 'json', 'pageSize': size, 'resultType': 'lite', 'sort': 'CITED desc'})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=60) as r: return json.load(r)
        except Exception as e:
            time.sleep(2 + 3 * attempt)
    raise RuntimeError('failed: ' + q[:80])

def main(symbols, names):
    done = set()
    if os.path.exists(OUT):
        for l in open(OUT): done.add(json.loads(l)['gene'])
    todo = [s for s in symbols if s not in done]
    print(len(symbols), 'symbols,', len(todo), 'to do', flush=True)
    lock = threading.Lock(); fh = open(OUT, 'a'); n = [0]
    def one(s):
        amb = ambiguous(s); name = names.get(s, '')
        q = lambda f: f'{f}:"{s}"' if not name else (f'{f}:"{name}"' if amb else f'{f}:"{s}" OR {f}:"{name}"')
        term = f'({q("TITLE")} OR {q("ABSTRACT")})'
        a = search(f'{term} AND {IMG}', 5); e = search(f'{term} AND {EM}', 1)
        top = [{'pmid': x.get('pmid', ''), 'pmcid': x.get('pmcid', ''), 'doi': x.get('doi', ''), 'title': x.get('title', '')[:200],
                'year': x.get('pubYear', ''), 'journal': x.get('journalTitle', ''), 'cited': x.get('citedByCount', 0)}
               for x in a.get('resultList', {}).get('result', [])]
        rec = {'gene': s, 'ambiguous': amb, 'name_used': name if amb else '', 'any': a.get('hitCount', 0), 'em': e.get('hitCount', 0), 'top': top}
        with lock:
            fh.write(json.dumps(rec) + '\n'); fh.flush(); n[0] += 1
            if n[0] % 500 == 0: print(n[0], 'done', time.strftime('%H:%M:%S'), flush=True)
        time.sleep(0.1)
    def safe(x):
        try: one(x)
        except Exception as e: print('skipped', x, e, flush=True)   # left out of the file, so a rerun picks it up
    with ThreadPoolExecutor(14) as ex: list(ex.map(safe, todo))
    print('finished', flush=True)

if __name__ == '__main__' and sys.argv[1:] != ['build']:
    syms = []
    for r in csv.DictReader(open('anchor_scores.tsv'), delimiter='\t'):
        g = r['gene_name']
        if g and not g.startswith('ENSG') and g not in syms: syms.append(g)
    names = {}
    hg = 'data/screens/raw/hgnc.tsv'
    if os.path.exists(hg):
        for r in csv.DictReader(open(hg), delimiter='\t'): names[r['symbol']] = r['name']
    main(syms, names)

def build():
    """Compact outputs: data/lit/imaging_counts.tsv for score.py, and lit/<letter>.json with the papers, loaded by the page on demand."""
    import collections
    rows = [json.loads(l) for l in open(OUT)]
    with open('data/lit/imaging_counts.tsv', 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t'); w.writerow(['gene', 'imaging_papers', 'em_papers', 'ambiguous_symbol', 'name_used'])
        for r in rows: w.writerow([r['gene'], r['any'], r['em'], 'yes' if r['ambiguous'] else '', r['name_used']])
    shards = collections.defaultdict(dict)
    for r in rows:
        if not r['top']: continue
        k = r['gene'][0].upper() if r['gene'][0].isalpha() else '0'
        shards[k][r['gene']] = {'n': r['any'], 'em': r['em'], 'name': r['name_used'],
                                'p': [[t['pmid'], t['pmcid'], t['year'], t['title'][:160], t['journal'][:60], t['cited']] for t in r['top']]}
    os.makedirs('lit', exist_ok=True)
    for k, d in shards.items(): json.dump(d, open(f'lit/{k}.json', 'w'), separators=(',', ':'))
    print('built', len(rows), 'genes,', len(shards), 'shards')

if __name__ == '__main__' and sys.argv[1:] == ['build']:
    build()
