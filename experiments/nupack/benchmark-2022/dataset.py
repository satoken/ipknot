#!/usr/bin/env python3
"""Reproduce the published 2022 single-sequence benchmark, without BPP caching.

Prepare official references; resume complete paired measurements; retain every
prediction and child resource measurement. No accuracy-dependent selection.
"""
import argparse, collections, concurrent.futures, csv, gzip, hashlib, json, os
from pathlib import Path
import multiprocessing, platform, random, statistics, subprocess, time, zipfile

MD5='01f491361b31b4520f5fd527d4ad5cbf'
SHA256='fd360f9d98820fb7f41b7cd3ea333416fc4e0ce12a1816171ed0fa12303d41a2'
MODES=['ilp','dd50','dd100','dd500']
COMMON=['-e','lpc','-t','auto,auto','-r','1','--beam-size','100']

def digest(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write_json(path, value):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value,separators=(',',':'),ensure_ascii=False)+'\n')
    temp.replace(path)
def load_index(root):
    with gzip.open(root/'references.json.gz','rt') as f: return json.load(f)
def read_bpseq(text):
    rows=[]
    for line in text.splitlines():
        if not line.strip() or line.startswith('#'): continue
        s=line.split()
        if len(s)!=3 or not s[0].isdigit() or not s[2].isdigit():
            raise ValueError('invalid BPSEQ row: '+line[:100])
        if int(s[0])!=len(rows)+1 or len(s[1])!=1: raise ValueError('invalid index/base')
        rows.append((s[1].upper().replace('T','U'),int(s[2])))
    if not rows: raise ValueError('empty BPSEQ')
    for i,(base,j) in enumerate(rows,1):
        if not 0<=j<=len(rows) or i==j or (j and rows[j-1][1]!=i):
            raise ValueError('asymmetric/out-of-range BPSEQ')
    return ''.join(b for b,j in rows),[(i,j) for i,(b,j) in enumerate(rows,1) if i<j]

def crossing_pairs(pairs):
    """Mark both participants of each i<k<j<l crossing in O(n log n)."""
    pairs=sorted(map(tuple,pairs))
    if not pairs: return set()
    n=max(j for i,j in pairs); tree=[0]*(n+1); result=set()
    def add(k):
        while k<=n: tree[k]+=1; k+=k&-k
    def count(k):
        v=0
        while k: v+=tree[k]; k-=k&-k
        return v
    for i,j in pairs:
        if count(j-1)-count(i): result.add((i,j))
        add(j)
    tree=[0]*(n+1)
    for i,j in reversed(pairs):
        k=j-1; v=0
        while k: v=max(v,tree[k]); k-=k&-k
        if v>j: result.add((i,j))
        k=i
        while k<=n: tree[k]=max(tree[k],j); k+=k&-k
    return result

def metrics(reference,prediction):
    ref=set(map(tuple,reference)); pred=set(map(tuple,prediction))
    tp=len(ref&pred); fp=len(pred-ref); fn=len(ref-pred)
    return {'tp':tp,'fp':fp,'fn':fn,'ppv':tp/len(pred) if pred else 0.0,
            'sen':tp/len(ref) if ref else 0.0,
            'f1':2*tp/(len(ref)+len(pred)) if ref or pred else 1.0}

def prepare(args):
    root=args.root; root.mkdir(parents=True,exist_ok=True)
    blob=args.archive.read_bytes()
    if hashlib.md5(blob).hexdigest()!=MD5 or hashlib.sha256(blob).hexdigest()!=SHA256:
        raise ValueError('official archive checksum mismatch')
    refs=[]; counts=collections.Counter(); alphabet=collections.Counter()
    with zipfile.ZipFile(args.archive) as z:
        for member in sorted(z.namelist()):
            if not member.endswith('.bpseq'): continue
            raw=z.read(member); seq,pairs=read_bpseq(raw.decode())
            group=member.split('/')[1]; dataset={'from_Rfam14.5':'Rfam14.5','from_bpRNA-1m':'bpRNA-1m'}[group]
            pk=bool(crossing_pairs(pairs)); length=len(seq)
            size='short' if length<=150 else 'medium' if length<=500 else 'long'
            identifier=hashlib.sha256(member.encode()).hexdigest()[:20]
            fasta=root/'inputs'/f'{identifier}.fa'; fasta.parent.mkdir(exist_ok=True)
            fasta.write_text('>'+Path(member).stem+'\n'+seq+'\n')
            refs.append({'id':identifier,'member':member,'dataset':dataset,'length':length,
                'size':size,'pk':pk,'sequence':seq,'pairs':pairs,'input_sha256':digest(fasta),
                'reference_sha256':hashlib.sha256(raw).hexdigest()})
            counts[(dataset,size,pk)]+=1; alphabet.update(seq)
    with gzip.open(root/'references.json.gz','wt') as f: json.dump(refs,f,separators=(',',':'))
    inventory={'doi':'10.5281/zenodo.4923158','url':'https://zenodo.org/records/4923158',
        'download_url':'https://zenodo.org/records/4923158/files/ipknot_dataset.zip?download=1',
        'archive_md5':MD5,'archive_sha256':SHA256,'total':len(refs),
        'minimum_length':min(r['length'] for r in refs),'maximum_length':max(r['length'] for r in refs),
        'alphabet':dict(alphabet),'counts':[{'dataset':k[0],'size':k[1],'pk':k[2],'n':v} for k,v in sorted(counts.items())],
        'paper_table1_single_total':12788,'published_archive_single_total':11843,
        'note':'Published archive lacks 945 Rfam single sequences relative to Table 1. No sequences were filtered here.'}
    assert len(refs)==11843 and len({r['id'] for r in refs})==11843
    write_json(root/'inventory.json',inventory); print(json.dumps(inventory,indent=2))

def init_worker(cpus):
    # Executor processes have unique 1-based identities. CPUs are homogeneous X925.
    import multiprocessing
    identity=multiprocessing.current_process()._identity[-1]
    os.sched_setaffinity(0,{cpus[(identity-1)%len(cpus)]})

def invoke(payload):
    entry,root,dd,ilp,launcher,modes,repetitions,timeout,phase=payload
    root=Path(root); dest=root/phase; dest.mkdir(exist_ok=True)
    completed=0
    for rep in range(repetitions):
        order=list(modes)
        random.Random(int(entry['id'],16)+rep*7919).shuffle(order)
        for mode in order:
            target=dest/f"{entry['id']}-{mode}-{rep}.json"
            if target.exists():
                old=json.loads(target.read_text())
                if old['input_sha256']!=entry['input_sha256']: raise ValueError('changed input')
                completed+=1; continue
            stem=target.with_suffix(''); prediction=Path(str(stem)+'.bpseq'); resources=Path(str(stem)+'.resources')
            binary=ilp if mode=='ilp' else dd
            options=['--decoder','ilp' if mode=='ilp' else 'dd',*COMMON]
            if mode!='ilp': options+=['--dd-max-iter',mode[2:],'--dd-patience','0']
            cmd=[binary,*options,'-B',str(prediction),str(root/'inputs'/f"{entry['id']}.fa")]
            env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1'}
            env.pop('IPKNOT_HIGHS_LOG',None)
            run=subprocess.run([launcher,str(resources),str(timeout),*cmd],capture_output=True,text=True,env=env)
            if not resources.exists(): raise RuntimeError('launcher failed: '+run.stderr)
            result=json.loads(resources.read_text())
            result.update(id=entry['id'],mode=mode,repetition=rep,phase=phase,
                input_sha256=entry['input_sha256'],cpu=next(iter(os.sched_getaffinity(0))),
                command=cmd,stdout_sha256=hashlib.sha256(run.stdout.encode()).hexdigest(),stderr=run.stderr)
            if result['exit_code']==0:
                seq,pairs=read_bpseq(prediction.read_text())
                if seq!=entry['sequence']: raise RuntimeError('prediction/input sequence mismatch: '+entry['member'])
                result.update(pairs=pairs,accuracy=metrics(entry['pairs'],pairs),
                    crossing_accuracy=metrics(crossing_pairs(entry['pairs']),crossing_pairs(pairs)))
            write_json(target,result)
            prediction.unlink(missing_ok=True); resources.unlink(missing_ok=True)
            completed+=1
    return completed

def run(args):
    refs=load_index(args.root); selected=select(refs,args.sample_per_stratum)
    cpus=[int(x) for x in args.cpus.split(',')]
    manifest={'phase':args.phase,'modes':args.modes,'n_sequences':len(selected),'repetitions':args.repetitions,
        'cpus':cpus,'workers':len(cpus),'common_options':COMMON,'dd_patience':0,
        'dd_beam':100,'dd_crossing_beam':100,'dd_witnesses':16,'dd_bounds':'default (off)',
        'pk_scores':'off','timeout_seconds':args.timeout,'dd_binary':str(args.dd),
        'ilp_binary':str(args.ilp),'dd_sha256':digest(args.dd),'ilp_sha256':digest(args.ilp),
        'launcher_sha256':digest(args.launcher),'script_sha256':digest(__file__),
        'architecture':platform.machine(),'uname':list(platform.uname()),
        'sample_per_stratum':args.sample_per_stratum,'sample_ids':[r['id'] for r in selected],
        'measurement':'C launcher CLOCK_MONOTONIC wall, wait4 child CPU and ru_maxrss (KiB); includes BPP/refinement; no traces',
        'scheduling':'One sequence owns one CPU for every mode; shuffled mode order determined only by input hash and repetition.',
        'started_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())}
    path=args.root/(args.phase+'-manifest.json')
    if path.exists():
        previous=json.loads(path.read_text())
        for k in ['modes','n_sequences','repetitions','cpus','dd_sha256','ilp_sha256','launcher_sha256','sample_ids','timeout_seconds']:
            if previous[k]!=manifest[k]: raise ValueError('resume manifest changed: '+k)
        manifest=previous
    else: write_json(path,manifest)
    # Expensive sequences first prevents a single long sequence from dominating the tail.
    selected.sort(key=lambda r:(-r['length'],r['id']))
    tasks=[(r,str(args.root),str(args.dd),str(args.ilp),str(args.launcher),args.modes,
        args.repetitions,args.timeout,args.phase) for r in selected]
    start=time.monotonic(); completed=0; last=start
    with concurrent.futures.ProcessPoolExecutor(max_workers=len(cpus),initializer=init_worker,initargs=(cpus,),mp_context=multiprocessing.get_context("fork")) as pool:
        futures={pool.submit(invoke,t):t[0]['id'] for t in tasks}
        for future in concurrent.futures.as_completed(futures):
            completed+=future.result(); now=time.monotonic()
            if now-last>30:
                print(json.dumps({'phase':args.phase,'measurements_completed':completed,
                    'measurements_total':len(tasks)*len(args.modes)*args.repetitions,'elapsed_seconds':round(now-start,1)}),flush=True)
                last=now
    manifest['completed_utc']=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
    manifest['invocation_elapsed_seconds']=time.monotonic()-start
    write_json(path,manifest); print('COMPLETE',completed,flush=True)

def select(refs,per):
    if not per: return list(refs)
    strata=collections.defaultdict(list)
    for r in refs: strata[(r['dataset'],r['size'],r['pk'])].append(r)
    out=[]
    for key,rows in sorted(strata.items()):
        rows.sort(key=lambda r:hashlib.sha256(('timing2022/'+r['member']).encode()).digest())
        out.extend(rows[:per])
    return out

def main():
    ap=argparse.ArgumentParser(__doc__); ap.add_argument('action',choices=['prepare','run'])
    ap.add_argument('--root',type=Path,required=True);ap.add_argument('--archive',type=Path)
    ap.add_argument('--dd',type=Path);ap.add_argument('--ilp',type=Path);ap.add_argument('--launcher',type=Path)
    ap.add_argument('--modes',nargs='+',default=MODES,choices=MODES)
    ap.add_argument('--cpus',default='19');ap.add_argument('--repetitions',type=int,default=1)
    ap.add_argument('--timeout',type=int,default=1800);ap.add_argument('--phase',default='full')
    ap.add_argument('--sample-per-stratum',type=int,default=0)
    args=ap.parse_args(); prepare(args) if args.action=='prepare' else run(args)
if __name__=='__main__': main()
