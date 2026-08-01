#!/usr/bin/env python3
"""
reproduce_all.py  —  Deterministic, dependency-light reproduction of EVERY
downstream table and figure in:

  "A design-aware reanalysis of ultra-low-biomass mesenteric-lymph 16S data
   reveals batch- and contamination-sensitive temporal associations"

Input : otu_taxon.xls  (Majorbio OTU-with-taxonomy table, 8,042 OTUs)
Output: ./repro_out/  (CSV per table + figure-source CSVs) and a printed log.
        Every number in the manuscript, supplementary tables S1-S6 and the
        sample-accounting table is generated here. Fixed RNG seeds make the
        output bit-for-bit reproducible.

Run   : python3 reproduce_all.py  path/to/otu_taxon.xls
Deps  : numpy, scipy, pandas, matplotlib  (no scikit-bio / statsmodels needed)
"""
import sys, re, os, warnings
import numpy as np, pandas as pd
from scipy.stats import fisher_exact, f as Fdist, false_discovery_control as fdrc
from scipy.spatial.distance import pdist, squareform, braycurtis
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

# Input OTU table: 1st CLI arg, else bundled data/otu_taxon.xls next to this script.
_here = os.path.dirname(os.path.abspath(__file__))
OTU_TAXON = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_here, "data", "otu_taxon.xls")
if not os.path.exists(OTU_TAXON):
    OTU_TAXON = "/mnt/user-data/uploads/合并/workflow_results/2.OTU_Analysis/OTU_Taxon/OtuTaxon_summary/otu_taxon.xls"
MIN_DEPTH = 5000
OUT = "./repro_out"; os.makedirs(OUT, exist_ok=True)
SEED_PERM = 1          # PERMANOVA restricted-permutation seed
SEED_RARE = 42         # rarefaction seed
NP = 1999              # permutations

LB0 = {"LB_1":"1","LB_2":"3","LB_3":"5","LB_4":"6","LB_5":"9"}
KITOME = {'acetobacter','ralstonia','pseudomonas','acinetobacter','phyllobacterium','bradyrhizobium',
 'sphingomonas','sphingobium','methylobacterium','stenotrophomonas','delftia','herbaspirillum','burkholderia',
 'burkholderia-caballeronia-paraburkholderia','comamonas','cutibacterium','propionibacterium','rhodococcus',
 'microbacterium','brevundimonas','chryseobacterium','pseudoxanthomonas','massilia','rhizobium','agrobacterium',
 'variovorax','novosphingobium','afipia','caulobacter','pelomonas','methylophilus','oxalobacter','solirubrobacter','devosia'}
GUT = {'lactobacillus','bacteroides','prevotella','roseburia','ruminococcus','lachnoclostridium','blautia',
 'faecalibacterium','romboutsia','turicibacter','colidextribacter','ucg-005'}

# ---------- load ----------
df = pd.read_csv(OTU_TAXON, sep="\t", index_col=0)
tax = df["taxonomy"]; C = df.drop(columns=["taxonomy"]).astype(float)
rank = lambda t,r:(re.search(r+r"__([^;]+)",str(t)).group(1) if re.search(r+r"__([^;]+)",str(t)) else "unclassified")
family = tax.map(lambda t:rank(t,"f")); genus = tax.map(lambda t:rank(t,"g")); depth = C.sum()

ctrl  = [c for c in C.columns if c.startswith("con_") or c=="water"]
fecal = [c for c in C.columns if c.endswith("_FB")]
balf  = [c for c in C.columns if c.endswith("_FP")]
lymph0 = list(LB0.keys())
later = [c for c in C.columns if re.match(r"(1|3|5|6|9)_(4|8|12|24)hLB", c)]
later_qc = [c for c in later if depth[c] >= MIN_DEPTH]
dropped  = [(c,int(depth[c])) for c in later if depth[c] < MIN_DEPTH]
an = np.array([re.match(r"(\d+)_",c).group(1) for c in later_qc])
tp = np.array([int(re.match(r"\d+_(\d+)h",c).group(1)) for c in later_qc])

bio = later_qc + fecal + balf                 # 18 + 15 + 15 = 48 biological samples
print(f"[design] biological samples in prevalence screen n={len(bio)} "
      f"(lymph 4-24h QC={len(later_qc)}, faecal={len(fecal)}, BALF={len(balf)}); "
      f"0h baseline n={len(lymph0)}; controls n={len(ctrl)}; dropped low-depth={dropped}")
assert len(bio) == 48, "biological-sample count must be 48"

# ---------- helpers ----------
def prevalence_flag(reference, biological):
    nc, nb = len(reference), len(biological); flagged={}
    for o in C.index:
        cp=(C.loc[o,reference]>0).sum(); bp=(C.loc[o,biological]>0).sum()
        p=fisher_exact([[cp,nc-cp],[bp,nb-bp]],alternative="greater")[1]
        if p<0.05: flagged[o]=(int(cp),int(bp),float(p))
    return flagged

def clr_Ftest(clean_counts, cols, an_, tp_):
    fam=clean_counts.groupby(family).sum()[cols]
    keep=fam[(fam>0).sum(1)>=4].index
    ra=fam.loc[keep].values+0.5
    clr=pd.DataFrame(np.log(ra)-np.log(ra).mean(0,keepdims=True),index=keep,columns=cols)
    n=len(cols)
    dm=lambda x:np.column_stack([(x==l).astype(float) for l in sorted(set(x))[1:]]) if len(set(x))>1 else np.empty((n,0))
    Xr=np.column_stack([np.ones(n),dm(an_)]); Xf=np.column_stack([np.ones(n),dm(an_),dm(tp_)])
    rss=lambda X,y:((y-X@np.linalg.lstsq(X,y,rcond=None)[0])**2).sum()
    d1=Xf.shape[1]-Xr.shape[1]; d2=n-Xf.shape[1]
    rows=[]
    for f in clr.index:
        y=clr.loc[f].values; F=((rss(Xr,y)-rss(Xf,y))/d1)/(rss(Xf,y)/d2)
        rows.append((f,F,float(Fdist.sf(F,d1,d2))))
    da=pd.DataFrame(rows,columns=["family","F","p"]).sort_values("p")
    da["BH_FDR"]=fdrc(da["p"].values); da["significant"]=da["BH_FDR"]<0.05
    return da.reset_index(drop=True)

def permanova(cols, clean_counts, an_, tp_, seed=SEED_PERM, restrict=True):
    """restrict=True: permute tp_ within animal strata (repeated-measures design,
    e.g. the 4-24h time-course). restrict=False: free permutation of group labels
    (independent-sample comparison, e.g. baseline vs controls)."""
    fam=clean_counts.groupby(family).sum()[cols]; rel=(fam/fam.sum()).T.values
    D=squareform(pdist(rel,"braycurtis"))
    def pF(g):
        N=len(g); tot=(D**2).sum()/(2*N); ssw=0
        for gg in set(g):
            idx=np.where(np.array(g)==gg)[0]
            if len(idx)>1: ssw+=(D[np.ix_(idx,idx)]**2).sum()/(2*len(idx))
        a=len(set(g)); return ((tot-ssw)/(a-1))/(ssw/(N-a))
    F0=pF(tp_); rng=np.random.RandomState(seed); cnt=0
    for _ in range(NP):
        if restrict:
            ptp=tp_.copy()
            for aa in set(an_):
                idx=np.where(an_==aa)[0]; ptp[idx]=rng.permutation(tp_[idx])
        else:
            ptp=rng.permutation(tp_)
        if pF(ptp)>=F0: cnt+=1
    N=len(tp_); tot=(D**2).sum()/(2*N); ssw=0
    for gg in set(tp_):
        idx=np.where(tp_==gg)[0]
        if len(idx)>1: ssw+=(D[np.ix_(idx,idx)]**2).sum()/(2*len(idx))
    R2=(tot-ssw)/tot
    return float(F0),float(R2),(cnt+1)/(NP+1)

def permdisp(groups_cols, clean_counts, labels, seed=SEED_PERM):
    allc=[c for cols in groups_cols for c in cols]
    fam=clean_counts.groupby(family).sum()[allc]; rel=(fam/fam.sum()).T.values
    D=squareform(pdist(rel,"braycurtis")); lab=np.array(labels)
    def dispstat(lab_):
        # distance to group centroid in PCoA space
        N=len(lab_); J=np.eye(N)-np.ones((N,N))/N
        val,vec=np.linalg.eigh(-0.5*J@(D**2)@J); pos=val>1e-8
        Y=vec[:,pos]*np.sqrt(val[pos])
        res=[]
        for g in set(lab_):
            idx=np.where(lab_==g)[0]; cen=Y[idx].mean(0)
            res.append(np.sqrt(((Y[idx]-cen)**2).sum(1)))
        gmeans=[r.mean() for r in res]; return abs(gmeans[0]-gmeans[1])
    F0=dispstat(lab); rng=np.random.RandomState(seed); cnt=0
    for _ in range(NP):
        if dispstat(rng.permutation(lab))>=F0: cnt+=1
    return (cnt+1)/(NP+1)

# ---------- Table 1 / S6-input: four decontamination rules ----------
prev = prevalence_flag(ctrl, bio)
prev_set=set(prev); kit=set(C.index[genus.str.lower().isin(KITOME)])
rules={"prevalence_only":prev_set,"kitome_only":kit,"union":prev_set|kit,"intersection":prev_set&kit}
drows=[]
for nm,fl in rules.items():
    r={"rule":nm,"n_OTU":len(fl)}
    for st,l in [(later_qc,"lymph"),(fecal,"fecal"),(balf,"BALF")]:
        tr=C[st].values.sum(); r[f"pct_{l}"]=round(100*C.loc[C.index.isin(fl),st].values.sum()/tr,2)
    drows.append(r)
decon=pd.DataFrame(drows); decon.to_csv(f"{OUT}/Table1_decontamination_rules.csv",index=False)
print("\n[Table 1] four decontamination rules\n",decon.to_string(index=False))

# ---------- S2: removed OTUs with Fisher p, BH-FDR, flags ----------
flag=rules["union"]
prev_p={o:prev[o][2] for o in prev}
bh=fdrc(np.array([prev_p[o] for o in prev_set])) if prev_set else np.array([])
prev_bh=dict(zip(list(prev_set),bh))
s2=[]
for o in sorted(flag):
    cp=(C.loc[o,ctrl]>0).sum(); bp=(C.loc[o,bio]>0).sum()
    s2.append({"OTU":o,"taxonomy":tax[o],"genus":genus[o],"family":family[o],
      "ctrl_prevalence":round(cp/len(ctrl),3),"bio_prevalence":round(bp/len(bio),3),
      "fisher_p":round(prev_p.get(o,np.nan),4) if o in prev_p else "",
      "fisher_BH_FDR":round(prev_bh[o],4) if o in prev_bh else "",
      "flagged_by_prevalence":o in prev_set,"flagged_by_catalogue":o in kit,
      "reads_all_lymph":int(C.loc[o,later_qc+lymph0].sum()),"reads_lymph_4_24h":int(C.loc[o,later_qc].sum())})
pd.DataFrame(s2).to_csv(f"{OUT}/S2_removed_OTUs.csv",index=False)
nbh=int((bh<0.05).sum()) if len(bh) else 0
print(f"[S2] union removed {len(flag)} OTUs; prevalence-flagged {len(prev_set)}, of which {nbh} survive BH-FDR<0.05")

# ---------- controls-vs-lymph-only sensitivity ----------
prev_lymphonly=set(prevalence_flag(ctrl, later_qc).keys())
tr=C[later_qc].values.sum()
pct_lo=100*C.loc[C.index.isin(prev_lymphonly),later_qc].values.sum()/tr
pct_union_lo=100*C.loc[C.index.isin(prev_lymphonly|kit),later_qc].values.sum()/tr
print(f"[sensitivity: controls vs lymph only] control-prevalence lymph removal={pct_lo:.2f}%; union={pct_union_lo:.2f}%")

# ---------- S1 / Table 2: primary DA (union, 4-24h) ----------
clean=C.drop(index=list(flag))
da=clr_Ftest(clean, later_qc, an, tp)
da.to_csv(f"{OUT}/S1_primary_families.csv",index=False)
print(f"\n[S1/Table 2] primary DA families={len(da)} FDR<0.05={int(da.significant.sum())}")
for t in ["norank_o__Lineage_IV","Chloroflexaceae","Bacteroidaceae","Enterococcaceae"]:
    r=da[da.family==t]
    if len(r): print(f"   {t}: F={r.F.values[0]:.2f} p={r.p.values[0]:.4f} FDR={r.BH_FDR.values[0]:.3f}")

# ---------- PERMANOVA + PERMDISP (primary) ----------
F0,R2,pp=permanova(later_qc, clean, an, tp)
print(f"[PERMANOVA 4-24h] pseudo-F={F0:.2f} R2={R2:.2f} p={pp:.3f}")

# ---------- baseline vs controls (Fig 3/4) ----------
def compo(cols):
    s=C[cols].sum(axis=1); tot=s.sum()
    return (round(100*s[genus.str.lower().isin(GUT)].sum()/tot,2),
            round(100*s[genus.str.lower().isin(KITOME)].sum()/tot,2))
g0,e0=compo(lymph0)
# family-level Bray-Curtis (consistent with PERMANOVA basis); exclude self-pairs
_bccols=lymph0+ctrl; _famrel=(C.groupby(family).sum()[_bccols]); _famrel=_famrel/_famrel.sum()
mbc=lambda A,B: float(np.mean([braycurtis(_famrel[a],_famrel[b]) for a in A for b in B if a!=b]))
F0b,R2b,ppb=permanova(lymph0+ctrl, C, np.array(["base"]*5+["ctrl"]*len(ctrl)),
                       np.array([0]*5+[1]*len(ctrl)), restrict=False)
disp_p=permdisp([lymph0,ctrl], C, ["base"]*5+["ctrl"]*len(ctrl))
print(f"[baseline] 0h gut={g0}% reagent={e0}%; PERMANOVA vs controls pseudo-F={F0b:.2f} R2={R2b:.2f} p={ppb:.2f}; "
      f"PERMDISP p={disp_p:.3f}; mean BC 0h->controls={mbc(lymph0,ctrl):.2f}, within-controls={mbc(ctrl,ctrl):.2f}")
pd.DataFrame({"group":["lymph_0h","controls","lymph_4_24h","fecal"],
   "gut_pct":[compo(x)[0] for x in [lymph0,ctrl,later_qc,fecal]],
   "reagent_pct":[compo(x)[1] for x in [lymph0,ctrl,later_qc,fecal]]}).to_csv(f"{OUT}/Fig3_baseline_composition.csv",index=False)

# ---------- leave-one-animal-out (Enterococcaceae) ----------
loo=[]
for hold in sorted(set(an)):
    keepidx=[i for i,a in enumerate(an) if a!=hold]
    cols=[later_qc[i] for i in keepidx]
    da_l=clr_Ftest(clean, cols, an[keepidx], tp[keepidx])
    r=da_l[da_l.family=="Enterococcaceae"]
    loo.append((hold,float(r.p.values[0]) if len(r) else np.nan))
print("[LOO Enterococcaceae] "+", ".join(f"drop {h}: p={p:.2f}" for h,p in loo)
      +f"  (range {min(p for _,p in loo):.2f}-{max(p for _,p in loo):.2f})")
pd.DataFrame(loo,columns=["held_out_animal","Enterococcaceae_p"]).to_csv(f"{OUT}/S_LOO_enterococcaceae.csv",index=False)

# ---------- Table 3 / S5: stepwise ablation ----------
def ablation_step(cols, an_, tp_, use_rm):
    if use_rm:
        da_s=clr_Ftest(clean, cols, an_, tp_)
    else:  # samples independent -> animal term dropped
        da_s=clr_Ftest(clean, cols, np.array(["all"]*len(cols)), tp_)
    ent=da_s[da_s.family=="Enterococcaceae"]
    return da_s, len(da_s), int(da_s.significant.sum()), (float(ent.p.values[0]) if len(ent) else np.nan)
# naive: 0h included, low-depth kept, independent
naive_cols=later+lymph0
naive_an=np.array([ (LB0[c] if c in LB0 else re.match(r"(\d+)_",c).group(1)) for c in naive_cols])
naive_tp=np.array([ (0 if c in LB0 else int(re.match(r"\d+_(\d+)h",c).group(1))) for c in naive_cols])
steps=[]
d0,n0,s0,e0p=ablation_step(naive_cols,naive_an,naive_tp,use_rm=False); steps.append(("0_naive",n0,s0,e0p))
# drop 0h (still independent, low-depth kept)
d1,n1,s1,e1p=ablation_step(later,np.array(["all"]*len(later)),
      np.array([int(re.match(r"\d+_(\d+)h",c).group(1)) for c in later]),use_rm=False); steps.append(("1_drop0h",n1,s1,e1p))
# depth QC (independent)
d2,n2,s2c,e2p=ablation_step(later_qc,np.array(["all"]*len(later_qc)),tp,use_rm=False); steps.append(("2_depthQC",n2,s2c,e2p))
# repeated measures (primary)
d3,n3,s3,e3p=ablation_step(later_qc,an,tp,use_rm=True); steps.append(("3_repeated_measures",n3,s3,e3p))
print("\n[Table 3] stepwise ablation (families_tested, FDR-sig, Enterococcaceae p)")
for nm,nt,sg,ep in steps: print(f"   {nm:22s} tested={nt} FDR-sig={sg} Ent.p={ep:.3f}")
s5=pd.concat([d.assign(analysis_step=nm) for d,nm in
   [(d0,"0_naive"),(d1,"1_drop0h"),(d2,"2_depthQC"),(d3,"3_repeated_measures")]])[
   ["analysis_step","family","F","p","BH_FDR","significant"]]
s5.to_csv(f"{OUT}/S5_ablation_full.csv",index=False)

# ---------- S6: four-rule statistics ----------
s6=[]
for nm,fl in rules.items():
    cl=C.drop(index=list(fl)); da_r=clr_Ftest(cl,later_qc,an,tp)
    ent=da_r[da_r.family=="Enterococcaceae"]
    F_,R_,p_=permanova(later_qc,cl,an,tp)
    # alpha on cleaned, rarefied
    def rrare(col,d=5000,seed=SEED_RARE):
        rng=np.random.RandomState(seed); v=col.values.astype(int)
        if v.sum()<d: return None
        idx=np.repeat(np.arange(len(v)),v); return np.bincount(rng.choice(idx,d,replace=False),minlength=len(v))
    sh_=lambda x:(-(x[x>0]/x.sum()*np.log(x[x>0]/x.sum())).sum())
    sh=np.array([sh_(rrare(cl[c])) for c in later_qc]); ob=np.array([(rrare(cl[c])>0).sum() for c in later_qc])
    def Ft_simple(y):
        n=len(later_qc); dm=lambda x:np.column_stack([(x==l).astype(float) for l in sorted(set(x))[1:]])
        Xr=np.column_stack([np.ones(n),dm(an)]); Xf=np.column_stack([np.ones(n),dm(an),dm(tp)])
        rss=lambda X:((y-X@np.linalg.lstsq(X,y,rcond=None)[0])**2).sum()
        d1=Xf.shape[1]-Xr.shape[1]; d2=n-Xf.shape[1]; F=((rss(Xr)-rss(Xf))/d1)/(rss(Xf)/d2); return float(Fdist.sf(F,d1,d2))
    s6.append({"decontam_rule":nm,"families_tested":len(da_r),"FDR_sig_families":int(da_r.significant.sum()),
      "Enterococcaceae_p":round(float(ent.p.values[0]),3),"PERMANOVA_pseudoF":round(F_,2),
      "PERMANOVA_R2":round(R_,3),"PERMANOVA_p":round(p_,3),"Shannon_p":round(Ft_simple(sh),3),"richness_p":round(Ft_simple(ob),3)})
pd.DataFrame(s6).to_csv(f"{OUT}/S6_four_rule_stats.csv",index=False)
print("\n[S6] four-rule statistics\n",pd.DataFrame(s6).to_string(index=False))

# ---------- Section 3.5 sensitivity: 0h + controls as uniform references ----------
ref2=ctrl+lymph0
prev2=set(prevalence_flag(ref2, bio).keys())
clean2=C.drop(index=list(prev2|kit))
tr=C[later_qc].values.sum()
rem1=100*C.loc[C.index.isin(prev_set),later_qc].values.sum()/tr
rem2=100*C.loc[C.index.isin(prev2),later_qc].values.sum()/tr
da_s=clr_Ftest(clean2, later_qc, an, tp); ent_s=da_s[da_s.family=="Enterococcaceae"]
Fs,Rs,ps=permanova(later_qc, clean2, an, tp)
print(f"\n[Section 3.5 sensitivity] prevalence-flagged {len(prev_set)}->{len(prev2)}; "
      f"control-supported lymph removal {rem1:.1f}%->{rem2:.1f}%; "
      f"FDR-sig={int(da_s.significant.sum())}/{len(da_s)}; Enterococcaceae p={float(ent_s.p.values[0]):.3f}; "
      f"PERMANOVA pseudo-F={Fs:.2f} p={ps:.3f}")

# ---------- S3 / sample accounting ----------
def batch_of(c):
    if c in LB0: return "baseline_recollection"
    if c.startswith("con_"): return "baseline_recollection"
    if c=="water": return "main_4_24h"
    return "main_4_24h"
def status_of(c):
    if c in LB0: return "included (0h baseline, descriptive + sensitivity only)"
    if c in later_qc: return "included (primary 4-24h)"
    if c in later: return f"excluded (low depth, {int(depth[c])} reads)"
    if c in ctrl: return "included (negative control, sensitivity only)"
    return "included (faecal/BALF, not in lymph time-course)"
acc=[]
for c in C.columns:
    typ=("lymph" if re.match(r"(1|3|5|6|9)_(0|4|8|12|24)h",c) or c in LB0 else
         "faecal" if c in fecal else "BALF" if c in balf else "control")
    aid=(LB0[c] if c in LB0 else re.match(r"(\d+)_",c).group(1) if re.match(r"\d+_",c) else "NA")
    th=(0 if c in LB0 else int(re.match(r"\d+_(\d+)h",c).group(1)) if re.match(r"\d+_\d+h",c) else "NA")
    acc.append({"sample":c,"specimen_type":typ,"animal_id":aid,"timepoint_h":th,
      "batch":batch_of(c),"read_depth":int(depth[c]),
      "control_type":("extraction_blank" if c.startswith("con_") else "water" if c=="water" else ""),
      "status":status_of(c)})
pd.DataFrame(acc).to_csv(f"{OUT}/S3_sample_accounting.csv",index=False)
print(f"\n[S3] sample accounting written: {len(acc)} rows")

# ---------- S4 genus catalogue ----------
pd.DataFrame({"contaminant_genus_catalogue":sorted(KITOME)}).to_csv(f"{OUT}/S4_genus_catalogue.csv",index=False)

# ---------- Figures 2-5 (Figure 1 is a schematic, not generated here) ----------
plt.rcParams.update({"font.size":13,"axes.spines.top":False,"axes.spines.right":False})
PAL={"lymph":"#4C72B0","fecal":"#55A868","BALF":"#C44E52","control":"#8C8C8C","baseline":"#9467BD"}
def _save(fig,name):
    fig.savefig(f"{OUT}/{name}.png",dpi=200,bbox_inches="tight")
    fig.savefig(f"{OUT}/{name}.pdf",bbox_inches="tight"); plt.close(fig)

# Figure 2: reads removed under four rules
fig,ax=plt.subplots(figsize=(7,4.2)); x=np.arange(len(decon)); w=0.26
for i,(l,c) in enumerate([("lymph",PAL["lymph"]),("fecal",PAL["fecal"]),("BALF",PAL["BALF"])]):
    ax.bar(x+(i-1)*w,decon[f"pct_{l}"],w,label=l,color=c)
ax.set_xticks(x); ax.set_xticklabels(decon["rule"],rotation=12); ax.set_ylabel("% reads removed")
ax.legend(frameon=False); _save(fig,"Figure2_decontam_rules")

# Figure 3: gut vs reagent read fractions by group
grp=["lymph 0h\n(LB_1-5)","controls","lymph 4-24h","fecal"]; G=[compo(lymph0),compo(ctrl),compo(later_qc),compo(fecal)]
fig,ax=plt.subplots(figsize=(7,4.2)); x=np.arange(len(grp))
ax.bar(x-0.2,[a for a,b in G],0.4,label="gut-associated genera %",color="#55A868")
ax.bar(x+0.2,[b for a,b in G],0.4,label="reagent/environmental genera %",color="#8C8C8C")
ax.set_xticks(x); ax.set_xticklabels(grp); ax.set_ylabel("% of reads"); ax.legend(frameon=False); _save(fig,"Figure3_baseline_composition")

# Figure 4: PCoA of all specimens + controls (union-decontaminated)
alls=later_qc+lymph0+fecal+balf+ctrl; famA=clean.groupby(family).sum()[alls]; relA=(famA/famA.sum()).T.values
Dm=squareform(pdist(relA,"braycurtis")); Jm=np.eye(len(alls))-np.ones((len(alls),)*2)/len(alls)
val,vec=np.linalg.eigh(-0.5*Jm@(Dm**2)@Jm); o=np.argsort(val)[::-1]; co=vec[:,o[:2]]*np.sqrt(np.abs(val[o[:2]]))
labs=(["lymph 4-24h"]*len(later_qc)+["lymph 0h (baseline)"]*len(lymph0)+["fecal"]*len(fecal)+["BALF"]*len(balf)+["control"]*len(ctrl))
cm={"lymph 4-24h":PAL["lymph"],"lymph 0h (baseline)":PAL["baseline"],"fecal":PAL["fecal"],"BALF":PAL["BALF"],"control":PAL["control"]}
fig,ax=plt.subplots(figsize=(7,5.2))
for lab in ["fecal","BALF","lymph 4-24h","lymph 0h (baseline)","control"]:
    idx=[i for i,l in enumerate(labs) if l==lab]; mk="*" if "0h" in lab else "o"; sz=300 if "0h" in lab else 60
    ax.scatter(co[idx,0],co[idx,1],s=sz,marker=mk,c=cm[lab],label=lab,edgecolor="w",linewidth=.6,alpha=.9)
tot=val[o].clip(min=0).sum(); ax.set_xlabel(f"PCo1 ({100*max(val[o[0]],0)/tot:.0f}%)"); ax.set_ylabel(f"PCo2 ({100*max(val[o[1]],0)/tot:.0f}%)")
ax.legend(frameon=False); _save(fig,"Figure4_pcoa")

# Figure 5: Enterococcaceae trajectories 4-24h
famrel=(clean.groupby(family).sum()/clean.groupby(family).sum().sum())
idxmap={c:i for i,c in enumerate(later_qc)}; tps=[4,8,12,24]
fig,ax=plt.subplots(figsize=(7,4.6))
for aa in sorted(set(an)):
    ys=[]
    for t in tps:
        cs=[c for c in later_qc if an[idxmap[c]]==aa and tp[idxmap[c]]==t]
        ys.append(famrel.loc["Enterococcaceae",cs[0]]*100 if cs and "Enterococcaceae" in famrel.index else np.nan)
    ax.plot(tps,ys,"-o",label=f"animal {aa}",linewidth=2,markersize=7,alpha=.85)
ax.set_xlabel("Hours post-CLP"); ax.set_ylabel("Enterococcaceae (% of sample)"); ax.set_xticks(tps)
ax.legend(frameon=False,ncol=2); _save(fig,"Figure5_enterococcaceae")
print("[figures] Figures 2-5 written (PNG+PDF)")

print("\n[done] all outputs in", OUT)
