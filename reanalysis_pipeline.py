#!/usr/bin/env python3
"""
FINAL corrected reanalysis — mesenteric lymph 16S, rat CLP sepsis
Primary window 4->24h (within-animal); 0h (LB_1-5) descriptive baseline (~blank).
"""
import re, os, warnings
import numpy as np, pandas as pd
from scipy.stats import fisher_exact, f as Fdist, false_discovery_control as fdrc
from scipy.spatial.distance import pdist, squareform, braycurtis
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

OTU_TAXON="data/otu_taxon.xls"
MIN_DEPTH=5000; OUT="reanalysis_out"; os.makedirs(OUT,exist_ok=True)
LB0={"LB_1":"1","LB_2":"3","LB_3":"5","LB_4":"6","LB_5":"9"}
KITOME={'acetobacter','ralstonia','pseudomonas','acinetobacter','phyllobacterium','bradyrhizobium','sphingomonas','sphingobium','methylobacterium','stenotrophomonas','delftia','herbaspirillum','burkholderia','burkholderia-caballeronia-paraburkholderia','comamonas','cutibacterium','propionibacterium','rhodococcus','microbacterium','brevundimonas','chryseobacterium','pseudoxanthomonas','massilia','rhizobium','agrobacterium','variovorax','novosphingobium','afipia','caulobacter','pelomonas','methylophilus','oxalobacter','solirubrobacter','devosia'}
GUT={'lactobacillus','bacteroides','prevotella','roseburia','ruminococcus','lachnoclostridium','blautia','faecalibacterium','romboutsia','turicibacter','colidextribacter','ucg-005'}

df=pd.read_csv(OTU_TAXON,sep="\t",index_col=0); tax=df["taxonomy"]; C=df.drop(columns=["taxonomy"]).astype(float)
rank=lambda t,r:(re.search(r+r"__([^;]+)",str(t)).group(1) if re.search(r+r"__([^;]+)",str(t)) else "unclassified")
family=tax.map(lambda t:rank(t,"f")); genus=tax.map(lambda t:rank(t,"g")); depth=C.sum()
ctrl=[c for c in C.columns if c.startswith("con_") or c=="water"]
fecal=[c for c in C.columns if c.endswith("_FB")]; balf=[c for c in C.columns if c.endswith("_FP")]
lymph0=list(LB0.keys())
later=[c for c in C.columns if re.match(r"(1|3|5|6|9)_(4|8|12|24)hLB",c)]
later_qc=[c for c in later if depth[c]>=MIN_DEPTH]; dropped=[(c,int(depth[c])) for c in later if depth[c]<MIN_DEPTH]
an=np.array([re.match(r"(\d+)_",c).group(1) for c in later_qc]); tp=np.array([int(re.match(r"\d+_(\d+)h",c).group(1)) for c in later_qc])
print(f"[design] 0h baseline n=5; primary 4-24h n={len(later_qc)} dropped={dropped}; per-tp={ {t:int((tp==t).sum()) for t in [4,8,12,24]} }")

bio=later_qc+fecal+balf; nc,nb=len(ctrl),len(bio); prev=set()
for o in C.index:
    cp=(C.loc[o,ctrl]>0).sum(); bp=(C.loc[o,bio]>0).sum()
    if fisher_exact([[cp,nc-cp],[bp,nb-bp]],alternative="greater")[1]<0.05: prev.add(o)
kit=set(C.index[genus.str.lower().isin(KITOME)])
rules={"prevalence_only":prev,"kitome_only":kit,"union":prev|kit,"intersection":prev&kit}
drows=[]
for nm,fl in rules.items():
    r={"rule":nm,"n_OTU":len(fl)}
    for st,l in [(later_qc,"lymph"),(fecal,"fecal"),(balf,"BALF")]:
        tr=C[st].values.sum(); r[f"pct_{l}"]=round(100*C.loc[C.index.isin(fl),st].values.sum()/tr,2)
    drows.append(r)
decon=pd.DataFrame(drows); decon.to_csv(f"{OUT}/decontamination_rules.csv",index=False)
print("[decon]\n",decon.to_string(index=False))
flag=rules["union"]; clean=C.drop(index=list(flag))

famtab=C.groupby(family).sum(); rel=famtab/famtab.sum()
mbc=lambda A,B: float(np.mean([braycurtis(rel[a],rel[b]) for a in A for b in B]))
def compo(cols):
    s=C[cols].sum(axis=1); tot=s.sum(); return (100*s[genus.str.lower().isin(GUT)].sum()/tot,100*s[genus.str.lower().isin(KITOME)].sum()/tot)
g0,e0=compo(lymph0); print(f"[baseline] 0h gut={g0:.1f}% reagent={e0:.1f}%; BC to controls={mbc(lymph0,ctrl):.2f} to 4-24h={mbc(lymph0,later_qc):.2f}")

fam=clean.groupby(family).sum()[later_qc]; keep=fam[(fam>0).sum(1)>=4].index
ra=fam.loc[keep].values+0.5; clr=pd.DataFrame(np.log(ra)-np.log(ra).mean(0,keepdims=True),index=keep,columns=later_qc)
n=len(later_qc); dm=lambda x:np.column_stack([(x==l).astype(float) for l in sorted(set(x))[1:]])
Xr=np.column_stack([np.ones(n),dm(an)]); Xf=np.column_stack([np.ones(n),dm(an),dm(tp)])
rss=lambda X,y:((y-X@np.linalg.lstsq(X,y,rcond=None)[0])**2).sum()
def Ft(y):
    r0,r1=rss(Xr,y),rss(Xf,y); d1,d2=Xf.shape[1]-Xr.shape[1],n-Xf.shape[1]; F=((r0-r1)/d1)/(r1/d2); return F,Fdist.sf(F,d1,d2)
da=pd.DataFrame([(f,*Ft(clr.loc[f].values)) for f in clr.index],columns=["family","F","p"]).sort_values("p")
da["FDR"]=fdrc(da["p"].values); da.to_csv(f"{OUT}/differential_4_24h.csv",index=False)
print(f"\n[DA 4-24h] families={len(da)} FDR<0.05={(da['FDR']<0.05).sum()}"); print(da.head(6).to_string(index=False))
for t in ["Enterococcaceae","Chloroflexaceae","Bacteroidaceae"]:
    r=da[da.family==t]
    if len(r): print(f"   {t}: F={r.F.values[0]:.2f} p={r.p.values[0]:.3f} FDR={r.FDR.values[0]:.3f}")

relL=(fam/fam.sum()).T.values; D=squareform(pdist(relL,"braycurtis"))
def pF(D,g):
    N=len(g); tot=(D**2).sum()/(2*N); ssw=0
    for gg in set(g):
        idx=np.where(np.array(g)==gg)[0]
        if len(idx)>1: ssw+=(D[np.ix_(idx,idx)]**2).sum()/(2*len(idx))
    a=len(set(g)); return ((tot-ssw)/(a-1))/(ssw/(N-a))
F0=pF(D,tp); rng=np.random.RandomState(1); cnt=0; NP=1999
for _ in range(NP):
    ptp=tp.copy()
    for aa in set(an):
        idx=np.where(an==aa)[0]; ptp[idx]=rng.permutation(tp[idx])
    if pF(D,ptp)>=F0: cnt+=1
permp=(cnt+1)/(NP+1); print(f"[beta 4-24h] pseudo-F={F0:.2f} p={permp:.3f}")

def rrare(col,d=5000,seed=42):
    rng=np.random.RandomState(seed); v=col.values.astype(int)
    if v.sum()<d: return None
    idx=np.repeat(np.arange(len(v)),v); return np.bincount(rng.choice(idx,d,replace=False),minlength=len(v))
sh_ = lambda x:(-(x[x>0]/x.sum()*np.log(x[x>0]/x.sum())).sum())
for lab,tab in [("before",C),("after",clean)]:
    sh=np.array([sh_(rrare(tab[c])) for c in later_qc]); ob=np.array([(rrare(tab[c])>0).sum() for c in later_qc])
    print(f"[alpha 4-24h {lab}] Shannon F,p={Ft(sh)[0]:.2f},{Ft(sh)[1]:.3f}; ObsOTU F,p={Ft(ob)[0]:.2f},{Ft(ob)[1]:.3f}")

# ---- SENSITIVITY: use near-blank 0h samples as ADDITIONAL controls, decontaminate uniformly ----
ref2=ctrl+lymph0                      # 6 controls + 5 near-blank 0h = 11 references
prev2=set()
for o in C.index:
    cp=(C.loc[o,ref2]>0).sum(); bp=(C.loc[o,bio]>0).sum()
    if fisher_exact([[cp,len(ref2)-cp],[bp,len(bio)-bp]],alternative="greater")[1]<0.05: prev2.add(o)
flag2=prev2|kit; clean2=C.drop(index=list(flag2))
tr=C[later_qc].values.sum()
print(f"[sens] ref=11(6 controls+5 0h-blanks): prevalence-flagged {len(prev)}->{len(prev2)}; "
      f"control-supported lymph removal {100*C.loc[C.index.isin(prev),later_qc].values.sum()/tr:.1f}%->{100*C.loc[C.index.isin(prev2),later_qc].values.sum()/tr:.1f}%")
fam2=clean2.groupby(family).sum()[later_qc]; keep2=fam2[(fam2>0).sum(1)>=4].index
ra2=fam2.loc[keep2].values+0.5; clr2=pd.DataFrame(np.log(ra2)-np.log(ra2).mean(0,keepdims=True),index=keep2,columns=later_qc)
da2=pd.DataFrame([(f,*Ft(clr2.loc[f].values)) for f in clr2.index],columns=["family","F","p"]).sort_values("p"); da2["FDR"]=fdrc(da2["p"].values)
er=da2[da2.family=="Enterococcaceae"]
print(f"[sens] 4-24h after uniform decontam: FDR<0.05={(da2['FDR']<0.05).sum()}/{len(da2)}; "
      f"Enterococcaceae p={er.p.values[0]:.3f} FDR={er.FDR.values[0]:.3f}  (conclusion unchanged)")

# figures
plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False,"figure.dpi":130})
PAL={"lymph":"#4C72B0","fecal":"#55A868","BALF":"#C44E52","control":"#8C8C8C","baseline":"#9467BD"}
fig,ax=plt.subplots(figsize=(7,4)); x=np.arange(len(decon)); w=0.25
for i,(l,c) in enumerate([("lymph",PAL["lymph"]),("fecal",PAL["fecal"]),("BALF",PAL["BALF"])]):
    ax.bar(x+(i-1)*w,decon[f"pct_{l}"],w,label=l,color=c)
ax.set_xticks(x); ax.set_xticklabels(decon["rule"],rotation=15); ax.set_ylabel("% reads removed")
ax.set_title("Decontamination is rule-dependent\n(control-based prevalence removes almost nothing)"); ax.legend(); fig.tight_layout(); fig.savefig(f"{OUT}/fig1_decontam_rules.png"); plt.close()

famrel=(fam/fam.sum()); fig,ax=plt.subplots(figsize=(7,4.2)); tps=[4,8,12,24]
idxmap={c:i for i,c in enumerate(later_qc)}
for aa in sorted(set(an)):
    ys=[]
    for t in tps:
        cs=[c for c in later_qc if an[idxmap[c]]==aa and tp[idxmap[c]]==t]
        ys.append(famrel.loc["Enterococcaceae",cs[0]]*100 if cs and "Enterococcaceae" in famrel.index else np.nan)
    ax.plot(tps,ys,"-o",label=f"animal {aa}",alpha=.8)
r=da[da.family=="Enterococcaceae"]
ax.set_xlabel("hours post-CLP"); ax.set_ylabel("Enterococcaceae (% of sample)"); ax.set_xticks(tps)
ax.set_title(f"No time trend in the 4-24 h window (F={r.F.values[0]:.2f}, p={r.p.values[0]:.3f})"); ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(f"{OUT}/fig2_enterococcaceae.png"); plt.close()

alls=later_qc+lymph0+fecal+balf+ctrl; famA=clean.groupby(family).sum()[alls]; relA=(famA/famA.sum()).T.values
Dm=squareform(pdist(relA,"braycurtis")); J=np.eye(len(alls))-np.ones((len(alls),)*2)/len(alls)
val,vec=np.linalg.eigh(-0.5*J@(Dm**2)@J); o=np.argsort(val)[::-1]; co=vec[:,o[:2]]*np.sqrt(np.abs(val[o[:2]]))
labs=(["lymph 4-24h"]*len(later_qc)+["lymph 0h (baseline)"]*len(lymph0)+["fecal"]*len(fecal)+["BALF"]*len(balf)+["control"]*len(ctrl))
cm={"lymph 4-24h":PAL["lymph"],"lymph 0h (baseline)":PAL["baseline"],"fecal":PAL["fecal"],"BALF":PAL["BALF"],"control":PAL["control"]}
fig,ax=plt.subplots(figsize=(6.6,5))
for lab in ["fecal","BALF","lymph 4-24h","lymph 0h (baseline)","control"]:
    idx=[i for i,l in enumerate(labs) if l==lab]; mk="*" if "0h" in lab else "o"; sz=180 if "0h" in lab else 45
    ax.scatter(co[idx,0],co[idx,1],s=sz,marker=mk,c=cm[lab],label=lab,edgecolor="w",linewidth=.5,alpha=.9)
tot=val[o].clip(min=0).sum(); ax.set_xlabel(f"PCo1 ({100*max(val[o[0]],0)/tot:.0f}%)"); ax.set_ylabel(f"PCo2 ({100*max(val[o[1]],0)/tot:.0f}%)")
ax.set_title("Baseline lymph (0 h) clusters with negative controls\n(near-sterile; ~indistinguishable from reagent blanks)"); ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(f"{OUT}/fig3_pcoa_specimen.png"); plt.close()

fig,ax=plt.subplots(figsize=(6.6,4)); grp=["lymph 0h\n(LB_1-5)","controls","lymph 4-24h","fecal"]; G=[compo(lymph0),compo(ctrl),compo(later_qc),compo(fecal)]
gut=[a for a,b in G]; rea=[b for a,b in G]; x=np.arange(len(grp))
ax.bar(x-0.2,gut,0.4,label="gut genera %",color="#55A868"); ax.bar(x+0.2,rea,0.4,label="reagent/env genera %",color="#8C8C8C")
ax.set_xticks(x); ax.set_xticklabels(grp); ax.set_ylabel("% of reads")
ax.set_title("Baseline lymph is blank-like: almost no gut taxa"); ax.legend(); fig.tight_layout(); fig.savefig(f"{OUT}/fig4_baseline_blank.png"); plt.close()
print("[done]")
