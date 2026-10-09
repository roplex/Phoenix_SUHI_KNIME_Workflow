# Builds Figure 7.11 (forest plot of rural-baseline scenarios) and Figure 7.12
# (per-scene baseline / urban LST / SUHI) from the outputs of
# rural_reference_zone_sensitivity_analysis.R.
# Reads  <PHX_ROOT>/Sensitivity/sensitivity_summary.csv
#        <PHX_ROOT>/Sensitivity/sensitivity_temporal_by_scene.csv
# Writes Figure_7_11_sensitivity_forest.png, Figure_7_12_sensitivity_temporal.png
#        into the same folder. Requires matplotlib. Usage: PHX_ROOT=/path python make_sensitivity_figures.py
import csv, os, datetime as dt
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
D=os.path.join(os.environ.get('PHX_ROOT','PhoenixData2'),'Sensitivity')+os.sep
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.edgecolor':'#52514e','axes.labelcolor':'#0b0b0b','text.color':'#0b0b0b','xtick.color':'#52514e','ytick.color':'#52514e'})
BLUE,ORG,AQ,MAG,GRY='#2a78d6','#eb6834','#1baf7a','#e87ba4','#52514e'
rows={ (r['family'],r['scenario']):r for r in csv.DictReader(open(D+'sensitivity_summary.csv'))}
def get(fam,pref):
    m=[r for (f,s),r in rows.items() if f==fam and s.startswith(pref)]; assert len(m)==1,(pref,m); return float(m[0]['baseline_C']), m[0]['n_pixels']
G=[]  # (group, color, [(label, value, n, hollow)])
G.append(('Reference',GRY,[('Chapter definition (±100 m, 60 km)',)+get('Primary','Chapter definition')+(False,),
  ('No elevation matching (full rural sample)',)+get('Primary','Same, no elevation')+(False,),
  ('Original municipal-boundary mask (as reported)',44.8,'–',False)]))
A=[('Matched ±%d m'%t,)+get('A_elevation','Matched +/-%d m'%t)+(t>200,) for t in (25,50,75,150,200,300,500)]
A+=[('Lower half-window (urban −100 m to urban)',)+get('A_elevation','Lower half')+(False,),
    ('Upper half-window (urban to urban +100 m)',)+get('A_elevation','Upper half')+(False,),
    ('Regression-adjusted, linear',)+get('A_elevation','Regression-adjusted, linear')+(False,),
    ('Regression-adjusted, quadratic',)+get('A_elevation','Regression-adjusted, quadratic')+(False,),
    ('Re-weighted to urban elevation distribution',)+get('A_elevation','Histogram')+(False,)]
G.append(('A. Elevation handling',BLUE,A))
B=[('Shrub+grass+crop+bare, ≥80% of cell',)+get('B_landcover','Shrub+grass+crop+bare, >=80%')+(False,),
   ('Excluding cropland',)+get('B_landcover','Excluding cropland')+(False,),
   ('Native desert only (shrub+grass)',)+get('B_landcover','Native desert')+(False,),
   ('All non-built, non-water classes',)+get('B_landcover','All non-built')+(False,),
   ('Built-up fraction ≤2%',)+get('B_landcover','Any land cover with built-up fraction <= 2%')+(False,),
   ('Built-up fraction ≤20%',)+get('B_landcover','Any land cover with built-up fraction <= 20%')+(False,)]
G.append(('B. Land-cover exclusions',ORG,B))
C=[('Inner exclusion %g km (outer 60 km)'%r,)+get('C_buffer','Inner exclusion %g km, outer radius 60 km'%r)+(False,) for r in (0,2.5,5,10,15,20)]
C+=[('Inner 0 km, outer 20 km',)+get('C_buffer','Inner exclusion 0 km, outer radius 20 km')+(False,),
    ('Inner 20 km, outer 30 km',)+get('C_buffer','Inner exclusion 20 km, outer radius 30 km')+(False,)]
G.append(('C. Buffer from urbanized area',AQ,C))
Dg=[('Leave out %s sector'%s,)+get('D_distribution','Leave out %s'%s)+(False,) for s in 'ESW']
Dg+=[('Sector-balanced mean',)+get('D_distribution','Sector-balanced')+(False,),
     ('Block bootstrap median (95% interval shaded)',)+get('D_distribution','Spatial block bootstrap')+(False,)]
G.append(('D. Retained-pixel distribution',MAG,Dg))
n_rows=sum(len(g[2])+1 for g in G)
fig,ax=plt.subplots(figsize=(7.4,0.225*n_rows+0.95))
from matplotlib.transforms import blended_transform_factory as btf
tr=btf(ax.transAxes,ax.transData)
ytl=[]
y=0; yt=[];yl=[]
ax.axvspan(47.60,48.32,color='#d9d8d2',alpha=.6,lw=0,zorder=0)
ax.axvline(47.96,color='#0b0b0b',ls='--',lw=.9,zorder=1)
for gname,col,items in G:
    ax.text(-0.575,y,gname,transform=tr,fontweight='bold',fontsize=9.2,va='center',ha='left',clip_on=False); y+=1
    for lab,v,n,hol in items:
        ax.plot([44,v],[y,y],color='#e6e5df',lw=.7,zorder=1)
        ax.plot(v,y,'o',ms=5.2,mfc='white' if hol else col,mec=col,mew=1.3,zorder=3)
        ax.text(-0.015,y,lab,transform=tr,va='center',ha='right',fontsize=8.6,color='#0b0b0b' if not hol else GRY,clip_on=False)
        ax.text(v+.09,y,'%.2f'%v,va='center',fontsize=7.6,color=GRY)
        y+=1
ax.set_ylim(y-.3,-.8); ax.set_xlim(44,49.6); ax.set_yticks([])
ax.set_xlabel('Rural baseline LST (°C)')
for s in ('top','right','left'): ax.spines[s].set_visible(False)
ax.grid(axis='x',color='#e6e5df',lw=.6); ax.set_axisbelow(True)
ax.legend(handles=[Line2D([],[],color='#0b0b0b',ls='--',lw=.9,label='Chapter baseline (47.96 °C)'),
  plt.Rectangle((0,0),1,1,fc='#d9d8d2',label='Bootstrap 95% interval (47.60–48.32 °C)'),
  Line2D([],[],marker='o',mfc='white',mec=GRY,ls='',label='Diagnostic only (tolerance no longer matches urban core)')],
  loc='upper right',bbox_to_anchor=(0.0,-0.06),ncol=1,frameon=False,fontsize=7)
fig.subplots_adjust(left=0.575,right=0.985,top=0.995,bottom=0.15)
fig.savefig(D+'Figure_7_11_sensitivity_forest.png',dpi=300)
# temporal
T=list(csv.DictReader(open(D+'sensitivity_temporal_by_scene.csv')))
d=[dt.date.fromisoformat(r['date']) for r in T]
ru=[float(r['rural_C']) for r in T]; ur=[float(r['urban_C']) for r in T]; su=[float(r['SUHI_C']) for r in T]
fig,(a1,a2)=plt.subplots(2,1,figsize=(6.6,4.9),sharex=True,gridspec_kw={'height_ratios':[1.25,1]})
a1.plot(d,ru,'-o',color=ORG,ms=4.5,lw=1.6,label='Rural baseline (elevation-matched cells)')
a1.plot(d,ur,'-o',color=BLUE,ms=4.5,lw=1.6,label='City-wide mean LST (sharpened, within city limit)')
a1.axhline(47.96,color='#0b0b0b',ls='--',lw=.9)
a1.text(d[0],47.96+.25,'composite baseline 47.96 °C',ha='left',fontsize=7,color=GRY)
a1.set_ylabel('LST (°C)'); a1.legend(loc='lower right',frameon=False,fontsize=7.2)
a2.axhline(0,color=GRY,lw=.8)
mean=sum(su)/len(su); a2.axhline(mean,color='#0b0b0b',ls='--',lw=.9)
a2.text(d[-1],mean-.35,'mean −%.2f °C'%abs(mean),ha='right',fontsize=7,color=GRY)
a2.bar(d,su,width=5,color=AQ)
a2.set_ylabel('City-wide SUHI (°C)\n(urban − rural, same scene)'); a2.set_xlabel('Start date of 8-day MOD11A2 scene, 2024')
for a in (a1,a2):
    for s in ('top','right'): a.spines[s].set_visible(False)
    a.grid(axis='y',color='#e6e5df',lw=.6); a.set_axisbelow(True)
import matplotlib.dates as md
a2.xaxis.set_major_locator(md.MonthLocator()); a2.xaxis.set_major_formatter(md.DateFormatter('%b'))
fig.tight_layout(); fig.savefig(D+'Figure_7_12_sensitivity_temporal.png',dpi=300)
print(len(d),min(su),max(su),mean)
