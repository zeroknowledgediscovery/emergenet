import pandas as pd
geof=pd.read_csv('finalseq.csv')

import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
import contextily as ctx
import geopandas
import geoplot as gplt
import pylab as plt
import pandas as pd
import matplotlib.cm as cm
import matplotlib as mpl
import numpy as np
import json

geo_DF=geopandas.GeoDataFrame(
    geof, crs="EPSG:4326",geometry=geopandas.points_from_xy(geof.longitude, geof.latitude))
df=geo_DF.to_crs('epsg:4326')
df_=df.to_crs('epsg:4326')
df__ = df.to_crs(epsg=3857) # reproject it in Web mercator
geo_DF__ = geo_DF.to_crs(epsg=3857) # reproject it in Web mercator

def getColor(x,cmap='jet',VMIN=.5,VMAX=1.0,alpha=None):
    if alpha is None:
        alpha=1
    else:
        alpha = (((x-VMIN)/(VMAX-VMIN)))
    if alpha == 1:
        alpha=.999
    norm = mpl.colors.Normalize(vmin=VMIN, vmax=VMAX)
    m = cm.ScalarMappable(norm=norm, cmap=cmap)
    ctup=np.array(list(m.to_rgba(x)))
    ctup[3]=alpha
    return tuple(ctup)

def getColor(x,cmap='jet',VMIN=.5,VMAX=1.0):
    norm = mpl.colors.Normalize(vmin=VMIN, vmax=VMAX)
    m = cm.ScalarMappable(norm=norm, cmap=cmap)
    return m.to_rgba(x)


    
def saveFIG(filename='tmp.pdf'):
    import pylab as plt
    plt.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, 
            hspace = 0, wspace = 0)
    plt.margins(0,0)
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    plt.savefig(filename,dpi=600, bbox_inches = 'tight',pad_inches = 0,transparent=True) 
    return

geo_DF__=geo_DF__.rename(columns={'IRATe':'emergence_risk'})
allriskystrains=pd.read_csv('../../../paper_data_v1/irat_qnet/results/animal_predictions/allriskystrains_collapsed.csv',index_col=0)
allriskystrains=allriskystrains[allriskystrains.emergence_risk>6]
nonrisky_geo_df=geo_DF__[geo_DF__.emergence_risk<6]


riskystrains=[x for x in geo_DF__.id if x in allriskystrains.index.values]
risky_geo_df=geo_DF__.set_index('id').loc[riskystrains,:].reset_index()
pf=risky_geo_df.copy()
npf=nonrisky_geo_df.copy()

def mssize2(x):
    return ((5**(x-2))) + (x>6)*300 + 20


def plotRisk2(df,ax,variable='emergence_risk',ALPHA=.2,COL=None,mssize=mssize2,colR=1,colG=.3,colB=.3,colalpha=.2,
             cmap='jet',VMIN=0,VMAX=None,markersize=20,markeredgecolor='w'):
    if VMAX is None:
        VMAX=1

    MS = lambda x: mssize2(x) 
    
    df.plot(
        ax=ax,
        markersize=markersize* MS(df[variable]),
        edgecolor=markeredgecolor,lw=.5,
        color=(colR,colG,colB,colalpha),#'k',#getColor(df[variable],cmap=cmap,VMIN=VMIN,VMAX=VMAX),
        #alpha=ALPHA
    )
    
    return ax    

cmap = cm.gist_rainbow
cmap = cm.gray
TVAR='emergence_risk'

fname1='./WB_Coastlines_10m/WB_Coastlines_10m.shp'
df1 = geopandas.read_file(fname1)
plt.style.use('dark_background')
plt.style.use('seaborn-whitegrid')
df1__ = df1.to_crs(epsg=3857) # reproject it in Web mercator
ax1 = df1__.plot(figsize=(20,20), alpha=.0, edgecolor='w')
VMIN=pf[TVAR].min()
VMAX=pf[TVAR].max()

ALPHA=1
MS=1.5*.5

plotRisk2(npf,ax=ax1,cmap=cmap,markersize=MS,ALPHA=ALPHA,mssize=mssize2,
         markeredgecolor=(0,0,0,1),variable=TVAR,VMIN=VMIN,colalpha=.025,colR=.1,colG=0.1,colB=.1,
         VMAX=VMAX)

plotRisk2(pf[pf.subtype=='H7N9'],ax=ax1,cmap=cmap,markersize=MS,ALPHA=ALPHA,mssize=mssize2,
         markeredgecolor='#777777',variable=TVAR,VMIN=VMIN,colalpha=.4,colR=0.3,colG=.7,colB=0.3,
         VMAX=VMAX)

plotRisk2(pf[pf.subtype=='H5N1'],ax=ax1,cmap=cmap,markersize=MS,ALPHA=ALPHA,mssize=mssize2,
         markeredgecolor='#446688',variable=TVAR,VMIN=VMIN,colalpha=.4,colR=0.1,colG=.4,colB=0.6,
         VMAX=VMAX)


plotRisk2(pf[pf.subtype=='H9N2'],ax=ax1,cmap=cmap,markersize=MS,ALPHA=ALPHA,mssize=mssize2,
         markeredgecolor='#BD890F',variable=TVAR,VMIN=VMIN,colalpha=1,colR=0.99,colG=.72,colB=0.07,
         VMAX=VMAX)

plotRisk2(pf[pf.subtype=='H1N1'],ax=ax1,cmap=cmap,markersize=MS,ALPHA=ALPHA,mssize=mssize2,
         markeredgecolor='#eeeeee',variable=TVAR,VMIN=VMIN,colalpha=.5,colR=1,colG=0,colB=0,
         VMAX=VMAX)

plotRisk2(pf[pf.subtype=='H3N2'],ax=ax1,cmap=cmap,markersize=MS,ALPHA=ALPHA,mssize=mssize2,
         markeredgecolor=(1,1,1,1),variable=TVAR,VMIN=VMIN,colalpha=.5,colR=0,colG=0,colB=1,
         VMAX=VMAX)
# plotRisk2(sf,ax=ax1,cmap=cmap,markersize=MS,ALPHA=ALPHA,mssize=mssize2,
#          markeredgecolor=(0,0,0,1),variable=TVAR,VMIN=VMIN,colalpha=.5,colR=0,colG=0,colB=1,
#          VMAX=VMAX)

#ctx.add_basemap(ax1,source=ctx.providers.Stamen.TonerLite,alpha=.5)
ctx.add_basemap(ax1,source=ctx.providers.CartoDB.Positron,alpha=.7)
#ctx.add_basemap(ax1,source=ctx.providers.Stamen.Toner)
#ctx.add_basemap(ax1,source=ctx.providers.Stamen.Watercolor,alpha=.5)
#ctx.add_basemap(ax1,source=ctx.providers.NASAGIBS.ViirsEarthAtNight2012,alpha=1)
#ax1.set_ylim(-.85e7,.985e7)
ax1.set_ylim(-.1e7,1.1e7)
ax1.set_xlim(-1.6e7,1.65e7)
ax1.set_axis_off()
saveFIG(filename='bionorad_upd_X.pdf')


