import pyRserve
conn = pyRserve.connect()
conn.r.ls()

import pandas as pd
mydata = pd.read_csv("data/OutputData.csv")

coords= mydata.copy()
coords.drop(coords.columns[0:3],axis=1,inplace=True)
conn.r.coords=coords.to_numpy()

size=mydata.centeroid.to_numpy()
conn.r.size=size

conn.eval('require(geomorph)')
#conn.eval('require(Morpho)')
conn.eval('arr=arrayspecs(coords,p=41,k=3)')
conn.voidEval('gdf=geomorph.data.frame(Size=size,Coords=arr)')
conn.voidEval('outlm=procD.lm(Coords ~ Size,data=gdf)')
model=conn.eval('outlm')
#conn.eval('remove(list=ls())')
conn.shutdown 
#to shutdown Rserve. Otherwise it will stay open and contain all variable from the current session).
