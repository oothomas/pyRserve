library(SlicerMorphR)
library(geomorph)
#log=file.choose()

log="~/Desktop/2025-02-25_13_40_06/analysis.json"
parsed=parser2(log)

dat=read.csv(paste(parsed$output.path, parsed$OutputData, sep="/"))
size=dat$centroid
coords=arrayspecs(dat[,4:ncol(dat)], p=parsed$no.LM, k=3)
gdf=geomorph.data.frame(Size=size, Coords=coords)


size.model = as.formula(Coords~Size)
outlm=procD.lm(size.model, data=gdf)
