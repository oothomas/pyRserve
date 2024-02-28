library(SlicerMorphR)
library(geomorph)
#log=file.choose()

log="./LMs/2023-12-29_10_41_10/analysis.log"
parsed=parser(log)

dat=read.csv(paste(parsed$output.path, parsed$OutputData, sep="/"))
size=dat$centeroid
coords=arrayspecs(dat[,4:ncol(dat)], p=parsed$no.LM, k=3)
gdf=geomorph.data.frame(Size=size, Coords=coords)


size.model = as.formula(Coords~Size)
outlm=procD.lm(size.model, data=gdf)


