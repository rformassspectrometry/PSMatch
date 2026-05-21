## This file explains how the spBoekweg and psmBoekweg data are created.
## In reality, it is a subset of the bulk Boekweg files from the `MsDataHub`
## package. See ?Boekweg2022 for more details.

library(PSMatch)
library(sager)
library(MsDataHub)
library(Spectra)

############
### bulk ###
############
f <- MsDataHub::OR11_20160122_PG_HeLa_CVB3_CT_A.mzML()
f <- "/home/guillaumedeflandre/Documents/drive_UCL/PHD/data/boekweg2022/OR11_20160122_PG_HeLa_CVB3_CT_A.mzML"
sp <- Spectra(f)

p <- MsDataHub::OR11_20160122_PG_HeLa_CVB3_CT_A.sage.tsv()

psms <- sagePSM(p)

# # remove duplicated spectrumId for same rank (for TMT_Erwinia)
# n <- which(duplicated(psms$spectrumID) & psms$rank == 1)
# psms <- psms[-n, ]

format(object.size(psms), units = "Mb") ## 44.7 Mb is too much

psmSubset <- sort(sample(NROW(psms), 10000))
psms <- psms[psmSubset, ]

format(object.size(psms), units = "Mb") ## 4.8 Mb is better

psmBoekweg <- psms

save(psmBoekweg, file = "../../data/psmBoekweg.rda")

psms <- filterPsmRank(psms)

head(psms$pkey <- paste0(basename(psms$filename), sub("^.+scan=", "::", psms$scannr)))
head(sp$pkey <- paste0(basename(sp$dataOrigin), sub("^.+scan=", "::", sp$spectrumId)))

spJoined <- joinSpectraData(sp, psms, by.x = "pkey")
format(object.size(spJoined), units = "Mb") ## 101.2 Mb: too big, take only a sample

ms2 <- spJoined[which(spJoined$msLevel == 2 & !is.na(spJoined$peptide))]
ms1 <- spJoined[which(spJoined$msLevel == 1)]

sampleMs2_prep <- sample(ms2, length(ms2) / 4)
sampleMs2 <- which(spJoined$scanIndex %in% sampleMs2_prep$scanIndex)
sampleMs1 <- which(spJoined$scanIndex %in% sampleMs2_prep$precScanNum)

spJoinedSampled <- spJoined[c(sampleMs1, sampleMs2)]
o <- order(spJoinedSampled$rtime, decreasing = FALSE)
spJoined <- spJoinedSampled[o] ## order by retention time

format(object.size(spJoined), units = "Mb") ## 1.6 Mb, much better

# spMatchedErwinia <- spJoined
# save(spMatchedErwinia, file = "../../data/spMatchedErwinia.rda")

## only subset original Spectra object but now linked to identifications
n <- which(sp$pkey %in% spJoined$pkey)
spBoekweg <- sp[n]

spBoekweg$dataOrigin <- basename(spBoekweg$dataOrigin)
spBoekweg <- setBackend(spBoekweg, MsBackendMemory())
save(spBoekweg, file =
    "/home/guillaumedeflandre/Documents/drive_UCL/PHD/rformassspectrometry/PSMatch-oriented/guideflandre/PSMatch/data/spBoekweg.rda")