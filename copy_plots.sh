#!/bin/bash


target="/import/kilauea-data/HTML_Monitor/figures/"

date=$(date -d -1days +%Y%m%d)
year=$(date -d -1days +%Y)

# copy Z spectra
cp /import/freenas-ffb-01-data/romy_plots/${year}/BJZ/spectra/BW_ROMY_10_BJZ_${date}.png ${target}BW_ROMY_10_BJZ.png

# copy UVW spectra
for r in U V W; do
    cp /import/freenas-ffb-01-data/romy_plots/${year}/BJ${r}/spectra/BW_ROMY__BJ${r}_${date}.png ${target}BW_ROMY__BJ${r}.png

done;

# copy sagnac spectra
path_to_images="/import/freenas-ffb-01-data/romy_plots/"
for r in Z U V; do
    # extend path
    path_to_spectra="${path_to_images}${year}/R${r}/spectra/"
    # find last modifed file
    filename=$(ls -tp ${path_to_spectra} | grep -v /$ | head -1)
    # copy images
    cp ${path_to_spectra}${filename} ${target}"html_sagnacspectra_R${r}.png"
done

# EOF
