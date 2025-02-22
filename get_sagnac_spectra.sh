#!/bin/bash

path_to_images="/import/freenas-ffb-01-data/romy_plots/"
path_to_html_figures="/import/kilauea-data/HTML_Monitor/figures/"

year=$(date --utc +%Y)

for r in Z U V W; do
        path_to_spectra="${path_to_images}${year}/R${r}/spectra/"
	file=$(ls -tp ${path_to_spectra} | grep -v /$ | tail -1)
	cp ${path_to_spectra}${file} ${path_to_html_figures}"html_sagnacspectra_R${r}.png"
done

#EOF
