#!/usr/bin/env bash

fluxed_spectrum=$1 # textfile listing fluxed spectra
output_spectrum=$2 # 1 for same detector/order, 2 for different
output_response=$3 # combined spectrum

rgsfmat << EOF
$fluxed_spectrum
$output_spectrum
$output_response
EOF
