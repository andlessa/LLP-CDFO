#!/bin/bash 
pdflatex SSDMv2-GaugeES.tex 
cd Diagrams 
find . -name "*.mp" -exec mpost {} \; 
cd .. 
pdflatex SSDMv2-GaugeES.tex 
echo "" 
echo "PDF for Model finished" 
echo "Thanks for using SARAH" 
echo "" 
