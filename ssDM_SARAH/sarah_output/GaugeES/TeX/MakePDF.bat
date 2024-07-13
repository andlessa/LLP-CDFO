pdflatex SSDMv2-GaugeES.tex 
cd Diagrams 
FOR %%I IN (*.mp) DO MPOST "%%I" 
cd .. 
pdflatex SSDMv2-GaugeES.tex 
