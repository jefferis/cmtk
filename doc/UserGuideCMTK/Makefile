#!gmake

NAME = UserGuideCMTK

BIBTEXDATABSENAME = ${NAME}


${NAME}.pdf: ${NAME}.dvi
	dvipdfm ${NAME}


${NAME}.bbl: ${BIBTEXDATABSENAME}.bib ${NAME}.aux
	bibtex ${NAME}


${NAME}.dvi: ${NAME}.tex  ${BIBTEXDATABSENAME}.bib
	latex ${NAME}.tex
	bibtex ${NAME}
	latex ${NAME}.tex
	latex ${NAME}.tex


${NAME}.aux: ${NAME}.tex  ${BIBTEXDATABSENAME}.bib
	latex ${NAME}.tex


clean:
	rm ${NAME}.dvi ${NAME}.ps ${NAME}.pdf ${NAME}.aux ${NAME}.log ${NAME}.bbl ${NAME}.blg ${NAME}.toc ${NAME}.out  ${NAME}.brf

