GTHMANUALTMP=gthmanual-tmp

GTHMANUALSRC=../../VERSION\
             ${GTHMANUALTMP}.tex\
             author.tex\
             LocalDefs.sty\
             speciestable.tex

ISPELLALL=author.ispell\
          gthmanual.ispell\
          optiontable.ispell\
          cutout.ispell\
          speciestable.ispell

all: ispellall\
     gthmanual.pdf

${GTHMANUALTMP}.tex: gthmanual.tex execcommand.pl
	#findemptyoutfiles.sh
	execcommand.pl $< > $@

gthmanual.pdf: ${GTHMANUALSRC}
	touch ${GTHMANUALTMP}.ind
	pdflatex ${GTHMANUALTMP}
	pdflatex ${GTHMANUALTMP}
	makeindex -q ${GTHMANUALTMP}.idx
	bibtex ${GTHMANUALTMP}
	pdflatex ${GTHMANUALTMP}
	bibtex ${GTHMANUALTMP}
	pdflatex ${GTHMANUALTMP}
	pdflatex ${GTHMANUALTMP}
	mv ${GTHMANUALTMP}.pdf $@

gthmanual.dvi: ${GTHMANUALSRC}
	touch ${GTHMANUALTMP}.ind
	latex ${GTHMANUALTMP}
	latex ${GTHMANUALTMP}
	makeindex -q ${GTHMANUALTMP}.idx
	bibtex ${GTHMANUALTMP}
	latex ${GTHMANUALTMP}
	bibtex ${GTHMANUALTMP}
	latex ${GTHMANUALTMP}
	latex ${GTHMANUALTMP}
	mv ${GTHMANUALTMP}.dvi $@

gthlicense.dvi: gthlicense.tex
	latex gthlicense
	latex gthlicense

gthlicense.pdf: gthlicense.tex
	pdflatex gthlicense
	pdflatex gthlicense

%.ispell: %.tex
	ispell -d american $<
	@touch $@

ispellall: ${ISPELLALL}

%.ps: %.dvi
	dvips $< -o $@

clean:
	rm -f *.toc *.ilg *.out *.idx *.ind *.log *.aux *.bbl *.blg *.lot *.lof
	rm -f ${GTHMANUALTMP}.tex
	rm -f output/*.out
	rm -f ceres.inter.gz new.inter.gz
	rm -f ceres_full.inter.gz.* 
	rm -f U89959.inter.gz.* 
	rm -f cdna
	rm -f *.scr*

cleanup: clean
	rm -f *.dvi *.ps *.ispell
