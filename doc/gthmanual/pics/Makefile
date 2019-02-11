all: dataflow.dvi

dataflow.inc: dataflow.pic
	pic -t dataflow.pic > $@

dataflow.dvi: dataflow.tex dataflow.inc dataflow.pic
	latex dataflow.tex
	latex dataflow.tex

dataflow.ps: dataflow.dvi
	dvips -o $@ $<

dataflow.pdf: dataflow.ps
	ps2pdf dataflow.ps
