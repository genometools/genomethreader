GTHLIBS:=lib/libgenomethreader.a\
         lib/libvmatch.a\
         lib/libvmengine.a\
         lib/libmkvtree.a\
         lib/libkurtz.a\
         lib/libkurtz-basic.a\
         lib/liblicensemanager.a\
         ../genometools/lib/libgenometools.a

all: bin/gth bin/gthconsensus

bin/gth: obj/gth.o $(GTHLIBS)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(CC) $^ -lm -ldl -o $@

bin/gthconsensus: obj/gthconsensus.o $(GTHLIBS)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(CC) $^ -lm -ldl -o $@
