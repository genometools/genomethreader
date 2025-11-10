#
# Copyright (c) 2006-2020 Gordon Gremme <gordon@gremme.org>
# Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#

INCLUDEOPT:=-I$(CURDIR)/../genometools/src \
            -I$(CURDIR)/../genometools/src/external/zlib-1.2.8 \
            -I$(CURDIR)/../genometools/src/external/bzip2-1.0.6 \
            -I$(CURDIR)/../genometools/src/external/expat-2.0.1/lib \
            -I$(CURDIR)/../genometools/src/external/lua-5.1.5/src \
            -I$(CURDIR)/src -I$(CURDIR)/obj \
# these variables are exported by the configuration script
CC:=gcc
CXX:=g++
EXP_CFLAGS:=$(CFLAGS)
EXP_CXXFLAGS:=$(CXXFLAGS)
EXP_CPPFLAGS:=$(CPPFLAGS)
EXP_LDLIBS:=$(LIBS) -lm

ifneq ($(fpic),no)
  FPIC:=-fPIC
endif

# ...while those starting with GTH_ are for internal purposes only
GTH_CFLAGS:=-g -Wall -Wunused-parameter -pipe $(FPIC) -Wpointer-arith -fno-stack-protector -Wno-unknown-pragmas
# mkvtree needs -DWITHLCP
# rnv needs -DUNISTD_H="<unistd.h>" -DEXPAT_H="<expat.h>" -DRNV_VERSION="\"1.7.8\""
EXT_FLAGS:= -DWITHLCP -DUNISTD_H="<unistd.h>" -DEXPAT_H="<expat.h>" \
            -DRNV_VERSION=\"1.7.8\"
EXP_CPPFLAGS+=-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 $(EXT_FLAGS)
GTH_CPPFLAGS:=$(INCLUDEOPT) -U_FORTIFY_SOURCE
GTH_CXXFLAGS:=-g -pipe
GTH_LDFLAGS:=-Llib
STEST_FLAGS:=
BUILDSTAMP:=$(shell date +'"%Y-%m-%d %H:%M:%S"')

# try to set RANLIB automatically
SYSTEM:=$(shell uname -s)
MACHINE:=$(shell uname -m)
ifeq ($(SYSTEM),Darwin)
  RANLIB:=ranlib
  SHARED:=-dynamiclib -undefined dynamic_lookup
  SHARED_OBJ_NAME_EXT:=.dylib
  ifeq ($(universal),yes)
    MACHINE:="Universal_Binary"
    GTH_CFLAGS+=-arch i386 -arch ppc -arch_errors_fatal
    GTH_LDFLAGS+=-arch i386 -arch ppc -arch_errors_fatal
  endif
  ifeq ($(ppc),yes)
    MACHINE:="Power_Macintosh"
    GTH_CFLAGS+=-arch ppc -arch_errors_fatal
    GTH_LDFLAGS+=-arch ppc -arch_errors_fatal
  endif
else
  SHARED_OBJ_NAME_EXT:=.so
  SHARED:=-shared
endif
ifeq ($(SYSTEM),Linux)
  ifneq ($(mingw),yes)
    EXP_LDLIBS+=-ldl -lpthread
  endif
endif
ifeq ($(SYSTEM),Windows)
  EXP_LDLIBS+=-liphlpapi
endif
ifeq ($(wrapmemcpy),yes)
  GTH_LDFLAGS+=-Wl,--wrap=memcpy
endif

ifneq ($(SYSTEM),Darwin)
ifneq ($(MACHINE),ARMv6_hf)
GTH_CFLAGS+=-Wno-error=misleading-indentation
endif
endif

# system specific stuff (concerning 64bit compilation)
ifeq ($(SYSTEM),Darwin)
  ifeq ($(64bit),yes)
    m64=yes
  else
    m32=yes
  endif
endif

ifeq ($(64bit),yes)
  ifneq ($(MACHINE),x86_64)
    m64=yes
  endif
else
  ifeq ($(MACHINE),x86_64)
    m32=yes
  endif
endif

ifeq ($(64bit),yes)
  BIT=64bit
else
  BIT=32bit
endif

ifeq ($(m32),yes)
  GTH_CFLAGS += -m32
  GTH_LDFLAGS += -m32
endif

ifeq ($(m64),yes)
  GTH_CFLAGS += -m64
  GTH_LDFLAGS += -m64
endif

# libraries for which we build replacements (that also appear in dependencies)
#EXP_LDLIBS+=-lz -lbz2 -lm -ldl
#OVERRIDELIBS:=../genometools/lib/libbz2.a

# compiled executables
GTMAIN_SRC:=src/gt.c src/gtr.c src/gtt.c src/interactive.c
GTMAIN_OBJ:=$(GTMAIN_SRC:%.c=obj/%.o)
GTMAIN_DEP:=$(GTMAIN_SRC:%.c=obj/%.d)

VSTREEDIR:=../vstree/src
GTHLIBS:=$(VSTREEDIR)/lib/$(CONFIGGUESS)/$(BIT)/libvmatch.a\
         $(VSTREEDIR)/lib/$(CONFIGGUESS)/$(BIT)/libvmengine.a\
         $(VSTREEDIR)/lib/$(CONFIGGUESS)/$(BIT)/libmkvtree.a\
         $(VSTREEDIR)/lib/$(CONFIGGUESS)/$(BIT)/libkurtz.a\
         $(VSTREEDIR)/lib/$(CONFIGGUESS)/$(BIT)/libkurtz-basic.a
GTH_CPPFLAGS+=-I$(CURDIR)/$(VSTREEDIR)/include \
             -I$(CURDIR)/$(VSTREEDIR)/Vmatch \

SERVER=stronghold
WWWBASEDIR=/var/www/htdocs

# process arguments
ifeq ($(assert),no)
  EXP_CPPFLAGS += -DNDEBUG
endif

ifneq ($(errorcheck),no)
  GTH_CFLAGS += -Werror
endif

ifneq ($(opt),no)
  GTH_CFLAGS += -O3
  GTH_CXXFLAGS += -O3
endif

ifeq ($(prof),yes)
  GTH_CFLAGS += -pg
  GTH_LDFLAGS += -pg
endif

ifdef gthtestdata
  STEST_FLAGS += -gthtestdata $(gthtestdata)
endif

ifeq ($(memcheck),yes)
  STEST_FLAGS += -memcheck
endif

RNV_DIR:=src/external/rnv-1.7.10
LIBRNV_SRC:=$(RNV_DIR)/rn.c $(RNV_DIR)/rnc.c $(RNV_DIR)/rnd.c $(RNV_DIR)/rnl.c \
            $(RNV_DIR)/rnv.c $(RNV_DIR)/rnx.c $(RNV_DIR)/drv.c \
            $(RNV_DIR)/ary.c $(RNV_DIR)/xsd.c $(RNV_DIR)/xsd_tm.c \
            $(RNV_DIR)/dxl.c $(RNV_DIR)/dsl.c $(RNV_DIR)/sc.c $(RNV_DIR)/u.c \
            $(RNV_DIR)/ht.c $(RNV_DIR)/er.c $(RNV_DIR)/xmlc.c $(RNV_DIR)/s.c \
            $(RNV_DIR)/m.c $(RNV_DIR)/rx.c
LIBRNV_OBJ:=$(LIBRNV_SRC:%.c=obj/%.o)
LIBRNV_DEP:=$(LIBRNV_SRC:%.c=obj/%.d)

RNVMAIN_SRC:=$(RNV_DIR)/xcl.c
RNVMAIN_OBJ:=$(RNVMAIN_SRC:%.c=obj/%.o)
RNVMAIN_DEP:=$(RNVMAIN_SRC:%.c=obj/%.d)

ifeq ($(wrapmemcpy),yes)
  RNVMAIN_OBJ+=../genometools/obj/src/memcpy.o
endif

LIBGENOMETHREADER_DIRS:= src/libgenomethreader src/gth

# the GenomeThreader library
LIBGENOMETHREADER_SRC:=$(foreach DIR,$(LIBGENOMETHREADER_DIRS),$(sort $(wildcard $(DIR)/*.c)))
LIBGENOMETHREADER_OBJ:=$(LIBGENOMETHREADER_SRC:%.c=obj/%.o)
LIBGENOMETHREADER_DEP:=$(LIBGENOMETHREADER_SRC:%.c=obj/%.d)

# set prefix for install target
prefix ?= /usr/local

# allow to set patch program
patch ?= patch

all: lib/libgenomethreader.a bin/gth bin/gthconsensus bin/gthbssmfileinfo \
     bin/gthbssmbuild bin/gthbssmprint bin/gthbssmrmsd bin/gthbssmtrain \
     bin/gthfilestat bin/gthmkbssmfiles bin/gthsplit bin/gthunit \
     bin/gthgetseq bin/align_dna bin/rnv

lib/libgenomethreader.a: obj/gth_config.h $(LIBGENOMETHREADER_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(AR) ru $@ $(LIBGENOMETHREADER_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

lib/librnv.a: $(LIBRNV_OBJ)
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(AR) ru $@ $(LIBRNV_OBJ)
ifdef RANLIB
	@$(RANLIB) $@
endif

define PROGRAM_template
$(1): $(2)
	@echo "[link $$(@F)]"
	@test -d $$(@D) || mkdir -p $$(@D)
	@$$(CC) $$(LDFLAGS) $$(GTH_LDFLAGS) $$(filter-out $$(OVERRIDELIBS),$$^) \
	  $$(filter-out $$(patsubst lib%.a,-l%,$$(notdir $$(OVERRIDELIBS))),\
	  $$(EXP_LDLIBS)) $$(OVERRIDELIBS) -o $$@
endef

$(eval $(call PROGRAM_template, bin/gth, obj/src/gth.o \
                                lib/libgenomethreader.a \
                                $(GTHLIBS) ../genometools/lib/libgenometools.a))

$(eval $(call PROGRAM_template, bin/gthconsensus, obj/src/gthconsensus.o \
                                lib/libgenomethreader.a \
                                $(GTHLIBS) ../genometools/lib/libgenometools.a))

$(eval $(call PROGRAM_template, bin/gthbssmfileinfo, obj/src/gthbssmfileinfo.o \
                                lib/libgenomethreader.a \
                                ../genometools/lib/libgenometools.a $(GTHLIBS)))

$(eval $(call PROGRAM_template, bin/gthbssmbuild, obj/src/gthbssmbuild.o \
                                lib/libgenomethreader.a \
                                ../genometools/lib/libgenometools.a $(GTHLIBS)))

$(eval $(call PROGRAM_template, bin/gthbssmprint, obj/src/gthbssmprint.o \
                                lib/libgenomethreader.a \
                                ../genometools/lib/libgenometools.a $(GTHLIBS)))

$(eval $(call PROGRAM_template, bin/gthbssmrmsd, obj/src/gthbssmrmsd.o \
                                lib/libgenomethreader.a \
                                ../genometools/lib/libgenometools.a $(GTHLIBS)))

$(eval $(call PROGRAM_template, bin/gthbssmtrain, obj/src/gthbssmtrain.o \
                                lib/libgenomethreader.a \
                                ../genometools/lib/libgenometools.a $(GTHLIBS)))

$(eval $(call PROGRAM_template, bin/gthfilestat, obj/src/gthfilestat.o \
                                lib/libgenomethreader.a \
                                $(GTHLIBS) ../genometools/lib/libgenometools.a))

$(eval $(call PROGRAM_template, bin/gthmkbssmfiles, obj/src/gthmkbssmfiles.o \
                                lib/libgenomethreader.a \
                                $(GTHLIBS) ../genometools/lib/libgenometools.a))

$(eval $(call PROGRAM_template, bin/gthsplit, obj/src/gthsplit.o \
                                lib/libgenomethreader.a \
                                $(GTHLIBS) ../genometools/lib/libgenometools.a))

$(eval $(call PROGRAM_template, bin/gthunit, obj/src/gthunit.o \
                                lib/libgenomethreader.a \
                                $(GTHLIBS) ../genometools/lib/libgenometools.a))

$(eval $(call PROGRAM_template, bin/gthgetseq, obj/src/gthgetseq.o \
                                lib/libgenomethreader.a \
                                ../genometools/lib/libgenometools.a $(GTHLIBS)))

$(eval $(call PROGRAM_template, bin/align_dna, obj/src/align_dna.o \
                                lib/libgenomethreader.a \
                                ../genometools/lib/libgenometools.a $(GTHLIBS)))

bin/rnv: $(RNVMAIN_OBJ) lib/librnv.a ../genometools/lib/libexpat.a
	@echo "[link $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(CC) $(EXP_LDFLAGS) $(GTH_LDFLAGS) $^ -o $@

obj/gth_config.h: VERSION
	@echo '[create $@]'
	@test -d $(@D) || mkdir -p $(@D)
	@(echo '#ifndef GTH_CONFIG_H' ;\
	echo '#define GTH_CONFIG_H' ;\
	echo '#define GTH_CC "'`$(CC) --version | head -n 1`\" ;\
	echo '#define GTH_CFLAGS "$(EXP_CFLAGS) $(GTH_CFLAGS)"' ;\
	echo '$(EXP_CPPFLAGS) $(GTH_CPPFLAGS)' | \
	sed -e 's/\([^\]\)"/\1\\"/g' -e 's/^"/\\"/g' -e 's/$$/"/' \
	    -e 's/^/#define GTH_CPPFLAGS "/'; ) > $@
	@echo '#define GTH_VERSION "'`cat VERSION`\" >> $@
	@cat VERSION | \
          sed 's/\([0-9]*\)\.[0-9]*\.[0-9]*/#define GTH_MAJOR_VERSION \1/' >> $@
	@cat VERSION | \
          sed 's/[0-9]*\.\([0-9]*\)\.[0-9]*/#define GTH_MINOR_VERSION \1/' >> $@
	@cat VERSION | \
          sed 's/[0-9]*\.[0-9]*\.\([0-9]*\)/#define GTH_MICRO_VERSION \1/' >> $@
	@echo '#endif' >> $@

obj/amalgamation.c: $(LIBGENOMETOOLS_PRESRC)
	@echo '[create $@]'
	@test -d $(@D) || mkdir -p $(@D)
	@scripts/create_amalgamation $(LIBGENOMETOOLS_PRESRC) > $@

define COMPILE_template
$(1): $(2)
	@echo "[compile $$(@F)]"
	@test -d $$(@D) || mkdir -p $$(@D)
	@$$(CC) -c $$< -o $$@ $$(EXP_CPPFLAGS) $$(GTH_CPPFLAGS) $$(EXP_CFLAGS) \
	  $$(GTH_CFLAGS) $(3)
	@$$(CC) -c $$< -o $$(@:.o=.d) $$(EXP_CPPFLAGS) $$(GTH_CPPFLAGS) \
        $(3) -MM -MP -MT $$@
endef

$(eval $(call COMPILE_template, obj/%.o, %.c))

obj/%.o: %.cxx
	@echo "[compile $(@F)]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(CXX) -c $< -o $@ $(EXP_CPPFLAGS) $(GTH_CPPFLAGS) \
	  $(EXP_CXXFLAGS) $(GTH_CXXFLAGS)
	@$(CXX) -c $< -o $(@:.o=.d) $(EXP_CPPFLAGS) $(GTH_CPPFLAGS) -MM -MP \
	  -MT $@

obj/%.o: %.cpp
	@echo "[compile $@]"
	@test -d $(@D) || mkdir -p $(@D)
	@$(CXX) -c $< -o $@ $(EXP_CPPFLAGS) $(GTH_CPPFLAGS) \
	  $(EXP_CXXFLAGS) $(GTH_CXXFLAGS)
	@$(CXX) -c $< -o $(@:.o=.d) $(EXP_CPPFLAGS) $(GTH_CPPFLAGS) -MM -MP \
	  -MT $@

src/libgenomethreader/gthversionfunc.c: obj/gth_config.h

# read dependencies
-include $(LIBGENOMETHREADER_DEP) \
         $(LIBRNV_DEP) \
       	 $(RNVMAIN_DEP)

.SUFFIXES:
.PHONY: dist libdist srcdist release gt install docs manuals installwww patch \
        push splint test train clean cleanup hmmer_get

VERSION:="`cat $(CURDIR)/VERSION`"
SYSTEMNAME:="$(SYSTEM)_$(MACHINE)"
DISTDIR:="$(CURDIR)/dist/$(SYSTEMNAME)"

GTHDISTBASENAME:="gth-$(VERSION)-$(SYSTEMNAME)-${BIT}"
GTHDISTDIR:="$(DISTDIR)/$(GTHDISTBASENAME)"

GTHLIBDISTBASENAME:="gthlibs-$(VERSION)-$(SYSTEMNAME)-${BIT}"
GTHLIBDISTDIR:="$(DISTDIR)/$(GTHLIBDISTBASENAME)"

STRIP:=strip

dist: all
	@echo "[build gth distribution]"
	@rm -rf $(GTHDISTDIR)
ifeq ($(SYSTEM),Windows)
	@rm -f $(DISTDIR)/$(GTHDISTBASENAME).zip
else
	@rm -f $(DISTDIR)/$(GTHDISTBASENAME).tar.gz
endif
	@mkdir -p $(GTHDISTDIR)
	@mkdir -p $(GTHDISTDIR)/bin/bssm
	@mkdir -p $(GTHDISTDIR)/bin/gthdata
	@mkdir -p $(GTHDISTDIR)/doc
	@mkdir -p $(GTHDISTDIR)/data
	@cp $(CURDIR)/gthdata/README $(GTHDISTDIR)
	@cp $(CURDIR)/CHANGELOG $(GTHDISTDIR)
ifeq ($(SYSTEM),Windows)
	@cp bin/gth $(GTHDISTDIR)/bin/gth.exe
	@$(STRIP) $(GTHDISTDIR)/bin/gth.exe
	@cp bin/gthconsensus $(GTHDISTDIR)/bin/gthconsensus.exe
	@$(STRIP) $(GTHDISTDIR)/bin/gthconsensus.exe
	@cp bin/gthbssmfileinfo $(GTHDISTDIR)/bin/gthbssmfileinfo.exe
	@$(STRIP) $(GTHDISTDIR)/bin/gthbssmfileinfo.exe
	@cp bin/gthbssmbuild $(GTHDISTDIR)/bin/gthbssmbuild.exe
	@$(STRIP) $(GTHDISTDIR)/bin/gthbssmbuild.exe
	@cp bin/gthbssmrmsd $(GTHDISTDIR)/bin/gthbssmrmsd.exe
	@$(STRIP) $(GTHDISTDIR)/bin/gthbssmrmsd.exe
	@cp bin/gthbssmtrain $(GTHDISTDIR)/bin/gthbssmtrain.exe
	@$(STRIP) $(GTHDISTDIR)/bin/gthbssmtrain.exe
	@cp bin/gthfilestat $(GTHDISTDIR)/bin/gthfilestat.exe
	@$(STRIP) $(GTHDISTDIR)/bin/gthfilestat.exe
	@cp bin/gthgetseq $(GTHDISTDIR)/bin/gthgetseq.exe
	@$(STRIP) $(GTHDISTDIR)/bin/gthgetseq.exe
	@cp bin/gthsplit $(GTHDISTDIR)/bin/gthsplit.exe
	@$(STRIP) $(GTHDISTDIR)/bin/gthsplit.exe
else
	@cp bin/gth $(GTHDISTDIR)/bin
	@$(STRIP) $(GTHDISTDIR)/bin/gth
	@cp bin/gthconsensus $(GTHDISTDIR)/bin
	@$(STRIP) $(GTHDISTDIR)/bin/gthconsensus
	@cp bin/gthbssmfileinfo $(GTHDISTDIR)/bin
	@$(STRIP) $(GTHDISTDIR)/bin/gthbssmfileinfo
	@cp bin/gthbssmbuild $(GTHDISTDIR)/bin
	@$(STRIP) $(GTHDISTDIR)/bin/gthbssmbuild
	@cp bin/gthbssmrmsd $(GTHDISTDIR)/bin
	@$(STRIP) $(GTHDISTDIR)/bin/gthbssmrmsd
	@cp bin/gthbssmtrain $(GTHDISTDIR)/bin
	@$(STRIP) $(GTHDISTDIR)/bin/gthbssmtrain
	@cp bin/gthfilestat $(GTHDISTDIR)/bin
	@$(STRIP) $(GTHDISTDIR)/bin/gthfilestat
	@cp bin/gthgetseq $(GTHDISTDIR)/bin
	@$(STRIP) $(GTHDISTDIR)/bin/gthgetseq
	@cp bin/gthsplit $(GTHDISTDIR)/bin
	@$(STRIP) $(GTHDISTDIR)/bin/gthsplit
	@cp $(CURDIR)/scripts/gthclean.sh $(GTHDISTDIR)/bin
	@cp $(CURDIR)/scripts/gthcleanrec.sh $(GTHDISTDIR)/bin
	@cp $(CURDIR)/scripts/gthsplit2dim.sh $(GTHDISTDIR)/bin
endif
	@cp $(CURDIR)/gthdata/bio/genomic_sequence_single_exon.fas $(GTHDISTDIR)/data
	@cp $(CURDIR)/gthdata/bio/cdna_sequence_single_exon.fas $(GTHDISTDIR)/data
	@cp bin/bssm/*.bssm $(GTHDISTDIR)/bin/bssm
	@cp $(CURDIR)/gthdata/gth/BLOSUM62 $(GTHDISTDIR)/bin/gthdata
	@cp $(CURDIR)/gthdata/gth/TransDNAX $(GTHDISTDIR)/bin/gthdata
	@cp $(CURDIR)/gthdata/gth/TransProt11 $(GTHDISTDIR)/bin/gthdata
	@cp $(CURDIR)/gthdata/gth/polyatail.dna $(GTHDISTDIR)/bin/gthdata
	@cp $(CURDIR)/doc/gthmanual/gthmanual.pdf $(GTHDISTDIR)/doc
ifeq ($(SYSTEM),Windows)
	@cd $(DISTDIR) && 7z a -tzip $(GTHDISTBASENAME).zip $(GTHDISTBASENAME)
	@echo "$(DISTDIR)/$(GTHDISTBASENAME).zip"
else
	@cd $(DISTDIR) && tar --numeric-owner -cf $(GTHDISTBASENAME).tar  $(GTHDISTBASENAME)
	@cd $(DISTDIR) && gzip -f -9 $(GTHDISTBASENAME).tar
	@echo "$(DISTDIR)/$(GTHDISTBASENAME).tar.gz"
endif

libdist: all
	@echo "[build gth library distribution]"
	@rm -rf $(GTHLIBDISTDIR)
	@rm -f $(DISTDIR)/$(GTHLIBDISTBASENAME).tar.gz
	@mkdir -p $(GTHLIBDISTDIR)
	@mkdir -p $(GTHLIBDISTDIR)/lib
	@mkdir -p $(GTHLIBDISTDIR)/obj
	@cp $(CURDIR)/gthdata/README.libdist $(GTHLIBDISTDIR)/README
	@cp $(CURDIR)/gthdata/Makefile $(GTHLIBDISTDIR)
	@cp lib/libgenomethreader.a $(GTHLIBDISTDIR)/lib
	@cp $(GTHLIBS) $(GTHLIBDISTDIR)/lib
	@cp obj/src/gth.o $(GTHLIBDISTDIR)/obj
	@cp obj/src/gthconsensus.o $(GTHLIBDISTDIR)/obj
	@cd $(DISTDIR) && tar --owner=root -cf $(GTHLIBDISTBASENAME).tar $(GTHLIBDISTBASENAME)
	@cd $(DISTDIR) && gzip -f -9 $(GTHLIBDISTBASENAME).tar
	@echo "$(DISTDIR)/$(GTHLIBDISTBASENAME).tar.gz"

srcdist:
	git archive --format=tar --prefix=genomethreader-`cat VERSION`/ HEAD | \
        gzip -9 > genomethreader-`cat VERSION`.tar.gz

release:
	git tag "v`cat VERSION`"
	git push --tags origin master

installwww:
# install genomethreader.org website
	rsync -rv --delete www/genomethreader.org/htdocs/ $(SERVER):$(WWWBASEDIR)/genomethreader.org

push:
	git push origin master

cflags:
	@echo ${GTH_CFLAGS}

RUBY:=ruby

test: all
	cd testsuite && env -i GTH_MEM_BOOKKEEPING=on MAKE=$(MAKE) PATH="$(PATH)" \
	  BSSMDIR=$(CURDIR)/bin/bssm GTHDATADIR=$(CURDIR)/gthdata/gth \
          CCACHE_DISABLE=yes HOME=$(HOME) GTHTESTDATADIR=$(GTHTESTDATADIR) \
          $(RUBY) -I. testsuite.rb \
          -testdata $(CURDIR)/testdata -bin $(CURDIR)/bin -cur $(CURDIR) \
          -gtruby $(CURDIR)/gtruby $(STEST_FLAGS)

clean:
	rm -rf obj lib
	rm -rf testsuite/stest_testsuite testsuite/stest_stest_tests

cleanup: clean
	rm -rf bin
	find gthdata -name '*.des' | xargs rm -f
	find gthdata -name '*.esq' | xargs rm -f
	find gthdata -name '*.md5' | xargs rm -f
	find gthdata -name '*.ois' | xargs rm -f
	find gthdata -name '*.sds' | xargs rm -f

obj/train: bin/gthbssmbuild
	@mkdir -p bin/bssm
	@echo "[train arabidopsis.bssm]"
	@bin/gthbssmbuild -bssmfile bin/bssm/arabidopsis.bssm \
                          -datapath gthdata/bssm/arabidopsis \
                          -gtdonor -gcdonor -agacceptor -gzip
	@echo "[train rice.bssm]"
	@bin/gthbssmbuild -bssmfile bin/bssm/rice.bssm \
                          -datapath gthdata/bssm/rice \
                          -gtdonor -gcdonor -agacceptor -gzip
	@echo "[train maize.bssm]"
	@bin/gthbssmbuild -bssmfile bin/bssm/maize.bssm \
                          -datapath gthdata/bssm/maize \
                          -gtdonor -gcdonor -agacceptor -gzip
	@echo "[train medicago.bssm]"
	@bin/gthbssmbuild -bssmfile bin/bssm/medicago.bssm \
                          -datapath gthdata/bssm/medicago \
                          -gtdonor -agacceptor -gzip
	@touch $@

obj/old_train: bin/gthmkbssmfiles obj/train
	@echo "[write old BSSMs]"
	@bin/gthmkbssmfiles bin/bssm
	@touch $@

obj/gthdata:
	@echo "[copy gthdata stuff to bin]"
	@mkdir -p bin/gthdata
	@cp $(CURDIR)/gthdata/gth/BLOSUM62 bin/gthdata
	@cp $(CURDIR)/gthdata/gth/TransDNAX bin/gthdata
	@cp $(CURDIR)/gthdata/gth/TransProt11 bin/gthdata
	@cp $(CURDIR)/gthdata/gth/polyatail.dna bin/gthdata
	@touch $@

train: obj/train obj/old_train obj/gthdata
