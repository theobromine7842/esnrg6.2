# Makefile for making ESNRG program 
#
ARCH   = `uname -s`
FOR    = gfortran 
FFLAGS = -O4 -align records -align dcommons \
         -extend_source -math_library accurate \
         -tune host -arch host -DARCH=$(ARCH) -u 
PROGRAM = esnrg6.2.exe
VERSION = 6.2

.SUFFIXES:
.SUFFIXES: .o .f .s .c 

#
SRCS= esnrg6.2.f numread.f assigngrp.f translate_newpbc.f id.f props.f \
      vdw.f nbond.f topread.f topread2.f born.f born2.f cs_interface.f \
      find_chg.f find_solv.f find_coi.f chgsolv.f chk_seg.f

OBJS= esnrg6.2.o numread.o assigngrp.o translate_newpbc.o id.o props.o \
      vdw.o nbond.o topread.o topread2.o born.o born2.o cs_interface.o \
      find_chg.o find_solv.o find_coi.o chgsolv.o chk_seg.o

INCS= dim.inc cmn.inc

.f.o:
	$(FOR) -c  $*.f 

esnrg:	$(OBJS) $(OBJP)
	$(FOR)  $(OBJS) $(OBJP) -o $(PROGRAM) 

dist:	$(SRCS) Makefile
	tar -cvzf $(PROGRAM)-$(VERSION).tar.gz Makefile $(INCS) $(SRCS) $(SRCPS)

clean:
	rm $(OBJS) $(OBJP) $(PROGRAM)

#
#
