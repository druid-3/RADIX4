################################################################################
################################################################################
##  developed by DRUID's labs.                                                ##
##  druid3@i.ua or odliudko@mail.ua                                           ##
##  http://dsp.la.net.ua                                                      ##
##                                                                            ##
##                 00  000   0  0  0  000   0   000                           ##
##                000  0  0  0  0  0  0  0     00        *******              ##
##               0100  000   0  0  0  0  0       00     * ****  *             ##
##              01100  0  0   000  0  000      000     *  *   *  *            ##
##             011 00                                  *  ****   *            ##
##            011  00     11   1      111               * *   * *             ##
##           00000000    1  1  111   11                  *******              ##
##          00000000011  1111  1  1    11                                     ##
##          11111111111  1  1  111   111                                      ##
##                                                                            ##
################################################################################
################################################################################
# Clean automated Makefile 4 C&C++ (prg) 
# distributed under GNU GPL v.3
WHERE_IS_GCC  = /usr/bin/

#INCLUDE=-I../src/include -I/usr/src/linux/include/

CC            = gcc
CPP           = g++
STRIP         = $(WHERE_IS_GCC)strip

CC_FLAGS      = -O0                                                            \
                -W                                                             \
                -Wall                                                          \
                -Winline                                                       \
                -Wmissing-noreturn                                             \
                -Wmissing-prototypes                                           \
                -Wnested-externs                                               \
                -Wconversion                                                   \
                -Wfloat-equal                                                  \
                -Wshadow                                                       \
                -Wsequence-point                                               \
                -pipe 

# <add, if you are developing a kernel module>
#KERNEL        = -DMODULE -D__KERNEL__ -DLINUX
#INCLUDE=-I/usr/include/linux/

LD_FLAGS      = -L./LIB/ $(ALIB) $(DLIB)
DLIB          = -lm
ALIB          =


PRG_NAME      = FFTkit_test

PRJDIRS       = ./OBJ                                                          \
                ./ASM                                                          \
                ./LIB

OBJFILES     := $(patsubst %.c,   ./OBJ/%.o,       $(wildcard *.c))
OBJFILES     += $(patsubst %.cpp, ./OBJ/%.o,       $(wildcard *.cpp))
OBJFILES     += $(patsubst %.s,   ./OBJ/%.o,       $(wildcard *.s))

C2ASM        := $(patsubst %.c,   ./ASM/%.S,       $(wildcard *.c))
CPP2ASM      := $(patsubst %.cpp, ./ASM/%.S,       $(wildcard *.cpp))


#  include <bla-bla-bla> error and stop if *.mk not available
# -include <bla-bla-bla> work if *.mk not available
-include ./options/options.mk

all: prj_dirs make_lib myprog show_src run_script

#PRG----------------------------------------------------------------------------
myprog: $(OBJFILES) 
	$(CC) -o $(PRG_NAME) $(OBJFILES) $(LD_FLAGS)

./OBJ/%.o: %.c
	$(CC) $(CC_FLAGS) -c -o $@ $<
./OBJ/%.o: %.cpp
	$(CPP) $(CC_FLAGS) -c -o $@ $<
./OBJ/%.o: %.s  
	$(CC) -o $@ -c $<
#-------------------------------------------------------------------------------

run: all
	./run

src2asm: c2asm cpp2asm

c2asm: $(C2ASM)

./ASM/%.S: %.c
	$(CC) -o $@ -S -fverbose-asm $<
	
cpp2asm: $(CPP2ASM)
	
./ASM/%.S: %.cpp
	$(CPP) -o $@ -S -fverbose-asm $<


prj_dirs:
	@$(foreach DIR, $(PRJDIRS),                                                    \
	if ! test -d $(DIR);                                                           \
		then                                                                   \
			mkdir $(DIR);                                                  \
			chmod 775 $(DIR);                                              \
		else                                                                   \
			echo ' $(DIR) directory already is available ';                \
	fi;                                                                            \
	)

run_script:
	@if ! test -f ./run;                                                           \
	then                                                                           \
		touch ./run;                                                           \
		chmod 775 ./run;                                                       \
		echo '#!/bin/bash'                                   >> ./run;         \
		echo 'echo "starting ..."'                           >> ./run;         \
		echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:./LIB' >> ./run;         \
		echo 'exec ./$(PRG_NAME) $$@'                        >> ./run;         \
	else                                                                           \
		echo ' run-file exists ';                                              \
	fi

#-------------------------------------------------------------------------------
# * !!! <attention> !!!
# <ALL the source code from this directory join to the project.>
# * !!! <attention> !!!
show_src:
	@echo "-------------------------------------------------------------------------------"
	@echo "src in project:"
	@echo
	@for i in `ls *.[c,h,s] *.[c,h]?*`;                                                     \
		do                                                                              \
			echo $${i};                                                             \
		done
	@echo "-------------------------------------------------------------------------------"
#-------------------------------------------------------------------------------
count:
	@wc -l *.[c,h,s] *.[c,h]?*

archive: clean
	@tar -zcvf ../$(PRG_NAME)__`date +%y.%m.%d`.tgz ./ --recursion -P

clean:
	@rm -f $(PRG_NAME) run
	@$(foreach DIR, ./ $(PRJDIRS), rm -f                                       \
	                                     $(DIR)/*.o                            \
	                                     $(DIR)/*.d                            \
	                                     $(DIR)/*.map                          \
	                                     $(DIR)/*.hex                          \
	                                     $(DIR)/*.listing                      \
	                                     $(DIR)/*.bak                          \
	                                     $(DIR)/*.out                          \
	                                     $(DIR)/*.swp                          \
	                                     $(DIR)/*~                             \
	                                     $(DIR)/core                           \
	                                     $(DIR)/*.tu                           \
	  )

striping:
	@$(STRIP) ./$(PRG_NAME)

doc++:
	@doc++ --all --html --tables-border --private --dir ./DOC *.[c,h,s] *.[c,h]?*
	
gnu_style:
	@astyle --style=gnu `ls *.[c,h]?*`
	@rm -f *.orig

#make_lib:
#	$(MAKE) -C ./LIB/
	
.PHONY: all
.PHONY: make_lib
.PHONY: myprog
.PHONY: run
.PHONY: src2asm
.PHONY: prj_dirs
.PHONY: run_script
.PHONY: show_src
.PHONY: count
.PHONY: archive
.PHONY: clean
.PHONY: striping
.PHONY: doc++
