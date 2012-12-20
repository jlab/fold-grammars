VIENNA_SOURCEDIR = /home/sjanssen/Desktop/Turner2004/ViennaRNA-2.0.5/
VIENNA_TAREGTDIR = Vienna/
FOLDGRAMMARSDIR = /home/sjanssen/Desktop/fold-grammars/
GAPC_BASEDIR = /stefan/

INCLUDES = -I $(GAPC_BASEDIR)/include/rtlib -I $(GAPC_BASEDIR)/include/librna -I $(GAPC_BASEDIR)/include -I ../
LIBDIRS = -Xlinker -rpath -Xlinker $(GAPC_BASEDIR)lib/ -L$(GAPC_BASEDIR)lib/

CPPFLAGS = -g3 -O0 -Wall -Wnon-virtual-dtor -Wno-parentheses -fmessage-length=0	#for debugging
#~ CPPFLAGS = -g -O3 -DNDEBUG -Wall -Wnon-virtual-dtor -Wno-parentheses -fmessage-length=0	#for production

all: copyVienna gapc compile
	
gapc:
	cd $(FOLDGRAMMARSDIR) && gapc -p "(alg_shape5 * (alg_mfe * alg_dotBracket)) suchthat range_shape_mfe_db" -t -o shape5mfedb.cc nodangle.gap
	cat $(FOLDGRAMMARSDIR)shape5mfedb.hh | sed "s/namespace gapc/namespace ns_shape5mfedb/" > ./shape5mfedb.hh
	#~ cp $(FOLDGRAMMARSDIR)shape5mfedb.hh .
	cp $(FOLDGRAMMARSDIR)shape5mfedb.cc .
	
	cd $(FOLDGRAMMARSDIR) && gapc -p "(alg_shape4 * (alg_mfe * alg_dotBracket)) suchthat range_shape_mfe_db" -t -o shape4mfedb.cc nodangle.gap
	cat $(FOLDGRAMMARSDIR)shape4mfedb.hh | sed "s/namespace gapc/namespace ns_shape4mfedb/" > ./shape4mfedb.hh
	#~ cp $(FOLDGRAMMARSDIR)shape4mfedb.hh .
	cp $(FOLDGRAMMARSDIR)shape4mfedb.cc .
	
	cd $(FOLDGRAMMARSDIR) && gapc -p "(alg_shape3 * (alg_mfe * alg_dotBracket)) suchthat range_shape_mfe_db" -t -o shape3mfedb.cc nodangle.gap
	cat $(FOLDGRAMMARSDIR)shape3mfedb.hh | sed "s/namespace gapc/namespace ns_shape3mfedb/" > ./shape3mfedb.hh
	#~ cp $(FOLDGRAMMARSDIR)shape3mfedb.hh .
	cp $(FOLDGRAMMARSDIR)shape3mfedb.cc .
	
	cd $(FOLDGRAMMARSDIR) && gapc -p "(alg_shape2 * (alg_mfe * alg_dotBracket)) suchthat range_shape_mfe_db" -t -o shape2mfedb.cc nodangle.gap
	cat $(FOLDGRAMMARSDIR)shape2mfedb.hh | sed "s/namespace gapc/namespace ns_shape2mfedb/" > ./shape2mfedb.hh
	#~ cp $(FOLDGRAMMARSDIR)shape2mfedb.hh .
	cp $(FOLDGRAMMARSDIR)shape2mfedb.cc .
	
	cd $(FOLDGRAMMARSDIR) && gapc -p "(alg_shape1 * (alg_mfe * alg_dotBracket)) suchthat range_shape_mfe_db" -t -o shape1mfedb.cc nodangle.gap
	cat $(FOLDGRAMMARSDIR)shape1mfedb.hh | sed "s/namespace gapc/namespace ns_shape1mfedb/" > ./shape1mfedb.hh
	#~ cp $(FOLDGRAMMARSDIR)shape1mfedb.hh .
	cp $(FOLDGRAMMARSDIR)shape1mfedb.cc .
	
	cp $(FOLDGRAMMARSDIR)pfunc_filter_foldrna.hh .
	cp $(FOLDGRAMMARSDIR)mfesubopt.hh .
	
compile:
	gengetopt < RNAshapes.ggo
	g++ -c RNAshapes_cmdl.c
	gcc -I $(VIENNA_TAREGTDIR)H/ -c $(VIENNA_TAREGTDIR)lib/utils.c -o utils.o
	g++ $(INCLUDES) $(CPPFLAGS) -c $(GAPC_BASEDIR)include/rtlib/string.cc -o string.o
	g++ $(INCLUDES) $(CPPFLAGS) -MMD -MP -o shape5mfedb.o -c shape5mfedb.cc
	g++ $(INCLUDES) $(CPPFLAGS) -MMD -MP -o shape4mfedb.o -c shape4mfedb.cc
	g++ $(INCLUDES) $(CPPFLAGS) -MMD -MP -o shape3mfedb.o -c shape3mfedb.cc
	g++ $(INCLUDES) $(CPPFLAGS) -MMD -MP -o shape2mfedb.o -c shape2mfedb.cc
	g++ $(INCLUDES) $(CPPFLAGS) -MMD -MP -o shape1mfedb.o -c shape1mfedb.cc
	g++ -I $(VIENNA_TAREGTDIR)H/ $(INCLUDES) $(CPPFLAGS) -MMD -MP -MF"RNAshapes.d" -MT"RNAshapes.d" -o RNAshapes.o -c RNAshapes.cc
	g++ -o RNAshapes RNAshapes_cmdl.o RNAshapes.o string.o shape5mfedb.o shape4mfedb.o shape3mfedb.o shape2mfedb.o shape1mfedb.o utils.o $(LIBDIRS) -lrna
	

	
copyVienna:
	test -d $(VIENNA_TAREGTDIR) || mkdir $(VIENNA_TAREGTDIR)
	cp $(VIENNA_SOURCEDIR)config.h $(VIENNA_TAREGTDIR)
	
	test -d $(VIENNA_TAREGTDIR)H/ || mkdir $(VIENNA_TAREGTDIR)H/
	cp $(VIENNA_SOURCEDIR)H/utils.h $(VIENNA_TAREGTDIR)H/
	
	test -d $(VIENNA_TAREGTDIR)lib/ || mkdir $(VIENNA_TAREGTDIR)lib/
	cp $(VIENNA_SOURCEDIR)lib/utils.c $(VIENNA_TAREGTDIR)lib/

	
clean:
	#~ rm -rf $(HDIR) $(CDIR)
	#~ rm -f config.h
	rm -f RNAshapes_cmdl*
	rm -f RNAshapes.o
	rm -f RNAshapes
	