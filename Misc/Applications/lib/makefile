PREFIX=$(HOME)/local/

MAKE=make
PERL=perl
BASEDIR=../../../

install:
	if [ ! -d "$(PREFIX)/lib" ]; then install -d $(PREFIX)/lib; fi;
	install -d $(PREFIX)/lib/foldGrammars
	for file in foldGrammars/*; do \
		install -m 644 $$file $(PREFIX)/lib/foldGrammars; \
	done;

clean:
	echo "Nothing to do"