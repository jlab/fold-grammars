PREFIX=$(HOME)/local/

MAKE=make
PERL=perl
BASEDIR=../../../

install:
	install -d $(PREFIX)/lib
	install -d $(PREFIX)/lib/foldGrammars
	for file in foldGrammars/*; do \
		install -m 644 $$file $(PREFIX)/lib/foldGrammars; \
	done;

clean:
	echo "Nothing to do"