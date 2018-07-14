PBCLIB=pbc/.libs/libpbc.a

all: test

test: bilinear_acc example
	./bilinear_acc < pbc/param/a.param
	./example < pbc/param/a.param

bilinear_acc: bilinear_acc.c $(PBCLIB)
	$(CC) -o $@ $< -Ipbc/include -Ipbc $(PBCLIB) -lgmp

example: example.modified.c $(PBCLIB)
	$(CC) -o $@ $< -Ipbc/include -Ipbc $(PBCLIB) -lgmp

$(PBCLIB): pbc/configure
	make -C pbc

pbc/configure:
	cd pbc && autoreconf --install && ./configure
