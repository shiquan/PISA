$(ZLIBDIR)/libz.a:
	+cd $(ZLIBDIR) && ./configure && $(MAKE) static
