all:
	./setup.sh

install:
	mkdir -p $(PREFIX)
	cp -rp ../westpa $(PREFIX)
	#bash $(PREFIX)/westpa/westpa.sh
