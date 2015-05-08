SHELL = /bin/sh
.SUFFIXES: .c .o

# Notes:
# When compiling in Windows, change the "kc.so" below to "kc.dll". The "-fPIC" can be optionally omitted
# When compiling in Mac OS X, change the "kc.so" below to "kc.dylib", and change "-shared" to "-dynamiclib"
# When compiling in Solaris, change "gcc" to "cc"



CC = gcc
LD = gcc
CFLAGS = `perl -MExtUtils::Embed -e ccopts` -fPIC
LDFLAGS = `perl -MExtUtils::Embed -e ldopts`
ARCH = `perl -MConfig -e 'print $$Config{archname}'`
VER = `perl -MConfig -e 'print $$Config{version}'`

kc.so: kc_wrap.o kc.o 
	$(LD) -shared -o $@ $^ $(LDFLAGS)
	mkdir -p $(VER)
	mkdir -p $(VER)/$(ARCH)/
	mkdir -p $(VER)/$(ARCH)/auto/
	mv $@ $(VER)/$(ARCH)/

clean:
	rm -f *.o *.so
