CLEANFILES = *.d *.so *_rdict.* *_rdict_* *.rootmap
MAINTAINERCLEANFILES = Makefile.in

# EXTRA_DIST = some extra files, e.g. test-data

lib_LTLIBRARIES = librad-spec-tools.la

librad_spec_tools_la_SOURCES = \
	Binning.cxx \
	HistAnalysis.cxx \
	SDPeak.cxx

librad_spec_tools_la_headers = \
	Binning.h \
	HistAnalysis.h \
	SDPeak.h

pkginclude_HEADERS = $(librad_spec_tools_la_headers)

# bin_PROGRAMS = rad-spec-tools

# noinst_PROGRAMS =

# rad_spec_tools_SOURCES = rad-spec-tools.cxx
# rad_spec_tools_LDADD = librad-spec-tools.la
# rad_spec_tools_LDFLAGS = -static


rootmapdir = $(libdir)
rootmap_DATA =


librad_spec_tools_la_SOURCES += librad-spec-tools_rdict.cxx
rootmap_DATA += .libs/librad-spec-tools.rootmap

if WITH_CLING

rootmap_DATA += .libs/librad@ROOTPCMDASH@spec@ROOTPCMDASH@tools_rdict.pcm

librad-spec-tools_rdict.cxx: $(librad_spec_tools_la_headers) rad-spec-tools_LinkDef.h
	$(ROOTCLING) -f $@.tmp -s librad-spec-tools@SHLIBEXT@ -rml librad-spec-tools@SHLIBEXT@ -rmf librad-spec-tools.rootmap.tmp -c $(CPPFLAGS) $(CXXFLAGS) -I$(includedir) $+
	@# Some magic to prefix header names with "$(PACKAGE)/", and only that, in dictionary and rootmap:
	$(GREP) -F -v '"'"`pwd`"'/",' $@.tmp | $(SED) 's|"\([^"]*/\)\?\([^/"]*[.]h\)",|"'$(PACKAGE)/'\2",| ; s|\\"\([^"]*/\)\?\([^/"]*[.]h\)\\"\\n"|\\"'$(PACKAGE)/'\2\\"\\n"|' > $@.tmp2
	$(SED) 's|\$$clingAutoload\$$\([^/""]\+/\)*|$$clingAutoload$$'$(PACKAGE)'/|; /.*DICTPAYLOAD(.*/,/.*)DICTPAYLOAD.*/ s|#include "\([^/"]\+/\)*\(.*\)"|#include <'$(PACKAGE)'/\2>|' $@.tmp2 > $@ && $(RM) $@.tmp $@.tmp2
	$(SED) 's|\(header \+\)\([^ ].*/\)\?\([^ ].*[.]h\)|\1'$(PACKAGE)/'\3|' librad-spec-tools.rootmap.tmp > librad-spec-tools.rootmap && $(RM) librad-spec-tools.rootmap.tmp

.libs/librad-spec-tools.rootmap: librad-spec-tools.la
	cp librad-spec-tools.rootmap librad@ROOTPCMDASH@spec@ROOTPCMDASH@tools_rdict.pcm .libs/

.libs/librad@ROOTPCMDASH@spec@ROOTPCMDASH@tools_rdict.pcm: .libs/librad-spec-tools.rootmap

else # not WITH_CLING

librad-spec-tools_rdict.cxx: $(librad_spec_tools_la_headers) rad-spec-tools_LinkDef.h
	$(ROOTCINT) -f $@ -c -p $(CPPFLAGS) `echo $(CXXFLAGS) | sed 's/\(^\|\s\)\(-std=[A-Za-z0-9+]\+\)\($$\|\s\)/\1\3/g'` $+

.libs/librad-spec-tools.rootmap: librad-spec-tools.la | $(librad_spec_tools_la_headers) rad-spec-tools_LinkDef.h
	$(RLIBMAP) -o $@ -l .libs/librad-spec-tools@SHLIBEXT@ -c $|

endif # WITH_CLING
