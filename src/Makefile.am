CLEANFILES = *.d *.so *-rdict.cxx *-rdict.h
MAINTAINERCLEANFILES = Makefile.in

# EXTRA_DIST = some extra files, e.g. test-data

lib_LTLIBRARIES = librad-spec-tools.la

librad_spec_tools_la_SOURCES =

librad_spec_tools_la_headers =

pkginclude_HEADERS = $(librad_spec_tools_la_headers)

# bin_PROGRAMS = rad-spec-tools

# noinst_PROGRAMS =

# rad_spec_tools_SOURCES = rad-spec-tools.cxx
# rad_spec_tools_LDADD = librad-spec-tools.la
# rad_spec_tools_LDFLAGS = -static


librad_spec_tools_la_SOURCES += librad-spec-tools-rdict.cxx

librad-spec-tools-rdict.cxx: $(librad_spec_tools_la_headers) rad-spec-tools-LinkDef.h
	$(ROOTCINT) -f $@ -c -p $(CPPFLAGS) $(CXXFLAGS) $+

rootmapdir = $(libdir)
rootmap_DATA = .libs/librad-spec-tools.rootmap

.libs/librad-spec-tools.rootmap: librad-spec-tools.la | $(librad_spec_tools_la_headers) rad-spec-tools-LinkDef.h
	$(RLIBMAP) -o $@ -l .libs/librad-spec-tools@SHLIBEXT@ -c $|
