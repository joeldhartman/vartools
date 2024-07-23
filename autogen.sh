#! /bin/bash
case `uname` in Darwin*) glibtoolize --force --copy --ltdl --recursive ;;
  *) libtoolize --force --copy --ltdl --recursive ;; esac
aclocal --force
autoheader --force
echo 'example_files = \' > examplefiles.mk
find EXAMPLES -type f -print | grep -v '\.svn' | sed 's/^/ /;$q;s/$/ \\/' >> examplefiles.mk
automake --force --copy --add-missing
autoconf --force
cp Makefile.in Makefile.in.bak
< Makefile.in.bak \
 sed -e 's/ACLOCAL = @ACLOCAL@/ACLOCAL = @ACLOCAL@ || echo/' \
     -e 's/AUTOCONF = @AUTOCONF@/AUTOCONF = @AUTOCONF@ || echo/' \
     -e 's/AUTOHEADER = @AUTOHEADER@/AUTOHEADER = @AUTOHEADER@ || echo/' \
     -e 's/AUTOMAKE = @AUTOMAKE@/AUTOMAKE = @AUTOMAKE@ || echo/' \
> Makefile.in
cd libltdl
cp Makefile.in Makefile.in.bak
< Makefile.in.bak \
 sed -e 's/ACLOCAL = @ACLOCAL@/ACLOCAL = @ACLOCAL@ || echo/' \
     -e 's/AUTOCONF = @AUTOCONF@/AUTOCONF = @AUTOCONF@ || echo/' \
     -e 's/AUTOHEADER = @AUTOHEADER@/AUTOHEADER = @AUTOHEADER@ || echo/' \
     -e 's/AUTOMAKE = @AUTOMAKE@/AUTOMAKE = @AUTOMAKE@ || echo/' \
> Makefile.in
cd ..
