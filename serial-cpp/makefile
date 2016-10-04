# USE CLANG COMPILER FOR MAC OS X MAVERICKS / YOSEMITE

CC=clang++

# DEFINE WARN FLAG FOR VERBOSE COMPILER WARNINGS

#WARN=-Weverything -Wno-global-constructors -Wno-documentation -Wno-missing-prototypes -Wno-c++98-compat -Wno-missing-noreturn -Wno-exit-time-destructors -Wno-tautological-compare -ferror-limit=1

# USE -DNDEBUG FLAG TO TURN ON DEBUGGING WITH ASSERT

CFLAGS=-O2 -stdlib=libc++ -std=c++11
#CFLAGS=-O3 -DNDEBUG -stdlib=libc++ -std=c++11

# BUILD CODE FOR NORMAL (DEFAULT) OR STRANGE (DEFINE STRANGE MACRO) CORRELATIONS

aklt-normal: aklt.cpp
	$(CC) -o $@ aklt.cpp $(CFLAGS) $(WARN)

aklt-strange: aklt.cpp
	$(CC) -o $@ aklt.cpp $(CFLAGS) $(WARN) -DSTRANGE

clean:
	@rm aklt aklt-full

