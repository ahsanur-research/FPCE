CC=gcc -std=gnu99 -Wall
CFLAGS=-O3
DEBUG=-g
PROF=#-pg

denseSG: ODES_master.o ODES_slave.o util.o odes.h denseSG.c
	$(CC) $(CFLAGS) $(DEBUG) $(PROF) denseSG.c ODES_master.o \
	-o denseSG  ODES_slave.o util.o -lpthread
        
ODES_master.o: ODES_master.c odes.h
	$(CC) -c $(CFLAGS) $(DEBUG) $(PROF) ODES_master.c
        
ODES_slave.o: ODES_slave.c odes.h
	$(CC) -c           $(DEBUG) $(PROF) ODES_slave.c
        
util.o: util.c odes.h
	$(CC) -c $(CFLAGS) $(DEBUG) $(PROF) util.c
        
clean: 
	rm *o

distclean:
	if [ -e denseSG ]; then rm denseSG; fi
	rm *o


# A breadth-first algorithm to find dense subgraphs with density >= 1/2
# Copyright (C) 2010 James Long
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
