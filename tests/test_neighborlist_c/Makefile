#
# CDDL HEADER START
#
# The contents of this file are subject to the terms of the Common Development
# and Distribution License Version 1.0 (the "License").
#
# You can obtain a copy of the license at
# http://www.opensource.org/licenses/CDDL-1.0.  See the License for the
# specific language governing permissions and limitations under the License.
#
# When distributing Covered Code, include this CDDL HEADER in each file and
# include the License file in a prominent location with the name LICENSE.CDDL.
# If applicable, add the following below this CDDL HEADER, with the fields
# enclosed by brackets "[]" replaced with your own identifying information:
#
# Portions Copyright (c) [yyyy] [name of copyright owner]. All rights reserved.
#
# CDDL HEADER END
#

#
# Copyright (c) 2012, Regents of the University of Minnesota.  All rights reserved.
#
# Contributors:
#    Ryan S. Elliott
#

#
# Release: This file is part of the openkim-api.git repository.
#


TEST_NAME := test_neighborlist_c 
FILES = ../../neighborlist.c ../../cvec.c test_neighborlist_c.c

ifndef KIM_API_DIR
   include $(KIM_DIR)KIM_API/Include.mk
else
   include $(KIM_API_DIR)Include.mk
endif


LOCALOBJ =

all: $(TEST_NAME)
#$(TEST_NAME) $(KIM_LIB_FILE)
$(TEST_NAME): 
	$(LINKCOMPILER) $(CPPFLAG) $(CPPLIBFLAG) $(FILES)  $< $(KIM_LIB) $(OUTPUTIN) $@

#$(TEST_NAME).o: $(LOCALOBJ) Makefile

Makefile: $(KIM_API_DIR)Include.mk
	@touch Makefile

clean:
	rm -f *.mod *.o $(TEST_NAME)
