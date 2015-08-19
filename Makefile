SUBDIRS = utility filetypes matrices blocked MMF 

.PHONY: all objects tests clean $(SUBDIRS)

objects:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir objects; \
	done

tests:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir tests; \
	done

clean:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean; \
	done

all: objects tests

# declare dependencies among subdirectories
filetypes: utility
matrices: filetypes
blocked: matrices
MMF: blocked

tests: objects

anew: clean all 
