subdirs := eval finder glearn glearn2 gspan
.PHONY: all $(subdirs)
all: $(subdirs)
$(subdirs):
	$(MAKE) -C $@
.PHONY: clean
clean:
	@for dir in $(subdirs); do (cd $$dir; $(MAKE) clean); done
	rm *~ *.out *.mat
