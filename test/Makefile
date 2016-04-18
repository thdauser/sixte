# -*- mode: Make -*-

scripts := $(wildcard *.csh)
 .PHONY: $(scripts)	


test:
	@echo " "
	@echo " *************************** "
	@echo "   testing the sixte tools   "
	@echo " *************************** "
	@echo " "
	make $(scripts)	
	make clean

$(scripts): 
	./$(basename $@).csh > log/$(basename $@).log
	@echo " "


make clean:
	rm -f output/*.fits
