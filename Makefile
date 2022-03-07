PREFIX = .
CCFLAGS = -std=c++11 -O3 -Wall


ifeq ($(OS),Windows_NT)
	CCFLAGS += -static -D AMD64
else
	UNAME_S := $(shell uname -s)
	ifeq ($(UNAME_S),Linux)
		CCFLAGS += -static-libstdc++ -D LINUX
	endif
	ifeq ($(UNAME_S),Darwin)
		CCFLAGS += -D OSX
	endif
endif


all:
	mkdir -p $(PREFIX)/bin
	mkdir -p $(PREFIX)/util

	$(CXX) src/data_norm.cpp -o $(PREFIX)/bin/data_norm $(CCFLAGS)
	$(CXX) src/data_filter.cpp -o $(PREFIX)/bin/data_filter $(CCFLAGS)
	$(CXX) src/network_build.cpp -o $(PREFIX)/bin/network_build -pthread $(CCFLAGS)
	$(CXX) src/module_identify.cpp -o $(PREFIX)/bin/module_identify -pthread $(CCFLAGS)
	$(CXX) src/annotate.cpp -o $(PREFIX)/bin/annotate -pthread $(CCFLAGS)
	$(CXX) src/rwr.cpp -o $(PREFIX)/bin/rwr $(CCFLAGS)

	$(CXX) src/data_stat.cpp -o $(PREFIX)/util/data_stat $(CCFLAGS)
	$(CXX) src/network_stat.cpp -o $(PREFIX)/util/network_stat $(CCFLAGS)
	$(CXX) src/network_merge.cpp -o $(PREFIX)/util/network_merge $(CCFLAGS)
	$(CXX) src/network_extract.cpp -o $(PREFIX)/util/network_extract $(CCFLAGS)
	$(CXX) src/network_shuffle.cpp -o $(PREFIX)/util/network_shuffle $(CCFLAGS)
	$(CXX) src/enrich.cpp -o $(PREFIX)/util/enrich $(CCFLAGS)
	$(CXX) src/calculate_accuracy.cpp -o $(PREFIX)/util/calculate_accuracy -pthread $(CCFLAGS)
	$(CXX) src/generate_expr_matrix_from_rsem.cpp -o $(PREFIX)/util/generate_expr_matrix_from_rsem $(CCFLAGS)
	$(CXX) src/generate_expr_matrix_from_stringtie.cpp -o $(PREFIX)/util/generate_expr_matrix_from_stringtie $(CCFLAGS)
	$(CXX) src/csv_to_tsv.cpp -o $(PREFIX)/util/csv_to_tsv $(CCFLAGS)
	$(CXX) src/tsv_to_csv.cpp -o $(PREFIX)/util/tsv_to_csv $(CCFLAGS)
