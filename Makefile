CXX = g++
LFLAGS = -std=c++17 -static-libstdc++ -O3 -Wall
MFLAGS = -std=c++17 -static-libgcc -O3 -Wall
WFLAGS = -std=c++17 -static -O3 -Wall


linux: 
	mkdir -p bin
	mkdir -p util

	$(CXX) src/data_norm.cpp -o bin/data_norm $(LFLAGS)
	$(CXX) src/data_filter.cpp -o bin/data_filter $(LFLAGS)
	$(CXX) src/network_build.cpp -o bin/network_build -pthread $(LFLAGS)
	$(CXX) src/module_identify.cpp -o bin/module_identify -pthread $(LFLAGS)
	$(CXX) src/annotate.cpp -o bin/annotate -pthread -lstdc++fs $(LFLAGS)
	$(CXX) src/rwr.cpp -o bin/rwr $(LFLAGS)

	$(CXX) src/data_stat.cpp -o util/data_stat $(LFLAGS)
	$(CXX) src/network_merge.cpp -o util/network_merge $(LFLAGS)
	$(CXX) src/enrich.cpp -o util/enrich $(LFLAGS)
	$(CXX) src/generate_expr_matrix_from_rsem.cpp -o util/generate_expr_matrix_from_rsem $(LFLAGS)
	$(CXX) src/generate_expr_matrix_from_stringtie.cpp -o util/generate_expr_matrix_from_stringtie $(LFLAGS)


macos: 
	mkdir -p bin
	mkdir -p util

	$(CXX) src/data_norm.cpp -o bin/data_norm $(MFLAGS)
	$(CXX) src/data_filter.cpp -o bin/data_filter $(MFLAGS)
	$(CXX) src/network_build.cpp -o bin/network_build -pthread $(MFLAGS)
	$(CXX) src/module_identify.cpp -o bin/module_identify -pthread $(MFLAGS)
	$(CXX) src/annotate.cpp -o bin/annotate -pthread -lstdc++fs $(MFLAGS)
	$(CXX) src/rwr.cpp -o bin/rwr $(MFLAGS)

	$(CXX) src/data_stat.cpp -o util/data_stat $(MFLAGS)
	$(CXX) src/network_merge.cpp -o util/network_merge $(MFLAGS)
	$(CXX) src/enrich.cpp -o util/enrich $(MFLAGS)
	$(CXX) src/generate_expr_matrix_from_rsem.cpp -o util/generate_expr_matrix_from_rsem $(LFLAGS)
	$(CXX) src/generate_expr_matrix_from_stringtie.cpp -o util/generate_expr_matrix_from_stringtie $(LFLAGS)


windows: 
	mkdir -p bin
	mkdir -p util

	$(CXX) src/data_norm.cpp -o bin/data_norm $(WFLAGS)
	$(CXX) src/data_filter.cpp -o bin/data_filter $(WFLAGS)
	$(CXX) src/network_build.cpp -o bin/network_build -pthread $(WFLAGS)
	$(CXX) src/module_identify.cpp -o bin/module_identify -pthread $(WFLAGS)
	$(CXX) src/annotate.cpp -o bin/annotate -pthread -lstdc++fs $(WFLAGS)
	$(CXX) src/rwr.cpp -o bin/rwr $(WFLAGS)

	$(CXX) src/data_stat.cpp -o util/data_stat $(WFLAGS)
	$(CXX) src/network_merge.cpp -o util/network_merge $(WFLAGS)
	$(CXX) src/enrich.cpp -o util/enrich $(WFLAGS)
	$(CXX) src/generate_expr_matrix_from_rsem.cpp -o util/generate_expr_matrix_from_rsem $(LFLAGS)
	$(CXX) src/generate_expr_matrix_from_stringtie.cpp -o util/generate_expr_matrix_from_stringtie $(LFLAGS)


clean:
	rm -rf bin
	rm -rf util
