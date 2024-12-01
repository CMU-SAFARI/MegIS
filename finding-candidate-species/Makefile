all:
	$(CXX) --std=c++17 -O3 -Wall -Werror -DNDEBUG -fopenmp -march=native -o intersection intersection.cpp progress_bar.cpp

profile:
	$(CXX) --std=c++17 -O2 -Wall -Werror -DNDEBUG -fopenmp -march=native -pg -g -o intersection_profile intersection.cpp progress_bar.cpp
