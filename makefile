CXX = g++
CXXFLAGS = -Wall

mwis3:
	$(CXX) $(CXXFLAGS) -std=c++11 -O3 -o mwis3 MWIS3.cpp Graph.cpp
	
clean:
	rm mwis3
