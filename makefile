SHARKFLAGS = -shared
CXX = gcc #-std=gnu++11

obj/%.o: src/%.cpp
	$(CXX) -fPIC -O2 $(SHARKFLAGS) -c -o $@ $<

LIBshark = lib/libshark.so

Oshark = \
	shark.o \
	polylogarithm/src/cpp/Li2.o \
	polylogarithm/src/cpp/Li3.o \
	polylogarithm/src/cpp/Li4.o

$(LIBshark): $(patsubst %.o,obj/%.o,$(Oshark))
	$(CXX) -o $@ $+ $(SHARKFLAGS)

clean:
	rm obj/*.o obj/polylogarithm/src/cpp/*.o $(LIBshark)
