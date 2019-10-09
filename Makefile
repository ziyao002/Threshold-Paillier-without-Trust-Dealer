LIBS:= ntl gmp
LIBFLAGS:=$(addprefix -l, $(LIBS));

main : main.cpp threshold_paillier_without_trust_dealer.cpp
	g++ -std=c++11 -lpthread -g -O3 $^ -o $@ $(LIBFLAGS)
