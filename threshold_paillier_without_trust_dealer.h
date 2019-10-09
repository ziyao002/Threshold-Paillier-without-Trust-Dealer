#include <NTL/ZZ.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/vector.h>


class Paillier_wt {
    public:
    Paillier_wt(const long keyLength);
    NTL::ZZ encrypt(NTL::ZZ& message); 
    NTL::ZZ partial_decrypt(NTL::ZZ& ciphertext, NTL::ZZ& fi); 
    NTL::ZZ combine_partial_decrypt(NTL::ZZ& c1, NTL::ZZ& c2, NTL::ZZ& c3);

	
	NTL::ZZ pick_pq(const long keyLength, long i);
	void coefficient_generation();
	void compute_tuple(NTL::ZZ& j, NTL::ZZ& PP, NTL::ZZ& pi, NTL::ZZ& qi, NTL::ZZ& pij, NTL::ZZ& ppij, NTL::ZZ& qij, NTL::ZZ& qqij, NTL::ZZ& hij, NTL::ZZ& hhij);
	bool biprimality_check(NTL::ZZ & N);
	NTL::ZZ distributed_RSA_modulus_generation(const long keyLength);
	
	NTL::ZZ PP;							// public
	NTL::ZZ p1, p2, p3, q1, q2, q3;		// pi and qi are private to Pi
	NTL::ZZ a1;							// private to Pi
	NTL::ZZ ppi, aa1;					// private to Pi
	NTL::ZZ b1;							// private to Pi
	NTL::ZZ qqi, bb1;					// private to Pi
	NTL::ZZ c0, c1, c2;					// private to Pi
	NTL::ZZ cc0, cc1, cc2;				// private to Pi
	
	NTL::ZZ N;							// public
	NTL::ZZ generator;					// public
	NTL::ZZ delta;						// public
	NTL::ZZ f1, f2, f3;					// fi is private to Pi
	NTL::ZZ theta;						// public
	NTL::ZZ ak1, ak2;					// public
	
	private:
	NTL::ZZ L_function(const NTL::ZZ& x) { return (x - 1) / N; }
};

