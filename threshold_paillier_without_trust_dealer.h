#include <NTL/ZZ.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/vector.h>


class Paillier_wt {
    public:
	
	NTL::ZZ K; 							// the security parameter such that (2**K)**-1 is negligible
	long t, l;
	long keyLength;
	
    Paillier_wt(const long input_keyLength);
    NTL::ZZ encrypt(NTL::ZZ& message); 
    NTL::ZZ partial_decrypt(NTL::ZZ& ciphertext, NTL::ZZ& fi); 
    NTL::ZZ combine_partial_decrypt(NTL::ZZ& c1, NTL::ZZ& c2, NTL::ZZ& c3);
	
	NTL::ZZ pick_pq(long i);
	void coefficient_generation();
	void compute_tuple(NTL::ZZ& j, NTL::ZZ& PP, NTL::ZZ& pi, NTL::ZZ& qi, NTL::ZZ& pij, NTL::ZZ& ppij, NTL::ZZ& qij, NTL::ZZ& qqij, NTL::ZZ& hij, NTL::ZZ& hhij);
	bool biprimality_check(NTL::ZZ & N);
	void Ri_sharing(NTL::ZZ& r1, NTL::ZZ& r2, NTL::ZZ& r3, NTL::ZZ& R1, NTL::ZZ& R2, NTL::ZZ& R3, NTL::ZZ& modulus);
	NTL::ZZ distributed_RSA_modulus_generation();
	
	void ZKP_gen_R(NTL::ZZ& c, NTL::ZZ& r1, NTL::ZZ& R1, NTL::ZZ& R2);
	NTL::ZZ ZKP_gen_cc();
	NTL::ZZ ZKP_comput_z(NTL::ZZ& r, NTL::ZZ& cc, NTL::ZZ& fi);
	bool ZKP_check(NTL::ZZ& c, NTL::ZZ& ci, NTL::ZZ& R1, NTL::ZZ& R2, NTL::ZZ& z, NTL::ZZ& cc, NTL::ZZ& vki);
	bool ZKP_for_partial_decryption(NTL::ZZ& cc_massage, NTL::ZZ& c1, NTL::ZZ& c2, NTL::ZZ& c3);
	
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
	
	NTL::ZZ vk;							// public
	NTL::ZZ vk1, vk2, vk3;				// public
	
	private:
	NTL::ZZ L_function(const NTL::ZZ& x) { return (x - 1) / N; }
};

