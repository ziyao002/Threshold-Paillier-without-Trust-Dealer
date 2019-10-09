#include <iostream>
#include <stdlib.h>
#include "threshold_paillier_without_trust_dealer.h"

using namespace std;
using namespace NTL;

/* Reference: Nishide, T., & Sakurai, K. (2010, August). Distributed paillier cryptosystem without trusted dealer. */

NTL::ZZ Gen_Coprime(const NTL::ZZ& n){
	 /* Coprime generation function. Generates a random coprime number of n.
	 *
     * Parameters
     * ==========
     * NTL::ZZ n : a prime number.
     *
     * Returns
     * =======
     * NTL:ZZ ret : a random coprime number of n.
     */
     NTL::ZZ ret;
    while (true) {
        ret = RandomBnd(n);
        if (NTL::GCD(ret, n) == 1) { return ret; }
    }
}

NTL::ZZ lcm(const NTL::ZZ& x, const NTL::ZZ& y){
	/* Least common multiple function. Computes the least common multiple of x and y.
	 *
     * Parameters
     * ==========
     * NTL::ZZ x, y: signed, arbitrary length integers.
     * 
     * Returns
     * =======
     * NTL:ZZ lcm : the least common multiple of x and y.
     */
	NTL::ZZ lcm;
	lcm = (x * y) / NTL::GCD(x,y);
	return lcm;
}

NTL::ZZ Paillier_wt::pick_pq(const long keyLength = 512, long i = 1){
	/* Pick pi, qi for Prime candidate generation function.
	 *
     * Parameters
     * ==========
     * long keyLength: the length of the key.
     * long i: the index of Pi.
     * 
     * Returns
     * =======
     * NTL:ZZ pq : prime candidate.
     */
	NTL::ZZ pq;
	pq = 2;
	if (i == 1){
		while (NTL::MulMod(pq, NTL::ZZ(1), NTL::ZZ(4)) != 3 || (pq - NTL::power(NTL::ZZ(2), keyLength - 1) <= 0)){
				pq = RandomBnd(NTL::power(NTL::ZZ(2), keyLength));
		}
	}
	else{
		while (NTL::MulMod(pq, NTL::ZZ(1), NTL::ZZ(4)) != 0 || (pq - NTL::power(NTL::ZZ(2), keyLength - 1) <= 0)){
				pq = RandomBnd(NTL::power(NTL::ZZ(2), keyLength));
			}
	}
	return pq;
}

void Paillier_wt::coefficient_generation(){
	/* Coefficient generation function for computing and verifying p * q = sum(pi) * sum(qi) in BGW protocol.
	 *
     * Parameters
     * ==========
     * f(x) = pi + a1 * x
     * ff(x) = ppi + aa1 * x
     * g(x) = qi + b1 * x
     * gg(x) = qqi + bb1 * x
     * h(x) = c0 + c1 * x + c2 * x ** 2
     * hh(x) = cc0 + cc1 * x + cc2 * x ** 2
     */
	a1	= RandomBnd(PP);
	aa1 = RandomBnd(PP);
	b1 	= RandomBnd(PP);
	bb1 = RandomBnd(PP);
	ppi = RandomBnd(PP);
	qqi = RandomBnd(PP);
	c0 	= 0;
	c1 	= RandomBnd(PP);
	c2 	= RandomBnd(PP);
	cc0 = RandomBnd(PP);
	cc1 = RandomBnd(PP);
	cc2 = RandomBnd(PP);

}

void Paillier_wt::compute_tuple(NTL::ZZ& j, NTL::ZZ& PP, NTL::ZZ& pi, NTL::ZZ& qi, NTL::ZZ& pij, NTL::ZZ& ppij, NTL::ZZ& qij, NTL::ZZ& qqij, NTL::ZZ& hij, NTL::ZZ& hhij){
	/* Compute the tuple to be sent to other parties.
	 *
	 * f(x) = pi + a1 * x
     * ff(x) = ppi + aa1 * x
     * g(x) = qi + b1 * x
     * gg(x) = qqi + bb1 * x
     * h(x) = c0 + c1 * x + c2 * x ** 2
     * hh(x) = cc0 + cc1 * x + cc2 * x ** 2
	 *
     * Parameters
     * ==========
     * NTL::ZZ j: the index of Pi.
     * NTL::ZZ PP: the large prime parties agreen on in advance.
     * NTL::ZZ pi: prime candidate generated in advance.
     * NTL::ZZ qi: prime candidate generated in advance.
     * NTL::ZZ pij: pij = f(j).
     * NTL::ZZ ppij: ppij = ff(j).
     * NTL::ZZ qij: qij = g(j).
     * NTL::ZZ qqij: qqij = gg(j).
     * NTL::ZZ hij: hij = h(j).
     * NTL::ZZ hhij: hhij = hh(j).
     */
	pij = (pi + a1 * j) % PP;
	ppij = (ppi + aa1 * j) % PP;
	qij = (qi + b1 * j) % PP;
	qqij = (qqi + bb1 * j) % PP;
	hij = (c0 + c1 * j + c2 * j * j) % PP;
	hhij = (cc0 + cc1 * j + cc2 * j * j) % PP;	
}

bool Paillier_wt::biprimality_check(NTL::ZZ & N){
	/* Biprimality check function.
	 *
     * Parameters
     * ==========
     * NTL::ZZ N : a large integer.
     *
     * Returns
     * =======
     * bool biprimality_check: biprimality_check = 1 if N is the product of two primes. 
     */
	NTL::ZZ gg;
	gg = Gen_Coprime(N);
	while (NTL::Jacobi(gg, N) != 1){
		gg = Gen_Coprime(N);
	}
	
	NTL::ZZ Q1, Q2, Q3, Q2_inv, Q3_inv;
	Q1 = NTL::PowerMod(gg, (N + 1 - p1 - q1)/4, N);
	Q2 = NTL::PowerMod(gg, (p2 + q2)/4, N);
	Q3 = NTL::PowerMod(gg, (p3 + q3)/4, N);		
	Q2_inv = NTL::InvMod(Q2, N);
	Q3_inv = NTL::InvMod(Q3, N);
	bool biprimality_check = (Q1 * Q2_inv * Q3_inv % N) == (NTL::ZZ(1) % N) || (Q1 * Q2_inv * Q3_inv % N) == (NTL::ZZ(-1) % N);		// biprimality check
	return biprimality_check;
}

NTL::ZZ Paillier_wt::distributed_RSA_modulus_generation(const long keyLength = 512){
	/* Distributed RSA modulus generation function under BGW protocol with biprimality check.
	 *
     * Parameters
     * ==========
     * long keyLength: the length of the key.
     */
	 
	int generation_done = 0;																// the flag of N passing the biprimality check.
	cout << "Distributed RSA modulus is being generated..." << endl;
	
	// Step 1: parties agree on a large prime PP
		PP = 0;
		long err = 80;
		while(PP <= NTL::power(3 * 3 * NTL::power(NTL::ZZ(2), keyLength - 1), 2)){
			PP = NTL::GenPrime_ZZ(5 + 2 * keyLength, err);									// PP > {n * (3 * 2^(k-1))}^2 > 2N (81 -> 2^7)
		}
		
	while (!generation_done){
		
		// Step 2: pick up pi, qi
		p1 = pick_pq(keyLength, 1);
		q1 = pick_pq(keyLength, 1);
		p2 = pick_pq(keyLength, 2);
		q2 = pick_pq(keyLength, 2);
		p3 = pick_pq(keyLength, 3);
		q3 = pick_pq(keyLength, 3);
	
		// Step 3: compute N = p * q using BGW protocol with zero-knowledge proof
		NTL::ZZ p11, pp11, q11, qq11, h11, hh11;
		NTL::ZZ p12, pp12, q12, qq12, h12, hh12;
		NTL::ZZ p13, pp13, q13, qq13, h13, hh13;
		NTL::ZZ p21, pp21, q21, qq21, h21, hh21;
		NTL::ZZ p22, pp22, q22, qq22, h22, hh22;
		NTL::ZZ p23, pp23, q23, qq23, h23, hh23;
		NTL::ZZ p31, pp31, q31, qq31, h31, hh31;
		NTL::ZZ p32, pp32, q32, qq32, h32, hh32;
		NTL::ZZ p33, pp33, q33, qq33, h33, hh33;
		NTL::ZZ j1, j2, j3;
		j1 = 1; 
		j2 = 2;
		j3 = 3;
		coefficient_generation();
		compute_tuple(j1, PP, p1, q1, p11, pp11, q11, qq11, h11, hh11); 				// 1 -> 1
		compute_tuple(j2, PP, p1, q1, p12, pp12, q12, qq12, h12, hh12); 				// 1 -> 2
		compute_tuple(j3, PP, p1, q1, p13, pp13, q13, qq13, h13, hh13); 				// 1 -> 3
		coefficient_generation();
		compute_tuple(j1, PP, p2, q2, p21, pp21, q21, qq21, h21, hh21); 				// 2 -> 1
		compute_tuple(j2, PP, p2, q2, p22, pp22, q22, qq22, h22, hh22); 				// 2 -> 2
		compute_tuple(j3, PP, p2, q2, p23, pp23, q23, qq23, h23, hh23); 				// 2 -> 3
		coefficient_generation();
		compute_tuple(j1, PP, p3, q3, p31, pp31, q31, qq31, h31, hh31); 				// 3 -> 1
		compute_tuple(j2, PP, p3, q3, p32, pp32, q32, qq32, h32, hh32); 				// 3 -> 2
		compute_tuple(j3, PP, p3, q3, p33, pp33, q33, qq33, h33, hh33); 				// 3 -> 3
		
		NTL::ZZ N1, N2, N3;																// partial modulus
		N1 = ((p11 + p21 + p31) * (q11 + q21 + q31) + (h11 + h21 + h31)) % PP;
		N2 = ((p12 + p22 + p32) * (q12 + q22 + q32) + (h12 + h22 + h32)) % PP;
		N3 = ((p13 + p23 + p33) * (q13 + q23 + q33) + (h13 + h23 + h33)) % PP;
		
		NTL::ZZ L1, L2, L3;
		L1 = (0 - 2) * (0 - 3) / ((1 - 2) * (1 - 3));									// Lagrange interpolation
		L2 = (0 - 1) * (0 - 3) / ((2 - 1) * (2 - 3));
		L3 = (0 - 1) * (0 - 2) / ((3 - 1) * (3 - 2));
		N = (N1 * L1 + N2 * L2 + N3 * L3) % PP;

		// Step 4: biprimality test
		if (biprimality_check(N)){generation_done = 1;}
	}
	cout << "RSA modulus generation is done!" << endl;
	return N;
}

Paillier_wt::Paillier_wt(const long keyLength = 512) {
	/* Distributed Paillier parameters generation function withour trust dealer.
	 *
     * Parameters
     * ==========
     * long keyLength: the length of the key.
     * 
     * =======
     * public key  = (modulus, generator, theta).
	 * private key = delta * phi * beta.
     */
	
	distributed_RSA_modulus_generation(keyLength);																// RSA modulus generation
	generator = N + 1;																							// generator = n + 1
	delta = 6;																									// 3 partirs = 3!
	
	NTL::ZZ K, beta1, beta2, beta3, r1, r2, r3;																	// parameters for secret-sharing secrect key over the integer
	NTL::ZZ phi, theta1, theta2, theta3;
	K = NTL::power(NTL::ZZ(2), 80);																				// the security parameter such that (2**K)**-1 is negligible
	beta1 =  RandomBnd(K * N);
	beta2 =  RandomBnd(K * N);
	beta3 =  RandomBnd(K * N);
	r1 = RandomBnd(K * K * N);
	r2 = RandomBnd(K * K * N);
	r3 = RandomBnd(K * K * N);
	
	phi = N + 1 - (p1 + q1) - (p2 + q2) - (p3 + q3);
	theta1 = delta * phi * beta1 + N * delta * r1;
	theta2 = delta * phi * beta2 + N * delta * r2;
	theta3 = delta * phi * beta3 + N * delta * r3;
	theta = (theta1 + theta2 + theta3) % N;
	
	ak1 = RandomBnd(N);																							// parties agree on ak1 and ak2
	ak2 = RandomBnd(N);
	f1 = -(delta * phi * beta1 + delta * phi * beta2 + delta * phi * beta3 + ak1 * 1 + ak2 * 1 * 1);			// the secret key is shared over the integer
	f2 = -(delta * phi * beta1 + delta * phi * beta2 + delta * phi * beta3 + ak1 * 2 + ak2 * 2 * 2);
	f3 = -(delta * phi * beta1 + delta * phi * beta2 + delta * phi * beta3 + ak1 * 3 + ak2 * 3 * 3);
}

NTL::ZZ Paillier_wt::encrypt(NTL::ZZ& message){
	/* Paillier encryption function. Takes in a message in F(modulus), and returns a message in F(modulus**2).
     *
     * Parameters
     * ==========
     * NTL::ZZ message : the message to be encrypted.
     *
     * Returns
     * =======
     * NTL:ZZ ciphertext : the encyrpted message.
     */
    NTL::ZZ random = Gen_Coprime(N);
    NTL::ZZ ciphertext = NTL::PowerMod(generator, message, N * N) * NTL::PowerMod(random, N, N * N);
    return ciphertext % (N * N);
}

NTL::ZZ Paillier_wt::partial_decrypt( NTL::ZZ& ciphertext, NTL::ZZ& fi){
	/* Paillier partial decryption function. Takes in a ciphertext in F(modulus**2), and returns a partial decryption in the same space.
     *
      * Parameters
     * ==========
     * NTL::ZZ cipertext : the encryption of the original message.
     * NTL::ZZ fi : the Pi's share of secret key, i.e., delta * phi * beta.
     *
     * Returns
     * =======
     * NTL::ZZ partial_decryption : The partial decryption of the original message.
     */
    NTL::ZZ partial_decryption = NTL::PowerMod(ciphertext, 2 * delta * fi, N * N);
    return partial_decryption;
}

NTL::ZZ Paillier_wt::combine_partial_decrypt(NTL::ZZ& c1, NTL::ZZ& c2, NTL::ZZ& c3){
	/* Combine the partial decryptions to obtain the decryption of the original ciphertext.
     *
     * Parameters
     * ==========
     * NTL::ZZ c1, c2, c3 : the partial decryptions.
     *
     * Returns
     * =======
     * NTL::ZZ M: the decryption of the original message.
     */
	NTL::ZZ lamda1, lamda2, lamda3;		// parameters for 3-party Paillier
	lamda1= 3;
	lamda2 = -3;
	lamda3 = 1;
	NTL::ZZ u1 = delta * lamda1;
	NTL::ZZ u2 = delta * lamda2;
	NTL::ZZ u3 = delta * lamda3;

    NTL::ZZ product_1 = NTL::PowerMod(c1, 2 * u1, N * N);
    NTL::ZZ product_2 = NTL::PowerMod(c2, 2 * u2, N * N);
    NTL::ZZ product_3 = NTL::PowerMod(c3, 2 * u3, N * N);
    NTL::ZZ product = NTL::MulMod(NTL::MulMod(product_1, product_2, N * N), product_3, N * N);
	NTL::ZZ Inv_temp = NTL::InvMod(-4 * delta * delta * theta % N, N);
	NTL::ZZ M = NTL::MulMod(L_function(product), Inv_temp, N);

    return M;
}























