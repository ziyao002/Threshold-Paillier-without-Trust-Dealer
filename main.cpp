#include <iostream>
#include <stdlib.h>
#include "threshold_paillier_without_trust_dealer.h"

using namespace std;
using namespace NTL;

int main()
{
	Paillier_wt paillier_wt(128);											// key generation and distribution
	
	ZZ massage;
	massage = 34614984;														// plaintext massage
		
    ZZ c = paillier_wt.encrypt(massage);									// ciphertext c = encryption(massage)
	
	ZZ c1 = paillier_wt.partial_decrypt(c, paillier_wt.f1);					// partial decryption, fi is the share of the secrey key
	ZZ c2 = paillier_wt.partial_decrypt(c, paillier_wt.f2);
	ZZ c3 = paillier_wt.partial_decrypt(c, paillier_wt.f3);
	
	ZZ dec_c = paillier_wt.combine_partial_decrypt(c1, c2, c3);				// combine partial decryptions
	
	if (dec_c == massage){
        cout << "Encryption and distributed decryption is successful" << endl;
    }
	
    return 0;
}