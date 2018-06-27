//build with:
//gcc -std=c11 min-bn-math-rsa-test.c

#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <string.h>
#include <stdint.h>

#define TEST_WORD_WIDTH 64

#if TEST_WORD_WIDTH==64
typedef uint64_t WORD;
#elif TEST_WORD_WIDTH==32
typedef uint32_t WORD;
#elif TEST_WORD_WIDTH==16
typedef uint16_t WORD;
#endif

#define BN_SAFE_WIDTH 1024
#include "min-bn-math.h"
#include "min-bn-io.h"

//RSA1024 test key (the key generation is not covered in this demo)
//d: private key exponent
//e: public key exponent
//n: modulus for both the public and private keys
//so public key is e and n
//private key is d

#if TEST_WORD_WIDTH==64
WORD rsa_d[16] = {
  0x50DA25003B88028B,0x98DDC2A578D21185,
  0x2AABD374681F29E8,0x73475B4383A133A1,
  0xC76CF85FA8BABBA3,0xA75C52F84DD16051,
  0x6EB10BA0C8B066E3,0xB996A1C8BBC76244,
  0x3E3BEF4D7CD48CB9,0x69C7F5258CFB26F8,
  0xCA743C5FF004E1CA,0x9559E413328345F0,
  0xFDB470EC894BB434,0xCAAC69D2E9716050,
  0x2EB86A9A1AF1CA32,0x3F5E23F45658D400
};

WORD rsa_e[16] = {
  0x0000000000000011,0x0000000000000000,
  0x0000000000000000,0x0000000000000000,
  0x0000000000000000,0x0000000000000000,
  0x0000000000000000,0x0000000000000000,
  0x0000000000000000,0x0000000000000000,
  0x0000000000000000,0x0000000000000000,
  0x0000000000000000,0x0000000000000000,
  0x0000000000000000,0x0000000000000000
};

WORD rsa_n[16] = {
  0x6ABFE1A178EDE91D,0xC7044B6186BCD342,
  0x1B4E9E7503FAE507,0xADFFF3421A09ADE8,
  0x386E1EB17B2A30A7,0xFD2660A5A62E4463,
  0x8AE9D64CA46282A9,0x2C32EA89267D2A65,
  0x1D736F4C991C9BB4,0x2FAFDE7410D9BE45,
  0xA050006E5A3DA2A0,0x2AB8C0F584DB925C,
  0x6D44FFF2628CE72E,0x89B7D2BA7475B557,
  0x90688F50B05FB6CB,0xC3DD295050B57800
};
#elif TEST_WORD_WIDTH==32
WORD rsa_d[32] = {
  0x3B88028B,0x50DA2500,0x78D21185,0x98DDC2A5,
  0x681F29E8,0x2AABD374,0x83A133A1,0x73475B43,
  0xA8BABBA3,0xC76CF85F,0x4DD16051,0xA75C52F8,
  0xC8B066E3,0x6EB10BA0,0xBBC76244,0xB996A1C8,
  0x7CD48CB9,0x3E3BEF4D,0x8CFB26F8,0x69C7F525,
  0xF004E1CA,0xCA743C5F,0x328345F0,0x9559E413,
  0x894BB434,0xFDB470EC,0xE9716050,0xCAAC69D2,
  0x1AF1CA32,0x2EB86A9A,0x5658D400,0x3F5E23F4
};

WORD rsa_e[32] = {
  0x00000011,0x00000000,0x00000000,0x00000000,
  0x00000000,0x00000000,0x00000000,0x00000000,
  0x00000000,0x00000000,0x00000000,0x00000000,
  0x00000000,0x00000000,0x00000000,0x00000000,
  0x00000000,0x00000000,0x00000000,0x00000000,
  0x00000000,0x00000000,0x00000000,0x00000000,
  0x00000000,0x00000000,0x00000000,0x00000000,
  0x00000000,0x00000000,0x00000000,0x00000000
};

WORD rsa_n[32] = {
  0x78EDE91D,0x6ABFE1A1,0x86BCD342,0xC7044B61,
  0x03FAE507,0x1B4E9E75,0x1A09ADE8,0xADFFF342,
  0x7B2A30A7,0x386E1EB1,0xA62E4463,0xFD2660A5,
  0xA46282A9,0x8AE9D64C,0x267D2A65,0x2C32EA89,
  0x991C9BB4,0x1D736F4C,0x10D9BE45,0x2FAFDE74,
  0x5A3DA2A0,0xA050006E,0x84DB925C,0x2AB8C0F5,
  0x628CE72E,0x6D44FFF2,0x7475B557,0x89B7D2BA,
  0xB05FB6CB,0x90688F50,0x50B57800,0xC3DD2950
};
#elif TEST_WORD_WIDTH==16
WORD rsa_d[64] = {
  0x028B,0x3B88,0x2500,0x50DA,0x1185,0x78D2,0xC2A5,0x98DD,
  0x29E8,0x681F,0xD374,0x2AAB,0x33A1,0x83A1,0x5B43,0x7347,
  0xBBA3,0xA8BA,0xF85F,0xC76C,0x6051,0x4DD1,0x52F8,0xA75C,
  0x66E3,0xC8B0,0x0BA0,0x6EB1,0x6244,0xBBC7,0xA1C8,0xB996,
  0x8CB9,0x7CD4,0xEF4D,0x3E3B,0x26F8,0x8CFB,0xF525,0x69C7,
  0xE1CA,0xF004,0x3C5F,0xCA74,0x45F0,0x3283,0xE413,0x9559,
  0xB434,0x894B,0x70EC,0xFDB4,0x6050,0xE971,0x69D2,0xCAAC,
  0xCA32,0x1AF1,0x6A9A,0x2EB8,0xD400,0x5658,0x23F4,0x3F5E,
};

WORD rsa_e[64] = {
  0x0011,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,
  0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,
  0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,
  0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,
  0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,
  0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,
  0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,
  0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000,0x0000
};

WORD rsa_n[64] = {
  0xE91D,0x78ED,0xE1A1,0x6ABF,0xD342,0x86BC,0x4B61,0xC704,
  0xE507,0x03FA,0x9E75,0x1B4E,0xADE8,0x1A09,0xF342,0xADFF,
  0x30A7,0x7B2A,0x1EB1,0x386E,0x4463,0xA62E,0x60A5,0xFD26,
  0x82A9,0xA462,0xD64C,0x8AE9,0x2A65,0x267D,0xEA89,0x2C32,
  0x9BB4,0x991C,0x6F4C,0x1D73,0xBE45,0x10D9,0xDE74,0x2FAF,
  0xA2A0,0x5A3D,0x006E,0xA050,0x925C,0x84DB,0xC0F5,0x2AB8,
  0xE72E,0x628C,0xFFF2,0x6D44,0xB557,0x7475,0xD2BA,0x89B7,
  0xB6CB,0xB05F,0x8F50,0x9068,0x7800,0x50B5,0x2950,0xC3DD
};
#endif



int main(void){
	printf("\nsizeof(WORD)=%lu\n\n",sizeof(WORD));
	WORD d[BN_WORDS];
	WORD e[BN_WORDS];
	WORD n[BN_WORDS];

	memset(d,0,sizeof(d));
	memset(e,0,sizeof(e));
	memset(n,0,sizeof(n));
	memcpy(d,rsa_d,sizeof(rsa_d));
	memcpy(e,rsa_e,sizeof(rsa_e));
	memcpy(n,rsa_n,sizeof(rsa_n));

	bn_dumphex("d=",d);
	bn_dumphex("e=",e);
	bn_dumphex("n=",n);

	WORD x[BN_WORDS];
	WORD y[BN_WORDS];
	WORD z[BN_WORDS];

	bn_mov0(x);
	strcpy((char*)&x,"Basic RSA enc/dec demo, use for educational purposes only !");

	printf("x=%s\n",(char*)&x);
	bn_dumpdec("x=",x);


	bn_modexp(y, x,e,n);

	bn_dumpdec("y=",y);

    //ref value computed using Python
    WORD expected_y[BN_WORDS];
    bn_mov0(expected_y);
    bn_movdec(expected_y,"79837308394276082832394099451962646557483853653027472358011015333284633012246912623845904086398177260311904569437181421386349546778813042583078163460695958969718198429293438416380710780491695216741741413468987522171899038568777375760893736736604396826856704709437643616981249394518637405444655334823570539370");

    if(bn_cmp(y,expected_y)){
        printf("\nERROR: incorrect encryption\n\n");
        exit(-1);
    }

	bn_modexp(z, y,d,n);

	bn_dumpdec("z=",z);
	printf("z=%s\n",(char*)&z);
	if(strcmp((char*)x,(char*)z)){
		printf("\nrsa_test fail\n");
		return -1;
	}
	printf("\nrsa_test pass\n");
	return 0;
}
