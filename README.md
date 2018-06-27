# min-bn-math
Minimalist big number math C library / RSA encryption demo.
Acknowledgement: This is based on some code found on a random website, apparently written by "z0mbie@i.am".

Test on windows/x86-64:

	C:\>gcc --version
	gcc (tdm64-1) 5.1.0
	Copyright (C) 2015 Free Software Foundation, Inc.
	This is free software; see the source for copying conditions.  There is NO
	warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

	C:\>gcc -std=c11 min-bn-math-test.c

	C:\>a.exe

	Constants
	BN_WORD_WIDTH     = (sizeof(WORD)*8) = 00000008
	BN_BIT            = ((56)+((sizeof(WORD)*8))) = 00000040
	BN_BYTES          = (((56)+((sizeof(WORD)*8)))/8) = 00000008
	BN_WORDS          = ((((56)+((sizeof(WORD)*8))))/((sizeof(WORD)*8))) = 00000008
	BN_WORDS_MSB_MASK = (((WORD)1)<<(((sizeof(WORD)*8))-1)) = 00000080

	Launching unit tests
	test_bn_mov0 pass
	test_bn_movc pass
	test_bn_copy pass
	test_bn_shl1 pass
	test_bn_add pass
	test_bn_sub pass
	test_bn_cmp pass
	test_bn_msb pass
	test_bn_mul pass
	test_bn_mod pass
	test_bn_div pass
	test_bn_modmul pass
	test_bn_modexp pass

	min_bn_math_self_test pass

	C:\>gcc -std=c11 min-bn-math-rsa-test.c

	C:\>a.exe

	sizeof(WORD)=8

	d=00000000000000003F5E23F45658D4002EB86A9A1AF1CA32CAAC69D2E9716050FDB470EC894BB4349559E413328345F0CA743C5FF004E1CA69C7F5258CFB26F83E3BEF4D7CD48CB9B996A1C8BBC762446EB10BA0C8B066E3A75C52F84DD16051C76CF8
	5FA8BABBA373475B4383A133A12AABD374681F29E898DDC2A578D2118550DA25003B88028B
	e=000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000
	00000000000000000000000000000000000000000000000000000000000000000000000011
	n=0000000000000000C3DD295050B5780090688F50B05FB6CB89B7D2BA7475B5576D44FFF2628CE72E2AB8C0F584DB925CA050006E5A3DA2A02FAFDE7410D9BE451D736F4C991C9BB42C32EA89267D2A658AE9D64CA46282A9FD2660A5A62E4463386E1E
	B17B2A30A7ADFFF3421A09ADE81B4E9E7503FAE507C7044B6186BCD3426ABFE1A178EDE91D
	x=Basic RSA enc/dec demo, use for educational purposes only !
	x=1577967905144362673266433177505423038122437959018968973900038225966870973606946826990427549786636453787275334088243397676835686379289079669058
	y=798373083942760828323940994519626465574838536530274723580110153332846330122469126238459040863981772603119045694371814213863495467788130425830781634606959589697181984292934384163807107804916952167417
	41413468987522171899038568777375760893736736604396826856704709437643616981249394518637405444655334823570539370
	z=1577967905144362673266433177505423038122437959018968973900038225966870973606946826990427549786636453787275334088243397676835686379289079669058
	z=Basic RSA enc/dec demo, use for educational purposes only !

	rsa_test pass

Confirmation using Python:

	C:\>python
	Python 3.6.5 (v3.6.5:f59c0932b4, Mar 28 2018, 17:00:18) [MSC v.1900 64 bit (AMD64)] on win32
	Type "help", "copyright", "credits" or "license" for more information.
	>>>d=44498337576932202591647937465965054805909571173905436721476387529105872932121462757382798827401460305961764413899350875562835523529416264032190619461235667960492995416894733958133247754192399449400041278095951041530225524273928470983356352923286280410417333809349064585688133716873643098654267641882774864523
	>>> e=17
	>>>n=137540316146881353465093624894801078490993219992071349866381561453599970881102703068274105466513604582063635461143448160830582527272741179735861914698364815342417850390232276308151352819576826952043882255852851892751378362922856527910065995914800713670770928829933567130814744455131846232076915939629570320669
	>>> x=1577967905144362673266433177505423038122437959018968973900038225966870973606946826990427549786636453787275334088243397676835686379289079669058
	>>> print("y=%d"%pow(x,e,n))
	y=79837308394276082832394099451962646557483853653027472358011015333284633012246912623845904086398177260311904569437181421386349546778813042583078163460695958969718198429293438416380710780491695216741741413468987522171899038568777375760893736736604396826856704709437643616981249394518637405444655334823570539370
	>>> print("z=%d"%pow(y,d,n))
	z=1577967905144362673266433177505423038122437959018968973900038225966870973606946826990427549786636453787275334088243397676835686379289079669058
	>>> x==z
	True
