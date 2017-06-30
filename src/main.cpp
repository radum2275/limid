/*
 ============================================================================
 Name        : main.cpp
 Author      : Radu Marinescu
 Version     :
 Copyright   : Copyright (c) IBM Corp. 2015
 Description :
 ============================================================================
 */

#include <valve.h>
#include "be.h"
#include "mbe.h"
#include "spu.h"

int main(void) {

	std::cout << VERSIONINFO << std::endl << COPYRIGHT << std::endl;

	merlin::limid gm;
	gm.read("/home/radu/git/limid/examples/car.uai");
	merlin::be s1(gm);
	s1.run();
	merlin::valve s2(gm);
	s2.run();

//	merlin::mbe m(gm);
//	m.run();
//	m.run2();

//	merlin::limid gm;
//	gm.read("/home/radu/git/limid/examples/chain.uai");
//	merlin::spu s(gm);
//	s.run();
//	s.brute_force();

	return 0;
}


