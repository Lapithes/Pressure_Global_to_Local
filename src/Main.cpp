#include "Header.h"
#include "GLmapper.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>
#include <chrono>

using String = std::string;

int main() {


	InitLog();
	Log("Press enter to run...");
	std::cin.get();

	GLmapper* mapper = new GLmapper(".\\input\\generalinfo.txt", ".\\input\\meshfiles.txt", ".\\input\\pressurefiles.txt", ".\\input\\mapping.txt");
	auto start = std::chrono::steady_clock::now();

	mapper->run();

	delete mapper;

	auto stop = std::chrono::steady_clock::now();
	auto seconds = std::chrono::duration_cast<std::chrono::seconds>(stop - start);

	std::cout << "Runtime: " << seconds.count() << "s" << std::endl;
	//////////////////////////
	std::cin.get();
}

