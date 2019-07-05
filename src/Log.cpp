#include "Header.h"
#include <iostream>
#include <string>


void InitLog() {
	Log("Initializing Log");
}

void Log(std::string message) {
	std::cout << message << std::endl;
}

