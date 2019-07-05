#pragma once

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include <array>
#include <vector>

using String = std::string;

class GLmapper {
private:
	const int MAX_ELEMS;
	const String GENERALINFO;
	const String MESHFILES;
	const String PRESSUREFILES;
	const String MAPPINGFILE;
	std::ifstream inFile1;

	String file_gelem, file_lelem, file_gnode, file_lnode;
	String file_gpressure, file_lpressure;


	String GNodes, GELEM, GPressure, LNodes, LELEM;

	std::vector<int> gn, ln, le;
	//std::vector<std::array<int, 5>> GE_N;
	//std::vector<std::array<int, 4>>	LE_N;
	//std::vector<std::array<int, 2>> GPressure_MAC;
	//std::vector<std::array<int, 5>> GXYZ;
	//std::vector<std::array<int, 3>> LXYZ;
	//std::vector<std::array<int, 5>> LXYZ_AVG;
	//std::vector<std::array<int, 5>> LPressure;

	const int COLS_ge_n;
	const int COLS_le_n;
	const int COLS_gpressure_mac;
	const int COLS_gxyz;
	const int COLS_lxyz;
	const int COLS_gxyz_avg;
	const int COLS_lxyz_avg;
	const int COLS_lpressure;
	const int HEADER_LENGTH;

	int* ge_n;
	int* le_n;
	double* gpressure_mac;
	double* gxyz;
	double* lxyz;
	double* gxyz_avg;
	double* lxyz_avg;
	String* lpressure;
	int m, N, o, p, q, r;
	int ilocal, jglobal;
	double offset_X, offset_Y, offset_Z;
	double outoffset_x, outoffset_y, outoffset_z;
	double radius, tolerance;


	void check(std::ifstream& inFile);

public:
	GLmapper(String s1, String s2, String s3, String s4);
	~GLmapper();
	void run();
	void inputdata(bool found);
	void mapping();
	//bool normal(double[] local_coords, double[] global_coords);
	int minDistance(double local_point[]);
	void pressure_cycle(int cycle);
	void output(int cycle, bool found);


};