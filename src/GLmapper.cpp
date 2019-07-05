#include "Header.h"
#include "GLmapper.h"
#include <iostream>
#include <iomanip>

using String = std::string;

GLmapper::GLmapper(String s1, String s2, String s3, String s4) : GENERALINFO(s1), MESHFILES(s2), PRESSUREFILES(s3), MAPPINGFILE(s4),
MAX_ELEMS(999999), COLS_ge_n(5), COLS_le_n(4), COLS_gpressure_mac(1), COLS_gxyz(5), COLS_lxyz(3), 
COLS_lxyz_avg(5), COLS_gxyz_avg(5), COLS_lpressure(5), HEADER_LENGTH(9)
{
	ge_n = new int[MAX_ELEMS * 5];
	le_n = new int[MAX_ELEMS * 4];
	gxyz = new double[MAX_ELEMS * 5];
	lxyz = new double[MAX_ELEMS * 3];
	gxyz_avg = new double[MAX_ELEMS * 5];
	lxyz_avg = new double[MAX_ELEMS * 5];
	lpressure = new String[MAX_ELEMS * 5];
	gpressure_mac = new double[MAX_ELEMS];

	for (int i = 0; i < MAX_ELEMS; i++) {
		gpressure_mac[i] = 999.0;
	}

}

GLmapper::~GLmapper() {
	delete[] ge_n;
	delete[] le_n;
	delete[] gpressure_mac;
	delete[] gxyz;
	delete[] lxyz;
	delete[] gxyz_avg;
	delete[] lxyz_avg;
	delete[] lpressure;
}

void GLmapper::run() {
	int i, j, k;
	int cycle;
	bool found = false;
	String line;
	

	inFile1.open(MAPPINGFILE);
	if (inFile1.fail()) {
		Log("An existing mapping file cannot be found. Proceeding to generate a new mapping.");
	}
	else {
		bool valid = false;
		String keyboardin;

		std::cout << "An existing mapping file \"" << MAPPINGFILE << "\" has been found." << std::endl;

		while (!valid) {
			Log("Would you like to use this mapping file? (y/n)");

			getline(std::cin, keyboardin);
			if (keyboardin == "y") {
				valid = true;
				found = true;
			}
			else if (keyboardin == "n") {
				std::cout << "\"" << MAPPINGFILE << "\" will be overwritten. Press enter to proceed." << std::endl;
				std::cin.get();
				valid = true;
			}
			else {
				Log("Invalid input, please try again.");
				Log("");
			}
		}
	}

	if (found) {
		String delimiter = "->";
		int pos = 0;
		int index = 0;
		while (getline(inFile1, line)) {
			pos = line.find(delimiter);
			lpressure[COLS_lpressure * index + 0] = line.substr(0, pos);
			line = line.substr(pos + delimiter.length());
			lpressure[COLS_lpressure * index + 1] = line.substr(0);
			index++;
		}
		o = index;
		inFile1.close();
	}
	else {
		inFile1.close();
		inputdata(found);
		Log("");
		Log("mapping global to local elements");
		mapping();
	}


	Log("");
	Log("writing pressure");

	cycle = 5;
	inFile1.open(PRESSUREFILES);
	check(inFile1);

	int index = 0; 
	while (getline(inFile1, line)) {
		index++;
	}
	inFile1.close();

	if (index == 0) {
		std::cerr << "No pressure files found." << std::endl;
		std::cin.get();
		exit(1);
	}
	else {
		std::cout << index << " pressure files found" << std::endl;
	}

	String g, l;
	inFile1.open(PRESSUREFILES);
	for (int i = 0; i < index; i++) {
		Log("");

		inFile1 >> g >> l;
		file_gpressure = g;
		file_lpressure = l;

		pressure_cycle(cycle);
		std::cout << file_gpressure << " read" << std::endl;
		output(cycle, found);
		std::cout << file_lpressure << " written" << std::endl;
		cycle++;
	}
	inFile1.close();
	Log("");
	Log("Done");
}

void GLmapper::inputdata(bool found) {
	int i, nstart, nstop, j, k;
	int ii, jj;

	Log("");
	//general info
	inFile1.open(GENERALINFO);
	check(inFile1);

	double number;
	double generalinfoarray[5];
	for (int i = 0; i < 5; i++) {
		inFile1 >> number;
		generalinfoarray[i] = number;
	}

	offset_X = generalinfoarray[0];
	offset_Y = generalinfoarray[1];
	offset_Z = generalinfoarray[2];
	tolerance = generalinfoarray[3];
	radius = generalinfoarray[4];

	inFile1.close();

	Log("general info read");


	//meshfiles
	inFile1.open(MESHFILES);
	//std::cin.get();
	check(inFile1);

	String word;
	String meshfilesarray[4];
	for(int i = 0; i < 4; i++) {
		inFile1 >> word;
		meshfilesarray[i] = word;
	}

	file_gelem = meshfilesarray[0];
	file_gnode = meshfilesarray[1];
	file_lelem = meshfilesarray[2];
	file_lnode = meshfilesarray[3];

	inFile1.close();

	Log("meshfiles read");

	//read from local and global element/node files
	//file_gelem
	inFile1.open(file_gelem);
	check(inFile1);

	String line;
	int index = 0; 
	while (getline(inFile1, line)) {
		ge_n[(COLS_ge_n)*index + 4] = std::stoi(line.substr(78, 6));
		ge_n[(COLS_ge_n)*index + 0] = std::stoi(line.substr(0, 6));
		ge_n[(COLS_ge_n)*index + 1] = std::stoi(line.substr(6, 6));
		ge_n[(COLS_ge_n)*index + 2] = std::stoi(line.substr(12, 6));
		ge_n[(COLS_ge_n)*index + 3] = std::stoi(line.substr(18, 6));
		index++;
	}
	N = index;
	inFile1.close();
	//nstart = ge_n[0 + 4];
	//nstop = ge_n[(MAX_ELEMS - 1) + 4];
	Log("file_gelem read");

	//file_gnode, global node file
	inFile1.open(file_gnode);
	check(inFile1);
	//Log("opened file_gnode");
	index = 0;
	while (getline(inFile1, line)) {
		//temp = node number
		int temp = std::stoi(line.substr(0, 8));
		//if ((temp >= nstart) && (temp <= nstop)) {
		gn.push_back(temp);
		gxyz[(COLS_gxyz) * (temp) + 0] = std::stod(line.substr(8, 20)) + offset_X;

		if (line.length() > 30) {
			gxyz[(COLS_gxyz) * (temp) + 1] = std::stod(line.substr(28, 20)) + offset_Y;
			if (line.length() > 50) {
				gxyz[(COLS_gxyz) * (temp) + 2] = std::stod(line.substr(48, 20)) + offset_Z;
			}
			else {
				gxyz[(COLS_gxyz) * (temp) + 2] = 0.0 + offset_Z;
			}
		}
		else {
			gxyz[(COLS_gxyz) * (temp) + 1] = 0.0 + offset_Y;
			gxyz[(COLS_gxyz) * (temp) + 2] = 0.0 + offset_Z;
		}
		//}
		index++;
	}

	m = index; //total number of nodes used to connect global elems
	inFile1.close();
	Log("file_gnode read");


	//file_lelem, local element file
	inFile1.open(file_lelem);
	check(inFile1);

	index = 0;
	while (getline(inFile1, line)) {
		le.push_back(std::stoi(line.substr(78, 6)));
		le_n[(COLS_le_n) * index + 0] = std::stoi(line.substr(0, 6));
		le_n[(COLS_le_n) * index + 1] = std::stoi(line.substr(6, 6));
		le_n[(COLS_le_n) * index + 2] = std::stoi(line.substr(12, 6));
		le_n[(COLS_le_n) * index + 3] = std::stoi(line.substr(18, 6));
		index++;
	}
	o = index;
	inFile1.close();
	Log("file_lelem read");


	//file_lnode, local node file
	inFile1.open(file_lnode);
	check(inFile1);

	index = 0;
	while (getline(inFile1, line)) {
		//temp = node number;
		int temp = std::stoi(line.substr(0, 8));
		ln.push_back(temp);
		lxyz[(COLS_lxyz) * (temp) + 0] = std::stod(line.substr(8, 20));

		if (line.length() > 30) {
			lxyz[(COLS_lxyz) * (temp) + 1] = std::stod(line.substr(28, 20));
			if (line.length() > 50) {
				lxyz[(COLS_lxyz) * (temp) + 2] = std::stod(line.substr(48, 20));
			}
			else {
				lxyz[(COLS_lxyz) * (temp) + 2] = 0.0;
			}
		}
		else {
			lxyz[(COLS_lxyz) * (temp) + 1] = 0.0;
			lxyz[(COLS_lxyz) * (temp) + 2] = 0.0;
		}
		index++;
	}

	p = index; //total number of nodes used to connect global elems
	inFile1.close();
	Log("file_lnode read");

}

void GLmapper::mapping() {
	int i, j, k, kk, jj, ii;
	const int temparray_SIZE = 999999;
	double xx[4];
	double yy[4];
	double zz[4];
	double point[3];
	double average[3];
	double *temp_lxyz_avg = new double[temparray_SIZE * 5];
	int inside;

	//find center of local elem
	for (int i = 0; i < o; i++) { // o is number of local elements
		lxyz_avg[COLS_lxyz_avg * i + 0] = 0;
		lxyz_avg[COLS_lxyz_avg * i + 1] = 0;
		lxyz_avg[COLS_lxyz_avg * i + 2] = 0;

		if (le_n[COLS_le_n * i + 2] == le_n[COLS_le_n * i + 3]) // if N3=N4, element is triangle
			kk = 3;
		else
			kk = 4;

		for (int k = 0; k < kk; k++) {
			j = le_n[COLS_le_n * i + k]; //node number in local element file is mapped to an index in the node coordinates array
			lxyz_avg[COLS_lxyz_avg * i + 0] = lxyz_avg[COLS_lxyz_avg * i + 0] + (lxyz[COLS_lxyz * j + 0] / kk);
			lxyz_avg[COLS_lxyz_avg * i + 1] = lxyz_avg[COLS_lxyz_avg * i + 1] + (lxyz[COLS_lxyz * j + 1] / kk);
			lxyz_avg[COLS_lxyz_avg * i + 2] = lxyz_avg[COLS_lxyz_avg * i + 2] + (lxyz[COLS_lxyz * j + 2] / kk);
		}
		lxyz_avg[COLS_lxyz_avg * i + 3] = le[i];
	}
	//find center of global elem
	for (int i = 0; i < N; i++) { //N is number of global elements
		gxyz_avg[COLS_gxyz_avg * i + 0] = 0;
		gxyz_avg[COLS_gxyz_avg * i + 1] = 0;
		gxyz_avg[COLS_gxyz_avg * i + 2] = 0;

		if (ge_n[COLS_ge_n * i + 2] == ge_n[COLS_ge_n * i + 3]) // if N3=N4, element is triangle
			kk = 3;
		else
			kk = 4;

		for (int k = 0; k < kk; k++) {
			j = ge_n[COLS_ge_n * i + k]; //node number in local element file is mapped to an index in the node coordinates array
			gxyz_avg[COLS_gxyz_avg * i + 0] = gxyz_avg[COLS_gxyz_avg * i + 0] + (gxyz[COLS_gxyz * j + 0] / kk);
			gxyz_avg[COLS_gxyz_avg * i + 1] = gxyz_avg[COLS_gxyz_avg * i + 1] + (gxyz[COLS_gxyz * j + 1] / kk);
			gxyz_avg[COLS_gxyz_avg * i + 2] = gxyz_avg[COLS_gxyz_avg * i + 2] + (gxyz[COLS_gxyz * j + 2] / kk);
		}
		gxyz_avg[COLS_gxyz_avg * i + 3] = ge_n[COLS_ge_n * i + 4]; //element number

	}
	
	//which global element contains the local element
	for (int ilocal = 0; ilocal < o; ilocal++) {
		point[0] = lxyz_avg[COLS_lxyz_avg * ilocal + 0];
		point[1] = lxyz_avg[COLS_lxyz_avg * ilocal + 1];
		point[2] = lxyz_avg[COLS_lxyz_avg * ilocal + 2];

		//vector calculation for later?

		int corresponding_global_elem = minDistance(point); //add temparray as parameter later?

		lpressure[COLS_lpressure * ilocal + 0] = std::to_string(le[ilocal]);
		lpressure[COLS_lpressure * ilocal + 1] = std::to_string(corresponding_global_elem);
		
	}

	delete[] temp_lxyz_avg;
}

int GLmapper::minDistance(double local_point[]) {
	int output = -1;
	double min = 999999999.0;
	double global_point[3];
	for (int iglobal = 0; iglobal < N; iglobal++) {
		global_point[0] = gxyz_avg[COLS_gxyz_avg * iglobal + 0];
		global_point[1] = gxyz_avg[COLS_gxyz_avg * iglobal + 1];
		global_point[2] = gxyz_avg[COLS_gxyz_avg * iglobal + 2];

		//calculate distance
		if (std::abs(local_point[0] - global_point[0]) < radius) {
			if (std::abs(local_point[1] - global_point[1]) < radius) {
				if (std::abs(local_point[2] - global_point[2]) < radius) {
					double delta_x = std::abs(local_point[0] - global_point[0]);
					double delta_y = std::abs(local_point[1] - global_point[1]);
					double delta_z = std::abs(local_point[2] - global_point[2]);
					double distance = std::abs(std::sqrt(std::pow(delta_x, 2.0) + std::pow(delta_y, 2.0) + std::pow(delta_z, 2.0)));
					if (distance < min) {
						min = distance;
						output = ge_n[COLS_ge_n * iglobal + 4];
					}
				}
			}
		}
	}
	return output;
}

void GLmapper::pressure_cycle(int cycle) {
	//const int temparray_SIZE = 999999;
	//double* temparray = new double[temparray_SIZE];
	int i, j, k, ii, kk;
	std::ifstream inFile2;

	String line;
	int index = 0;

	inFile2.open(file_gpressure);
	check(inFile2);
	String temp = "GLOC";
	while (getline(inFile2, line)&&index<HEADER_LENGTH) {
		if (line.substr(0, 5).compare(temp) == 0) {
			outoffset_x = std::stod(line.substr(6, 13)) + offset_X;
			outoffset_y = std::stod(line.substr(20, 13)) + offset_Y;
			outoffset_z = std::stod(line.substr(34, 13)) + offset_Z;
			index = 20;
		}
		index++;
	}
	inFile2.close();

	//read data from global pressure file
	inFile2.open(file_gpressure);
	ii = 0;
	k = 0;
	const int startingPos = 25;
	temp = "SFE,";
	String delimiter = ",";
	while (getline(inFile2, line)) {
		if (line.substr(0, 4).compare(temp) == 0) {
			double pressures[4];
			int gelement_index = std::stoi(line.substr(4, 6));

			line = line.substr(startingPos);
			for (int i = 0; i < 4; i++) {
				int newPos = line.find(delimiter);
				pressures[i] = std::stod(line.substr(0, newPos));
				if(i<3)
					line = line.substr(newPos + delimiter.length());
			}

			double average = (pressures[0] + pressures[1] + pressures[2] + pressures[3]) / 4;
			gpressure_mac[gelement_index] = average;
			ii++;
		}
		k++;
	}

	q = ii;
	inFile2.close();
		
	//global element number as index in array
	//int find1 = 0;
	//String empty = "";
	//for (int ilocal = 0; ilocal < o; ilocal++) {
	//	if (lpressure[COLS_lpressure * ilocal + 1].compare(empty)!=0) {
	//		double temp = gpressure_mac[std::stoi(lpressure[COLS_lpressure * ilocal + 1])];
	//		lpressure[COLS_lpressure * ilocal + 2] = (temp == -1.0) ? "Null" : temp;
	//	}
	//	else {
	//		lpressure[COLS_lpressure * ilocal + 2] = "Null";
	//	}
	//}

	//delete[] temparray;
}
void GLmapper::output(int cycle, bool found){
	int i, j, k;
	int output;
	output = cycle + 150;

	std::ofstream outFile;
	outFile.open(file_lpressure);

	if (!found) {
		std::ofstream mapFile;
		mapFile.open(MAPPINGFILE);
		for (int i = 0; i < o; i++) {
			mapFile << lpressure[COLS_lpressure * i + 0] << "->" << lpressure[COLS_lpressure * i + 1] << std::endl;
		}
		mapFile.close();
	}


	//outFile << " !Dynamic Imaginary Pressures for T=12.5000, Heading=0.000" << std::endl;
	//outFile << "!Units: Metric" << std::endl;
	//outFile << "!Wave amplitude (m):1" << std::endl;
	//outFile << "!Gravity (mm/sec^2) : 9810" << std::endl;
	//outFile << "!Water Density (kg/mm^3) : 1.025e-9" << std::endl;
	//outFile << "/solu" << std::endl;
	//outFile << "ACEL,0.00000e+000,0.00000e+000,0.00000e+000" << std::endl;
	//outFile << "FDELE,ALL,ALL" << std::endl;
	//outFile << "SFEDELE,ALL,ALL,ALL" << std::endl;
	//outFile << "CSYS,0" << std::endl;
	//outFile << "ACEL,1.398e+002,3.664e-003,-6.26621e-002," << std::endl;
	//outFile << "CGLOC," << outoffset_x << "," << outoffset_y << "," << outoffset_z << std::endl;
	//outFile << "DCGOMG,-2.338e-007,1.626e-004,-1.89675e-005," << std::endl;
	//outFile << "ALLSEL $CSYS,0" << std::endl;
	//outFile << "! Element Pressures " << std::endl;

	outFile << "\\solu" << std::endl;
	outFile << "! Element Pressures " << std::endl;
	
	for (int i = 0; i < o; i++) {
		String str = "Null";
		double result = gpressure_mac[std::stoi(lpressure[COLS_lpressure * i + 1])];
		if (result == 999.0) {
			outFile << "!SFE, " << lpressure[COLS_lpressure * i + 0] << " ,2,PRES,0, " << "Null" << std::endl;
		}
		else {
			outFile << "SFE, " << lpressure[COLS_lpressure * i + 0] << " ,2,PRES,0, " << std::scientific << std::setprecision(4) << result << std::endl;
		}
	}

	outFile.close();

}

void GLmapper::check(std::ifstream& inFile) {
	if (inFile.fail()) {
		std::cerr << "Error Opening File" << std::endl;
		std::cin.get();
		exit(1);
	}
}


