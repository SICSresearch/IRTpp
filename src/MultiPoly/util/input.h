/*
 * input.h
 *
 *  Created on: 14/04/2016
 *      Author: Milder
 */

#ifndef UTIL_INPUT_H_
#define UTIL_INPUT_H_

#include <fstream>
#include <cstring>
#include <string>
#include <cstdlib>
#include <iostream>
#include <algorithm>
#include "matrix.h"

namespace irtpp {

template<class T>
class input {
	private:
		/**
		 * Delimiter to split each file line
		 */
		char delimiter;

	public:

		input();
		input(char);
		virtual ~input();

		/**
		 * Imports matrixes from a csv
		 * */
		bool importData(std::string filename, matrix<T>& m);
		bool importData(std::string filename, std::vector<T>& m);

		/**
		 * Gets the delimiter used for inputting
		 * */
		char get_delimiter() const;

		/**
		 * Sets the delimiter for inputting text matrices
		 * */
		void set_delimiter(char);
};

template<class T>
input<T>::input() {
	delimiter = ',';
}

template<class T>
input<T>::input(char delimiter) {
	this->delimiter = delimiter;
}

template<class T>
input<T>::~input() {
}

template<class T>
bool input<T>::importData(std::string filename, matrix<T>& m) {
	std::string line;
	std::ifstream file(filename.c_str());
	if (file.is_open()) {
		while (std::getline(file, line)) {
			std::vector<T> splitted;
			std::string::const_iterator start = line.begin();
			std::string::const_iterator end = line.end();
			std::string::const_iterator next = std::find(start, end, delimiter);
			while (next != end) {
				std::string to_add(start, next);
				if ( !to_add.empty() )
					splitted.push_back(strtold(to_add.c_str(), NULL));
				start = next + 1;
				next = std::find(start, end, delimiter);
			}
			std::string to_add(start, next);
			if ( !to_add.empty() )
				splitted.push_back(strtold(to_add.c_str(), NULL));
			m.add_row(&splitted);
		}
		file.close();
		return true;
	}
	return false;
}


template<class T>
bool input<T>::importData(std::string filename, std::vector<T>& m) {
	std::string line;
	std::ifstream file(filename.c_str());
	if (file.is_open()) {
		while (std::getline(file, line)) {
			std::vector<T> splitted;
			std::string::const_iterator start = line.begin();
			std::string::const_iterator end = line.end();
			std::string::const_iterator next = std::find(start, end, delimiter);
			while (next != end) {
				std::string to_add(start, next);
				if ( !to_add.empty() )
					splitted.push_back(strtold(to_add.c_str(), NULL));
				start = next + 1;
				next = std::find(start, end, delimiter);
			}
			std::string to_add(start, next);
			if ( !to_add.empty() )
				splitted.push_back(strtold(to_add.c_str(), NULL));
			m.insert(m.end(), splitted.begin(), splitted.end());
		}
		file.close();
		return true;
	}
	return false;
}

template<class T>
char input<T>::get_delimiter() const {
	return delimiter;
}

template<class T>
void input<T>::set_delimiter(char del) {
	this->delimiter = delimiter;
}

} /* namespace irtpp */

#endif /* UTIL_INPUT_H_ */
