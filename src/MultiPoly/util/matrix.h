/*
 * matrix.h
 *
 *  Created on: 14/04/2016
 *      Author: Milder
 */

#ifndef UTIL_MATRIX_H_
#define UTIL_MATRIX_H_

#include <vector>
#include <iostream>

namespace irtpp {

/**
 *	Matrix class.
 *	Matrix class used to store information.
 * */

template<class T>
class matrix {
	private:
		std::vector<std::vector<T> > data;

	public:
		matrix();
		/**
		 * Creates a new matrix with size rows * cols
		 * */
		matrix(int rows, int cols);

		/**
		 * Creates a new matrix specifying only the number of rows
		 * rows can be added after using add_row(std::vector) method
		 * */
		matrix(int cols);

		/**
		 * Adds a new row to the current matrix
		 * */
		void add_row(std::vector<T> *new_row);

		/**
		 * Adds a new row to the current matrix
		 * */
		void add_row(std::vector<T> new_row);
		void add_row(T* new_row, int size);
		/**
		 * Adds an empty row
		 * */
		void add_row();

		/**
		 * Adds a empty row with specified size
		 * */
		void add_row(int);

		/**
		 * Clean the number of the specific row
		 * */
		void reset_row(int);

		void reset();

		/**
		 * Adds an element to the last row
		 * */
		void add_element(T e);
		virtual ~matrix();

		/**
		 * Returns the data
		 * */
		std::vector<std::vector<T> > get_data();

		/**
		 * Returns the i-th row in the matrix
		 * */
		std::vector<T> get_row(int i);

		/**
		 * Returns the i-th row in the matrix
		 * */
		std::vector<T>* get_pointer_row(int i);

		/**
		 * Returns number of rows of the matrix
		 * */
		int rows();

		/**
		 * Returns number of columns of a specific row
		 * */
		int columns(int row);

		/**
		 *  Accessing operator for an element
		 * */
		inline T &operator()(const int row, const int col);
};

template<class T>
matrix<T>::matrix(int rows) {
	data = std::vector<std::vector<T> >(rows);
}

template<class T>
matrix<T>::matrix() {
}

template<class T>
matrix<T>::~matrix() {
}

template<class T>
matrix<T>::matrix(int rows, int cols) {
	data = std::vector<std::vector<T> >(rows, std::vector<T>(cols));
}

template<class T>
void matrix<T>::add_row(std::vector<T> *new_row) {
	data.push_back(*new_row);
}

template<class T>
void matrix<T>::add_row(std::vector<T> new_row) {
	data.push_back(new_row);
}

template<class T>
void matrix<T>::add_row(T* new_row, int size) {
	data.push_back(std::vector<T>(new_row, new_row + size));
}

template<class T>
void matrix<T>::add_row() {
	data.push_back(std::vector<T>());
}

template<class T>
void matrix<T>::reset() {
	for ( int i = 0; i < data.size(); ++i )
		for ( int j = 0; j < data[i].size(); ++j )
			data[i][j] = 0;
}


template<class T>
void matrix<T>::add_row(int size) {
	data.push_back(std::vector<T>(size));
}

template<class T>
void matrix<T>::reset_row(int idx) {
	for ( unsigned int i = 0; i < data[idx].size(); ++i )
		data[idx][i] = 0;
}

template<class T>
void matrix<T>::add_element(T e) {
	data.back().push_back(e);
}

template<class T>
int matrix<T>::rows() {
	return data.size();
}

template<class T>
int matrix<T>::columns(int row) {
	if ( (unsigned)row >= data.size() )
		return -1;
	return data[row].size();
}

template<class T>
std::vector<std::vector<T> > matrix<T>::get_data() {
	return data;
}

template<class T>
std::vector<T> matrix<T>::get_row(int i) {
	return data[i];
}

template<class T>
std::vector<T>* matrix<T>::get_pointer_row(int i) {
	return &data[i];
}

template<class T>
inline T &matrix<T>::operator()(const int r, const int c) {
	return data.at(r).at(c);
}

template<class T>
std::ostream& operator<<(std::ostream &out, matrix<T> &m) {
	for (int i = 0; i < m.rows(); ++i) {
		for (int j = 0; j < m.columns(i); j++)
			out << (int)m(i, j) << '\t';
		out << '\n';
	}
	return (out);
}

} /* namespace irtpp */

#endif /* UTIL_MATRIX_H_ */
