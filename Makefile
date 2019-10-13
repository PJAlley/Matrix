CXX = clang++
CFLAGS = -std=c++17 -Wall -g
CCFLAGS = -c -std=c++17

all: matrixTest

matrixTest: Matrix.o SquareMatrix.o MatrixRow.o MatrixColumn.o RowVector.o ColumnVector.o QRMatrix.o LUMatrix.o LUPMatrix.o Eigen.o Hessen.o SVD.o Schur.o MatrixTest.o
	$(CXX) $(CFLAGS) Matrix.o SquareMatrix.o MatrixRow.o MatrixColumn.o RowVector.o ColumnVector.o QRMatrix.o LUMatrix.o LUPMatrix.o Eigen.o Hessen.o SVD.o Schur.o MatrixTest.o -o matrixTest.exe

charPolyTest: Matrix.o SquareMatrix.o MatrixRow.o MatrixColumn.o LUMatrix.o RowVector.o ColumnVector.o Schur.o charPolyTest.o
	$(CXX) $(CFLAGS) Matrix.o SquareMatrix.o MatrixRow.o MatrixColumn.o LUMatrix.o RowVector.o ColumnVector.o Schur.o charPolyTest.o -o charPolyTest.exe

vectorTest: Matrix.o ColumnVector.o MatrixRow.o RowVector.o vectorTest.o
	$(CXX) $(CFLAGS) Matrix.o ColumnVector.o MatrixRow.o RowVector.o vectorTest.o -o vectorTest.exe

eigenTest: Matrix.o SquareMatrix.o MatrixRow.o RowVector.o ColumnVector.o LUMatrix.o Eigen.o QRMatrix.o Hessen.o Schur.o eigenTest.o
	$(CXX) $(CFLAGS) Matrix.o SquareMatrix.o MatrixRow.o RowVector.o ColumnVector.o LUMatrix.o QRMatrix.o Eigen.o Hessen.o Schur.o eigenTest.o -o eigenTest.exe

decompTest: Matrix.o SquareMatrix.o MatrixRow.o RowVector.o ColumnVector.o QRMatrix.o LUMatrix.o SVD.o decompTest.o
	$(CXX) $(CFLAGS) Matrix.o SquareMatrix.o MatrixRow.o RowVector.o ColumnVector.o QRMatrix.o LUMatrix.o SVD.o decompTest.o -o decompTest.exe

lupTest: Matrix.o SquareMatrix.o MatrixRow.o ColumnVector.o RowVector.o LUPMatrix.o LUMatrix.o lupTest.o
	$(CXX) $(CFLAGS) Matrix.o SquareMatrix.o MatrixRow.o ColumnVector.o RowVector.o LUPMatrix.o LUMatrix.o lupTest.o -o lupTest.exe

luTest: Matrix.o SquareMatrix.o MatrixRow.o ColumnVector.o RowVector.o LUMatrix.o LUTest.o
	$(CXX) $(CFLAGS) Matrix.o SquareMatrix.o MatrixRow.o ColumnVector.o RowVector.o LUMatrix.o LUTest.o -o luTest.exe

svdTest: Matrix.o SquareMatrix.o MatrixRow.o ColumnVector.o RowVector.o LUMatrix.o SVD.o SVDTest.o
	$(CXX) $(CFLAGS) Matrix.o SquareMatrix.o MatrixRow.o ColumnVector.o RowVector.o LUMatrix.o SVD.o SVDTest.o -o svdTest.exe

qrTest: Matrix.o SquareMatrix.o MatrixRow.o ColumnVector.o RowVector.o QRMatrix.o QRTest.o
	$(CXX) $(CFLAGS) Matrix.o SquareMatrix.o MatrixRow.o ColumnVector.o RowVector.o LUMatrix.o QRMatrix.o QRTest.o -o qrTest.exe

rrefTest: Matrix.o SquareMatrix.o MatrixRow.o MatrixColumn.o ColumnVector.o RowVector.o rrefTest.o
	$(CXX) $(CFLAGS) Matrix.o SquareMatrix.o MatrixRow.o MatrixColumn.o ColumnVector.o RowVector.o rrefTest.o -o rrefTest.exe

rrefTest.o: rrefTest.cpp
	$(CXX) $(CCFLAGS) rrefTest.cpp

vectorTest.o: vectorTest.cpp
	$(CXX) $(CCFLAGS) vectorTest.cpp

charPolyTest.o: charPolyTest.cpp
	$(CXX) $(CCFLAGS) charPolyTest.cpp

MatrixTest.o: MatrixTest.cpp
	$(CXX) $(CCFLAGS) MatrixTest.cpp

eigenTest.o: eigenTest.cpp
	$(CXX) $(CCFLAGS) eigenTest.cpp

decompTest.o: decompTest.cpp
	$(CXX) $(CCFLAGS) decompTest.cpp

lupTest.o: lupTest.cpp
	$(CXX) $(CCFLAGS) lupTest.cpp

LUTest.o: LUTest.cpp
	$(CXX) $(CCFLAGS) LUTest.cpp

SVDTest.o: SVDTest.cpp
	$(CXX) $(CCFLAGS) SVDTest.cpp

QRTest.o: QRTest.cpp
	$(CXX) $(CCFLAGS) QRTest.cpp

Schur.o: Schur.cpp
	$(CXX) $(CCFLAGS) Schur.cpp

Hessen.o: Hessen.cpp
	$(CXX) $(CCFLAGS) Hessen.cpp

Eigen.o: Eigen.cpp
	$(CXX) $(CCFLAGS) Eigen.cpp

QRMatrix.o: QRMatrix.cpp
	$(CXX) $(CCFLAGS) QRMatrix.cpp

LUPMatrix.o: LUPMatrix.cpp
	$(CXX) $(CCFLAGS) LUPMatrix.cpp

LUMatrix.o: LUMatrix.cpp
	$(CXX) $(CCFLAGS) LUMatrix.cpp

SVD.o: SVD.cpp
	$(CXX) $(CCFLAGS) SVD.cpp

ColumnVector.o: ColumnVector.cpp
	$(CXX) $(CCFLAGS) ColumnVector.cpp

RowVector.o: RowVector.cpp
	$(CXX) $(CCFLAGS) RowVector.cpp

MatrixColumn.o: MatrixColumn.cpp
	$(CXX) $(CCFLAGS) MatrixColumn.cpp

MatrixRow.o: MatrixRow.cpp
	$(CXX) $(CCFLAGS) MatrixRow.cpp

SquareMatrix.o: SquareMatrix.cpp
	$(CXX) $(CCFLAGS) SquareMatrix.cpp

Matrix.o: Matrix.cpp
	$(CXX) $(CCFLAGS) Matrix.cpp

clean:
	rm *.o