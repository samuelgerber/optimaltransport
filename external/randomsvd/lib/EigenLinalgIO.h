//templated C-Style IO functions
#ifndef LINALGIO_H
#define LINALGIO_H

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <list>
#include <stdio.h>
#include <stdlib.h>

#include <Eigen/Dense>

namespace EigenLinalg{

template <typename TPrecision>
class LinalgIO{


  public:

  
  //Read Vector from header file
  static Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> readVector(const std::string &filename){
    std::ifstream hdr;
    hdr.open(filename.c_str());
    
    std::string token;
    getline(hdr, token);
    if(token.compare("DenseVector") != 0){
      throw "Not a vector header file";
    }

    getline(hdr, token, ' ');
    getline(hdr, token);
    int size = atoi(token.c_str());

    getline(hdr, token, ' ');
    getline(hdr, token);
    int elemSize = atoi(token.c_str());
    if(elemSize != sizeof(TPrecision)){
      throw "element size not equal to size of Vector precision";
    }

    getline(hdr, token, ' ');
    getline(hdr, token);

    Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> vector(size);
    size_t start = filename.find_last_of("/\\");
    std::stringstream ss;
    if(start != std::string::npos){
      ss << filename.substr(0, start+1);
    }
    ss << token;


    readVector(ss.str(), vector);

    return vector;
  };



  //Reads a  N Vector from binary file
  static bool readVector(const std::string &filename, Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> &vector){
    std::ifstream file;
    file.open(filename.c_str(), std::ios::binary);
    file.read((char*)vector.data(), sizeof(TPrecision)*vector.size());

    if(file.fail()){
      std::cout << "Reading failed" << std::endl;
      // get length of file:
      file.seekg (0, std::ios::beg);
      file.seekg (0, std::ios::end);
      long length = file.tellg();
      file.seekg (0, std::ios::beg);
      std::cout << length << std::endl;
      return false;  
    }  
    file.close();
    return true;
  };



  static void writeVector(const std::string &filename, 
      Eigen::Ref< Eigen::Matrix<TPrecision, Eigen::Dynamic, 1> > vector,
      bool writeHeader = true){
    
    std::ofstream file;
    file.open(filename.c_str(), std::ios::binary);
    char *data = (char*)vector.data();
    file.write(data, vector.size()*sizeof(TPrecision));
    file.close();

    if(writeHeader){
      size_t start = filename.find_last_of("/\\");
      std::string localFile;
      if(start != std::string::npos){
        localFile = filename.substr(start+1);
      }
      else{
        localFile = filename;
      };
      std::stringstream ss;
      ss << filename << ".hdr";

      std::ofstream hdr;
      hdr.open(ss.str().c_str());
      hdr << "DenseVector" << std::endl;
      hdr << "Size: " << vector.size() << std::endl;
      hdr << "ElementSize: " << sizeof(TPrecision) << std::endl;
      hdr << "DataFile: " << localFile << std::endl;
      hdr.close();
    }
  };



  //Read matrix from header file
  static Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> readMatrix(const std::string &filename){
    std::ifstream hdr;
    hdr.open(filename.c_str());
    
    std::string token;
    getline(hdr, token);
    if(token.compare("DenseMatrix") != 0){
      throw "Not a matrix header file";
    }

    getline(hdr, token, ' ');
    getline(hdr, token, ' ');
    int m = atoi(token.c_str());
    getline(hdr, token, ' ');
    getline(hdr, token);
    int n = atoi(token.c_str());

    getline(hdr, token, ' ');
    getline(hdr, token);
    int elemSize = atoi(token.c_str());
    if(elemSize != sizeof(TPrecision)){
      throw "Element size not equal to size of matrix precision";
    }

    getline(hdr, token, ' ');
    getline(hdr, token);
    bool rowMajor = atoi(token.c_str()) != 0;

    getline(hdr, token, ' ');
    getline(hdr, token);

    Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> matrix(m, n);

    size_t start = filename.find_last_of("/\\");
    std::stringstream ss;
    if(start != std::string::npos ){
      ss << filename.substr(0, start+1);
    }
    ss << token;

    readMatrix(ss.str(), matrix, rowMajor);

    return matrix;
  }



  //Reads a M x N Matrix from binary file
  static bool readMatrix(const std::string &filename, 
     Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic>  &matrix,
      bool rowMajor = false){
    std::ifstream file;
    file.open(filename.c_str(), std::ios::binary);
    file.read((char*)matrix.data(), sizeof(TPrecision)*matrix.rows()*matrix.cols());

    if(file.fail()){
      std::cout << "Reading binary matrix file failed" << std::endl;
      std::cout << matrix.rows() << " x " << matrix.cols() << std::endl;
      std::cout << sizeof(TPrecision) << std::endl;
      // get length of file
      file.seekg (0, std::ios::beg);
      file.seekg (0, std::ios::end);
      long length = file.tellg();
      file.seekg (0, std::ios::beg);
      std::cout << length << std::endl;
      return false;     
    }  
    
    file.close();

    //convert to column major if necessary
    if(rowMajor){
    }

    return true;
  };




  static void writeMatrix(const std::string &filename, 
      Eigen::Ref< Eigen::Matrix<TPrecision, Eigen::Dynamic, Eigen::Dynamic> >  matrix,
      bool writeHeader = true){
    std::ofstream file;
    file.open(filename.c_str(), std::ios::binary);
    file.write((char *)matrix.data(), matrix.rows()*matrix.cols()*sizeof(TPrecision));
    file.close();
    
    if(writeHeader){
      size_t start = filename.find_last_of("/\\");
      std::string localFile;
      if(start ==std::string::npos){
        localFile = filename;
      }
      else{
        localFile = filename.substr(start+1);
      }
      
      std::stringstream ss;
      ss << filename << ".hdr";

      std::ofstream hdr;
      hdr.open(ss.str().c_str());

      hdr << "DenseMatrix" << std::endl;
      hdr << "Size: " << matrix.rows() << " x " << matrix.cols() << std::endl;
      hdr << "ElementSize: " << sizeof(TPrecision) << std::endl;
      hdr << "RowMajor: " << false << std::endl;
      hdr << "DataFile: " << localFile << std::endl;
      hdr.close();
    }
  };


};

}

#endif
