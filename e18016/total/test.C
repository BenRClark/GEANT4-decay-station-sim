#include <iostream>
#include <vector>
#include <array>
#include "Riostream.h"


//defining a truncator function
bool truncator(int i){
if(i == 2){
    return false;
  }
  else{
    return true;
  }
}

//mapping proof of concept
void test(){

  double origin[4][2] = {0,1,2,3,4,5,6,7};
  int map[4] = {3,1,2,0};
  double aligned[4][2] = {};
  
  int mapper;
  for(int i = 0; i < 4; i++){
    for(int j = 0; j < 2; j++){
	mapper = map[i];
	//cout << mapper << endl;
	aligned[i][j] = origin[mapper][j];
	cout << "aligned[" << i << "][" << j << "]: " << aligned[i][j] << " " << "origin[" << i << "][" << j << "]: " << origin[i][j] << endl;
    }
  }


  //testing to see how to truncate two dimnesional array
  //let's say we want to remove the 2 'i' index from the origin array
 
  int mapping[3] = {0,1,3};
  int mapper1 = 0;
  int truncated[3][2] = {};
  int truncated1[3][2] = {};
  //method 1
  for(int i = 0; i < 2; i++){
    for(int j = 0; j < 4; j++){
      mapper1 = mapping[j];
      truncated[j][i] = origin[mapper1][i];

    }
  }
  //method 2
  for(int i = 0; i < 2; i++){
    int counter1 = 0;
    for(int j = 0; j < 4; j++){
      if(truncator(j)){
      truncated1[counter1][i] = origin[j][i];
      counter1++;
      }
    }
  }

  //printing out truncated
  for(int i = 0; i < 3; i++){
    for(int j = 0; j < 2; j++){
      cout << "truncated[" << i << "][" << j << "]: " << truncated[i][j] << endl;
      cout << "truncated1[" << i << "][" << j << "]: " << truncated1[i][j] << endl;
    }
  }


}


