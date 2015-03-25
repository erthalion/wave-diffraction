#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

int main(int argc,char **argv){
	
	if(argc!=3){
		cout<<"usage: process tec_file.dat grid_size\n"<<endl;
		return 1;
	}
	
	string tec_file,gnu_file,cur_line;
	int size=atoi(argv[2]),counter=0;
	ifstream tec;
	ofstream gnu;

	tec_file=argv[1];
	gnu_file="converted_"+tec_file;
	
	cout<<tec_file<<endl;
	cout<<gnu_file<<endl;
	
	tec.open(tec_file.data());
	gnu.open(gnu_file.data());

	if(!tec.is_open())
		cout<<"couldn't open tec file!\n"<<endl;
	if(!gnu.is_open())
		cout<<"couldn't open gnu file!\n"<<endl;

	while(tec.good()){
		getline(tec,cur_line);

		if(cur_line.substr(0,4)!="ZONE"){
			counter++;
			gnu<<cur_line<<endl;
			if(counter%size==0){
				gnu<<"\n";
			}

			if(counter%(size*size)==0){
				gnu<<"\n";
			}
		}
	}

	tec.close();
	gnu.close();

	return 0;
}
